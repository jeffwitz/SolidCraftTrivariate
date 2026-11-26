#!/usr/bin/env node

/**
 * Browser smoke test for the static engineering gallery.
 * Uses only Node built-ins and a local Chrome/Chromium executable.
 *
 * Usage:
 *   CHROME_BIN=/usr/bin/google-chrome node tests/references/gallery_browser.mjs [output-dir]
 *   node tests/references/gallery_browser.mjs --check-expressions
 */

import { spawn } from 'node:child_process';
import { createServer } from 'node:http';
import { createHash } from 'node:crypto';
import { inflateSync } from 'node:zlib';
import { mkdtemp, mkdir, readFile, rm, writeFile } from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const REPO_ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const SITE_ROOT = path.resolve(process.env.GALLERY_SITE_ROOT || path.join(REPO_ROOT, 'docs/site'));
const OUTPUT_DIR = path.resolve(process.argv[2] || path.join(REPO_ROOT, 'build/gallery-qa'));
const CHROME_BIN = process.env.CHROME_BIN || '/usr/bin/google-chrome';
const EXPECTED_ENTRY_COUNT = Number(process.env.GALLERY_EXPECTED_ENTRY_COUNT || 68);
const HOST = '127.0.0.1';
const MIME_TYPES = {
  '.css': 'text/css; charset=utf-8',
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.step': 'application/step',
};

const injectedExpressions = {
  galleryReady: `(() => {
    const status = document.querySelector('#canvas-state');
    const select = document.querySelector('#experiment-select');
    return Boolean(status?.hidden && select?.options.length > 0 && !select.disabled);
  })()`,
  bodySummary: 'document.body.innerText.slice(0, 1600)',
  publicUiLocale: `(() => ({
    lang: document.documentElement.lang,
    heading: document.querySelector('#page-title').textContent,
    canvasLabel: document.querySelector('#model-canvas').getAttribute('aria-label'),
    resetLabel: document.querySelector('#reset-view').textContent.trim(),
    experimentLabel: document.querySelector('label[for=experiment-select]').textContent,
  }))()`,
  canvasRegion: `(() => {
    const rect = document.querySelector('#model-canvas').getBoundingClientRect();
    return { x: rect.x, y: rect.y, width: rect.width, height: rect.height, scale: 1 };
  })()`,
  resources: `fetch(document.getElementById('catalog-url').href).then(response => response.json()).then(async experiments => {
    const comparisonResponse = await fetch('./comparison-data.json');
    const comparisonData = comparisonResponse.ok ? await comparisonResponse.json() : [];
    const comparisonIds = new Set(comparisonData.map(item => item.id));
    return Promise.all(experiments.map(async experiment => {
      const assetExtensions = ['.json', '.step', '.vtk', '.vtu', '.stl', '.obj', '.glb', '.gltf', '.png', '.jpg', '.jpeg'];
      const assets = Object.entries(experiment)
        .filter(([key, value]) => typeof value === 'string'
          && assetExtensions.some(extension => value.split('?')[0].split('#')[0].toLowerCase().endsWith(extension))
          && !value.startsWith('http://') && !value.startsWith('https://'))
        .map(([key, value]) => ({ key, path: value }));
    const modelResponse = await fetch(experiment.model);
    const modelBytes = await modelResponse.arrayBuffer();
    const model = JSON.parse(new TextDecoder().decode(modelBytes));
    const modelDigest = await crypto.subtle.digest('SHA-256', modelBytes);
    const modelSha256 = [...new Uint8Array(modelDigest)].map(value => value.toString(16).padStart(2, '0')).join('');
    const auditResponse = experiment.audit ? await fetch(experiment.audit) : null;
    const audit = auditResponse?.ok ? await auditResponse.json() : null;
    const representationResponse = experiment.representation ? await fetch(experiment.representation) : null;
    const representation = representationResponse?.ok ? await representationResponse.json() : null;
    const exactPatchesValid = !representation || representation.patches.every(patch =>
      Array.isArray(patch.degrees) && patch.degrees.length === 3
      && Array.isArray(patch.knots) && patch.knots.length === 3
      && Array.isArray(patch.control_shape) && patch.control_shape.length === 3
      && patch.knots.every((knots, axis) => Array.isArray(knots)
        && knots.length === patch.control_shape[axis] + patch.degrees[axis] + 1
        && knots.every((value, index) => Number.isFinite(value) && (index === 0 || value >= knots[index - 1])))
      && Array.isArray(patch.controls_homogeneous)
      && patch.controls_homogeneous.length === patch.control_shape.reduce((product, value) => product * value, 1)
      && patch.controls_homogeneous.every(control => Array.isArray(control) && control.length === 4
        && control.every(Number.isFinite) && control[3] > 0));
    const exactExteriorFaceKeys = representation?.boundary_faces
      ? representation.boundary_faces.map(face => (face.source_block || face.block) + ':' + face.side).sort()
      : (model.boundary_faces || model.topology?.exterior || []).map(face => Array.isArray(face)
        ? face[0] + ':' + face[1] : (face.source_block || face.block) + ':' + face.side).sort();
    const curvedBlocks = (model.blocks || []).filter(block => {
      const degrees = block.degrees || model.degrees;
      if (!Array.isArray(degrees) || !Array.isArray(block.controls || block.xyzw)) return false;
      const dimensions = block.shape || degrees.map(degree => degree + 1);
      const controls = block.controls || block.xyzw;
      for (let axis = 0; axis < 3; axis++) {
        if ((degrees[axis] || 1) < 2) continue;
        const shape = dimensions;
        const at = (i, j, k) => controls[(i * shape[1] + j) * shape[2] + k];
        const other = [0, 1, 2].filter(index => index !== axis);
        for (let a = 0; a < shape[other[0]]; a++) for (let b = 0; b < shape[other[1]]; b++) {
          const indices = [0, 0, 0]; indices[other[0]] = a; indices[other[1]] = b;
          const row = Array.from({ length: shape[axis] }, (_, index) => {
            indices[axis] = index; const h = at(...indices); return h.slice(0, 3).map(value => value / h[3]);
          });
          const first = row[0], last = row.at(-1);
          const chord = last.map((value, component) => value - first[component]);
          const length = Math.hypot(...chord) || 1;
          if (row.slice(1, -1).some(point => Math.hypot(...[
            (point[1] - first[1]) * chord[2] - (point[2] - first[2]) * chord[1],
            (point[2] - first[2]) * chord[0] - (point[0] - first[0]) * chord[2],
            (point[0] - first[0]) * chord[1] - (point[1] - first[1]) * chord[0],
          ]) / length > 1e-7)) return true;
        }
      }
      return false;
    });
    return {
      id: experiment.id,
      title: experiment.title || experiment.id,
      description: experiment.description || '',
      metrics: experiment.metrics || [],
      modelBlockCount: (model.blocks || []).length,
      modelInterfaceCount: (model.topology?.interfaces || []).length,
      exactExteriorFaceCount: representation?.boundary_faces
        ? representation.boundary_faces.length
        : (model.boundary_faces || model.topology?.exterior || []).length,
      exactExteriorFaceKeys,
      representationPatchCount: (representation?.patches || model.nurbs_patches || []).length,
      representationCellCount: new Set((representation?.patches || model.nurbs_patches || []).flatMap(patch => patch.cells || patch.source_blocks || [])).size,
      representationCellIds: [...new Set((representation?.patches || []).flatMap(patch => patch.cells || patch.source_blocks || []))].sort(),
      exactExtractionCount: representation?.exact_extraction_count ?? null,
      exactInterfaceCount: representation?.exact_interface_count ?? null,
      nativeNurbsPatchCount: (model.nurbs_patches || []).length,
      hasExactRepresentation: Boolean(representation),
      representationSourceMatches: !representation || representation.source_model_sha256 === modelSha256,
      exactPatchesValid,
      hasComparison: comparisonIds.has(experiment.id),
      modelAudit: audit?.model_audit || null,
      kernelRebuild: audit?.kernel_rebuild || null,
      representation: experiment.representation || null,
      section: experiment.section || null,
      links: {
        model: Boolean(experiment.model && experiment.model !== experiment.source),
        representation: Boolean(experiment.representation),
        source: Boolean(experiment.source),
        step: Boolean(experiment.step),
        audit: Boolean(experiment.audit),
        figure: Boolean(experiment.figure),
      },
      assets: await Promise.all(assets.map(async asset => {
        const response = await fetch(asset.path);
        const bytes = await response.arrayBuffer();
        return { ...asset, status: response.status, bytes: bytes.byteLength };
      })),
      curvedBlockCount: curvedBlocks.length,
      profileCount: (audit?.g1_oval_arcs || []).filter(item => Array.isArray(item.quadratic_arcs_mm)
        && item.quadratic_arcs_mm.length === 4).length,
    };}));})`,
  selectExperiment: id => `(() => {
    const select = document.querySelector('#experiment-select');
    select.value = ${JSON.stringify(id)};
    select.dispatchEvent(new Event('change'));
  })()`,
  blocksMode: "document.querySelector('[data-mode=blocks]').click()",
  zCut: fractionPercent => `(() => {
    const axis = document.querySelector('#clip-axis');
    axis.value = 'z';
    axis.dispatchEvent(new Event('change'));
    const fraction = document.querySelector('#clip-fraction');
    fraction.value = ${JSON.stringify(String(fractionPercent))};
    fraction.dispatchEvent(new Event('input'));
  })()`,
  s3SectionButton: "document.querySelector('#s3-face-pocket').click()",
  selectedClip: `(() => ({
    axis: document.querySelector('#clip-axis').value,
    fraction: Number(document.querySelector('#clip-fraction').value),
  }))()`,
  profileState: `(() => ({
    hidden: document.querySelector('#profile-card').hidden,
    count: document.querySelectorAll('#profile-card .profile-path').length,
    labels: [...document.querySelectorAll('#profile-card .profile-label')].map(item => item.textContent),
  }))()`,
  surfaceMode: "document.querySelector('[data-mode=surface]').click()",
  selectedId: "document.querySelector('#experiment-select').value",
  selectedState: `(() => ({
    id: document.querySelector('#experiment-select').value,
    title: document.querySelector('#view-title').textContent,
    stateHidden: document.querySelector('#canvas-state').hidden,
    stateText: document.querySelector('#canvas-state').textContent,
    boundaryFaceCount: document.querySelector('#model-canvas').dataset.boundaryFaceCount,
    exactExteriorFaceCount: document.querySelector('#model-canvas').dataset.exactExteriorFaceCount,
    excludedInternalFaceCount: document.querySelector('#model-canvas').dataset.excludedInternalFaceCount,
    boundarySource: document.querySelector('#model-canvas').dataset.boundarySource,
    boundaryFaceKeys: JSON.parse(document.querySelector('#model-canvas').dataset.boundaryFaceKeys || '[]'),
    metricText: document.querySelector('#metrics-list').textContent,
    metrics: Object.fromEntries([...document.querySelectorAll('#metrics-list > div')].map(row => [row.querySelector('dt')?.textContent, row.querySelector('dd')?.textContent])),
    representationLinkHidden: document.querySelector('#representation-link').hidden,
    comparisonLinkHidden: document.querySelector('.comparison-link').hidden,
    comparisonStatus: document.querySelector('.comparison-status').textContent,
    patchModeHidden: document.querySelector('[data-mode=patches]').hidden,
    patchModeLabel: document.querySelector('[data-mode=patches]').textContent,
    canvas: (() => {
      const rect = document.querySelector('#model-canvas').getBoundingClientRect();
      return { width: rect.width, height: rect.height };
    })(),
  }))()`,
  visibleVariantLinks: `(() => Object.fromEntries([
    ['model', '#model-link'], ['representation', '#representation-link'], ['source', '#source-link'], ['step', '#step-link'], ['audit', '#audit-link'],
    ['figure', '#figure-link'],
  ].map(([name, selector]) => [name, !document.querySelector(selector).hidden])))()`,
  modeControls: `(() => ({
    surface: !document.querySelector('[data-mode=surface]').disabled,
    blocks: !document.querySelector('[data-mode=blocks]').disabled,
    clip: !document.querySelector('#clip-axis').disabled,
  }))()`,
  mobileMetrics: `(() => ({
    viewportWidth: innerWidth,
    documentWidth: document.documentElement.scrollWidth,
    canvas: (() => {
      const rect = document.querySelector('#model-canvas').getBoundingClientRect();
      return { x: rect.x, y: rect.y, width: rect.width, height: rect.height };
    })(),
  }))()`,
};

function validateInjectedExpressions() {
  const expressions = [
    ['galleryReady', injectedExpressions.galleryReady],
    ['bodySummary', injectedExpressions.bodySummary],
    ['publicUiLocale', injectedExpressions.publicUiLocale],
    ['canvasRegion', injectedExpressions.canvasRegion],
    ['resources', injectedExpressions.resources],
    ['selectExperiment', injectedExpressions.selectExperiment('s3-face-pocket')],
    ['selectExperimentEscapedId', injectedExpressions.selectExperiment("quote'\\test")],
    ['blocksMode', injectedExpressions.blocksMode],
    ['zCut25', injectedExpressions.zCut(25)],
    ['zCut50', injectedExpressions.zCut(50)],
    ['selectedClip', injectedExpressions.selectedClip],
    ['profileState', injectedExpressions.profileState],
    ['surfaceMode', injectedExpressions.surfaceMode],
    ['selectedId', injectedExpressions.selectedId],
    ['selectedState', injectedExpressions.selectedState],
    ['visibleVariantLinks', injectedExpressions.visibleVariantLinks],
    ['modeControls', injectedExpressions.modeControls],
    ['mobileMetrics', injectedExpressions.mobileMetrics],
  ];
  for (const [name, expression] of expressions) {
    try {
      new Function(`return (${expression});`);
    } catch (error) {
      throw new SyntaxError(`Invalid Runtime.evaluate expression (${name}): ${error.message}`, { cause: error });
    }
  }
}

validateInjectedExpressions();

if (process.argv[2] === '--check-expressions') {
  console.log('Validated all injected Runtime.evaluate expressions.');
  process.exit(0);
}

await mkdir(OUTPUT_DIR, { recursive: true });

function sendFile(res, filePath, bytes) {
  res.writeHead(200, {
    'Content-Length': bytes.length,
    'Content-Type': MIME_TYPES[path.extname(filePath).toLowerCase()] || 'application/octet-stream',
    'Cache-Control': 'no-store',
  });
  res.end(bytes);
}

function inspectPng(bytes) {
  const signature = Buffer.from([137, 80, 78, 71, 13, 10, 26, 10]);
  if (!bytes.subarray(0, 8).equals(signature)) throw new Error('Chrome returned a non-PNG screenshot');
  let width = 0;
  let height = 0;
  let colorType = 0;
  let bitDepth = 0;
  let interlace = 0;
  const imageData = [];
  for (let offset = 8; offset < bytes.length;) {
    const length = bytes.readUInt32BE(offset);
    const type = bytes.toString('ascii', offset + 4, offset + 8);
    const data = bytes.subarray(offset + 8, offset + 8 + length);
    if (type === 'IHDR') {
      width = data.readUInt32BE(0);
      height = data.readUInt32BE(4);
      bitDepth = data[8];
      colorType = data[9];
      interlace = data[12];
    } else if (type === 'IDAT') imageData.push(data);
    else if (type === 'IEND') break;
    offset += length + 12;
  }
  const channels = colorType === 6 ? 4 : colorType === 2 ? 3 : 0;
  if (!width || !height || bitDepth !== 8 || !channels || interlace !== 0) {
    throw new Error(`Unsupported Chrome PNG format: ${width}x${height}, depth ${bitDepth}, type ${colorType}, interlace ${interlace}`);
  }
  const stride = width * channels;
  const packed = inflateSync(Buffer.concat(imageData));
  const pixels = Buffer.alloc(stride * height);
  const paeth = (a, b, c) => {
    const p = a + b - c;
    const pa = Math.abs(p - a);
    const pb = Math.abs(p - b);
    const pc = Math.abs(p - c);
    return pa <= pb && pa <= pc ? a : pb <= pc ? b : c;
  };
  for (let y = 0; y < height; y++) {
    const filter = packed[y * (stride + 1)];
    const source = y * (stride + 1) + 1;
    const row = y * stride;
    for (let x = 0; x < stride; x++) {
      const raw = packed[source + x];
      const left = x >= channels ? pixels[row + x - channels] : 0;
      const up = y > 0 ? pixels[row + x - stride] : 0;
      const upperLeft = y > 0 && x >= channels ? pixels[row + x - stride - channels] : 0;
      const predictor = filter === 0 ? 0 : filter === 1 ? left : filter === 2 ? up
        : filter === 3 ? Math.floor((left + up) / 2) : filter === 4 ? paeth(left, up, upperLeft) : null;
      if (predictor === null) throw new Error(`Unsupported PNG row filter ${filter}`);
      pixels[row + x] = (raw + predictor) & 255;
    }
  }
  const background = [240, 235, 222];
  let paintedPixels = 0;
  // Exclude the axis legend, help text and canvas edges.
  for (let index = Math.ceil(height * 0.15) * stride; index < Math.floor(height * 0.85) * stride; index += channels) {
    if (Math.max(...background.map((value, channel) => Math.abs(pixels[index + channel] - value))) > 7) {
      paintedPixels++;
    }
  }
  return { width, height, paintedPixels, pixelCount: width * height };
}

const siteServer = createServer(async (req, res) => {
  try {
    if (req.method !== 'GET' && req.method !== 'HEAD') {
      res.writeHead(405).end();
      return;
    }
    const pathname = decodeURIComponent(new URL(req.url, `http://${HOST}`).pathname);
    if (pathname === '/prior-visit.html') {
      res.writeHead(200, { 'Content-Type': 'text/html' }).end('<title>Old gallery</title>');
      return;
    }
    if (new URL(req.url, `http://${HOST}`).pathname === '/experiments.json' && !new URL(req.url, `http://${HOST}`).search) {
      res.writeHead(200, { 'Content-Type': 'application/json', 'Cache-Control': 'public, max-age=3600' })
        .end(JSON.stringify([{ id: 'g2-cross-plate', title: 'G2 plate with crossed bores', model: 'retired.json' }]));
      return;
    }
    const relativePath = pathname === '/' ? 'index.html' : pathname.slice(1);
    const filePath = path.resolve(SITE_ROOT, relativePath);
    if (!filePath.startsWith(`${SITE_ROOT}${path.sep}`)) throw new Error('Path escapes site root');
    const bytes = await readFile(filePath);
    if (req.method === 'HEAD') {
      res.writeHead(200, { 'Content-Length': bytes.length });
      res.end();
    } else {
      sendFile(res, filePath, bytes);
    }
  } catch {
    res.writeHead(404).end('Not found');
  }
});

await new Promise((resolve, reject) => {
  siteServer.once('error', reject);
  siteServer.listen(0, HOST, resolve);
});

const serverUrl = `http://${HOST}:${siteServer.address().port}/`;
const profileDir = await mkdtemp(path.join(os.tmpdir(), 'solidcraft-gallery-chrome-'));
const child = spawn(CHROME_BIN, [
  '--no-sandbox',
  '--disable-dev-shm-usage',
  '--disable-breakpad',
  '--no-first-run',
  '--headless=new',
  '--use-angle=swiftshader',
  '--enable-unsafe-swiftshader',
  '--remote-debugging-pipe',
  `--user-data-dir=${profileDir}`,
  'about:blank',
], { stdio: ['ignore', 'ignore', 'inherit', 'pipe', 'pipe'] });

const pendingCommands = new Map();
const runtimeErrors = [];
let commandId = 0;
let protocolBuffer = '';
let chromeExited = false;
let chromeExitCode = null;

child.stdio[4].on('data', chunk => {
  protocolBuffer += chunk.toString('utf8');
  let delimiter;
  while ((delimiter = protocolBuffer.indexOf('\0')) >= 0) {
    const message = JSON.parse(protocolBuffer.slice(0, delimiter));
    protocolBuffer = protocolBuffer.slice(delimiter + 1);
    if (message.id && pendingCommands.has(message.id)) {
      const [resolve, reject] = pendingCommands.get(message.id);
      pendingCommands.delete(message.id);
      if (message.error) reject(new Error(JSON.stringify(message.error)));
      else resolve(message.result);
    }
    if (message.method === 'Runtime.exceptionThrown') {
      const details = message.params.exceptionDetails;
      runtimeErrors.push(details.exception?.description || details.text || 'Unspecified JavaScript exception');
    }
  }
});

child.stdio[4].on('error', error => {
  for (const [, reject] of pendingCommands.values()) reject(error);
  pendingCommands.clear();
});
child.on('error', error => { chromeExited = true; chromeExitCode = error.code || error.message; });
child.on('exit', code => {
  chromeExited = true;
  chromeExitCode = code;
  for (const [, reject] of pendingCommands.values()) reject(new Error(`Chrome exited (${code})`));
  pendingCommands.clear();
});

function cdp(method, params = {}, sessionId) {
  if (chromeExited) return Promise.reject(new Error(`Chrome is unavailable (${chromeExitCode})`));
  return new Promise((resolve, reject) => {
    const id = ++commandId;
    pendingCommands.set(id, [resolve, reject]);
    child.stdio[3].write(`${JSON.stringify({ id, method, params, sessionId })}\0`, error => {
      if (!error) return;
      pendingCommands.delete(id);
      reject(error);
    });
  });
}

const pause = milliseconds => new Promise(resolve => setTimeout(resolve, milliseconds));
const watchdog = setTimeout(() => child.kill('SIGKILL'), 180_000);
const chromeStopped = new Promise(resolve => child.once('exit', resolve));

async function runBrowserChecks() {
  const { targetId } = await cdp('Target.createTarget', { url: 'about:blank' });
  const { sessionId } = await cdp('Target.attachToTarget', { targetId, flatten: true });
  const command = (method, params) => cdp(method, params, sessionId);
  const evaluate = async expression => {
    const response = await command('Runtime.evaluate', {
      expression,
      returnByValue: true,
      awaitPromise: true,
    });
    if (response.exceptionDetails) {
      throw new Error(response.exceptionDetails.exception?.description || response.exceptionDetails.text);
    }
    return response.result.value;
  };
  const waitForGallery = async expectedId => {
    for (let attempt = 0; attempt < 300; attempt++) {
      const ready = await evaluate(injectedExpressions.galleryReady);
      const selected = expectedId
        ? await evaluate(injectedExpressions.selectedId) === expectedId
        : true;
      if (ready && selected) {
        await pause(100);
        return;
      }
      await pause(200);
    }
    const summary = await evaluate(injectedExpressions.bodySummary);
    throw new Error(`Gallery did not become ready: ${summary}`);
  };
  const capture = async (name, region, analyze = false, persist = true) => {
    const result = await command('Page.captureScreenshot', {
      format: 'png',
      ...(region ? { clip: region } : {}),
    });
    const bytes = Buffer.from(result.data, 'base64');
    if (persist) await writeFile(path.join(OUTPUT_DIR, `${name}.png`), bytes);
    const hash = createHash('sha256').update(bytes).digest('hex');
    return analyze ? { hash, ...inspectPng(bytes) } : hash;
  };
  const canvasRegion = async () => evaluate(injectedExpressions.canvasRegion);
  const selectExperiment = async id => {
    await evaluate(injectedExpressions.selectExperiment(id));
    await waitForGallery(id);
    const state = await evaluate(injectedExpressions.selectedState);
    if (state.id !== id || !state.stateHidden || state.canvas.width < 1 || state.canvas.height < 1) {
      throw new Error(`Experiment ${id} did not recover to a rendered canvas: ${JSON.stringify(state)}`);
    }
    return state;
  };

  await command('Runtime.enable', {});
  await command('Page.enable', {});
  await command('Emulation.setDeviceMetricsOverride', {
    width: 1400, height: 960, deviceScaleFactor: 1, mobile: false,
  });
  await command('Page.navigate', { url: serverUrl + 'prior-visit.html' });
  for (let i = 0; i < 100 && !(await evaluate("document.title === 'Old gallery'")); i++) await pause(50);
  const oldId = await evaluate("fetch('./experiments.json').then(r => r.json()).then(rows => rows[0].id)");
  if (oldId !== 'g2-cross-plate') throw new Error('Historical catalog was not primed');
  const landingId = await evaluate("fetch('./experiments.json?catalog=qa').then(r => r.json()).then(rows => rows[0].id)");
  await command('Page.navigate', { url: serverUrl + '?model=retired-curved-catalog-id' });
  await waitForGallery(landingId);
  const publicUiLocale = await evaluate(injectedExpressions.publicUiLocale);
  if (publicUiLocale.lang !== 'en' || publicUiLocale.heading !== 'Sketch partitions in 3D'
    || !publicUiLocale.canvasLabel.startsWith('Interactive 3D model')
    || !publicUiLocale.resetLabel.endsWith('Reset view') || publicUiLocale.experimentLabel !== 'Experiment') {
    throw new Error(`Public gallery UI is not consistently in English: ${JSON.stringify(publicUiLocale)}`);
  }
  if (await evaluate("new URL(location.href).searchParams.get('model')") !== landingId) {
    throw new Error('Retired model permalink was not corrected');
  }

  const resources = await evaluate(injectedExpressions.resources);
  const historicalIds = [
    "plate_blind_pocket",
    "rectangle_internal_rectangle",
    "rectangle_offset_pocket",
    "c_section_two_holes",
    "l_bracket_polygon_holes_and_pocket",
    "c_section_lower_hole_x_plus_1",
    "fieldcraft_x40_y80",
    "fieldcraft_x50_y74",
    "fieldcraft_x51_y74",
    "fieldcraft_x64_y64",
    "fieldcraft_x80_y40",
  ];
  const ids = resources.map(item => item.id);
  const polyIds = ids.filter(id => id.startsWith('bench_poly_'));
  const bracketIds = ids.filter(id => id.startsWith('bench-bracket-'));
  const bracketCases = [
    'bench-bracket-original', 'bench-bracket-in-plane', 'bench-bracket-east-stack',
    'bench-bracket-north-stack', 'bench-bracket-east-offset-stack', 'bench-bracket-diagonal',
    'bench-bracket-three-levels', 'bench-bracket-both-stacked', 'bench-bracket-axial-z',
    'bench-bracket-unequal-stack', 'bench-bracket-unequal-transverse-bores',
  ];
  if (resources.length !== EXPECTED_ENTRY_COUNT || new Set(ids).size !== ids.length
    || polyIds.length !== 46 || bracketIds.length !== 11
    || historicalIds.some(id => !ids.includes(id))
    || !ids.includes('bench_poly_projection_diamonds')
    || bracketCases.some(id => !ids.includes(id))
    || ids[0] !== 'bench-bracket-east-stack') {
    throw new Error(`Unexpected gallery catalog: count=${resources.length}, poly=${polyIds.length}, brackets=${bracketIds.length}, first=${ids[0]}`);
  }
  if (resources.some(experiment => experiment.assets.length === 0
    || experiment.assets.some(asset => asset.status !== 200 || asset.bytes < 1))) {
    throw new Error(`Gallery resource requests failed: ${JSON.stringify(resources)}`);
  }
  const exactRepresentations = resources.filter(experiment => experiment.hasExactRepresentation);
  if (exactRepresentations.length !== 29 || exactRepresentations.some(experiment =>
    !experiment.representationSourceMatches
    || !experiment.exactPatchesValid
    || experiment.exactExtractionCount !== experiment.modelBlockCount
    || experiment.representationCellCount !== experiment.modelBlockCount
    || experiment.representationCellIds.length !== experiment.modelBlockCount
    || experiment.exactInterfaceCount !== experiment.modelInterfaceCount
  )) {
    throw new Error(`Exact NURBS supplements do not match their native models: ${JSON.stringify(exactRepresentations.map(({ id, representationSourceMatches, exactPatchesValid, exactExtractionCount, representationCellCount, modelBlockCount, exactInterfaceCount, modelInterfaceCount }) => ({ id, representationSourceMatches, exactPatchesValid, exactExtractionCount, representationCellCount, modelBlockCount, exactInterfaceCount, modelInterfaceCount })))}`);
  }
  if (resources.some(experiment => !experiment.modelAudit
    || experiment.modelAudit.status !== 'valid'
    || experiment.modelAudit.blocks !== experiment.modelBlockCount
    || experiment.modelAudit.shared_faces !== experiment.modelInterfaceCount
    || !experiment.kernelRebuild?.kernel_binary_sha256
    || !experiment.kernelRebuild?.recipe_sha256
    || !experiment.kernelRebuild?.model_sha256
    || !experiment.kernelRebuild?.step_sha256)) {
    throw new Error(`A gallery model is missing its current Rust audit/rebuild record: ${JSON.stringify(resources.filter(experiment => !experiment.modelAudit || experiment.modelAudit.blocks !== experiment.modelBlockCount || experiment.modelAudit.shared_faces !== experiment.modelInterfaceCount || !experiment.kernelRebuild?.kernel_binary_sha256).map(({ id, modelBlockCount, modelInterfaceCount, modelAudit, kernelRebuild }) => ({ id, modelBlockCount, modelInterfaceCount, modelAudit, kernelRebuild })))}`);
  }
  if (resources.some(experiment => experiment.assets.some(asset => !/[?&]v=[^&#]+/.test(asset.path)))) {
    throw new Error('A model, recipe, audit, STEP, or figure URL is missing the current cache version');
  }
  const legacyCurvedIds = new Set(historicalIds.slice(0, 6));
  if (resources.some(experiment => legacyCurvedIds.has(experiment.id) && experiment.curvedBlockCount < 1)) {
    throw new Error(`Historical curved geometry is missing: ${JSON.stringify(resources.filter(({ id }) => legacyCurvedIds.has(id)).map(({ id, curvedBlockCount }) => ({ id, curvedBlockCount })))}`);
  }
  if (resources.some(experiment => legacyCurvedIds.has(experiment.id) && experiment.profileCount < 1)) {
    throw new Error(`A curved-profile source audit has no four-arc profile data: ${JSON.stringify(resources.map(({ id, profileCount }) => ({ id, profileCount })))}`);
  }
  const version = await cdp('Browser.getVersion');
  const selectedExperiment = resources[0];

  await selectExperiment(selectedExperiment.id);
  const desktopHash = await capture('desktop');
  await evaluate(injectedExpressions.surfaceMode);
  await pause(500);
  const surfaceHash = await capture('surface', await canvasRegion());
  await evaluate(injectedExpressions.blocksMode);
  await pause(700);
  const blocksHash = await capture('blocks', await canvasRegion());
  const cutPercent = 80;
  await evaluate(injectedExpressions.zCut(cutPercent));
  await pause(1200);
  const blocksCutHash = await capture('blocks-cut', await canvasRegion());
  await evaluate(injectedExpressions.surfaceMode);
  await pause(500);
  const cutHash = await capture('cut', await canvasRegion());
  if (new Set([surfaceHash, blocksHash, blocksCutHash, cutHash]).size !== 4) {
    throw new Error('Surface, blocks, and cut screenshots were not visually distinct');
  }

  await evaluate("document.querySelector('#reset-view').click()");
  const experimentChecks = [];
  const modelScreens = [];
  const curvedSectionHashes = [];
  const screenshotIds = new Set([
    resources[0].id,
    'bench_poly_projection_diamonds',
    'bench-bracket-east-stack',
    'bench-bracket-unequal-transverse-bores',
  ]);
  for (const experiment of resources) {
    await command('Page.navigate', { url: `${serverUrl}?model=${encodeURIComponent(experiment.id)}` });
    await waitForGallery(experiment.id);
    const state = await evaluate(injectedExpressions.selectedState);
    if (state.id !== experiment.id || !state.stateHidden) {
      throw new Error(`Direct model URL ${experiment.id} did not render: ${JSON.stringify(state)}`);
    }
    if (experiment.representationPatchCount > 0 || experiment.nativeNurbsPatchCount > 0) {
      const merged = await evaluate(`({
        active: document.querySelector('[data-mode="patches"]').getAttribute('aria-pressed') === 'true',
        compare: new URL(document.querySelector('.comparison-link').href).searchParams.get('model'),
        count: document.querySelector('#metrics-list').textContent.includes('NURBS volumes')
      })`);
      if (!merged.active || !merged.count || (experiment.hasComparison && merged.compare !== experiment.id)) {
        throw new Error('Missing default NURBS view for ' + experiment.id + ': ' + JSON.stringify(merged));
      }
      if (experiment.id === 'bench_poly_seeded_convex_14') {
        const macro = await capture('seeded-14-volumes', await canvasRegion());
        await evaluate(injectedExpressions.blocksMode);
        await pause(100);
        if (macro === await capture('seeded-14-cells', await canvasRegion())) throw new Error('Seeded NURBS volumes and cells appear identical');
      }
    }
    const visibleLinks = await evaluate(injectedExpressions.visibleVariantLinks);
    const controls = await evaluate(injectedExpressions.modeControls);
    const isBoundary = experiment.representation === 'boundary';
    if (isBoundary && (controls.blocks || controls.clip)) {
      throw new Error(`Boundary experiment ${experiment.id} exposes unsupported block/cut controls: ${JSON.stringify(controls)}`);
    }
    for (const key of Object.keys(experiment.links)) {
      if (visibleLinks[key] !== experiment.links[key]) {
        throw new Error(`Experiment ${experiment.id} link ${key} visibility mismatch: expected ${experiment.links[key]}, got ${visibleLinks[key]}`);
      }
    }
    if (state.representationLinkHidden !== !experiment.links.representation) {
      throw new Error(`Experiment ${experiment.id} representation link visibility mismatch`);
    }
    if (state.comparisonLinkHidden !== !experiment.hasComparison) {
      throw new Error(`Experiment ${experiment.id} comparison availability is misreported`);
    }
    if (!experiment.hasComparison && !state.comparisonStatus.includes('No paired partition comparison')) {
      throw new Error(`Experiment ${experiment.id} does not explain the missing paired comparison`);
    }
    const hasPatchView = experiment.representationPatchCount > 0 || experiment.nativeNurbsPatchCount > 0;
    if (state.patchModeHidden === hasPatchView) throw new Error(`Experiment ${experiment.id} NURBS cover mode visibility is wrong`);
    if (experiment.representationPatchCount > 0) {
      const metrics = state.metrics;
      if (Number(metrics['NURBS volumes']) !== experiment.representationPatchCount
        || Number(metrics['Bézier cells']) !== experiment.representationCellCount) {
        throw new Error(`Experiment ${experiment.id} shows incorrect exact cover metrics: ${JSON.stringify(metrics)}`);
      }
    }
    if (state.title !== experiment.title) {
      throw new Error(`Experiment ${experiment.id} selected, but title is ${JSON.stringify(state.title)}`);
    }
    await evaluate(injectedExpressions.surfaceMode);
    await pause(100);
    const boundaryState = await evaluate(injectedExpressions.selectedState);
    if (experiment.exactExteriorFaceCount > 0) {
      const visibleCount = Number(boundaryState.boundaryFaceCount);
      const exactCount = Number(boundaryState.exactExteriorFaceCount);
      const excludedCount = Number(boundaryState.excludedInternalFaceCount);
      const totalFaces = experiment.modelBlockCount * 6;
      if (visibleCount !== experiment.exactExteriorFaceCount || exactCount !== visibleCount
        || excludedCount !== totalFaces - visibleCount
        || boundaryState.boundarySource === 'interfaces'
        || JSON.stringify([...boundaryState.boundaryFaceKeys].sort()) !== JSON.stringify(experiment.exactExteriorFaceKeys)) {
        throw new Error(`Renderer did not draw only the recorded exterior faces for ${experiment.id}: ${JSON.stringify(boundaryState)}`);
      }
    }
    if (hasPatchView) await evaluate("document.querySelector('[data-mode=patches]').click()");
    const profile = await evaluate(injectedExpressions.profileState);
    if (experiment.profileCount > 0 && (profile.hidden || profile.count !== Math.min(2, experiment.profileCount))) {
      throw new Error(`Curved profile vignette did not show audited arcs for ${experiment.id}: ${JSON.stringify({ profile, profileCount: experiment.profileCount })}`);
    }
    if (experiment.profileCount === 0 && !profile.hidden) {
      throw new Error(`Unexpected curved-profile vignette without audited arcs for ${experiment.id}: ${JSON.stringify(profile)}`);
    }
    if (experiment.id === 'plate_blind_pocket' && !await evaluate("document.querySelector('.profile-label').textContent.endsWith('20 mm')")) {
      throw new Error('The dimension must measure the profile, including its interior extrema');
    }
    const section = experiment.section || { axis: 'z', fraction: 0.5 };
    const activeClip = await evaluate(injectedExpressions.selectedClip);
    if (activeClip.axis !== section.axis || Math.abs(activeClip.fraction / 100 - section.fraction) > 0.011) {
      throw new Error(`Experiment ${experiment.id} did not apply its prescribed curved section: ${JSON.stringify({ activeClip, section })}`);
    }
    const modelScreen = await capture(
      `model-${experiment.id.replace(/[^a-z0-9_-]/gi, '-')}`,
      await canvasRegion(), true, screenshotIds.has(experiment.id),
    );
    if (modelScreen.paintedPixels < 250) {
      throw new Error(`Experiment ${experiment.id} canvas appears blank: ${JSON.stringify(modelScreen)}`);
    }
    modelScreens.push({ id: experiment.id, representation: experiment.representation || 'trivariate', ...modelScreen, visibleLinks, controls, state });
    if (screenshotIds.has(experiment.id)) {
      const cutScreen = await capture(`curved-section-${experiment.id.replace(/[^a-z0-9_-]/gi, '-')}`, await canvasRegion(), true, true);
      if (cutScreen.paintedPixels < 250) throw new Error(`Curved section ${experiment.id} rendered blank: ${JSON.stringify(cutScreen)}`);
      if (experiment.profileCount > 0) curvedSectionHashes.push(cutScreen.hash);
    }
    experimentChecks.push(experiment.id);
  }
  if (new Set(modelScreens.map(item => item.hash)).size < 2) {
    throw new Error('Canvas did not change across experiment selections; model rendering may be stale');
  }
  if (new Set(curvedSectionHashes).size !== curvedSectionHashes.length) {
    throw new Error(`Curved focal sections are not visually distinct: ${JSON.stringify(curvedSectionHashes)}`);
  }

  // Re-select the first entry after traversing the full catalog to verify recovery from
  // repeated model replacement and to leave a known model active for the mobile pass.
  await command('Emulation.setDeviceMetricsOverride', {
    width: 390, height: 844, deviceScaleFactor: 1, mobile: true,
  });
  let mobileMetrics;
  const mobileScreens = [];
  for (const experiment of resources) {
    await selectExperiment(experiment.id);
    await pause(100);
    mobileMetrics = await evaluate(injectedExpressions.mobileMetrics);
    if (mobileMetrics.documentWidth > mobileMetrics.viewportWidth + 1
      || mobileMetrics.canvas.width < 1 || mobileMetrics.canvas.height < 1) {
      throw new Error(`Mobile layout failed for ${experiment.id}: ${JSON.stringify(mobileMetrics)}`);
    }
    const mobileScreen = await capture(`mobile-${experiment.id.replace(/[^a-z0-9_-]/gi, '-')}`,
      await canvasRegion(), true, screenshotIds.has(experiment.id));
    if (mobileScreen.paintedPixels < 250) {
      throw new Error(`Mobile WebGL canvas for ${experiment.id} appears blank: ${JSON.stringify(mobileScreen)}`);
    }
    if (screenshotIds.has(experiment.id)) mobileScreens.push({ id: experiment.id, ...mobileScreen });
  }

  await command('Emulation.setDeviceMetricsOverride', {width: 1400, height: 960, deviceScaleFactor: 1, mobile: false});
  await command('Page.navigate', {url: serverUrl + 'compare.html?model=bench_poly_through_octagon_high_aspect'});
  for (let i = 0; i < 200 && !await evaluate("document.querySelector('#assemble-stats')?.textContent.includes('4 NURBS volumes')"); i++) await pause(50);
  if (!await evaluate("document.querySelector('#sketch')?.options.length === 46 && document.querySelector('#optimize-stats').textContent.includes('8 Bézier cells')")) {
    throw new Error('Octagon comparison did not load');
  }
  const regions = await evaluate(`['assemble','optimize'].map(id => {
    const r=document.getElementById(id).getBoundingClientRect();
    return {x:r.left,y:r.top,width:r.width,height:r.height,scale:1};
  })`);
  const views = async () => {
    await pause(120);
    const hashes = [];
    for (const r of regions) hashes.push((await capture('compare-view', r, true, false)).hash);
    return hashes;
  };
  const moveCamera = async key => {
    await evaluate(`document.getElementById('assemble').dispatchEvent(new KeyboardEvent('keydown', ${JSON.stringify(key)}))`);
  };
  let before = await views();
  for (const event of [{key:'ArrowLeft'},{key:'ArrowUp',shiftKey:true},{key:'+'}]) {
    await moveCamera(event);
    const after = await views();
    if (after.some((v,i)=>v===before[i])) throw new Error('Synchronized camera failed: '+JSON.stringify(event));
    before = after;
  }
  await evaluate("document.getElementById('sync').click()");
  await moveCamera({key:'ArrowRight'});
  const unlocked = await views();
  if (unlocked[0]===before[0] || unlocked[1]!==before[1]) throw new Error('Independent cameras failed');
  await evaluate("document.getElementById('sync').click()");
  if ((await views())[1] === unlocked[1]) throw new Error('Relinking views failed');
  await evaluate("document.getElementById('sketch').value='bench_poly_blind_two_levels';document.getElementById('sketch').dispatchEvent(new Event('change'))");
  if (!await evaluate("document.getElementById('optimize-stats').textContent.includes('24 NURBS volumes · 96 Bézier cells')")) throw new Error('Wrong two-depth pocket');
  const mergedSingleton = resources.find(item => item.hasComparison && item.hasExactRepresentation
    && /^bench_poly_seeded_convex_00$/.test(item.id));
  if (mergedSingleton) {
    await command('Page.navigate', {url: `${serverUrl}compare.html?model=${encodeURIComponent(mergedSingleton.id)}`});
    for (let i = 0; i < 200 && await evaluate("document.querySelector('#right-title')?.textContent !== 'Published exact NURBS cover'"); i++) await pause(50);
    const coverState = await evaluate(`({
      title: document.querySelector('#right-title').textContent,
      stats: document.querySelector('#optimize-stats').textContent,
      variant: document.querySelector('#right-variant').value,
      patchMode: document.querySelector('#colour').value,
      canvas: (() => { const r=document.getElementById('optimize').getBoundingClientRect(); return {x:r.left,y:r.top,width:r.width,height:r.height,scale:1}; })(),
    })`);
    if (coverState.title !== 'Published exact NURBS cover' || coverState.variant !== 'cover'
      || !coverState.stats.includes(`${mergedSingleton.representationPatchCount} merged NURBS patches`)
      || !coverState.stats.includes(`${mergedSingleton.representationCellCount} Bézier cells`)) {
      throw new Error(`Supplemented singleton comparison did not load its exact cover: ${JSON.stringify(coverState)}`);
    }
    const coverPixels = await capture('singleton-cover', coverState.canvas, true, true);
    await evaluate("document.querySelector('#colour').value='cell';document.querySelector('#colour').dispatchEvent(new Event('change'))");
    await pause(200);
    const cellPixels = await capture('singleton-cover-cells', coverState.canvas, true, true);
    if (coverPixels.hash === cellPixels.hash) throw new Error('Merged NURBS cover and Bézier cells look identical');
    await evaluate("document.querySelector('#cut').click();document.querySelector('#fraction').value='62';document.querySelector('#fraction').dispatchEvent(new Event('input'))");
    await pause(300);
    const cutPixels = await capture('singleton-cover-section', coverState.canvas, true, true);
    if (cutPixels.hash === cellPixels.hash) throw new Error('Exact NURBS cover section did not change the rendered boundary');
    await evaluate("document.querySelector('#right-variant').value='optimize';document.querySelector('#right-variant').dispatchEvent(new Event('change'))");
    await pause(150);
    if (!await evaluate("document.querySelector('#right-title').textContent === 'Native optimized'")) {
      throw new Error('Native optimized comparison candidate is not selectable');
    }
  }
  const noPair = resources.find(item => !item.hasComparison);
  if (noPair) {
    await command('Page.navigate', {url: `${serverUrl}compare.html?model=${encodeURIComponent(noPair.id)}`});
    const noPairStatus = await evaluate("({text:document.querySelector('#status').textContent,hidden:document.querySelector('.comparison').hidden})");
    if (!noPairStatus.hidden || !noPairStatus.text.includes('No paired partition comparison')) {
      throw new Error(`Unpaired model was given an unrelated comparison: ${JSON.stringify(noPairStatus)}`);
    }
  }
  await pause(150);
  await capture('nurbs-comparison-desktop', null, false, true);

  if (runtimeErrors.length) throw new Error(`Browser JavaScript errors: ${runtimeErrors.join('\n')}`);
  const result = {
    url: serverUrl,
    chrome: { product: version.product, userAgent: version.userAgent },
    experiments: experimentChecks,
    modelScreens,
    mobileScreens,
    selectedExperiment: selectedExperiment.id,
    cutPercent,
    screenshotHashes: { desktop: desktopHash, surface: surfaceHash, blocks: blocksHash, blocksCut: blocksCutHash, cut: cutHash },
    mobile: mobileMetrics,
    runtimeErrors,
  };
  await writeFile(path.join(OUTPUT_DIR, 'result.json'), `${JSON.stringify(result, null, 2)}\n`);
  return result;
}

try {
  const result = await runBrowserChecks();
  console.log(JSON.stringify(result, null, 2));
} catch (error) {
  console.error(error.stack || String(error));
  process.exitCode = 1;
} finally {
  clearTimeout(watchdog);
  if (!chromeExited) {
    child.kill('SIGTERM');
    await Promise.race([chromeStopped, pause(5_000)]);
    if (!chromeExited) {
      child.kill('SIGKILL');
      await Promise.race([chromeStopped, pause(2_000)]);
    }
  }
  await new Promise(resolve => siteServer.close(resolve));
  await rm(profileDir, { recursive: true, force: true, maxRetries: 8, retryDelay: 250 });
}
