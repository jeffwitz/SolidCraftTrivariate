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
  canvasRegion: `(() => {
    const rect = document.querySelector('#model-canvas').getBoundingClientRect();
    return { x: rect.x, y: rect.y, width: rect.width, height: rect.height, scale: 1 };
  })()`,
  resources: `fetch(document.getElementById('catalog-url').href).then(response => response.json()).then(async experiments =>
    Promise.all(experiments.map(async experiment => {
      const assetExtensions = ['.json', '.step', '.vtk', '.vtu', '.stl', '.obj', '.glb', '.gltf', '.png', '.jpg', '.jpeg'];
      const assets = Object.entries(experiment)
        .filter(([key, value]) => typeof value === 'string'
          && assetExtensions.some(extension => value.split('?')[0].split('#')[0].toLowerCase().endsWith(extension))
          && !value.startsWith('http://') && !value.startsWith('https://'))
        .map(([key, value]) => ({ key, path: value }));
    const modelResponse = await fetch(experiment.model);
    const model = await modelResponse.json();
    const auditResponse = experiment.audit ? await fetch(experiment.audit) : null;
    const audit = auditResponse?.ok ? await auditResponse.json() : null;
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
      representation: experiment.representation || null,
      section: experiment.section || null,
      links: {
        model: Boolean(experiment.model && experiment.model !== experiment.source),
        source: Boolean(experiment.source),
        step: Boolean(experiment.step),
        audit: Boolean(experiment.audit),
        figure: Boolean(experiment.figure),
      },
      assets: await Promise.all(assets.map(async asset => {
        const response = await fetch(asset.path, { method: 'HEAD' });
        return { ...asset, status: response.status };
      })),
      curvedBlockCount: curvedBlocks.length,
      profileCount: (audit?.g1_oval_arcs || []).filter(item => Array.isArray(item.quadratic_arcs_mm)
        && item.quadratic_arcs_mm.length === 4).length,
    };})))`,
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
    canvas: (() => {
      const rect = document.querySelector('#model-canvas').getBoundingClientRect();
      return { width: rect.width, height: rect.height };
    })(),
  }))()`,
  visibleVariantLinks: `(() => Object.fromEntries([
    ['model', '#model-link'], ['source', '#source-link'], ['step', '#step-link'], ['audit', '#audit-link'],
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
        .end(JSON.stringify([{ id: 'g2-cross-plate', title: 'Plaque G2 à perçages croisés', model: 'retired.json' }]));
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
  const capture = async (name, region, analyze = false) => {
    const result = await command('Page.captureScreenshot', {
      format: 'png',
      ...(region ? { clip: region } : {}),
    });
    const bytes = Buffer.from(result.data, 'base64');
    await writeFile(path.join(OUTPUT_DIR, `${name}.png`), bytes);
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
  await command('Page.navigate', { url: serverUrl + '?model=retired-curved-catalog-id' });
  await waitForGallery('plate_blind_pocket');
  if (await evaluate("new URL(location.href).searchParams.get('model')") !== 'plate_blind_pocket') {
    throw new Error('Retired model permalink was not corrected');
  }

  const resources = await evaluate(injectedExpressions.resources);
  const expected = [
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
  if (JSON.stringify(resources.map(item => item.id)) !== JSON.stringify(expected)) {
    throw new Error('Unexpected qualified gallery catalog');
  }
  if (resources.some(experiment => experiment.assets.length === 0
    || experiment.assets.some(asset => asset.status !== 200))) {
    throw new Error(`Gallery resource requests failed: ${JSON.stringify(resources)}`);
  }
  if (resources.some(experiment => experiment.assets.some(asset => !/[?&]v=[^&#]+/.test(asset.path)))) {
    throw new Error('A model, recipe, audit, STEP, or figure URL is missing the current cache version');
  }
  const legacyCurvedIds = new Set(expected.slice(0, 6));
  if (resources.some(experiment => experiment.curvedBlockCount < 1)) {
    throw new Error(`Expected actual non-linear geometry in every model: ${JSON.stringify(resources.map(({ id, curvedBlockCount }) => ({ id, curvedBlockCount })))}`);
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
  for (const experiment of resources) {
    await command('Page.navigate', { url: `${serverUrl}?model=${encodeURIComponent(experiment.id)}` });
    await waitForGallery(experiment.id);
    const state = await evaluate(injectedExpressions.selectedState);
    if (state.id !== experiment.id || !state.stateHidden) {
      throw new Error(`Direct model URL ${experiment.id} did not render: ${JSON.stringify(state)}`);
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
    if (state.title !== experiment.title) {
      throw new Error(`Experiment ${experiment.id} selected, but title is ${JSON.stringify(state.title)}`);
    }
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
    const modelScreen = await capture(`model-${experiment.id.replace(/[^a-z0-9_-]/gi, '-')}`, await canvasRegion(), true);
    if (modelScreen.paintedPixels < 250) {
      throw new Error(`Experiment ${experiment.id} canvas appears blank: ${JSON.stringify(modelScreen)}`);
    }
    modelScreens.push({ id: experiment.id, representation: experiment.representation || 'trivariate', ...modelScreen, visibleLinks, controls, state });
    const cutScreen = await capture(`curved-section-${experiment.id.replace(/[^a-z0-9_-]/gi, '-')}`, await canvasRegion(), true);
    if (cutScreen.paintedPixels < 250) throw new Error(`Curved section ${experiment.id} rendered blank: ${JSON.stringify(cutScreen)}`);
    if (experiment.profileCount > 0) curvedSectionHashes.push(cutScreen.hash);
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
  await command('Page.navigate', { url: serverUrl + '?model=c_section_lower_hole_x_plus_1' });
  await waitForGallery('c_section_lower_hole_x_plus_1');
  await selectExperiment(selectedExperiment.id);
  await command('Emulation.setDeviceMetricsOverride', {
    width: 390, height: 844, deviceScaleFactor: 1, mobile: true,
  });
  await command('Page.navigate', { url: `${serverUrl}?model=${encodeURIComponent(selectedExperiment.id)}` });
  await waitForGallery(selectedExperiment.id);
  await pause(500);
  const mobileMetrics = await evaluate(injectedExpressions.mobileMetrics);
  if (mobileMetrics.documentWidth > mobileMetrics.viewportWidth + 1
    || mobileMetrics.canvas.width < 1 || mobileMetrics.canvas.height < 1) {
    throw new Error(`Mobile layout failed: ${JSON.stringify(mobileMetrics)}`);
  }
  const mobileHash = await capture('mobile');
  await capture('mobile-canvas', await canvasRegion());

  if (runtimeErrors.length) throw new Error(`Browser JavaScript errors: ${runtimeErrors.join('\n')}`);
  const result = {
    url: serverUrl,
    chrome: { product: version.product, userAgent: version.userAgent },
    experiments: experimentChecks,
    modelScreens,
    selectedExperiment: selectedExperiment.id,
    cutPercent,
    screenshotHashes: { desktop: desktopHash, surface: surfaceHash, blocks: blocksHash, blocksCut: blocksCutHash, cut: cutHash, mobile: mobileHash },
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
