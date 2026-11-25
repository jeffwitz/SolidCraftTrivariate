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
import { mkdtemp, mkdir, readFile, rm, writeFile } from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const REPO_ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const SITE_ROOT = path.join(REPO_ROOT, 'docs/site');
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
  resources: `fetch('./experiments.json').then(response => response.json()).then(async experiments =>
    Promise.all(experiments.map(async experiment => ({
      id: experiment.id,
      title: experiment.title || experiment.id,
      status: await Promise.all(['model', 'source', 'step'].map(async key => {
        const response = await fetch(experiment[key]);
        return response.status;
      })),
    }))))`,
  selectExperiment: id => `(() => {
    const select = document.querySelector('#experiment-select');
    select.value = ${JSON.stringify(id)};
    select.dispatchEvent(new Event('change'));
  })()`,
  blocksMode: "document.querySelector('[data-mode=blocks]').click()",
  zCut: `(() => {
    const axis = document.querySelector('#clip-axis');
    axis.value = 'z';
    axis.dispatchEvent(new Event('change'));
    const fraction = document.querySelector('#clip-fraction');
    fraction.value = '50';
    fraction.dispatchEvent(new Event('input'));
  })()`,
  selectedId: "document.querySelector('#experiment-select').value",
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
    ['zCut', injectedExpressions.zCut],
    ['selectedId', injectedExpressions.selectedId],
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

const siteServer = createServer(async (req, res) => {
  try {
    if (req.method !== 'GET' && req.method !== 'HEAD') {
      res.writeHead(405).end();
      return;
    }
    const pathname = decodeURIComponent(new URL(req.url, `http://${HOST}`).pathname);
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
  const waitForGallery = async () => {
    for (let attempt = 0; attempt < 150; attempt++) {
      const ready = await evaluate(injectedExpressions.galleryReady);
      if (ready) return;
      await pause(200);
    }
    const summary = await evaluate(injectedExpressions.bodySummary);
    throw new Error(`Gallery did not become ready: ${summary}`);
  };
  const capture = async (name, region) => {
    const result = await command('Page.captureScreenshot', {
      format: 'png',
      ...(region ? { clip: region } : {}),
    });
    const bytes = Buffer.from(result.data, 'base64');
    await writeFile(path.join(OUTPUT_DIR, `${name}.png`), bytes);
    return createHash('sha256').update(bytes).digest('hex');
  };
  const canvasRegion = async () => evaluate(injectedExpressions.canvasRegion);
  const selectExperiment = async id => {
    await evaluate(injectedExpressions.selectExperiment(id));
    await waitForGallery();
  };

  await command('Runtime.enable', {});
  await command('Page.enable', {});
  await command('Emulation.setDeviceMetricsOverride', {
    width: 1400, height: 960, deviceScaleFactor: 1, mobile: false,
  });
  await command('Page.navigate', { url: serverUrl });
  await waitForGallery();

  const resources = await evaluate(injectedExpressions.resources);
  if (resources.some(experiment => experiment.status.some(status => status !== 200))) {
    throw new Error(`Gallery resource requests failed: ${JSON.stringify(resources)}`);
  }

  const version = await cdp('Browser.getVersion');
  const s3Experiment = resources.find(experiment => /(?:^|[-_])s3(?:[-_]|$)|face[-_]?pocket/i.test(`${experiment.id} ${experiment.title}`));
  const selectedExperiment = s3Experiment || resources[0];
  if (!selectedExperiment) throw new Error('No usable gallery experiments were found');

  await selectExperiment(selectedExperiment.id);
  const desktopHash = await capture('desktop');
  const surfaceHash = await capture('surface', await canvasRegion());
  await evaluate(injectedExpressions.blocksMode);
  await pause(700);
  const blocksHash = await capture('blocks', await canvasRegion());
  await evaluate(injectedExpressions.zCut);
  await pause(1200);
  const cutHash = await capture('cut', await canvasRegion());
  if (new Set([surfaceHash, blocksHash, cutHash]).size !== 3) {
    throw new Error('Surface, blocks, and cut screenshots were not visually distinct');
  }

  const experimentChecks = [];
  for (const experiment of resources) {
    await selectExperiment(experiment.id);
    const selectedId = await evaluate(injectedExpressions.selectedId);
    if (selectedId !== experiment.id) throw new Error(`Failed to select ${experiment.id}`);
    experimentChecks.push(experiment.id);
  }

  await selectExperiment(selectedExperiment.id);
  await command('Emulation.setDeviceMetricsOverride', {
    width: 390, height: 844, deviceScaleFactor: 1, mobile: true,
  });
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
    selectedS3Experiment: s3Experiment?.id || null,
    screenshotHashes: { desktop: desktopHash, surface: surfaceHash, blocks: blocksHash, cut: cutHash, mobile: mobileHash },
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
