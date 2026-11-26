import { createViewer } from './renderer.js?v=509b5502f556';

const byId = id => document.getElementById(id);
const status = byId('box-cylinder-status');
const viewers = {};
let syncing = false;
let view = 'volume';
let mode = 'blocks';
let tallViewer;
let tallManifest;
let tallMode = 'bezier';
let tallView = 'volume';
let tallLoadSerial = 0;
let previousTallState = 'union';
let tallCameraInitialized = false;

async function readModel(name) {
  const response = await fetch(`./box-cylinder/${name}.model.json`);
  if (!response.ok) throw new Error(`Could not load ${name} model (${response.status}).`);
  return response.json();
}

function setCamera(next) {
  syncing = true;
  for (const viewer of Object.values(viewers)) viewer.setCamera(next);
  syncing = false;
}

function setMode(next) {
  mode = next;
  for (const viewer of Object.values(viewers)) viewer.setMode(next);
  byId('mode-surface').setAttribute('aria-pressed', String(next === 'surface'));
  byId('mode-blocks').setAttribute('aria-pressed', String(next === 'blocks'));
  byId('mode-note').textContent = next === 'blocks'
    ? 'Blocks mode colors the four chart volumes and exposes their interfaces.'
    : 'Surface mode draws the recorded exterior faces of each model.';
}

function setView(next) {
  view = next;
  const current = viewers.result.getCamera();
  setCamera({
    ...current,
    yaw: next === 'top' ? 0 : 0.72,
    pitch: next === 'top' ? Math.PI / 2 - 1e-7 : 0.5,
  });
  byId('view-top').setAttribute('aria-pressed', String(next === 'top'));
  byId('view-volume').setAttribute('aria-pressed', String(next === 'volume'));
}

function setTallView(next) {
  tallView = next;
  const current = tallViewer.getCamera();
  tallViewer.setCamera({
    ...current,
    yaw: next === 'top' ? 0 : 0.72,
    pitch: next === 'top' ? Math.PI / 2 - 1e-7 : 0.5,
  });
  byId('tall-view-top').setAttribute('aria-pressed', String(next === 'top'));
  byId('tall-view-volume').setAttribute('aria-pressed', String(next === 'volume'));
}

function setTallMode(next) {
  tallMode = next;
  if (tallViewer && tallViewer.getSurfaceStats() !== null) {
    tallViewer.setMode(next === 'bezier' ? 'blocks' : 'patches');
  }
  for (const mode of ['bezier', 'nurbs', 'nurbs2', 'nurbs3']) {
    const button = byId(`tall-mode-${mode}`);
    button.setAttribute('aria-pressed', String(next === mode));
  }
  byId('tall-model-view').dataset.mode = next;
  const notes = {
    bezier: 'Showing the actual Bézier cell maps from the selected trivariate model.',
    nurbs: 'Showing the selected state’s exact NURBS patch map.',
    nurbs2: 'One square core and one exact annular NURBS volume. The two angular end traces of the annulus share field coefficients.',
    nurbs3: 'Showing three exact NURBS volumes with two half-ring patches and their shared patch boundaries.',
  };
  byId('tall-mode-note').textContent = notes[next];
  const isCylinder = byId('tall-state').value === 'cylinder';
  byId('tall-iga-scope').hidden = !isCylinder;
  byId('tall-iga-file').hidden = !isCylinder;
  byId('tall-torsion-file').hidden = !isCylinder;
}

async function readTallJson(path) {
  const response = await fetch(path);
  if (!response.ok) throw new Error(`Could not load ${path} (${response.status}).`);
  return response.json();
}

async function readTallModel(path) {
  const response = await fetch(path);
  if (!response.ok) throw new Error(`Could not load ${path} (${response.status}).`);
  const bytes = await response.arrayBuffer();
  const digest = await crypto.subtle.digest('SHA-256', bytes);
  const sha256 = [...new Uint8Array(digest)].map(byte => byte.toString(16).padStart(2, '0')).join('');
  return { model: JSON.parse(new TextDecoder().decode(bytes)), sha256 };
}

async function loadTallState() {
  const serial = ++tallLoadSerial;
  const stateId = byId('tall-state').value;
  const state = tallManifest.states.find(item => item.id === stateId);
  const modelPath = `./box-cylinder/tall/${stateId}.model.json`;
  const patchName = stateId === 'cylinder' && tallMode === 'nurbs2' ? 'cylinder-2'
    : stateId === 'cylinder' && tallMode === 'nurbs3' ? 'cylinder-3' : stateId;
  const patchesPath = `./box-cylinder/tall/${patchName}.patches.json`;
  const tallStatus = byId('tall-status');
  tallStatus.dataset.kind = 'loading';
  tallStatus.textContent = `Loading ${stateId} model and patch map…`;
  try {
    const [{ model, sha256 }, patches] = await Promise.all([readTallModel(modelPath), readTallJson(patchesPath)]);
    if (serial !== tallLoadSerial) return;
    if (!state || model.blocks?.length !== state.bezier_cells
      || (stateId !== 'cylinder' && patches.patches?.length !== state.nurbs_volumes)
      || (stateId === 'cylinder' && tallMode === 'nurbs2' && patches.patches?.length !== 2)
      || (stateId === 'cylinder' && tallMode === 'nurbs3' && patches.patches?.length !== 3)) {
      throw new Error(`The ${stateId} model, patch map and manifest counts disagree.`);
    }
    if (patches.schema !== 'gallery_exact_patch_map_v1') throw new Error(`The ${stateId} patch map has an unsupported schema.`);
    if (patches.source_model && patches.source_model !== `${stateId}.model.json`) {
      throw new Error(`The ${patchName} patch map references a different source model.`);
    }
    if (patches.source_model_sha256 && patches.source_model_sha256 !== sha256) {
      throw new Error(`The ${patchName} patch map does not match its source model fingerprint.`);
    }
    const savedCamera = tallCameraInitialized ? tallViewer.getCamera() : null;
    tallViewer.setModel(model, patches);
    tallViewer.setMode(tallMode === 'bezier' ? 'blocks' : 'patches');
    const camera = tallViewer.getCamera();
    tallViewer.setCamera(savedCamera || {
      ...camera, yaw: tallView === 'top' ? 0 : 0.72,
      pitch: tallView === 'top' ? Math.PI / 2 - 1e-7 : 0.5,
    });
    tallCameraInitialized = true;
    byId('tall-state-title').textContent = stateId.replaceAll('-', ' ').replace(/\b\w/g, char => char.toUpperCase());
    byId('tall-bezier-count').textContent = String(state.bezier_cells);
    byId('tall-nurbs-count').textContent = String(patches.patches.length);
    byId('tall-component-count').textContent = String(state.components);
    byId('tall-volume').textContent = `${Number(state.volume_mm3).toLocaleString('en-US', { maximumFractionDigits: 6 })} mm³`;
    const constant = state.volume_constant_mm3;
    const piCoefficient = state.volume_pi_coefficient_mm3;
    const exactVolume = piCoefficient === 0 ? String(constant)
      : constant === 0 ? `${piCoefficient}π`
      : `${constant} ${piCoefficient < 0 ? '−' : '+'} ${Math.abs(piCoefficient)}π`;
    byId('tall-volume-exact').textContent = `${exactVolume} mm³`;
    byId('tall-audit-file').href = `./box-cylinder/tall/${stateId}.audit.json`;
    byId('tall-model-file').href = modelPath;
    byId('tall-patches-file').href = patchesPath;
    byId('tall-patches-file').textContent = stateId === 'cylinder' && tallMode === 'nurbs2' ? 'NURBS 2 patch map'
      : stateId === 'cylinder' && tallMode === 'nurbs3' ? 'NURBS 3 patch map' : 'NURBS patch map';
    byId('tall-recipe-file').href = `./box-cylinder/tall/${stateId}.recipe.json`;
    byId('tall-model-view').dataset.state = stateId;
    byId('tall-model-view').dataset.mode = tallMode;
    byId('tall-contact-note').hidden = stateId !== 'xor';
    tallStatus.textContent = `Loaded ${stateId} · ${state.bezier_cells} Bézier cells · ${patches.patches.length} NURBS volumes in this map.`;
    tallStatus.dataset.kind = 'ready';
  } catch (error) {
    if (serial !== tallLoadSerial) return;
    tallStatus.textContent = `Could not load ${stateId}: ${error.message}`;
    tallStatus.dataset.kind = 'error';
  }
}


async function start() {
  viewers.result = createViewer(byId('result-view'), camera => {
    if (!syncing && viewers.classic) viewers.classic.setCamera(camera);
  });
  viewers.classic = createViewer(byId('classic-view'), camera => {
    if (!syncing && viewers.result) viewers.result.setCamera(camera);
  });
  status.textContent = 'Loading the two recorded models…';
  status.dataset.kind = 'loading';
  const [result, classic] = await Promise.all([readModel('result'), readModel('classic')]);
  for (const [name, model] of [['result', result], ['classic', classic]]) {
    if (!Array.isArray(model.blocks) || model.blocks.length !== 4
      || (model.topology?.interfaces || []).length !== 4
      || (model.boundary_faces || model.topology?.exterior || []).length !== 16) {
      throw new Error(`The ${name} model does not match the four-chart reference counts.`);
    }
  }
  const sameRationalChart = (first, second) => {
    if (JSON.stringify(first.degrees) !== JSON.stringify(second.degrees)
      || first.controls.length !== second.controls.length) return false;
    const integer = value => Number.isSafeInteger(value) ? BigInt(value) : null;
    const firstBaseWeight = integer(first.controls[0][3]);
    const secondBaseWeight = integer(second.controls[0][3]);
    if (!firstBaseWeight || !secondBaseWeight) return false;
    return first.controls.every((a, index) => {
      const b = second.controls[index];
      const aw = integer(a[3]);
      const bw = integer(b[3]);
      if (!aw || !bw) return false;
      const sameRelativeWeight = aw * secondBaseWeight === bw * firstBaseWeight;
      const samePoint = [0, 1, 2].every(axis => {
        const av = integer(a[axis]);
        const bv = integer(b[axis]);
        return av !== null && bv !== null && av * bw === bv * aw;
      });
      return sameRelativeWeight && samePoint;
    });
  };
  const unmatchedClassicCharts = [...classic.blocks];
  for (const block of result.blocks) {
    const match = unmatchedClassicCharts.findIndex(candidate => sameRationalChart(block, candidate));
    if (match < 0) throw new Error('The result and classical replay do not have pointwise-identical chart maps.');
    unmatchedClassicCharts.splice(match, 1);
  }
  if (unmatchedClassicCharts.length) {
    throw new Error('The result and classical replay do not have pointwise-identical chart maps.');
  }
  viewers.result.setModel(result);
  viewers.classic.setModel(classic);
  setMode(mode);
  setCamera(viewers.result.getCamera());
  byId('result-count').textContent = `${result.blocks.length} charts · 4 interfaces · 16 faces`;
  byId('classic-count').textContent = `${classic.blocks.length} charts · 4 interfaces · 16 faces`;
  status.textContent = 'Both trivariate models loaded · 4 charts each · 4 shared interfaces · 16 exterior faces.';
  status.dataset.kind = 'ready';
  byId('view-top').addEventListener('click', () => setView('top'));
  byId('view-volume').addEventListener('click', () => setView('volume'));
  byId('mode-surface').addEventListener('click', () => setMode('surface'));
  byId('mode-blocks').addEventListener('click', () => setMode('blocks'));
  byId('reset-views').addEventListener('click', () => {
    viewers.result.reset();
    viewers.classic.reset();
    if (view === 'top') setView('top');
  });

  tallViewer = createViewer(byId('tall-model-view'));
  tallViewer.setClip('x', 0.5);
  const manifestResponse = await fetch('./box-cylinder/tall/manifest.json');
  if (!manifestResponse.ok) throw new Error('Could not load the tall-cylinder state manifest.');
  tallManifest = await manifestResponse.json();
  const expectedStateIds = ['box', 'cylinder', 'union', 'intersection', 'a-minus-b', 'b-minus-a', 'xor'];
  if (!Array.isArray(tallManifest.states) || expectedStateIds.some(id => !tallManifest.states.some(state => state.id === id))) {
    throw new Error('The tall-cylinder manifest does not contain all seven Boolean states.');
  }
  const stateFromUrl = new URLSearchParams(location.search).get('state');
  if (expectedStateIds.includes(stateFromUrl)) byId('tall-state').value = stateFromUrl;
  previousTallState = byId('tall-state').value;
  if (previousTallState === 'cylinder') tallMode = 'nurbs2';
  const updateTallModeOptions = () => {
    const cylinder = byId('tall-state').value === 'cylinder';
    byId('tall-mode-nurbs').hidden = cylinder;
    byId('tall-mode-nurbs2').hidden = !cylinder;
    byId('tall-mode-nurbs3').hidden = !cylinder;
    setTallMode(tallMode);
  };
  updateTallModeOptions();
  byId('tall-state').addEventListener('change', () => {
    const nextState = byId('tall-state').value;
    if (nextState === 'cylinder' && previousTallState !== 'cylinder') tallMode = 'nurbs2';
    else if (nextState !== 'cylinder' && previousTallState === 'cylinder') tallMode = 'bezier';
    previousTallState = nextState;
    updateTallModeOptions();
    loadTallState();
  });
  byId('tall-mode-bezier').addEventListener('click', () => { setTallMode('bezier'); loadTallState(); });
  byId('tall-mode-nurbs').addEventListener('click', () => { setTallMode('nurbs'); loadTallState(); });
  byId('tall-mode-nurbs2').addEventListener('click', () => { setTallMode('nurbs2'); loadTallState(); });
  byId('tall-mode-nurbs3').addEventListener('click', () => { setTallMode('nurbs3'); loadTallState(); });
  byId('tall-view-top').addEventListener('click', () => setTallView('top'));
  byId('tall-view-volume').addEventListener('click', () => setTallView('volume'));
  const setTallSection = () => {
    const axis = byId('tall-section-axis').value || null;
    const fraction = Number(byId('tall-section-position').value) / 100;
    byId('tall-section-value').textContent = `${byId('tall-section-position').value}%`;
    byId('tall-model-view').dataset.sectionAxis = axis || 'off';
    byId('tall-model-view').dataset.sectionFraction = String(fraction);
    tallViewer.setClip(axis, fraction);
  };
  byId('tall-section-axis').addEventListener('change', setTallSection);
  byId('tall-section-position').addEventListener('input', setTallSection);
  byId('tall-reset').addEventListener('click', () => {
    tallViewer.reset();
    if (tallView === 'top') setTallView('top');
  });
  setTallSection();
  await loadTallState();
}

start().catch(error => {
  status.textContent = `Could not load the box-minus-cylinder comparison: ${error.message}`;
  status.dataset.kind = 'error';
});
