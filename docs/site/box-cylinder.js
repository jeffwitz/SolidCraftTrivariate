import { createViewer } from './renderer.js?v=758bdf44b4a0';

const byId = id => document.getElementById(id);
const status = byId('box-cylinder-status');
const viewers = {};
let syncing = false;
let view = 'volume';
let mode = 'blocks';

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
}

start().catch(error => {
  status.textContent = `Could not load the box-minus-cylinder comparison: ${error.message}`;
  status.dataset.kind = 'error';
});
