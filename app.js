const byId = (id) => document.getElementById(id);
const ui = {
  canvas: byId('model-canvas'), state: byId('canvas-state'), wrap: byId('canvas-wrap'),
  selector: byId('experiment-select'), description: byId('experiment-description'),
  viewTitle: byId('view-title'), metrics: byId('metrics-list'), model: byId('model-link'), source: byId('source-link'), step: byId('step-link'), hex8: byId('hex8-link'), audit: byId('audit-link'), figure: byId('figure-link'),
  axis: byId('clip-axis'), fraction: byId('clip-fraction'), value: byId('clip-value'), reset: byId('reset-view'),
  sectionButton: byId('s3-face-pocket'),
  fillet: byId('fillet-step'), chamfer: byId('chamfer-step'),
};

let experiments = [];
let viewer = null;
let mode = 'surface';
let selectionVersion = 0;
let modelReady = false;
let viewerError = 'Le moteur de visualisation n’est pas disponible.';

function showState(message, kind = 'info') {
  ui.state.textContent = message;
  ui.state.dataset.kind = kind;
  ui.state.hidden = false;
}

function hideState() { ui.state.hidden = true; }

function safeLocalUrl(path) {
  if (typeof path !== 'string' || !path.trim()) return null;
  try {
    const url = new URL(path, document.baseURI);
    if (url.origin !== location.origin) return null;
    return url.href;
  } catch { return null; }
}

function safeLinkUrl(path) {
  if (typeof path !== 'string' || !path.trim()) return null;
  try {
    const url = new URL(path, document.baseURI);
    const isLocal = url.origin === location.origin;
    const isProjectGithub = url.protocol === 'https:' && url.hostname === 'github.com' && url.pathname.startsWith('/jeffwitz/SolidCraftTrivariate/');
    return isLocal || isProjectGithub ? url.href : null;
  } catch { return null; }
}

function setLink(link, path) {
  const url = safeLinkUrl(path);
  if (!url) { link.hidden = true; link.removeAttribute('href'); return; }
  link.href = url;
  if (new URL(url).origin !== location.origin) { link.target = '_blank'; link.rel = 'noopener noreferrer'; }
  else { link.removeAttribute('target'); link.removeAttribute('rel'); }
  link.hidden = false;
}

function renderMetrics(metrics) {
  ui.metrics.replaceChildren();
  if (!Array.isArray(metrics) || metrics.length === 0) {
    const empty = document.createElement('div');
    empty.className = 'metric-empty';
    empty.textContent = 'Aucune mesure publiée pour cette expérience.';
    ui.metrics.append(empty);
    return;
  }
  for (const metric of metrics) {
    if (!metric || typeof metric.label !== 'string' || metric.value == null) continue;
    const row = document.createElement('div');
    const label = document.createElement('dt');
    const value = document.createElement('dd');
    label.textContent = metric.label;
    value.textContent = String(metric.value);
    row.append(label, value);
    ui.metrics.append(row);
  }
}

async function readModel(experiment) {
  const url = safeLocalUrl(experiment.model);
  if (!url) throw new Error('Chemin du modèle absent ou hors du site.');
  const response = await fetch(url);
  if (!response.ok) throw new Error(`Le modèle n’a pas pu être chargé (${response.status}).`);
  return response.json();
}

async function selectExperiment(id) {
  const experiment = experiments.find((item) => item.id === id);
  if (!experiment) return;
  const version = ++selectionVersion;
  modelReady = false;
  ui.selector.value = experiment.id;
  ui.viewTitle.textContent = experiment.title || experiment.id;
  ui.description.textContent = experiment.description || '';
  renderMetrics(experiment.metrics);
  const boundary = experiment.representation === 'boundary';
  const blockButton = document.querySelector('[data-mode="blocks"]');
  blockButton.disabled = boundary;
  ui.axis.disabled = boundary;
  if (boundary) {
    ui.axis.value = '';
    if (mode === 'blocks') document.querySelector('[data-mode="surface"]').click();
  }
  applyClip();
  ui.source.firstElementChild.textContent = experiment.sourceLabel || 'Recette source';
  ui.step.firstElementChild.textContent = experiment.stepLabel || 'Modèle STEP';
  ui.model.firstElementChild.textContent = boundary ? 'Surface triangulée (JSON)' : 'Modèle trivarié (JSON)';
  setLink(ui.model, experiment.model === experiment.source ? null : experiment.model);
  setLink(ui.source, experiment.source);
  setLink(ui.step, experiment.step);
  setLink(ui.hex8, experiment.hex8);
  setLink(ui.audit, experiment.audit);
  setLink(ui.figure, experiment.figure);
  setLink(ui.fillet, experiment.filletStep);
  setLink(ui.chamfer, experiment.chamferStep);
  const section = experiment.section;
  const hasSection = section && ['x', 'y', 'z'].includes(section.axis) && Number.isFinite(section.fraction) && section.fraction >= 0 && section.fraction <= 1;
  ui.sectionButton.hidden = !hasSection;
  ui.sectionButton.textContent = 'Voir la poche en coupe';
  if (location.hash === '#s3-face-pocket' && hasSection) requestAnimationFrame(() => ui.sectionButton.scrollIntoView({ block: 'center' }));
  if (!viewer) {
    showState(`${viewerError} Les mesures et les liens de source restent accessibles.`, 'error');
    return;
  }
  showState('Chargement du modèle…');
  try {
    const model = await readModel(experiment);
    if (version !== selectionVersion) return;
    viewer.setModel(model);
    viewer.setMode(mode);
    modelReady = true;
    applyClip();
    hideState();
  } catch (error) {
    if (version !== selectionVersion) return;
    showState(`${error.message} Vérifiez que le site est servi en HTTP et que le modèle est publié.`, 'error');
  }
}

function applyClip() {
  const axis = ui.axis.value || null;
  const fraction = Number(ui.fraction.value) / 100;
  ui.fraction.disabled = !axis;
  ui.value.textContent = `${ui.fraction.value} %`;
  if (modelReady) viewer.setClip(axis, fraction);
}

function initControls() {
  ui.selector.addEventListener('change', () => {
    const url = new URL(location.href);
    url.searchParams.set('model', ui.selector.value);
    url.hash = '';
    history.replaceState(null, '', url);
    selectExperiment(ui.selector.value);
  });
  ui.axis.addEventListener('change', applyClip);
  ui.fraction.addEventListener('input', applyClip);
  ui.reset.addEventListener('click', () => {
    ui.axis.value = '';
    ui.fraction.value = '50';
    applyClip();
    viewer?.reset();
  });
  ui.sectionButton.addEventListener('click', () => {
    const section = activeExperiment()?.section;
    if (!section || !['x', 'y', 'z'].includes(section.axis) || !Number.isFinite(section.fraction)) return;
    ui.axis.value = section.axis;
    ui.fraction.value = String(Math.round(section.fraction * 100));
    applyClip();
    ui.canvas.focus();
  });
  document.querySelectorAll('[data-mode]').forEach((button) => button.addEventListener('click', () => {
    mode = button.dataset.mode;
    document.querySelectorAll('[data-mode]').forEach((candidate) => {
      const selected = candidate === button;
      candidate.classList.toggle('is-active', selected);
      candidate.setAttribute('aria-pressed', String(selected));
    });
    if (modelReady) viewer.setMode(mode);
  }));
}

function activeExperiment() { return experiments.find((item) => item.id === ui.selector.value); }

async function start() {
  initControls();
  try {
    const response = await fetch(new URL('./experiments.json', document.baseURI));
    if (!response.ok) throw new Error(`Le catalogue d’expériences est indisponible (${response.status}).`);
    const data = await response.json();
    if (!Array.isArray(data) || data.length === 0) throw new Error('Le catalogue ne contient aucune expérience.');
    experiments = data.filter((item) => item && typeof item.id === 'string' && typeof item.model === 'string');
    if (experiments.length === 0) throw new Error('Aucune expérience ne contient de modèle exploitable.');
    ui.selector.replaceChildren();
    const groups = new Map();
    for (const experiment of experiments) {
      const option = document.createElement('option');
      option.value = experiment.id;
      option.textContent = experiment.title || experiment.id;
      const group = experiment.group || 'Modèles';
      if (!groups.has(group)) {
        const container = document.createElement('optgroup');
        container.label = group;
        groups.set(group, container);
        ui.selector.append(container);
      }
      groups.get(group).append(option);
    }
    ui.selector.disabled = false;
    try {
      const { createViewer } = await import('./renderer.js');
      viewer = createViewer(ui.canvas);
      if (!viewer || typeof viewer.setModel !== 'function') throw new Error('L’API du moteur de visualisation est incomplète.');
      hideState();
    } catch (error) {
      viewerError = /webgl|context|gpu/i.test(error.message) ? 'WebGL n’est pas disponible dans ce navigateur.' : `Moteur WebGL indisponible : ${error.message}`;
      showState(`${viewerError} Les mesures et les liens de source restent accessibles.`, 'error');
    }
    const requested = new URL(location.href).searchParams.get('model');
    const initial = experiments.find(entry => entry.id === requested) || experiments[0];
    await selectExperiment(initial.id);
  } catch (error) {
    showState(`${error.message} Le catalogue attendu est ./experiments.json.`, 'error');
    ui.viewTitle.textContent = 'Aucun modèle chargé';
    ui.description.textContent = '';
    renderMetrics([]);
  }
}

start();
