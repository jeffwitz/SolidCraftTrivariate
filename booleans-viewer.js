import { createViewer } from './renderer.js?v=758bdf44b4a0';

const byId = (id) => document.getElementById(id);
const sectionSelect = byId('boolean-section');
const operationSelect = byId('boolean-operation');
const status = byId('boolean-status');
const titles = { union: 'Union · A ∪ B', intersection: 'Intersection · A ∩ B', difference: 'Difference · A − B' };
const viewers = {
  a: createViewer(byId('operand-a')),
  b: createViewer(byId('operand-b')),
  result: createViewer(byId('boolean-result')),
};
for (const viewer of Object.values(viewers)) viewer.setClip(null, 0.5);

const slug = (section, role) => `booleans-${section}-${role}`;
const url = (section, role, ext) => `./${slug(section, role)}.${ext}`;

async function readJson(path) {
  const response = await fetch(path, { cache: 'force-cache' });
  if (!response.ok) throw new Error(`Could not load ${path} (${response.status}).`);
  return response.json();
}

function renderAudit(element, audit) {
  element.textContent = `${audit.method.replaceAll('_', ' ')} · ${audit.blocks} cells · ${audit.components} connected component${audit.components === 1 ? '' : 's'}`;
}

function setArtifactLinks(section, operation) {
  const box = byId('boolean-artifacts');
  box.replaceChildren();
  const text = document.createElement('span');
  text.textContent = 'Selected result files';
  box.append(text);
  for (const [label, role, ext] of [
    ['Recipe JSON', operation, 'recipe.json'],
    ['Native model JSON', operation, 'model.json'],
    ['Audit JSON', operation, 'audit.json'],
  ]) {
    const link = document.createElement('a');
    link.href = url(section, role, ext);
    link.textContent = label;
    link.download = '';
    box.append(link);
  }
}

async function loadSelection() {
  const section = sectionSelect.value;
  const operation = operationSelect.value;
  status.textContent = 'Loading and auditing native models…';
  status.dataset.kind = 'loading';
  try {
    const [aModel, aAudit, bModel, bAudit, resultModel, resultAudit] = await Promise.all([
      readJson(url(section, 'a', 'model.json')),
      readJson(url(section, 'a', 'audit.json')),
      readJson(url(section, 'b', 'model.json')),
      readJson(url(section, 'b', 'audit.json')),
      readJson(url(section, operation, 'model.json')),
      readJson(url(section, operation, 'audit.json')),
    ]);
    viewers.a.setModel(aModel);
    viewers.b.setModel(bModel);
    viewers.result.setModel(resultModel);
    byId('audit-a').textContent = `10 mm axial interval · ${aAudit.blocks} cells · volume ${aAudit.ideal_volume_mm3} mm³`;
    byId('audit-b').textContent = `10 mm axial interval · ${bAudit.blocks} cells · volume ${bAudit.ideal_volume_mm3} mm³`;
    renderAudit(byId('result-audit'), resultAudit);
    byId('result-title').textContent = titles[operation];
    byId('result-name').textContent = titles[operation];
    byId('result-volume').textContent = `Volume ${Number(resultAudit.ideal_volume_mm3).toLocaleString('en-US')} mm³`;
    setArtifactLinks(section, operation);
    status.textContent = 'Audited model loaded · operands and result are native trivariate models.';
    status.dataset.kind = 'ready';
  } catch (error) {
    status.textContent = `Model unavailable: ${error.message}`;
    status.dataset.kind = 'error';
  }
}

sectionSelect.addEventListener('change', loadSelection);
operationSelect.addEventListener('change', loadSelection);
byId('section-toggle').addEventListener('change', (event) => {
  viewers.result.setClip(event.target.checked ? 'z' : null, 0.5);
});
byId('reset-views').addEventListener('click', () => Object.values(viewers).forEach(viewer => viewer.reset()));
loadSelection();
