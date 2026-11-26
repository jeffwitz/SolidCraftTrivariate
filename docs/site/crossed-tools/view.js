import {createViewer} from '../renderer.js?v=758bdf44b4a0';

const el = id => document.getElementById(id);
const operation = el('operation');
const rise = el('rise');
const status = el('model-state');
let version = 0;

async function start() {
  const viewer = createViewer(el('crossed-model'));
  const response = await fetch('./crossed-tools/verification.txt');
  if (!response.ok) throw Error('The verification record could not be loaded.');
  const rows = (await response.text()).trim().split(/\r?\n/).map(line => line.split(/\s+/));
  async function update() {
    const request = ++version;
    const row = rows.find(([height, op]) => height === rise.value && op === operation.value);
    const path = './crossed-tools/' + operation.value + '-' + rise.value + '.json';
    status.hidden = false;
    status.textContent = 'Loading model…';
    try {
      const r = await fetch(path);
      if (!r.ok || !row) throw Error('The result or verification record could not be loaded.');
      const model = await r.json();
      if (request !== version) return;
      if (model.operation !== operation.value || model.blocks.length !== Number(row[2])) throw Error('The model and verification record disagree.');
      viewer.setModel(model);
      viewer.setMode('blocks');
      el('result-title').textContent = operation.selectedOptions[0].text;
      el('counts').textContent = row[2] + ' volume charts · ' + Number(row[4]).toLocaleString('en-US', {maximumFractionDigits:6}) + ' mm³';
      el('model-file').href = path;
      el('recipe-file').href = './crossed-tools/rise-' + rise.value + '.recipe.json';
      status.hidden = true;
    } catch (error) { if (request === version) status.textContent = error.message; }
  }
  operation.onchange = rise.onchange = update;
  el('reset').onclick = () => viewer.reset();
  await update();
}
start().catch(error => { status.hidden = false; status.textContent = error.message; });
