import { createViewer } from './renderer.js?v=758bdf44b4a0';

const canvas = document.querySelector('#edge-model-canvas');
const partition = document.querySelector('#edge-partition-select');
const select = document.querySelector('#edge-state-select');
const status = document.querySelector('#edge-view-state');
const reset = document.querySelector('#edge-reset-view');

try {
  const viewer = createViewer(canvas);
  let request = 0;

  async function showState() {
    const family = partition.value;
    const key = select.value;
    const volumes = family === 'four' ? 4 : 10;
    const cells = family === 'four' ? 12 : 128;
    const prefix = family === 'four' ? 'edges-four' : 'edges';
    const url = `./${prefix}-${key}.json`;
    const ticket = ++request;
    status.textContent = `Loading the ${volumes}-volume ${key} annulus…`;
    partition.disabled = true;
    select.disabled = true;
    try {
      const response = await fetch(url);
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const model = await response.json();
      if (ticket !== request) return;
      viewer.setModel(model);
      status.textContent = `Loaded ${family === 'four' ? 'four-sector' : 'prepared'} ${key} annulus · ${volumes} volumes · ${cells} Bézier cells.`;
      canvas.dataset.ready = `${family}-${key}`;
    } catch (error) {
      if (ticket !== request) return;
      status.textContent = `Could not load the recorded ${family}-sector ${key} geometry: ${error.message}`;
    } finally {
      if (ticket === request) {
        partition.disabled = false;
        select.disabled = false;
      }
    }
  }

  select.addEventListener('change', showState);
  partition.addEventListener('change', showState);
  reset.addEventListener('click', () => {
    canvas.focus();
    canvas.dispatchEvent(new KeyboardEvent('keydown', { key: 'Home', bubbles: true }));
  });
  showState();
} catch (error) {
  status.textContent = `Interactive WebGL is unavailable: ${error.message}`;
  partition.disabled = true;
  select.disabled = true;
  reset.disabled = true;
}
