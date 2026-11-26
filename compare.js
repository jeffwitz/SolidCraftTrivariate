import { createViewer } from "./renderer.js?v=e5860fed8412";

const byId = id => document.getElementById(id);
const modes = ["assemble", "optimize"];
const viewers = {};
let active = "assemble";
let selected;
let cases;

function colour(index) {
  const hue = (index * 0.61803398875 + 0.04) % 1;
  return [0, 4, 2].map(phase => {
    const x = (hue * 6 + phase) % 6;
    return 0.61 - 0.31 * Math.max(-1, Math.min(x - 1, 3 - x, 1));
  });
}

function asModel(pack) {
  const blocks = pack.blocks.map(([patch, ...points], index) => ({
    id: String(index), degrees: [1,1,1],
    controls: points.map(i => pack.points[i]),
    displayColor: colour(byId("colour").value === "patch" ? patch : index),
  }));
  const interfaces = pack.interfaces.map(([a,f,b,g]) => ({first:[String(a),f],second:[String(b),g]}));
  return {schema:"sketch_layer_partition_v1",blocks,topology:{interfaces}};
}

function draw(preserveCamera) {
  for (const mode of modes) {
    const view = viewers[mode];
    const camera = preserveCamera && view.getCamera();
    const data = selected.modes[mode];
    const stats = byId(mode + "-stats");
    stats.textContent = data.patches + " NURBS volumes · " + data.cells + " Bézier cells · κ₃ " + data.kappa_3d.toFixed(2);
    view.setModel(asModel(selected.geometry[mode] ?? selected.geometry.assemble));
    view.setMode("blocks");
    if (camera) view.setCamera(camera);
  }
  align();
}

function align() {
  if (byId("sync").checked) viewers[active === "assemble" ? "optimize" : "assemble"].setCamera(viewers[active].getCamera());
}

function choose() {
  selected = cases.find(item => item.id === byId("sketch").value);
  const url = new URL(location.href);
  url.searchParams.set("model", selected.id);
  history.replaceState(null, "", url);
  draw(false);
}

async function start() {
  const response = await fetch("./comparison-data.json?v=bc4679ffa680");
  if (!response.ok) throw Error("Could not load the comparison models.");
  cases = await response.json();
  for (const mode of modes) {
    viewers[mode] = createViewer(byId(mode), camera => {
      active = mode;
      const peer = viewers[mode === "assemble" ? "optimize" : "assemble"];
      if (peer && byId("sync").checked) peer.setCamera(camera);
    });
  }
  const demos = ["bench_poly_through_octagon_high_aspect","bench_poly_blind_two_levels","bench_poly_projection_diamonds"];
  const rank = id => demos.includes(id) ? demos.indexOf(id) : demos.length;
  cases.sort((a,b) => rank(a.id)-rank(b.id) || a.title.localeCompare(b.title));
  for (const item of cases) byId("sketch").add(new Option(item.title, item.id));
  const requested = new URL(location.href).searchParams.get("model");
  byId("sketch").value = cases.some(item => item.id === requested) ? requested : demos[0];
  byId("sketch").onchange = choose;
  byId("colour").onchange = () => draw(true);
  byId("sync").onchange = align;
  byId("tool").onchange = () => Object.values(viewers).forEach(v => v.setDragTool(byId("tool").value));
  byId("reset").onclick = () => Object.values(viewers).forEach(v => v.reset());
  const clip = () => {
    byId("fraction").disabled = !byId("cut").checked;
    byId("height").textContent = byId("fraction").value + "%";
    Object.values(viewers).forEach(v => v.setClip(byId("cut").checked ? "z" : null, Number(byId("fraction").value)/100));
  };
  byId("cut").onchange = clip;
  byId("fraction").oninput = clip;
  choose();
  byId("status").textContent = "";
}
start().catch(error => { byId("status").textContent = error.message; });
