import { createViewer, partitionColour } from "./renderer.js?v=758bdf44b4a0";

const byId = id => document.getElementById(id);
const modes = ["assemble", "optimize"];
const viewers = {};
let active = "assemble";
let selected;
let cases;
let experimentById = new Map();
let rightVariant = "optimize";
let renderSerial = 0;
const supplementedModels = new Map();

function asModel(pack) {
  const patchRecords = new Map();
  const blocks = pack.blocks.map(([patch, ...points], index) => ({
    id: String(index), degrees: [1,1,1],
    controls: points.map(i => pack.points[i]),
    displayColor: partitionColour(index),
  }));
  pack.blocks.forEach(([patch], index) => {
    const id = String(patch);
    if (!patchRecords.has(id)) patchRecords.set(id, { id, cells: [] });
    patchRecords.get(id).cells.push(String(index));
  });
  const interfaces = pack.interfaces.map(([a,f,b,g]) => ({first:[String(a),f],second:[String(b),g]}));
  return {schema:"sketch_layer_partition_v1",blocks,nurbs_patches:[...patchRecords.values()],topology:{interfaces}};
}

async function loadSupplementedModel(entry) {
  if (!supplementedModels.has(entry.id)) {
    const promise = Promise.all([fetch(entry.model), fetch(entry.representation)])
      .then(async ([modelResponse, representationResponse]) => {
        if (!modelResponse.ok || !representationResponse.ok) throw Error(`Could not load the exact NURBS cover for ${entry.id}.`);
        const bytes = await modelResponse.arrayBuffer();
        const model = JSON.parse(new TextDecoder().decode(bytes));
        const representation = await representationResponse.json();
        const digest = await crypto.subtle.digest('SHA-256', bytes);
        const sourceHash = [...new Uint8Array(digest)].map(value => value.toString(16).padStart(2, '0')).join('');
        if (representation.source_model_sha256 && representation.source_model_sha256 !== sourceHash) {
          throw Error(`The NURBS cover does not match the published native model for ${entry.id}.`);
        }
        return [model, representation];
      });
    supplementedModels.set(entry.id, promise);
  }
  return supplementedModels.get(entry.id);
}

async function draw(preserveCamera) {
  if (!selected) return;
  const serial = ++renderSerial;
  const colourMode = byId("colour").value === "patch" ? "patches" : "blocks";
  const leftCamera = preserveCamera && viewers.assemble.getCamera();
  const rightCamera = preserveCamera && viewers.optimize.getCamera();
  const assembled = selected.modes.assemble;
  viewers.assemble.setModel(asModel(selected.geometry.assemble));
  viewers.assemble.setMode(colourMode);
  byId("assemble-stats").textContent = assembled.patches + " NURBS volumes · " + assembled.cells + " Bézier cells · κ₃ " + assembled.kappa_3d.toFixed(2);

  if (rightVariant === "cover") {
    const entry = experimentById.get(selected.id);
    if (!entry?.representation) throw Error(`No exact NURBS cover is published for ${selected.id}.`);
    const [model, representation] = await loadSupplementedModel(entry);
    if (serial !== renderSerial) return;
    viewers.optimize.setModel(model, representation);
    viewers.optimize.setMode(colourMode);
    const cells = new Set(representation.patches.flatMap(patch => patch.cells || patch.source_blocks || []));
    byId("right-title").textContent = "Published exact NURBS cover";
    byId("right-subtitle").textContent = "Exact patch extraction and joins";
    byId("optimize-stats").textContent = representation.patches.length + " merged NURBS patches · " + cells.size + " Bézier cells · paired candidate κ₃ assembled/optimized " + selected.modes.assemble.kappa_3d.toFixed(2) + "/" + selected.modes.optimize.kappa_3d.toFixed(2);
  } else {
    const optimized = selected.modes.optimize;
    viewers.optimize.setModel(asModel(selected.geometry.optimize ?? selected.geometry.assemble));
    viewers.optimize.setMode(colourMode);
    byId("right-title").textContent = "Native optimized";
    byId("right-subtitle").textContent = "Reselected separators";
    byId("optimize-stats").textContent = optimized.patches + " NURBS volumes · " + optimized.cells + " Bézier cells · κ₃ " + optimized.kappa_3d.toFixed(2);
  }
  if (leftCamera) viewers.assemble.setCamera(leftCamera);
  if (rightCamera) viewers.optimize.setCamera(rightCamera);
  align();
}

function align() {
  if (byId("sync").checked) viewers[active === "assemble" ? "optimize" : "assemble"].setCamera(viewers[active].getCamera());
}

async function choose() {
  selected = cases.find(item => item.id === byId("sketch").value);
  const url = new URL(location.href);
  url.searchParams.set("model", selected.id);
  history.replaceState(null, "", url);
  const entry = experimentById.get(selected.id);
  const hasCover = Boolean(entry?.representation);
  byId("right-variant-control").hidden = !hasCover;
  rightVariant = hasCover ? "cover" : "optimize";
  byId("right-variant").value = rightVariant;
  await draw(false);
}

async function start() {
  const response = await fetch("./comparison-data.json?v=bc4679ffa680");
  if (!response.ok) throw Error("Could not load the comparison models.");
  cases = await response.json();
  const requested = new URL(location.href).searchParams.get("model");
  if (requested && !cases.some(item => item.id === requested)) {
    byId("status").textContent = `No paired partition comparison is published for “${requested}”.`;
    byId("sketch").disabled = true;
    document.querySelector(".comparison").hidden = true;
    byId("right-variant-control").hidden = true;
    return;
  }
  const catalogResponse = await fetch(byId("compare-catalog-url").href);
  if (catalogResponse.ok) {
    const catalog = await catalogResponse.json();
    experimentById = new Map(catalog.map(item => [item.id, item]));
  }
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
  byId("sketch").value = requested || demos[0];
  byId("sketch").onchange = choose;
  byId("colour").onchange = () => draw(true);
  byId("right-variant").onchange = () => { rightVariant = byId("right-variant").value; draw(true); };
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
  await choose();
  byId("status").textContent = "";
}
start().catch(error => { byId("status").textContent = error.message; });
