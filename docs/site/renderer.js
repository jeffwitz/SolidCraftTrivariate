/* Native WebGL renderer for rational trivariate model JSON. */
import { surfaceMesh } from './triangle-surface.js?v=509f11476ad1';

const VERTEX_SHADER = `
attribute vec3 aPosition;
attribute vec3 aNormal;
attribute vec3 aColor;
uniform mat4 uMvp;
uniform mat4 uModel;
uniform vec4 uClip;
uniform float uUseClip;
varying vec3 vNormal;
varying vec3 vColor;
varying float vSide;
void main() {
  vec4 world = vec4(aPosition, 1.0);
  vSide = dot(uClip.xyz, world.xyz) - uClip.w;
  gl_Position = uMvp * uModel * world;
  vNormal = mat3(uModel) * aNormal;
  vColor = aColor;
}`;

const FRAGMENT_SHADER = `
precision mediump float;
varying vec3 vNormal;
varying vec3 vColor;
varying float vSide;
uniform float uUseClip;
void main() {
  if (uUseClip > 0.5 && vSide > 0.00001) discard;
  vec3 n = normalize(vNormal);
  if (!gl_FrontFacing) n = -n;
  vec3 lightDirection = normalize(vec3(-0.45, 0.72, 0.8));
  float diffuse = max(dot(n, lightDirection), 0.0);
  float lighting = 0.48 + 0.52 * diffuse;
  gl_FragColor = vec4(vColor * lighting, 1.0);
}`;

const SURFACE_COLOR = [0.78, 0.67, 0.52];
const SECTION_COLOR = [0.20, 0.48, 0.46];
const DEFAULT_CAMERA = { yaw: 0.72, pitch: 0.5, distance: 2, panX: 0, panY: 0 };
const CAMERA_FILL = 0.85;
const GRID = 12;
const CAP_GRID = 12;
const FACE_SIDES = [
  ['u0', 0, 0], ['u1', 0, 1],
  ['v0', 1, 0], ['v1', 1, 1],
  ['w0', 2, 0], ['w1', 2, 1],
];
const CUBE_CORNERS = [
  [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
];
const CUBE_TETRAHEDRA = [
  [0, 5, 1, 6], [0, 1, 2, 6], [0, 2, 3, 6],
  [0, 3, 7, 6], [0, 7, 4, 6], [0, 4, 5, 6],
];
const TETRA_EDGES = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];

const clamp = (value, low, high) => Math.max(low, Math.min(high, value));

function normalize(vector) {
  const length = Math.hypot(...vector) || 1;
  return vector.map(value => value / length);
}

function cross(a, b) {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}

function subtract(a, b) {
  return a.map((value, index) => value - b[index]);
}

function multiplyMatrices(a, b) {
  const result = new Float32Array(16);
  for (let column = 0; column < 4; column++) {
    for (let row = 0; row < 4; row++) {
      result[column * 4 + row] = a[row] * b[column * 4]
        + a[4 + row] * b[column * 4 + 1]
        + a[8 + row] * b[column * 4 + 2]
        + a[12 + row] * b[column * 4 + 3];
    }
  }
  return result;
}

function perspective(fieldOfView, aspect, near, far) {
  const f = 1 / Math.tan(fieldOfView / 2);
  const result = new Float32Array(16);
  result[0] = f / aspect;
  result[5] = f;
  result[10] = (far + near) / (near - far);
  result[11] = -1;
  result[14] = (2 * far * near) / (near - far);
  return result;
}

function lookAt(eye, target, up) {
  const z = normalize(subtract(eye, target));
  const x = normalize(cross(up, z));
  const y = cross(z, x);
  return new Float32Array([
    x[0], y[0], z[0], 0,
    x[1], y[1], z[1], 0,
    x[2], y[2], z[2], 0,
    -x.reduce((sum, value, i) => sum + value * eye[i], 0),
    -y.reduce((sum, value, i) => sum + value * eye[i], 0),
    -z.reduce((sum, value, i) => sum + value * eye[i], 0),
    1,
  ]);
}

function basis(degree, t) {
  if (degree === 1) return [1 - t, t];
  if (degree === 2) return [(1 - t) ** 2, 2 * t * (1 - t), t ** 2];
  if (degree === 3) return [(1 - t) ** 3, 3 * t * (1 - t) ** 2, 3 * t ** 2 * (1 - t), t ** 3];
  throw new TypeError(`Unsupported Bézier degree ${degree}`);
}

function splineBasisVector(degree, parameter, knots, count) {
  const low = knots[degree];
  const high = knots[count];
  const value = low + (high - low) * parameter;
  const memo = new Map();
  function recurse(i, p) {
    const key = `${i}:${p}`;
    if (memo.has(key)) return memo.get(key);
    let result;
    if (p === 0) {
      result = (value >= knots[i] && value < knots[i + 1])
        || (parameter === 1 && i === count - 1 && value === high) ? 1 : 0;
    } else {
      const leftDenominator = knots[i + p] - knots[i];
      const rightDenominator = knots[i + p + 1] - knots[i + 1];
      const left = leftDenominator ? ((value - knots[i]) / leftDenominator) * recurse(i, p - 1) : 0;
      const right = rightDenominator ? ((knots[i + p + 1] - value) / rightDenominator) * recurse(i + 1, p - 1) : 0;
      result = left + right;
    }
    memo.set(key, result);
    return result;
  }
  return Array.from({ length: count }, (_, index) => recurse(index, degree));
}

export function partitionColour(index) {
  const hue = (index * 0.61803398875 + 0.04) % 1;
  return [0, 4, 2].map(phase => {
    const x = (hue * 6 + phase) % 6;
    return 0.61 - 0.31 * Math.max(-1, Math.min(x - 1, 3 - x, 1));
  });
}

function partitionColor(block, mode) {
  return mode === 'patches' && block.patchColor ? block.patchColor : (block.displayColor ?? blockColor(block.id));
}

function blockColor(identifier) {
  let hash = 2166136261;
  for (let index = 0; index < identifier.length; index++) {
    hash ^= identifier.charCodeAt(index);
    hash = Math.imul(hash, 16777619);
  }
  const hue = (hash >>> 0) % 360;
  const chroma = 0.28;
  const x = chroma * (1 - Math.abs((hue / 60) % 2 - 1));
  const m = 0.58;
  let rgb;
  if (hue < 60) rgb = [chroma, x, 0];
  else if (hue < 120) rgb = [x, chroma, 0];
  else if (hue < 180) rgb = [0, chroma, x];
  else if (hue < 240) rgb = [0, x, chroma];
  else if (hue < 300) rgb = [x, 0, chroma];
  else rgb = [chroma, 0, x];
  return rgb.map(value => value + m);
}

function validateModel(model) {
  if (model && Array.isArray(model.vertices) && Array.isArray(model.triangles)) {
    if (model.vertices.length % 3 || model.triangles.length % 3
      || model.vertices.some(value => !Number.isFinite(value))
      || model.triangles.some(index => !Number.isInteger(index) || index < 0 || index >= model.vertices.length / 3)) {
      throw new TypeError('Triangle surface has invalid vertices or triangle indices');
    }
    return;
  }
  model = normalizeTrivariateModel(model);
  if (!model || typeof model !== 'object' || !Array.isArray(model.blocks) || model.blocks.length === 0) {
    throw new TypeError('Expected a trivariate model JSON object with a non-empty blocks array');
  }
  const identifiers = new Set();
  for (const block of model.blocks) {
    if (!block || typeof block.id !== 'string' || identifiers.has(block.id)) {
      throw new TypeError('Each block must have a unique string id');
    }
    identifiers.add(block.id);
    if (!Array.isArray(block.degrees) || block.degrees.length !== 3
      || block.degrees.some(degree => !Number.isInteger(degree) || degree < 1 || degree > 3)) {
      throw new TypeError(`Block ${block.id} has unsupported degrees; supported degrees are 1 through 3`);
    }
    const dimensions = block.shape || block.degrees.map(value => value + 1);
    const expected = dimensions.reduce((product, value) => product * value, 1);
    if (!Array.isArray(block.controls) || block.controls.length !== expected) {
      throw new TypeError(`Block ${block.id} has an invalid control count`);
    }
    if (block.knots && (!Array.isArray(block.knots) || block.knots.length !== 3
      || block.knots.some((knots, axis) => !Array.isArray(knots)
        || knots.length !== dimensions[axis] + block.degrees[axis] + 1
        || knots.some((value, index) => !Number.isFinite(value) || (index > 0 && value < knots[index - 1]))
        || !(knots[block.degrees[axis]] < knots[dimensions[axis]])))) {
      throw new TypeError(`Block ${block.id} has invalid knot vectors`);
    }
    for (const control of block.controls) {
      if (!Array.isArray(control) || control.length !== 4
        || control.some(value => !Number.isFinite(value)) || control[3] <= 0) {
        throw new TypeError(`Block ${block.id} has invalid homogeneous controls`);
      }
    }
  }
  const blockIds = new Set(model.blocks.map(block => block.id));
  const exterior = model.boundary_faces || model.topology?.boundary_faces || model.topology?.exterior;
  if (exterior !== undefined && exterior !== null) {
    if (!Array.isArray(exterior)) throw new TypeError('Certified exterior faces must be an array');
    const faceKeys = new Set();
    for (const face of exterior) {
      const block = Array.isArray(face) ? face[0] : face?.source_block ?? face?.block;
      const side = Array.isArray(face) ? face[1] : face?.side;
      const key = `${block}:${side}`;
      if (typeof block !== 'string' || !blockIds.has(block)
        || !['u0', 'u1', 'v0', 'v1', 'w0', 'w1'].includes(side) || faceKeys.has(key)) {
        throw new TypeError('Certified exterior contains an invalid or repeated block face');
      }
      if (!Array.isArray(face) && face.status && face.status !== 'extracted_exactly_from_model') {
        throw new TypeError('Exterior face has an unsupported extraction status');
      }
      faceKeys.add(key);
    }
    const expected = model.blocks.length * 6 - (model.topology?.interfaces || []).length * 2;
    if (model.topology?.interfaces && faceKeys.size !== expected) {
      throw new TypeError(`Exterior face count ${faceKeys.size} does not close the interface graph (${expected} expected)`);
    }
  }
  for (const item of model.topology?.interfaces || []) {
    for (const side of [item.first, item.second]) {
      if (!Array.isArray(side) || side.length !== 2 || !blockIds.has(side[0])
        || !['u0', 'u1', 'v0', 'v1', 'w0', 'w1'].includes(side[1])) {
        throw new TypeError('Model contains an invalid block interface');
      }
    }
  }
}

function normalizeTrivariateModel(model) {
  if (!model || typeof model !== 'object') return model;
  if (Array.isArray(model.blocks) && model.blocks.every(block => Array.isArray(block.controls))) return model;
  if (Array.isArray(model.blocks) && Array.isArray(model.degrees) && Array.isArray(model.knots)) {
    return { ...model, blocks: model.blocks.map(block => ({
      ...block, degrees: model.degrees, knots: model.knots,
      shape: model.knots.map((knots, axis) => knots.length - model.degrees[axis] - 1),
      controls: block.xyzw,
    })) };
  }
  if (Array.isArray(model.charts)) {
    return { ...model, blocks: model.charts.map(chart => ({
      ...chart, controls: chart.xyzw,
    })) };
  }
  return model;
}

function isTriangleModel(model) {
  return Boolean(model && Array.isArray(model.vertices) && Array.isArray(model.triangles));
}

function parameterStations(block, axis, minimumDivisions = GRID) {
  if (!block.knots) return Array.from({ length: minimumDivisions + 1 }, (_, index) => index / minimumDivisions);
  const degree = block.degrees[axis];
  const count = (block.shape || block.degrees.map(value => value + 1))[axis];
  const knots = block.knots[axis];
  const low = knots[degree];
  const high = knots[count];
  const spans = [...new Set(knots.filter(value => value >= low && value <= high))].sort((a, b) => a - b);
  const minimumPerSpan = axis < 2 ? 16 : 6;
  const divisionsPerSpan = Math.max(minimumPerSpan, Math.ceil(minimumDivisions / Math.max(1, spans.length - 1)));
  const stations = [0];
  for (let span = 0; span < spans.length - 1; span++) {
    for (let step = 1; step <= divisionsPerSpan; step++) {
      stations.push((spans[span] + (spans[span + 1] - spans[span]) * step / divisionsPerSpan - low) / (high - low));
    }
  }
  return stations;
}

function boundsOf(model) {
  const low = [Infinity, Infinity, Infinity];
  const high = [-Infinity, -Infinity, -Infinity];
  function include(point) {
    for (let axis = 0; axis < 3; axis++) {
      low[axis] = Math.min(low[axis], point[axis]);
      high[axis] = Math.max(high[axis], point[axis]);
    }
  }
  for (const block of model.blocks) {
    if (block.knots) {
      const evaluate = makeEvaluator(block);
      const axes = [0, 1, 2].map(axis => parameterStations(block, axis, GRID));
      for (let fixedAxis = 0; fixedAxis < 3; fixedAxis++) {
        const otherAxes = [0, 1, 2].filter(axis => axis !== fixedAxis);
        for (const side of [0, 1]) {
          for (const first of axes[otherAxes[0]]) for (const second of axes[otherAxes[1]]) {
            const parameters = [0, 0, 0];
            parameters[fixedAxis] = side;
            parameters[otherAxes[0]] = first;
            parameters[otherAxes[1]] = second;
            include(evaluate(...parameters));
          }
        }
      }
      continue;
    }
    for (const control of block.controls) {
      include(control.slice(0, 3).map(value => value / control[3]));
    }
  }
  return [low, high];
}

function makeHomogeneousEvaluator(block) {
  const degree = block.degrees;
  const dimensions = block.shape || degree.map(value => value + 1);
  const controls = block.controls;
  return (u, v, w) => {
    const parameters = [u, v, w];
    const weights = degree.map((value, axis) => block.knots
      ? splineBasisVector(value, parameters[axis], block.knots[axis], dimensions[axis])
      : basis(value, parameters[axis]));
    const [bu, bv, bw] = weights;
    const h = [0, 0, 0, 0];
    let index = 0;
    for (let i = 0; i < dimensions[0]; i++) {
      for (let j = 0; j < dimensions[1]; j++) {
        for (let k = 0; k < dimensions[2]; k++) {
          const weight = bu[i] * bv[j] * bw[k];
          const control = controls[index++];
          for (let component = 0; component < 4; component++) h[component] += weight * control[component];
        }
      }
    }
    if (!(h[3] > 0) || !Number.isFinite(h[3])) throw new TypeError(`Block ${block.id} evaluates with an invalid weight`);
    return h;
  };
}

function makeEvaluator(block) {
  const evaluateHomogeneous = makeHomogeneousEvaluator(block);
  return (...parameters) => {
    const h = evaluateHomogeneous(...parameters);
    return [h[0] / h[3], h[1] / h[3], h[2] / h[3]];
  };
}

function triangleNormal(a, b, c) {
  return normalize(cross(subtract(b, a), subtract(c, a)));
}

function appendVertex(output, point, normal, color) {
  output.push(point[0], point[1], point[2], normal[0], normal[1], normal[2], color[0], color[1], color[2]);
}

function appendTriangle(output, a, b, c, color) {
  const normal = triangleNormal(a, b, c);
  appendVertex(output, a, normal, color);
  appendVertex(output, b, normal, color);
  appendVertex(output, c, normal, color);
}

function buildBlockGrid(evaluate, divisions, block = null) {
  const axes = block?.knots
    ? [0, 1, 2].map(axis => parameterStations(block, axis, divisions))
    : [0, 1, 2].map(() => Array.from({ length: divisions + 1 }, (_, index) => index / divisions));
  const sizes = axes.map(axis => axis.length);
  const grid = new Array(sizes[0] * sizes[1] * sizes[2]);
  const index = (i, j, k) => (i * sizes[1] + j) * sizes[2] + k;
  for (let i = 0; i < sizes[0]; i++) {
    for (let j = 0; j < sizes[1]; j++) {
      for (let k = 0; k < sizes[2]; k++) {
        grid[index(i, j, k)] = evaluate(axes[0][i], axes[1][j], axes[2][k]);
      }
    }
  }
  return { grid, index, divisions: sizes.map(size => size - 1), axes };
}

// Polygonal sketch layers extrude convex quadrilaterals by translation.
// Their faces and plane sections are straight-sided, including sheared layers.
function meshDivisions(model, block, standard) {
  return model.schema === 'sketch_layer_partition_v1' && !block.knots
    && block.degrees.every(degree => degree === 1) ? 1 : standard;
}

function buildSurfaceMeshes(model, mode) {
  const interfaces = model.topology?.interfaces || inferChartInterfaces(model);
  const blockIds = new Set(model.blocks.map(block => block.id));
  const recordedExterior = model.boundary_faces || model.topology?.boundary_faces || model.topology?.exterior || null;
  const exteriorKeys = recordedExterior && new Set(recordedExterior.map(face => {
    const side = Array.isArray(face) ? face : [face?.source_block ?? face?.block, face?.side];
    if (!Array.isArray(side) || side.length < 2 || !blockIds.has(String(side[0]))
      || !['u0', 'u1', 'v0', 'v1', 'w0', 'w1'].includes(side[1])) {
      throw new TypeError('Recorded exterior contains an invalid block face');
    }
    return `${side[0]}:${side[1]}`;
  }));
  const hiddenSides = new Map();
  const samePatchSides = new Map();
  const duplicateSides = new Map();
  const patchForCell = new Map();
  for (const [patchIndex, patch] of (model.nurbs_patches || []).entries()) {
    for (const cell of patch.cells || []) patchForCell.set(cell, patch.id ?? String(patchIndex));
  }
  for (const item of interfaces) {
    const [duplicateBlock, duplicateSide] = item.second;
    const [firstBlock, firstSide] = item.first;
    for (const [id, side] of [item.first, item.second]) {
      if (!hiddenSides.has(id)) hiddenSides.set(id, new Set());
      hiddenSides.get(id).add(side);
    }
    if (mode === 'patches' && patchForCell.has(firstBlock) && patchForCell.get(firstBlock) === patchForCell.get(duplicateBlock)) {
      for (const [id, side] of [item.first, item.second]) {
        if (!samePatchSides.has(id)) samePatchSides.set(id, new Set());
        samePatchSides.get(id).add(side);
      }
    }
    if (!duplicateSides.has(duplicateBlock)) duplicateSides.set(duplicateBlock, new Set());
    duplicateSides.get(duplicateBlock).add(duplicateSide);
  }

  const meshes = [];
  let renderedBoundaryFaces = 0;
  const renderedFaceKeys = [];
  for (const block of model.blocks) {
    const evaluate = makeEvaluator(block);
    const divisions = meshDivisions(model, block, GRID);
    const axes = [0, 1, 2].map(axis => parameterStations(block, axis, divisions));
    const color = mode === 'surface' ? SURFACE_COLOR : partitionColor(block, mode);
    const skip = mode === 'surface'
      ? (exteriorKeys
        ? new Set(FACE_SIDES.map(([name]) => name).filter(name => !exteriorKeys.has(`${block.id}:${name}`)))
        : hiddenSides.get(block.id) || new Set())
      : mode === 'patches'
        ? new Set([
          ...(samePatchSides.get(block.id) || []),
          ...(duplicateSides.get(block.id) || []),
        ])
      : duplicateSides.get(block.id) || new Set();
    const output = [];

    for (const [name, axis, side] of FACE_SIDES) {
      if (skip.has(name)) continue;
      if (mode === 'surface') {
        renderedBoundaryFaces++;
        renderedFaceKeys.push(`${block.id}:${name}`);
      }
      const otherAxes = [0, 1, 2].filter(value => value !== axis);
      const first = axes[otherAxes[0]];
      const second = axes[otherAxes[1]];
      const faceGrid = first.map(a => second.map(b => {
        const parameters = [0, 0, 0];
        parameters[axis] = side;
        parameters[otherAxes[0]] = a;
        parameters[otherAxes[1]] = b;
        return evaluate(...parameters);
      }));
      for (let i = 0; i < first.length - 1; i++) {
        for (let j = 0; j < second.length - 1; j++) {
          let corners = [faceGrid[i][j], faceGrid[i + 1][j], faceGrid[i + 1][j + 1], faceGrid[i][j + 1]];

          const reverse = (axis === 0 && side === 0) || (axis === 1 && side === 1) || (axis === 2 && side === 0);
          if (reverse) corners = [corners[0], corners[3], corners[2], corners[1]];
          appendTriangle(output, corners[0], corners[1], corners[2], color);
          appendTriangle(output, corners[0], corners[2], corners[3], color);
        }
      }
    }
    if (output.length) meshes.push({ id: block.id, values: new Float32Array(output), count: output.length / 9 });
  }
  if (mode === 'patches' && model.nurbs_representation?.representation_kind === 'exact_disk_cylinder_tensor_nurbs_cover') {
    const representation = model.nurbs_representation;
    const seamRecords = [
      ...(representation.self_seams || []),
      ...(representation.interfaces || []),
    ].filter(item => item.status === 'exact_match');
    const limits = boundsOf(model);
    const center = limits[0].map((value, axis) => (value + limits[1][axis]) / 2);
    const extent = Math.max(...limits[0].map((value, axis) => limits[1][axis] - value)) || 1;
    const halfWidth = extent * 0.0016;
    const offset = extent * 0.00018;
    const seamColor = [0.10, 0.18, 0.15];
    const lineValues = [];
    const emitted = new Set();
    let boundaryCurveCount = 0;
    const appendBoundary = endpoint => {
      const blockId = endpoint.patch_id;
      const sideName = endpoint.patch_side;
      const block = model.blocks.find(item => item.id === blockId);
      if (!block || typeof sideName !== 'string') return;
      const evaluate = makeEvaluator(block);
      const axis = 'uvw'.indexOf(sideName[0]);
      if (axis < 0 || (sideName[1] !== '0' && sideName[1] !== '1')) return;
      const other = [0, 1, 2].filter(value => value !== axis);
      const intervals = endpoint.parameter_intervals || [[0, 1], [0, 1], [0, 1]];
      const domains = block.knots.map((knots, direction) => {
        const degree = block.degrees[direction];
        const count = block.shape[direction];
        return [knots[degree], knots[count]];
      });
      const normalizeParameter = (value, direction) => (value - domains[direction][0]) / (domains[direction][1] - domains[direction][0]);
      const curves = [];
      for (const freeAxis of other) for (const fixed of [0, 1]) {
        const variableAxis = other.find(value => value !== freeAxis);
        const [start, end] = intervals[freeAxis];
        const fixedValue = intervals[variableAxis][fixed];
        if (!domains[variableAxis].includes(fixedValue)) continue;
        const divisions = 32;
        const curve = [];
        for (let step = 0; step <= divisions; step++) {
          const parameters = [0, 0, 0];
          parameters[axis] = Number(sideName[1]);
          parameters[variableAxis] = normalizeParameter(fixedValue, variableAxis);
          parameters[freeAxis] = normalizeParameter(start + (end - start) * step / divisions, freeAxis);
          curve.push(evaluate(...parameters));
        }
        curves.push(curve);
      }
      for (const curve of curves) {
        const key = curve.map(point => point.map(value => Math.round(value * 1e8)).join(',')).join(';');
        const reverseKey = [...curve].reverse().map(point => point.map(value => Math.round(value * 1e8)).join(',')).join(';');
        if (emitted.has(key) || emitted.has(reverseKey)) continue;
        emitted.add(key);
        boundaryCurveCount++;
        for (let index = 0; index < curve.length - 1; index++) {
          const first = curve[index];
          const second = curve[index + 1];
          const tangent = normalize(subtract(second, first));
          const radial = normalize(subtract(first, center));
          let widthDirection = normalize(cross(tangent, radial));
          if (Math.hypot(...widthDirection) < 0.5) widthDirection = normalize(cross(tangent, [0, 0, 1]));
          const shift = point => point.map((value, axisIndex) => value + radial[axisIndex] * offset);
          const p0 = shift(first).map((value, axisIndex) => value + widthDirection[axisIndex] * halfWidth);
          const p1 = shift(first).map((value, axisIndex) => value - widthDirection[axisIndex] * halfWidth);
          const p2 = shift(second).map((value, axisIndex) => value - widthDirection[axisIndex] * halfWidth);
          const p3 = shift(second).map((value, axisIndex) => value + widthDirection[axisIndex] * halfWidth);
          appendTriangle(lineValues, p0, p1, p2, seamColor);
          appendTriangle(lineValues, p0, p2, p3, seamColor);
        }
      }
    };
    for (const seam of seamRecords) for (const endpoint of [seam.first, seam.second]) {
      appendBoundary(endpoint);
    }
    if (lineValues.length) meshes.push({ id: 'exact-patch-boundaries', values: new Float32Array(lineValues), count: lineValues.length / 9 });
    meshes.boundaryCurveCount = boundaryCurveCount;
  }
  meshes.boundaryStats = {
    renderedFaces: mode === 'surface' ? renderedBoundaryFaces : null,
    exteriorFaceCount: exteriorKeys?.size ?? null,
    excludedFaceCount: mode === 'surface' ? model.blocks.length * FACE_SIDES.length - renderedBoundaryFaces : null,
    source: exteriorKeys ? (model.boundary_faces || model.topology?.boundary_faces ? 'boundary_faces' : 'topology.exterior') : 'interfaces',
    faceKeys: renderedFaceKeys,
    patchBoundaryCurves: meshes.boundaryCurveCount ?? 0,
  };
  return meshes;
}

function facePoint(evaluate, sideName, a, b) {
  const axis = 'uvw'.indexOf(sideName[0]);
  const parameters = [0.5, 0.5, 0.5];
  parameters[axis] = Number(sideName[1]);
  const other = [0, 1, 2].filter(index => index !== axis);
  parameters[other[0]] = a;
  parameters[other[1]] = b;
  return evaluate(...parameters);
}

function inferChartInterfaces(model) {
  if (!model.blocks?.some(block => block.shape)) return [];
  const sample = [0, 0.5, 1];
  const transforms = [
    (a, b) => [a, b], (a, b) => [1 - a, b], (a, b) => [a, 1 - b], (a, b) => [1 - a, 1 - b],
    (a, b) => [b, a], (a, b) => [1 - b, a], (a, b) => [b, 1 - a], (a, b) => [1 - b, 1 - a],
  ];
  const sides = ['u0', 'u1', 'v0', 'v1', 'w0', 'w1'];
  const descriptors = model.blocks.map(block => ({ block, evaluate: makeEvaluator(block) }));
  const sampleIndex = value => sample.indexOf(value);
  const maps = transforms.map(transform => sample.flatMap(a => sample.map(b => {
    const [x, y] = transform(a, b);
    return sampleIndex(x) * sample.length + sampleIndex(y);
  })));
  for (const descriptor of descriptors) {
    descriptor.faceSamples = Object.fromEntries(sides.map(side => [
      side, sample.flatMap(a => sample.map(b => facePoint(descriptor.evaluate, side, a, b))),
    ]));
  }
  const interfaces = [];
  const claimed = new Set();
  const scale = Math.max(1, ...boundsOf(model).flat().map(Math.abs));
  const close = 2e-9 * scale;
  for (let firstIndex = 0; firstIndex < descriptors.length; firstIndex++) {
    for (let firstSideIndex = 0; firstSideIndex < sides.length; firstSideIndex++) {
      const first = descriptors[firstIndex];
      const firstSide = sides[firstSideIndex];
      const key = `${first.block.id}:${firstSide}`;
      if (claimed.has(key)) continue;
      const firstPoints = first.faceSamples[firstSide];
      let found = null;
      for (let secondIndex = firstIndex; secondIndex < descriptors.length && !found; secondIndex++) {
        for (let secondSideIndex = 0; secondSideIndex < sides.length && !found; secondSideIndex++) {
          if (secondIndex === firstIndex && secondSideIndex <= firstSideIndex) continue;
          const second = descriptors[secondIndex];
          const secondSide = sides[secondSideIndex];
          if (first.block.id === second.block.id) continue;
          const secondKey = `${second.block.id}:${secondSide}`;
          if (claimed.has(secondKey)) continue;
          for (const mapping of maps) {
            let matches = true;
            outer: for (let pointIndex = 0; pointIndex < firstPoints.length; pointIndex++) {
              const point = firstPoints[pointIndex];
              const other = second.faceSamples[secondSide][mapping[pointIndex]];
              if (Math.hypot(...subtract(point, other)) > close) { matches = false; break outer; }
            }
            if (matches) { found = { second, secondSide, secondKey }; break; }
          }
        }
      }
      if (found) {
        claimed.add(key);
        claimed.add(found.secondKey);
        interfaces.push({ first: [first.block.id, firstSide], second: [found.second.block.id, found.secondSide] });
      }
    }
  }
  return interfaces;
}

function intersectTetrahedron(points, axis, plane) {
  const intersections = [];
  for (const [first, second] of TETRA_EDGES) {
    const a = points[first];
    const b = points[second];
    const da = a[axis] - plane;
    const db = b[axis] - plane;
    if (!((da <= 0 && db >= 0) || (da >= 0 && db <= 0))) continue;
    const denominator = da - db;
    if (Math.abs(denominator) < 1e-14) continue;
    const t = da / denominator;
    if (t < 0 || t > 1) continue;
    const point = a.map((value, index) => value + t * (b[index] - value));
    if (!intersections.some(other => Math.hypot(...subtract(point, other)) < 1e-8)) intersections.push(point);
  }
  return intersections;
}

function coplanarBlockFace(evaluate, axis, plane, tolerance) {
  const samples = [0, 0.25, 0.5, 0.75, 1];
  for (let parameterAxis = 0; parameterAxis < 3; parameterAxis++) {
    for (let side = 0; side < 2; side++) {
      let onPlane = true;
      for (const first of samples) for (const second of samples) {
        const parameters = [0.5, 0.5, 0.5];
        parameters[parameterAxis] = side;
        const otherAxes = [0, 1, 2].filter(index => index !== parameterAxis);
        parameters[otherAxes[0]] = first;
        parameters[otherAxes[1]] = second;
        if (Math.abs(evaluate(...parameters)[axis] - plane) > tolerance) onPlane = false;
      }
      if (!onPlane) continue;

      const parameters = [0.5, 0.5, 0.5];
      parameters[parameterAxis] = side === 0 ? 1e-4 : 1 - 1e-4;
      const inwardValue = evaluate(...parameters)[axis];
      return {
        sideName: `${'uvw'[parameterAxis]}${side}`,
        retainedSide: inwardValue <= plane + tolerance,
      };
    }
  }
  return null;
}

function zSliceLevels(block) {
  const degreeW = block.degrees[2];
  if (degreeW !== 1) return null;
  const close = (a, b) => Math.abs(a - b) <= 1e-12 * Math.max(1, Math.abs(a), Math.abs(b));
  const dimensions = block.shape || block.degrees.map(value => value + 1);
  if (dimensions[2] !== 2) return null;
  const countV = dimensions[1];
  let z0 = null;
  let z1 = null;
  for (let i = 0; i < dimensions[0]; i++) {
    for (let j = 0; j < dimensions[1]; j++) {
      const first = block.controls[(i * countV + j) * 2];
      const last = block.controls[(i * countV + j) * 2 + 1];
      const firstZ = first[2] / first[3];
      const lastZ = last[2] / last[3];
      if (z0 === null) {
        z0 = firstZ;
        z1 = lastZ;
      } else if (!close(z0, firstZ) || !close(z1, lastZ)) {
        return null;
      }
    }
  }
  return z1 > z0 ? { z0, z1 } : null;
}

function buildCapMeshes(model, axis, fraction, bounds, mode) {
  if (axis === null) return [];
  const axisIndex = { x: 0, y: 1, z: 2 }[axis];
  const plane = bounds[0][axisIndex] + (bounds[1][axisIndex] - bounds[0][axisIndex]) * fraction;
  const tolerance = Math.max(1, bounds[1][axisIndex] - bounds[0][axisIndex]) * 1e-9;
  const zSlices = axis === 'z' ? model.blocks.map(zSliceLevels) : null;
  const interfaces = new Map();
  const duplicateSides = new Map();
  for (const item of model.topology?.interfaces || []) {
    for (const side of [item.first, item.second]) {
      if (!interfaces.has(side[0])) interfaces.set(side[0], new Set());
      interfaces.get(side[0]).add(side[1]);
    }
    const [duplicateBlock, duplicateSide] = item.second;
    if (!duplicateSides.has(duplicateBlock)) duplicateSides.set(duplicateBlock, new Set());
    duplicateSides.get(duplicateBlock).add(duplicateSide);
  }
  const meshes = [];
  for (let blockIndex = 0; blockIndex < model.blocks.length; blockIndex++) {
    const block = model.blocks[blockIndex];
    const evaluate = makeEvaluator(block);
    const output = [];
    const coplanar = coplanarBlockFace(evaluate, axisIndex, plane, tolerance);
    if (coplanar) {
      const interfaceFace = interfaces.get(block.id)?.has(coplanar.sideName) || false;
      const duplicateFace = duplicateSides.get(block.id)?.has(coplanar.sideName) || false;
      const faceVisible = mode !== 'surface' ? !duplicateFace : !interfaceFace;
      if (faceVisible || !coplanar.retainedSide) continue;
      // The retained side needs one cap where its shared interface is hidden.
    }
    const color = mode !== 'surface' ? partitionColor(block, mode) : SECTION_COLOR;
    if (zSlices?.[blockIndex]) {
      const { z0, z1 } = zSlices[blockIndex];
      const w = (plane - z0) / (z1 - z0);
      if (w < -1e-12 || w > 1 + 1e-12) continue;
      const divisions = meshDivisions(model, block, GRID);
      const axisU = parameterStations(block, 0, divisions);
      const axisV = parameterStations(block, 1, divisions);
      const evaluateHomogeneous = makeHomogeneousEvaluator(block);
      const points = axisU.map(u => axisV.map(v => {
        const first = evaluateHomogeneous(u, v, 0);
        const last = evaluateHomogeneous(u, v, 1);
        const denominator = (last[2] - first[2]) - plane * (last[3] - first[3]);
        if (Math.abs(denominator) < 1e-14) throw new TypeError(`Block ${block.id} has a degenerate z-section map`);
        const parameterW = clamp((plane * first[3] - first[2]) / denominator, 0, 1);
        return evaluate(u, v, parameterW);
      }));
      for (let i = 0; i < axisU.length - 1; i++) {
        for (let j = 0; j < axisV.length - 1; j++) {
          const p00 = points[i][j];
          const p10 = points[i + 1][j];
          const p11 = points[i + 1][j + 1];
          const p01 = points[i][j + 1];
          appendTriangle(output, p00, p10, p11, color);
          appendTriangle(output, p00, p11, p01, color);
        }
      }
      if (output.length) meshes.push({ id: `${block.id}:cap`, values: new Float32Array(output), count: output.length / 9 });
      continue;
    }

    const { grid, index, divisions } = buildBlockGrid(evaluate, meshDivisions(model, block, CAP_GRID), block);
    const point = (u, v, w) => grid[index(u, v, w)];
    for (let i = 0; i < divisions[0]; i++) {
      for (let j = 0; j < divisions[1]; j++) {
        for (let k = 0; k < divisions[2]; k++) {
          const base = [i, j, k];
          const corners = CUBE_CORNERS.map(offset => point(
            base[0] + offset[0], base[1] + offset[1], base[2] + offset[2],
          ));
          for (const tetrahedron of CUBE_TETRAHEDRA) {
            const tetraPoints = tetrahedron.map(corner => corners[corner]);
            const polygon = intersectTetrahedron(tetraPoints, axisIndex, plane);
            if (polygon.length < 3) continue;
            const u = (axisIndex + 1) % 3;
            const v = (axisIndex + 2) % 3;
            const center = [0, 0, 0];
            for (const vertex of polygon) for (let component = 0; component < 3; component++) center[component] += vertex[component] / polygon.length;
            polygon.sort((a, b) => Math.atan2(a[v] - center[v], a[u] - center[u]) - Math.atan2(b[v] - center[v], b[u] - center[u]));
            for (let triangle = 1; triangle < polygon.length - 1; triangle++) {
              appendTriangle(output, polygon[0], polygon[triangle], polygon[triangle + 1], color);
            }
          }
        }
      }
    }
    if (output.length) meshes.push({ id: `${block.id}:cap`, values: new Float32Array(output), count: output.length / 9 });
  }
  return meshes;
}

function createProgram(gl) {
  function compile(type, source) {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(shader));
    return shader;
  }
  const program = gl.createProgram();
  gl.attachShader(program, compile(gl.VERTEX_SHADER, VERTEX_SHADER));
  gl.attachShader(program, compile(gl.FRAGMENT_SHADER, FRAGMENT_SHADER));
  gl.linkProgram(program);
  if (!gl.getProgramParameter(program, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(program));
  return program;
}

/** Create an interactive viewer for an accepted trivariate model JSON value. */
export function createViewer(canvas, onCameraChange = () => {}) {
  const gl = canvas.getContext('webgl', { antialias: true, alpha: false });
  if (!gl) throw new Error('WebGL is unavailable');
  const program = createProgram(gl);
  const locations = {
    position: gl.getAttribLocation(program, 'aPosition'),
    normal: gl.getAttribLocation(program, 'aNormal'),
    color: gl.getAttribLocation(program, 'aColor'),
    mvp: gl.getUniformLocation(program, 'uMvp'),
    model: gl.getUniformLocation(program, 'uModel'),
    clip: gl.getUniformLocation(program, 'uClip'),
    useClip: gl.getUniformLocation(program, 'uUseClip'),
  };

  let model = null;
  let patchVisualizationModel = null;
  let bounds = null;
  let mode = 'surface';
  let clipAxis = null;
  let clipFraction = 0.5;
  let center = [0, 0, 0];
  let scale = 1;
  let camera = { ...DEFAULT_CAMERA };
  let cameraAdjusted = false;
  let lastCanvasSize = null;
  let drag = null;
  let dragTool = 'orbit';
  let frame = 0;
  let surfaceMeshes = [];
  let capMeshes = [];
  let triangleMesh = null;
  let surfaceStats = null;
  const activePointers = new Map();

  if (typeof ResizeObserver === 'function') {
    new ResizeObserver(scheduleDraw).observe(canvas);
  }

  function deleteMeshes(meshes) {
    for (const mesh of meshes) if (mesh.buffer) gl.deleteBuffer(mesh.buffer);
  }

  function uploadMeshes(meshes) {
    return meshes.map(mesh => {
      const buffer = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
      gl.bufferData(gl.ARRAY_BUFFER, mesh.values, gl.STATIC_DRAW);
      return { buffer, count: mesh.count };
    });
  }

  function scheduleDraw() {
    if (!frame) frame = requestAnimationFrame(draw);
  }

  function draw() {
    frame = 0;
    const pixelRatio = Math.min(globalThis.devicePixelRatio || 1, 2);
    const width = Math.max(1, Math.floor(canvas.clientWidth * pixelRatio));
    const height = Math.max(1, Math.floor(canvas.clientHeight * pixelRatio));
    if (canvas.width !== width || canvas.height !== height) {
      canvas.width = width;
      canvas.height = height;
    }
    if (!lastCanvasSize || lastCanvasSize[0] !== width || lastCanvasSize[1] !== height) {
      lastCanvasSize = [width, height];
      if (model && !cameraAdjusted) camera.distance = fittedDistance(width, height);
    }
    gl.viewport(0, 0, width, height);
    gl.enable(gl.DEPTH_TEST);
    gl.disable(gl.CULL_FACE);
    gl.clearColor(0.94, 0.92, 0.87, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    if (!model) return;

    gl.useProgram(program);
    const cosPitch = Math.cos(camera.pitch);
    const eye = [
      camera.distance * Math.sin(camera.yaw) * cosPitch,
      camera.distance * Math.cos(camera.yaw) * cosPitch,
      camera.distance * Math.sin(camera.pitch),
    ];
    const view = lookAt(eye, [0, 0, 0], [0, 0, 1]);
    view[12] += camera.panX;
    view[13] += camera.panY;
    const projection = perspective(Math.PI / 4, width / height, 0.01, 100);
    const mvp = multiplyMatrices(projection, view);
    const modelMatrix = new Float32Array([
      1 / scale, 0, 0, 0,
      0, 1 / scale, 0, 0,
      0, 0, 1 / scale, 0,
      -center[0] / scale, -center[1] / scale, -center[2] / scale, 1,
    ]);
    gl.uniformMatrix4fv(locations.mvp, false, mvp);
    gl.uniformMatrix4fv(locations.model, false, modelMatrix);
    const axisIndex = clipAxis === null ? 0 : { x: 0, y: 1, z: 2 }[clipAxis];
    const plane = clipAxis === null ? 0 : bounds[0][axisIndex] + (bounds[1][axisIndex] - bounds[0][axisIndex]) * clipFraction;
    const clipNormal = [0, 0, 0];
    if (clipAxis !== null) clipNormal[axisIndex] = 1;
    gl.uniform4f(locations.clip, clipNormal[0], clipNormal[1], clipNormal[2], plane);
    gl.uniform1f(locations.useClip, clipAxis === null ? 0 : 1);

    for (const mesh of [...surfaceMeshes, ...capMeshes]) {
      gl.bindBuffer(gl.ARRAY_BUFFER, mesh.buffer);
      const stride = 9 * Float32Array.BYTES_PER_ELEMENT;
      gl.enableVertexAttribArray(locations.position);
      gl.vertexAttribPointer(locations.position, 3, gl.FLOAT, false, stride, 0);
      gl.enableVertexAttribArray(locations.normal);
      gl.vertexAttribPointer(locations.normal, 3, gl.FLOAT, false, stride, 3 * Float32Array.BYTES_PER_ELEMENT);
      gl.enableVertexAttribArray(locations.color);
      gl.vertexAttribPointer(locations.color, 3, gl.FLOAT, false, stride, 6 * Float32Array.BYTES_PER_ELEMENT);
      gl.drawArrays(gl.TRIANGLES, 0, mesh.count);
    }
  }

  function rebuildSurfaces() {
    deleteMeshes(surfaceMeshes);
    surfaceStats = null;
    if (isTriangleModel(model)) {
      surfaceMeshes = uploadMeshes([{ values: triangleMesh.values, count: triangleMesh.count }]);
    } else {
      const visibleModel = mode === 'patches' && patchVisualizationModel ? patchVisualizationModel : model;
      const generated = buildSurfaceMeshes(visibleModel, mode);
      surfaceStats = generated.boundaryStats;
      surfaceMeshes = uploadMeshes(generated);
    }
    canvas.dataset.boundaryFaceCount = surfaceStats?.renderedFaces == null ? '' : String(surfaceStats.renderedFaces);
    canvas.dataset.exactExteriorFaceCount = surfaceStats?.exteriorFaceCount == null ? '' : String(surfaceStats.exteriorFaceCount);
    canvas.dataset.excludedInternalFaceCount = surfaceStats?.excludedFaceCount == null ? '' : String(surfaceStats.excludedFaceCount);
    canvas.dataset.boundarySource = surfaceStats?.source || (isTriangleModel(model) ? 'triangulated-surface' : '');
    canvas.dataset.boundaryFaceKeys = JSON.stringify(surfaceStats?.faceKeys || []);
    canvas.dataset.patchBoundaryCurves = String(surfaceStats?.patchBoundaryCurves ?? 0);
    canvas.dataset.patchRenderSource = patchVisualizationModel ? 'exact-patch-net' : '';
    const multispanAxes = patchVisualizationModel ? [0, 1, 2].filter(axis => patchVisualizationModel.blocks.some(block => {
      const degree = block.degrees[axis];
      const count = block.shape[axis];
      return new Set(block.knots[axis].filter(value => value >= block.knots[axis][degree]
        && value <= block.knots[axis][count])).size > 2;
    })) : [];
    canvas.dataset.patchMultispanAxes = JSON.stringify(multispanAxes);
  }

  function rebuildCaps() {
    deleteMeshes(capMeshes);
    const visibleModel = mode === 'patches' && patchVisualizationModel ? patchVisualizationModel : model;
    capMeshes = uploadMeshes(isTriangleModel(visibleModel) ? [] : buildCapMeshes(visibleModel, clipAxis, clipFraction, bounds, mode));
  }

  function fittedDistance(width, height) {
    if (!bounds) return DEFAULT_CAMERA.distance;
    const direction = [
      Math.sin(camera.yaw) * Math.cos(camera.pitch),
      Math.cos(camera.yaw) * Math.cos(camera.pitch),
      Math.sin(camera.pitch),
    ];
    const right = normalize([-Math.cos(camera.yaw), Math.sin(camera.yaw), 0]);
    const up = normalize(cross(direction, right));
    const aspect = Math.max(1e-6, width / height);
    const tangentY = Math.tan(Math.PI / 8);
    const tangentX = tangentY * aspect;
    let requiredDistance = 0;
    for (let corner = 0; corner < 8; corner++) {
      const relative = [0, 1, 2].map(axis => (
        ((corner >> axis) & 1 ? bounds[1][axis] : bounds[0][axis]) - center[axis]
      ) / scale);
      const depth = relative.reduce((sum, coordinate, axis) => sum + coordinate * direction[axis], 0);
      const horizontal = Math.abs(relative.reduce((sum, coordinate, axis) => sum + coordinate * right[axis], 0));
      const vertical = Math.abs(relative.reduce((sum, coordinate, axis) => sum + coordinate * up[axis], 0));
      requiredDistance = Math.max(
        requiredDistance,
        depth + horizontal / (tangentX * CAMERA_FILL),
        depth + vertical / (tangentY * CAMERA_FILL),
      );
    }
    return Math.max(0.85, requiredDistance * 1.04);
  }

  canvas.style.touchAction = 'none';
  function notifyCamera() {
    cameraAdjusted = true;
    scheduleDraw();
    onCameraChange({ ...camera });
  }

  function pan(dx, dy) {
    const units = 2 * camera.distance * Math.tan(Math.PI / 8) / Math.max(canvas.clientHeight, 1);
    camera.panX += dx * units;
    camera.panY -= dy * units;
  }

  canvas.addEventListener('contextmenu', event => event.preventDefault());
  canvas.addEventListener('pointerdown', event => {
    canvas.setPointerCapture(event.pointerId);
    activePointers.set(event.pointerId, { x: event.clientX, y: event.clientY });
    if (activePointers.size === 1) {
      drag = event.shiftKey || event.button === 1 || event.button === 2 ? 'pan' : dragTool;
    }
  });
  canvas.addEventListener('pointermove', event => {
    const previous = activePointers.get(event.pointerId);
    if (!previous) return;
    const dx = event.clientX - previous.x, dy = event.clientY - previous.y;
    const other = [...activePointers.entries()].find(([id]) => id !== event.pointerId)?.[1];
    activePointers.set(event.pointerId, { x: event.clientX, y: event.clientY });
    if (other) {
      const before = Math.hypot(previous.x - other.x, previous.y - other.y);
      const after = Math.hypot(event.clientX - other.x, event.clientY - other.y);
      if (before > 0 && after > 0) camera.distance = clamp(camera.distance * before / after, 0.85, 9);
      pan(dx / 2, dy / 2);
    } else if (drag === 'pan' || event.shiftKey) {
      pan(dx, dy);
    } else {
      camera.yaw += dx * 0.008;
      camera.pitch = clamp(camera.pitch + dy * 0.008, -1.45, 1.45);
    }
    notifyCamera();
  });
  function releasePointer(event) {
    activePointers.delete(event.pointerId);
    if (!activePointers.size) drag = null;
  }
  for (const event of ['pointerup','pointercancel','lostpointercapture']) canvas.addEventListener(event, releasePointer);
  canvas.addEventListener('wheel', event => {
    event.preventDefault();
    camera.distance = clamp(camera.distance * Math.exp(event.deltaY * 0.001), 0.85, 9);
    notifyCamera();
  }, { passive: false });
  canvas.addEventListener('keydown', event => {
    const orbitStep = event.shiftKey ? 0.14 : 0.06;
    if (event.key === 'Home') {
      camera = { ...DEFAULT_CAMERA };
      camera.distance = fittedDistance(canvas.clientWidth, canvas.clientHeight);
      cameraAdjusted = false;
    } else if (event.shiftKey && event.key.startsWith('Arrow')) {
      const directions = {ArrowLeft: [-20,0], ArrowRight: [20,0], ArrowUp: [0,-20], ArrowDown: [0,20]};
      pan(...directions[event.key]);
    } else if (event.key === 'ArrowLeft') camera.yaw -= orbitStep;
    else if (event.key === 'ArrowRight') camera.yaw += orbitStep;
    else if (event.key === 'ArrowUp') camera.pitch = clamp(camera.pitch - orbitStep, -1.45, 1.45);
    else if (event.key === 'ArrowDown') camera.pitch = clamp(camera.pitch + orbitStep, -1.45, 1.45);
    else if (event.key === '+' || event.key === '=' || event.code === 'NumpadAdd') { camera.distance = clamp(camera.distance * 0.9, 0.85, 9); cameraAdjusted = true; }
    else if (event.key === '-' || event.code === 'NumpadSubtract') { camera.distance = clamp(camera.distance / 0.9, 0.85, 9); cameraAdjusted = true; }
    else return;
    if (event.key !== 'Home') cameraAdjusted = true;
    event.preventDefault();
    scheduleDraw();
    onCameraChange({ ...camera });
  });

  return {
    setDragTool(tool) {
      if (!['orbit', 'pan'].includes(tool)) throw new TypeError('Unknown drag tool');
      dragTool = tool;
    },
    getCamera() { return { ...camera }; },
    getSurfaceStats() { return surfaceStats ? { ...surfaceStats } : null; },
    setCamera(next) {
      camera = { ...next };
      cameraAdjusted = true;
      scheduleDraw();
    },
    setModel(value, representation = null) {
      validateModel(value);
      model = isTriangleModel(value) ? value : normalizeTrivariateModel(value);
      patchVisualizationModel = null;
      if (!isTriangleModel(model) && !model.topology && model.blocks?.some(block => block.shape)) {
        model = { ...model, topology: { interfaces: inferChartInterfaces(model) } };
      }
      if (Array.isArray(model.nurbs_patches)) {
        const colors = new Map();
        model.nurbs_patches.forEach((patch, index) => {
          const color = partitionColour(index);
          patch.cells.forEach(id => colors.set(id, color));
        });
        model = {...model, blocks: model.blocks.map(block => ({...block, patchColor: colors.get(block.id)}))};
      }
      if (representation) {
        if (!Array.isArray(representation.patches) || !Array.isArray(representation.boundary_faces)) {
          throw new TypeError('Exact NURBS representation needs patches and boundary_faces arrays');
        }
        const blockIds = new Set(model.blocks.map(block => block.id));
        const patchIds = new Set();
        const assignedBlocks = new Set();
        const patchForBlock = new Map();
        const patchRecords = representation.patches.map(patch => {
          if (!patch || typeof patch.id !== 'string' || patchIds.has(patch.id)) {
            throw new TypeError('Exact NURBS representation has a missing or repeated patch id');
          }
          patchIds.add(patch.id);
          const cells = patch.source_blocks || patch.cells;
          if (!Array.isArray(cells) || cells.length === 0 || cells.some(id => !blockIds.has(id))) {
            throw new TypeError(`NURBS patch ${patch.id} has an invalid source cell list`);
          }
          for (const id of cells) {
            if (assignedBlocks.has(id)) throw new TypeError(`Source block ${id} appears in more than one NURBS patch`);
            assignedBlocks.add(id);
            patchForBlock.set(id, patch.id);
          }
          return { id: patch.id, cells };
        });
        if (assignedBlocks.size !== blockIds.size) throw new TypeError('NURBS patches do not cover every native model block');
        for (const face of representation.boundary_faces) {
          const sourceBlock = face?.source_block ?? face?.block;
          if (!face || !blockIds.has(sourceBlock) || !patchIds.has(face.patch_id)
            || !['u0', 'u1', 'v0', 'v1', 'w0', 'w1'].includes(face.side)
            || face.status !== 'extracted_exactly_from_model') {
            throw new TypeError('Exact NURBS representation has an invalid exterior face reference');
          }
          if (patchForBlock.get(sourceBlock) !== face.patch_id) {
            throw new TypeError('Exterior face is assigned to the wrong NURBS patch');
          }
        }
        const nativeInterfaces = model.topology?.interfaces || [];
        const sourceFaceKey = (block, side) => `${block}:${side}`;
        const nativeExterior = new Set((model.topology?.exterior || []).map(face => sourceFaceKey(face[0], face[1])));
        const supplementExterior = new Set(representation.boundary_faces.map(face => sourceFaceKey(face.source_block ?? face.block, face.side)));
        if (nativeExterior.size && (nativeExterior.size !== supplementExterior.size
          || [...nativeExterior].some(key => !supplementExterior.has(key)))) {
          throw new TypeError('Exact NURBS boundary does not match the native model exterior faces');
        }
        const interfaceSides = new Set(nativeInterfaces.flatMap(item => [
          sourceFaceKey(item.first[0], item.first[1]), sourceFaceKey(item.second[0], item.second[1]),
        ]));
        if ([...supplementExterior].some(key => interfaceSides.has(key))
          || supplementExterior.size !== blockIds.size * 6 - nativeInterfaces.length * 2) {
          throw new TypeError('Exact NURBS boundary does not close the native cell interface graph');
        }
        const endpoint = value => Array.isArray(value)
          ? { patch_id: value[0], source_block: value[0], side: value[1] } : value;
        const sourcePairKey = (first, second) => [sourceFaceKey(first.source_block, first.side), sourceFaceKey(second.source_block, second.side)].sort().join('|');
        const sourceInterfacePairs = new Set(nativeInterfaces.map(item => sourcePairKey(
          { source_block: item.first[0], side: item.first[1] },
          { source_block: item.second[0], side: item.second[1] },
        )));
        const expectedPatchInterfaces = new Set();
        for (const item of nativeInterfaces) {
          if (patchForBlock.get(item.first[0]) !== patchForBlock.get(item.second[0])) {
            expectedPatchInterfaces.add(sourcePairKey(
              { source_block: item.first[0], side: item.first[1] },
              { source_block: item.second[0], side: item.second[1] },
            ));
          }
        }
        const observedPatchInterfaces = new Set();
        for (const item of representation.interfaces || []) {
          if (item.status !== 'exact_match') {
            throw new TypeError('NURBS representation contains an unjoined patch interface');
          }
          for (const value of [item.first, item.second]) {
            const side = endpoint(value);
            if (!side || !patchIds.has(side.patch_id) || !blockIds.has(side.source_block)
              || !['u0', 'u1', 'v0', 'v1', 'w0', 'w1'].includes(side.side)
              || patchForBlock.get(side.source_block) !== side.patch_id) {
              throw new TypeError('NURBS representation contains an invalid patch interface');
            }
          }
          const first = endpoint(item.first);
          const second = endpoint(item.second);
          const pair = sourcePairKey(first, second);
          if (!sourceInterfacePairs.has(pair) || patchForBlock.get(first.source_block) === patchForBlock.get(second.source_block)) {
            throw new TypeError('NURBS patch interface does not match a cross-patch native cell interface');
          }
          observedPatchInterfaces.add(pair);
        }
        if (observedPatchInterfaces.size !== (representation.interfaces || []).length
          || observedPatchInterfaces.size !== expectedPatchInterfaces.size
          || [...expectedPatchInterfaces].some(key => !observedPatchInterfaces.has(key))) {
          throw new TypeError('NURBS patch interfaces do not cover all cross-patch cell joins');
        }
        const boundaryFaces = representation.boundary_faces.map(face => ({
          ...face, block: face.source_block ?? face.block,
        }));
        model = { ...model, boundary_faces: boundaryFaces, nurbs_patches: patchRecords, nurbs_representation: representation };
        const colors = new Map();
        patchRecords.forEach((patch, index) => patch.cells.forEach(id => colors.set(id, partitionColour(index))));
        model = { ...model, blocks: model.blocks.map(block => ({ ...block, patchColor: colors.get(block.id) })) };
        if (representation.representation_kind === 'exact_disk_cylinder_tensor_nurbs_cover') {
          const patchBlocks = representation.patches.map((patch, index) => ({
            id: patch.id,
            degrees: patch.degrees,
            shape: patch.control_shape,
            knots: patch.knots,
            controls: patch.controls_homogeneous,
            displayColor: partitionColour(index),
          }));
          const patchEndpoint = endpoint => Array.isArray(endpoint)
            ? endpoint : [endpoint.patch_id, endpoint.patch_side];
          const patchInterfaces = [...(representation.interfaces || []), ...(representation.self_seams || [])]
            .filter(item => item.status === 'exact_match')
            .map(item => ({ first: patchEndpoint(item.first), second: patchEndpoint(item.second) }));
          patchVisualizationModel = {
            blocks: patchBlocks,
            topology: { interfaces: patchInterfaces },
            nurbs_patches: patchBlocks.map(block => ({ id: block.id, cells: [block.id] })),
            nurbs_representation: representation,
          };
          validateModel(patchVisualizationModel);
        }
      }
      triangleMesh = isTriangleModel(model) ? surfaceMesh(model) : null;
      bounds = triangleMesh ? triangleMesh.bounds : boundsOf(model);
      center = bounds[0].map((coordinate, axis) => (coordinate + bounds[1][axis]) / 2);
      scale = Math.max(...bounds[0].map((coordinate, axis) => bounds[1][axis] - coordinate)) || 1;
      rebuildSurfaces();
      rebuildCaps();
      camera = { ...DEFAULT_CAMERA };
      camera.distance = fittedDistance(canvas.clientWidth, canvas.clientHeight);
      cameraAdjusted = false;
      scheduleDraw();
    },
    setMode(value) {
      if (!['surface', 'blocks', 'patches'].includes(value)) throw new TypeError('Unknown display mode');
      if (isTriangleModel(model) && value !== 'surface') throw new TypeError('Triangle surfaces have no logical block color mode');
      if (mode === value) return;
      mode = value;
      if (model) {
        rebuildSurfaces();
        rebuildCaps();
      }
      scheduleDraw();
    },
    setClip(axis, fraction) {
      if (axis !== null && !['x', 'y', 'z'].includes(axis)) throw new TypeError('Clip axis must be x, y, z, or null');
      if (!Number.isFinite(fraction) || fraction < 0 || fraction > 1) throw new RangeError('Clip fraction must be in [0,1]');
      if (clipAxis === axis && clipFraction === fraction) return;
      clipAxis = axis;
      clipFraction = fraction;
      if (model) rebuildCaps();
      scheduleDraw();
    },
    reset() {
      camera = { ...DEFAULT_CAMERA };
      camera.distance = fittedDistance(canvas.clientWidth, canvas.clientHeight);
      cameraAdjusted = false;
      scheduleDraw();
      onCameraChange({ ...camera });
    },
  };
}
