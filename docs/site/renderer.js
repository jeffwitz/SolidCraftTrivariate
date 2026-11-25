/* Native WebGL renderer for rational trivariate model JSON. */

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
const DEFAULT_CAMERA = { yaw: 0.75, pitch: 0.52, distance: 2.7 };
const GRID = 12;
const CAP_GRID = 6;
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
    const expected = block.degrees.reduce((product, degree) => product * (degree + 1), 1);
    if (!Array.isArray(block.controls) || block.controls.length !== expected) {
      throw new TypeError(`Block ${block.id} has an invalid control count`);
    }
    for (const control of block.controls) {
      if (!Array.isArray(control) || control.length !== 4
        || control.some(value => !Number.isFinite(value)) || control[3] <= 0) {
        throw new TypeError(`Block ${block.id} has invalid homogeneous controls`);
      }
    }
  }
  const blockIds = new Set(model.blocks.map(block => block.id));
  for (const item of model.topology?.interfaces || []) {
    for (const side of [item.first, item.second]) {
      if (!Array.isArray(side) || side.length !== 2 || !blockIds.has(side[0])
        || !['u0', 'u1', 'v0', 'v1', 'w0', 'w1'].includes(side[1])) {
        throw new TypeError('Model contains an invalid block interface');
      }
    }
  }
}

function boundsOf(model) {
  const low = [Infinity, Infinity, Infinity];
  const high = [-Infinity, -Infinity, -Infinity];
  for (const block of model.blocks) {
    for (const control of block.controls) {
      for (let axis = 0; axis < 3; axis++) {
        const coordinate = control[axis] / control[3];
        low[axis] = Math.min(low[axis], coordinate);
        high[axis] = Math.max(high[axis], coordinate);
      }
    }
  }
  return [low, high];
}

function makeEvaluator(block) {
  const degree = block.degrees;
  const dimensions = degree.map(value => value + 1);
  const controls = block.controls;
  return (u, v, w) => {
    const bu = basis(degree[0], u);
    const bv = basis(degree[1], v);
    const bw = basis(degree[2], w);
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

function buildBlockGrid(evaluate, divisions) {
  const side = divisions + 1;
  const grid = new Array(side * side * side);
  const index = (i, j, k) => (i * side + j) * side + k;
  for (let i = 0; i <= divisions; i++) {
    for (let j = 0; j <= divisions; j++) {
      for (let k = 0; k <= divisions; k++) {
        grid[index(i, j, k)] = evaluate(i / divisions, j / divisions, k / divisions);
      }
    }
  }
  return { grid, index, divisions };
}

function buildSurfaceMeshes(model, mode) {
  const hiddenSides = new Map();
  for (const item of model.topology?.interfaces || []) {
    for (const [id, side] of [item.first, item.second]) {
      if (!hiddenSides.has(id)) hiddenSides.set(id, new Set());
      hiddenSides.get(id).add(side);
    }
  }

  const meshes = [];
  for (const block of model.blocks) {
    const evaluate = makeEvaluator(block);
    const { grid, index, divisions } = buildBlockGrid(evaluate, GRID);
    const color = mode === 'surface' ? SURFACE_COLOR : blockColor(block.id);
    const skip = mode === 'surface' ? hiddenSides.get(block.id) || new Set() : new Set();
    const output = [];
    const point = (u, v, w) => grid[index(u, v, w)];

    for (const [name, axis, side] of FACE_SIDES) {
      if (skip.has(name)) continue;
      for (let i = 0; i < divisions; i++) {
        for (let j = 0; j < divisions; j++) {
          let corners;
          if (axis === 0) corners = [
            point(side * divisions, i, j), point(side * divisions, i + 1, j),
            point(side * divisions, i + 1, j + 1), point(side * divisions, i, j + 1),
          ];
          else if (axis === 1) corners = [
            point(i, side * divisions, j), point(i + 1, side * divisions, j),
            point(i + 1, side * divisions, j + 1), point(i, side * divisions, j + 1),
          ];
          else corners = [
            point(i, j, side * divisions), point(i + 1, j, side * divisions),
            point(i + 1, j + 1, side * divisions), point(i, j + 1, side * divisions),
          ];

          const reverse = (axis === 0 && side === 0) || (axis === 1 && side === 1) || (axis === 2 && side === 0);
          if (reverse) corners = [corners[0], corners[3], corners[2], corners[1]];
          appendTriangle(output, corners[0], corners[1], corners[2], color);
          appendTriangle(output, corners[0], corners[2], corners[3], color);
        }
      }
    }
    if (output.length) meshes.push({ id: block.id, values: new Float32Array(output), count: output.length / 9 });
  }
  return meshes;
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

function buildCapMeshes(model, axis, fraction, bounds) {
  if (axis === null) return [];
  const axisIndex = { x: 0, y: 1, z: 2 }[axis];
  const plane = bounds[0][axisIndex] + (bounds[1][axisIndex] - bounds[0][axisIndex]) * fraction;
  const meshes = [];
  for (const block of model.blocks) {
    const evaluate = makeEvaluator(block);
    const { grid, index, divisions } = buildBlockGrid(evaluate, CAP_GRID);
    const output = [];
    const point = (u, v, w) => grid[index(u, v, w)];
    for (let i = 0; i < divisions; i++) {
      for (let j = 0; j < divisions; j++) {
        for (let k = 0; k < divisions; k++) {
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
              appendTriangle(output, polygon[0], polygon[triangle], polygon[triangle + 1], SURFACE_COLOR);
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
export function createViewer(canvas) {
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
  let bounds = null;
  let mode = 'surface';
  let clipAxis = null;
  let clipFraction = 0.5;
  let center = [0, 0, 0];
  let scale = 1;
  let camera = { ...DEFAULT_CAMERA };
  let drag = null;
  let frame = 0;
  let surfaceMeshes = [];
  let capMeshes = [];
  const activePointers = new Map();
  let pinchDistance = null;

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
      camera.distance * Math.sin(camera.pitch),
      camera.distance * Math.cos(camera.yaw) * cosPitch,
    ];
    const view = lookAt(eye, [0, 0, 0], [0, 1, 0]);
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
    surfaceMeshes = uploadMeshes(buildSurfaceMeshes(model, mode));
  }

  function rebuildCaps() {
    deleteMeshes(capMeshes);
    capMeshes = uploadMeshes(buildCapMeshes(model, clipAxis, clipFraction, bounds));
  }

  canvas.style.touchAction = 'none';
  canvas.addEventListener('pointerdown', event => {
    activePointers.set(event.pointerId, { x: event.clientX, y: event.clientY });
    canvas.setPointerCapture(event.pointerId);
    if (activePointers.size === 1) {
      drag = { x: event.clientX, y: event.clientY, yaw: camera.yaw, pitch: camera.pitch };
    } else if (activePointers.size === 2) {
      const [first, second] = [...activePointers.values()];
      pinchDistance = Math.hypot(first.x - second.x, first.y - second.y);
      drag = null;
    }
  });
  canvas.addEventListener('pointermove', event => {
    const previous = activePointers.get(event.pointerId);
    if (!previous) return;
    activePointers.set(event.pointerId, { x: event.clientX, y: event.clientY });
    if (activePointers.size >= 2) {
      const [first, second] = [...activePointers.values()];
      const nextDistance = Math.hypot(first.x - second.x, first.y - second.y);
      if (pinchDistance > 0 && nextDistance > 0) camera.distance = clamp(camera.distance * pinchDistance / nextDistance, 1.1, 9);
      pinchDistance = nextDistance;
      scheduleDraw();
      return;
    }
    if (!drag) return;
    camera.yaw = drag.yaw + (event.clientX - drag.x) * 0.008;
    camera.pitch = clamp(drag.pitch + (event.clientY - drag.y) * 0.008, -1.45, 1.45);
    scheduleDraw();
  });
  function releasePointer(event) {
    activePointers.delete(event.pointerId);
    pinchDistance = null;
    if (activePointers.size === 1) {
      const [remaining] = activePointers.values();
      drag = { x: remaining.x, y: remaining.y, yaw: camera.yaw, pitch: camera.pitch };
    } else {
      drag = null;
    }
  }
  canvas.addEventListener('pointerup', releasePointer);
  canvas.addEventListener('pointercancel', releasePointer);
  canvas.addEventListener('wheel', event => {
    event.preventDefault();
    camera.distance = clamp(camera.distance * Math.exp(event.deltaY * 0.001), 1.1, 9);
    scheduleDraw();
  }, { passive: false });
  canvas.addEventListener('keydown', event => {
    const orbitStep = event.shiftKey ? 0.14 : 0.06;
    if (event.key === 'Home') {
      camera = { ...DEFAULT_CAMERA };
    } else if (event.key === 'ArrowLeft') camera.yaw -= orbitStep;
    else if (event.key === 'ArrowRight') camera.yaw += orbitStep;
    else if (event.key === 'ArrowUp') camera.pitch = clamp(camera.pitch - orbitStep, -1.45, 1.45);
    else if (event.key === 'ArrowDown') camera.pitch = clamp(camera.pitch + orbitStep, -1.45, 1.45);
    else if (event.key === '+' || event.key === '=' || event.code === 'NumpadAdd') camera.distance = clamp(camera.distance * 0.9, 1.1, 9);
    else if (event.key === '-' || event.code === 'NumpadSubtract') camera.distance = clamp(camera.distance / 0.9, 1.1, 9);
    else return;
    event.preventDefault();
    scheduleDraw();
  });

  return {
    setModel(value) {
      validateModel(value);
      model = value;
      bounds = boundsOf(model);
      center = bounds[0].map((coordinate, axis) => (coordinate + bounds[1][axis]) / 2);
      scale = Math.max(...bounds[0].map((coordinate, axis) => bounds[1][axis] - coordinate)) || 1;
      rebuildSurfaces();
      rebuildCaps();
      camera = { ...DEFAULT_CAMERA };
      scheduleDraw();
    },
    setMode(value) {
      if (value !== 'surface' && value !== 'blocks') throw new TypeError("Mode must be 'surface' or 'blocks'");
      if (mode === value) return;
      mode = value;
      if (model) rebuildSurfaces();
      scheduleDraw();
    },
    setClip(axis, fraction) {
      if (axis !== null && !['x', 'y', 'z'].includes(axis)) throw new TypeError('Clip axis must be x, y, z, or null');
      if (!Number.isFinite(fraction) || fraction < 0 || fraction > 1) throw new RangeError('Clip fraction must be in [0,1]');
      clipAxis = axis;
      clipFraction = fraction;
      if (model) rebuildCaps();
      scheduleDraw();
    },
    reset() {
      camera = { ...DEFAULT_CAMERA };
      scheduleDraw();
    },
  };
}
