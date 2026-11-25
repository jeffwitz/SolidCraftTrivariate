// Prepare the tessellation of an independent STEP boundary for the gallery.
export function surfaceMesh({ vertices, triangles }) {
  if (!Array.isArray(vertices) || !vertices.length || vertices.length % 3
    || !vertices.every(Number.isFinite) || !Array.isArray(triangles)
    || !triangles.length || triangles.length % 3) {
    throw new TypeError('The surface triangulation is invalid.');
  }
  const values = new Float32Array(triangles.length * 9);
  const bounds = [[Infinity, Infinity, Infinity], [-Infinity, -Infinity, -Infinity]];
  for (let i = 0; i < vertices.length; i++) {
    const axis = i % 3;
    bounds[0][axis] = Math.min(bounds[0][axis], vertices[i]);
    bounds[1][axis] = Math.max(bounds[1][axis], vertices[i]);
  }
  for (let i = 0; i < triangles.length; i += 3) {
    const [a, b, c] = triangles.slice(i, i + 3).map(index => {
      if (!Number.isSafeInteger(index) || index < 0 || 3 * index >= vertices.length) {
        throw new RangeError('Triangle index is outside the model.');
      }
      return vertices.slice(index * 3, index * 3 + 3);
    });
    const u = b.map((x, j) => x - a[j]);
    const v = c.map((x, j) => x - a[j]);
    const normal = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]];
    const length = Math.hypot(...normal) || 1;
    for (const [j, point] of [a,b,c].entries()) {
      values.set([...point, ...normal.map(x => x / length), 0.77, 0.67, 0.47], (i+j)*9);
    }
  }
  return { values, count: triangles.length, bounds };
}
