import { Program } from './gl.mjs';
import * as gm from './gm.mjs';

// Global variables for mesh data
window._meshData = null;
window._meshProgram = null;

// Parse STL binary data into triangle mesh
export function parseSTL(buffer) {
  const view = new DataView(buffer);
  const triangles = [];
  const littleEndian = true;
  const numTriangles = view.getUint32(80, littleEndian);
  let offset = 84;
  for (let i = 0; i < numTriangles; i++) {
    // Skip normal (12 bytes)
    offset += 12;
    const a = [view.getFloat32(offset, littleEndian), view.getFloat32(offset + 4, littleEndian), view.getFloat32(offset + 8, littleEndian)];
    offset += 12;
    const b = [view.getFloat32(offset, littleEndian), view.getFloat32(offset + 4, littleEndian), view.getFloat32(offset + 8, littleEndian)];
    offset += 12;
    const c = [view.getFloat32(offset, littleEndian), view.getFloat32(offset + 4, littleEndian), view.getFloat32(offset + 8, littleEndian)];
    offset += 12;
    // Skip attribute byte count (2 bytes)
    offset += 2;
    triangles.push({ a3: a, b3: b, c3: c });
  }
  return triangles;
}

// Setup GL mesh from STL buffer
export function setupGLMeshFromSTL(arrayBuffer) {
  const triangles = parseSTL(arrayBuffer);

  let positions = [];
  let normals = [];

  for (const tri of triangles) {
    const { a3, b3, c3 } = tri;
    const u = gm.sub(b3, a3);
    const v = gm.sub(c3, a3);
    const n = gm.normalize(gm.cross(u, v));

    positions.push(...a3, ...b3, ...c3);
    normals.push(...n, ...n, ...n);
  }

  const gl = window.gl;

  const meshData = {
    positions: new Float32Array(positions),
    normals: new Float32Array(normals),
    count: positions.length / 3,
  };

  meshData.positionBuffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, meshData.positionBuffer);
  gl.bufferData(gl.ARRAY_BUFFER, meshData.positions, gl.STATIC_DRAW);

  meshData.normalBuffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, meshData.normalBuffer);
  gl.bufferData(gl.ARRAY_BUFFER, meshData.normals, gl.STATIC_DRAW);

  window._meshData = meshData;
  
  // Create the shader program if it doesn't exist yet
  if (!window._meshProgram) {
    const vsSource = `
      attribute vec3 aPosition;
      attribute vec3 aNormal;
      uniform mat4 uProjection;
      uniform mat4 uView;
      varying vec3 vNormal;
      void main() {
        gl_Position = uProjection * uView * vec4(aPosition, 1.0);
        vNormal = aNormal;
      }
    `;
    
    const fsSource = `
      precision mediump float;
      varying vec3 vNormal;
      void main() {
        vec3 light = normalize(vec3(1.0, 1.0, 1.0));
        float diffuse = max(dot(normalize(vNormal), light), 0.5);
        gl_FragColor = vec4(vec3(1.0, 1.0, 1.0) * diffuse, 0.5);  // vec4(vec3(R, G, B), Aplha)
      }
    `;
    
    // Create the program using the existing Program class
    window._meshProgram = new Program(gl, vsSource, fsSource);
  }
  
  // Request a redraw
  window.requestRedraw();
}

// Draw the mesh
export function drawGLMesh() {
  if (!window._meshData || !window._meshProgram) return;

  const gl = window.gl;
  const meshData = window._meshData;
  const program = window._meshProgram;
  
  gl.useProgram(program.program);

  // Setup attribute pointers
  if (program.attribLocations.aPosition !== undefined) {
    gl.bindBuffer(gl.ARRAY_BUFFER, meshData.positionBuffer);
    gl.enableVertexAttribArray(program.attribLocations.aPosition);
    gl.vertexAttribPointer(program.attribLocations.aPosition, 3, gl.FLOAT, false, 0, 0);
  }

  if (program.attribLocations.aNormal !== undefined) {
    gl.bindBuffer(gl.ARRAY_BUFFER, meshData.normalBuffer);
    gl.enableVertexAttribArray(program.attribLocations.aNormal);
    gl.vertexAttribPointer(program.attribLocations.aNormal, 3, gl.FLOAT, false, 0, 0);
  }

  // Set uniforms
  gl.uniformMatrix4fv(program.uniformLocations.uProjection, false, window.cameraProjectionMatrix());
  gl.uniformMatrix4fv(program.uniformLocations.uView, false, window.cameraViewMatrix());

  // Enable blending for transparency
  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  
  // Draw the triangles
  gl.drawArrays(gl.TRIANGLES, 0, meshData.count);
  
  // Restore state
  gl.disable(gl.BLEND);
}

export function cleanupGLMesh() {
  if (window._meshData) {
    const gl = window.gl;
    
    // delete buffers
    if (window._meshData.positionBuffer) {
      gl.deleteBuffer(window._meshData.positionBuffer);
    }
    
    if (window._meshData.normalBuffer) {
      gl.deleteBuffer(window._meshData.normalBuffer);
    }
    
    if (window._meshData.colorBuffer) {
      gl.deleteBuffer(window._meshData.colorBuffer);
    }
    
    // initialize meshData to null
    window._meshData = null;
  }
}