import { createProgram } from './webgl/shader';
import vertexSrc from './shaders/fullscreen.vert.glsl?raw';
import fragmentSrc from './shaders/blackhole.frag.glsl?raw';
import gridVS from './shaders/grid.vert.glsl?raw';
import gridFS from './shaders/grid.frag.glsl?raw';
import { OrbitCamera } from './orbitCamera';

const canvas = document.getElementById('glcanvas') as HTMLCanvasElement;
const errDiv = document.getElementById('err') as HTMLDivElement | null;
const glCtx = canvas.getContext('webgl2', { antialias: false, preserveDrawingBuffer: false });
if (!glCtx) {
  if (errDiv) errDiv.textContent = 'WebGL2 not supported or blocked. Try a different browser/device.';
  throw new Error('WebGL2 not supported');
}
const gl = glCtx as WebGL2RenderingContext;
const maxStepsSlider = document.getElementById('maxSteps') as HTMLInputElement;
const maxStepsVal = document.getElementById('maxStepsVal') as HTMLSpanElement;
let maxSteps = 20000;
const obj1Toggle = document.getElementById('obj1Toggle') as HTMLInputElement;
const obj2Toggle = document.getElementById('obj2Toggle') as HTMLInputElement;
let obj1Enabled = true;
let obj2Enabled = true;
if (maxStepsSlider && maxStepsVal) {
  const parseSteps = () => {
    maxSteps = parseInt(maxStepsSlider.value, 10);
    // show compact thousand suffix
    const k = Math.round(maxSteps / 1000);
    maxStepsVal.textContent = `${k}k`;
  };
  maxStepsSlider.addEventListener('input', parseSteps);
  parseSteps();
}
if (obj1Toggle) {
  obj1Toggle.addEventListener('change', () => { obj1Enabled = obj1Toggle.checked; });
  obj1Enabled = obj1Toggle.checked;
}
if (obj2Toggle) {
  obj2Toggle.addEventListener('change', () => { obj2Enabled = obj2Toggle.checked; });
  obj2Enabled = obj2Toggle.checked;
}

// Low internal resolution for pixelated look; CSS scales to window
let internalWidth = 200;
let internalHeight = 150;

function resizeCanvas() {
  // Keep internal low-res; only adjust when window is too small
  const w = Math.max(100, Math.floor(window.innerWidth / 8));
  const h = Math.max(75, Math.floor(window.innerHeight / 8));
  internalWidth = Math.min(320, w);
  internalHeight = Math.min(240, h);
  canvas.width = internalWidth;
  canvas.height = internalHeight;
  gl.viewport(0, 0, internalWidth, internalHeight);
}
resizeCanvas();
window.addEventListener('resize', resizeCanvas);

const program = createProgram(gl, vertexSrc, fragmentSrc);
gl.useProgram(program);

// Grid pipeline
const gridProgram = createProgram(gl, gridVS, gridFS);
const gridVAO = gl.createVertexArray();
const gridVBO = gl.createBuffer();
const gridEBO = gl.createBuffer();
let gridIndexCount = 0;

function generateGrid() {
  const gridSize = 25;
  const spacing = 1.0; // world units (rs = 1)
  const gridYOffset = -10.0; // push grid below the black hole (z-up, negative is down)
  const vertices: number[] = [];
  const indices: number[] = [];

  for (let z = 0; z <= gridSize; z++) {
    for (let x = 0; x <= gridSize; x++) {
      const worldX = (x - gridSize / 2) * spacing;
      const worldY = (z - gridSize / 2) * spacing;

      // Warp y using a Schwarzschild-inspired embedding surface
      const rs = 1.0;
      const r = Math.hypot(worldX, worldY);
      let y = 0;
      if (r > rs) {
        y = gridYOffset + 2.0 * Math.sqrt(Math.max(rs * (r - rs), 0.0));
      } else {
        // keep center continuous and cupped upward toward the BH
        y = gridYOffset + 2.0 * Math.sqrt(rs * rs);
      }
      // z-up: height goes into Z; plane is X-Y
      vertices.push(worldX, worldY, y);
    }
  }

  for (let z = 0; z < gridSize; z++) {
    for (let x = 0; x < gridSize; x++) {
      const i = z * (gridSize + 1) + x;
      indices.push(i, i + 1);
      indices.push(i, i + gridSize + 1);
    }
  }

  gl.bindVertexArray(gridVAO);
  gl.bindBuffer(gl.ARRAY_BUFFER, gridVBO);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.DYNAMIC_DRAW);
  gl.enableVertexAttribArray(0);
  gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 12, 0);

  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, gridEBO);
  gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(indices), gl.STATIC_DRAW);
  gridIndexCount = indices.length;
}

generateGrid();

// Fullscreen quad
const vao = gl.createVertexArray();
gl.bindVertexArray(vao);
const quad = new Float32Array([
  -1,  1, 0, 1,
  -1, -1, 0, 0,
   1, -1, 1, 0,
  -1,  1, 0, 1,
   1, -1, 1, 0,
   1,  1, 1, 1,
]);
const vbo = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, vbo);
gl.bufferData(gl.ARRAY_BUFFER, quad, gl.STATIC_DRAW);
gl.enableVertexAttribArray(0);
gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 16, 0);
gl.enableVertexAttribArray(1);
gl.vertexAttribPointer(1, 2, gl.FLOAT, false, 16, 8);

// Uniform locations
const loc = {
  camPos: gl.getUniformLocation(program, 'u_camPos'),
  camRight: gl.getUniformLocation(program, 'u_camRight'),
  camUp: gl.getUniformLocation(program, 'u_camUp'),
  camForward: gl.getUniformLocation(program, 'u_camForward'),
  tanHalfFov: gl.getUniformLocation(program, 'u_tanHalfFov'),
  aspect: gl.getUniformLocation(program, 'u_aspect'),
  moving: gl.getUniformLocation(program, 'u_camMoving'),
  resolution: gl.getUniformLocation(program, 'u_resolution'),
  rs: gl.getUniformLocation(program, 'u_rs'),
  disk: {
    r1: gl.getUniformLocation(program, 'u_disk_r1'),
    r2: gl.getUniformLocation(program, 'u_disk_r2'),
  },
  highSteps: gl.getUniformLocation(program, 'u_highSteps'),
};

const camera = new OrbitCamera({
  radius: 6.34194,
  minRadius: 1.0,
  maxRadius: 60.0,
  orbitSpeed: 0.012,
  zoomSpeed: 1.08,
});
camera.attach(canvas);

let highSteps = false;
window.addEventListener('keydown', (e) => {
  if (e.key.toLowerCase() === 'g') highSteps = !highSteps;
});

function render() {
  gl.useProgram(program);
  gl.clearColor(0, 0, 0, 1);
  gl.clear(gl.COLOR_BUFFER_BIT);

  // Camera basis
  const pos = camera.position;
  // Keep both raytracer and grid looking toward +Z axis (BH at origin)
  const target: [number, number, number] = [0, 0, 0];
  const forward = normalize(subtract(target, pos));
  // Ensure world up is +Z to match shader's polar axis and grid warp
  const upWorld: [number, number, number] = [0, 0, 1];
  // Right-handed basis: right = forward × upWorld, up = right × forward
  const right = normalize(cross(forward, upWorld));
  const up = normalize(cross(right, forward));

  // Uniforms
  gl.uniform3fv(loc.camPos, pos);
  gl.uniform3fv(loc.camRight, right);
  gl.uniform3fv(loc.camUp, up);
  gl.uniform3fv(loc.camForward, forward);
  gl.uniform1f(loc.tanHalfFov, Math.tan((60 * Math.PI / 180) * 0.5));
  gl.uniform1f(loc.aspect, internalWidth / internalHeight);
  gl.uniform1i(loc.moving, camera.moving ? 1 : 0);
  gl.uniform2f(loc.resolution, internalWidth, internalHeight);
  gl.uniform1i(loc.highSteps, highSteps ? 1 : 0);
  // Ensure a sane escape radius (render distance) so rays don't bail immediately
  const escapeRLoc = gl.getUniformLocation(program, 'u_escapeR');
  if (escapeRLoc) gl.uniform1f(escapeRLoc, 50000.0);
  const maxStepsLoc = gl.getUniformLocation(program, 'u_maxSteps');
  if (maxStepsLoc) gl.uniform1i(maxStepsLoc, maxSteps);
  // time and simple orbiting object (for frame of reference)
  const t = performance.now() * 0.001;
  const timeLoc = gl.getUniformLocation(program, 'u_time');
  if (timeLoc) gl.uniform1f(timeLoc, t);
  const radius = 20.0;
  const objPos: [number, number, number] = [radius * Math.cos(t*0.2), radius * Math.sin(t*0.2), 0.0];
  const objR = 1.5;
  const objLoc = gl.getUniformLocation(program, 'u_objPosRadius');
  if (objLoc) gl.uniform4f(objLoc, objPos[0], objPos[1], objPos[2], objR);
  const objColorLoc = gl.getUniformLocation(program, 'u_objColor');
  if (objColorLoc) gl.uniform3f(objColorLoc, 1.0, 0.9, 0.6);
  const obj1EnabledLoc = gl.getUniformLocation(program, 'u_obj1Enabled');
  if (obj1EnabledLoc) gl.uniform1i(obj1EnabledLoc, obj1Enabled ? 1 : 0);
  // Second star closer & larger
  const radius2 = 14.0; // farther out
  const obj2Pos: [number, number, number] = [radius2 * Math.cos(t*0.08), radius2 * Math.sin(t*0.08), 0.0];
  const obj2R = 3.0;
  const obj2Loc = gl.getUniformLocation(program, 'u_obj2PosRadius');
  if (obj2Loc) gl.uniform4f(obj2Loc, obj2Pos[0], obj2Pos[1], obj2Pos[2], obj2R);
  const obj2ColorLoc = gl.getUniformLocation(program, 'u_obj2Color');
  if (obj2ColorLoc) gl.uniform3f(obj2ColorLoc, 1.0, 0.85, 0.4); // warmer yellow-orange
  const obj2EnabledLoc = gl.getUniformLocation(program, 'u_obj2Enabled');
  if (obj2EnabledLoc) gl.uniform1i(obj2EnabledLoc, obj2Enabled ? 1 : 0);
  // Background dotted plane behind BH
  const bgZLoc = gl.getUniformLocation(program, 'u_bgZ');
  const bgPeriodLoc = gl.getUniformLocation(program, 'u_bgDotPeriod');
  const bgSizeLoc = gl.getUniformLocation(program, 'u_bgDotSize');
  const bgColorLoc = gl.getUniformLocation(program, 'u_bgColor');
  if (bgZLoc) gl.uniform1f(bgZLoc, -200.0);
  if (bgPeriodLoc) gl.uniform1f(bgPeriodLoc, 8.0);
  if (bgSizeLoc) gl.uniform1f(bgSizeLoc, 0.17);
  if (bgColorLoc) gl.uniform3f(bgColorLoc, 0.9, 0.95, 1.0);
  // Cloud: a static purple-ish sphere off to one side
  const cloudPos = [30.0, 15.0, -50.0] as [number, number, number];
  const cloudR = 12.0;
  const cloudPosLoc = gl.getUniformLocation(program, 'u_cloudPosRadius');
  const cloudColorLoc = gl.getUniformLocation(program, 'u_cloudColor');
  if (cloudPosLoc) gl.uniform4f(cloudPosLoc, cloudPos[0], cloudPos[1], cloudPos[2], cloudR);
  if (cloudColorLoc) gl.uniform3f(cloudColorLoc, 0.7, 0.4, 0.9);

  // Units: scale so Schwarzschild radius = 1.0
  // Keep camera well outside event horizon in scaled units
  gl.uniform1f(loc.rs, 1.0);
  gl.uniform1f(loc.disk.r1, 2.2);
  gl.uniform1f(loc.disk.r2, 5.2);

  // Draw grid first (feed same camera basis to grid shader)
  gl.useProgram(gridProgram);
  const u_color = gl.getUniformLocation(gridProgram, 'u_color');
  gl.uniform3fv(gl.getUniformLocation(gridProgram, 'u_camPos'), pos);
  gl.uniform3fv(gl.getUniformLocation(gridProgram, 'u_camRight'), right);
  gl.uniform3fv(gl.getUniformLocation(gridProgram, 'u_camUp'), up);
  gl.uniform3fv(gl.getUniformLocation(gridProgram, 'u_camForward'), forward);
  gl.uniform1f(gl.getUniformLocation(gridProgram, 'u_aspect'), internalWidth / internalHeight);
  gl.uniform1f(gl.getUniformLocation(gridProgram, 'u_tanHalfFov'), Math.tan((60 * Math.PI / 180) * 0.5));
  gl.uniform4f(u_color, 0.6, 0.6, 0.6, 0.85);

  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  gl.disable(gl.DEPTH_TEST);
  gl.bindVertexArray(gridVAO);
  gl.drawElements(gl.LINES, gridIndexCount, gl.UNSIGNED_INT, 0);
  gl.enable(gl.DEPTH_TEST);
  gl.disable(gl.BLEND);

  // Draw raytraced quad over grid using alpha to composite BH/disk
  gl.useProgram(program);
  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  gl.disable(gl.DEPTH_TEST);
  gl.bindVertexArray(vao);
  gl.drawArrays(gl.TRIANGLES, 0, 6);
  gl.enable(gl.DEPTH_TEST);
  gl.disable(gl.BLEND);

  camera.endFrame();
  requestAnimationFrame(render);
}

requestAnimationFrame(render);

function subtract(a: [number, number, number], b: [number, number, number]): [number, number, number] {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}
function cross(a: [number, number, number], b: [number, number, number]): [number, number, number] {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}
function normalize(a: [number, number, number]): [number, number, number] {
  const l = Math.hypot(a[0], a[1], a[2]) || 1;
  return [a[0] / l, a[1] / l, a[2] / l];
}

// Minimal mat4 utilities (column-major)
function computeViewProj(eye: [number, number, number], center: [number, number, number], up: [number, number, number], fovy: number, aspect: number, near: number, far: number) {
  // Standard right-handed lookAt (column-major)
  const z = normalize([eye[0]-center[0], eye[1]-center[1], eye[2]-center[2]]); // camera forward
  const x = normalize(cross(up, z));
  const y = cross(z, x);

  const view = new Float32Array([
    x[0], x[1], x[2], 0,
    y[0], y[1], y[2], 0,
    z[0], z[1], z[2], 0,
    -dot(x, eye), -dot(y, eye), -dot(z, eye), 1,
  ]);

  const t = Math.tan(fovy / 2);
  // Column-major projection matrix
  const proj = new Float32Array([
    1/(aspect*t), 0, 0, 0,
    0, 1/t, 0, 0,
    0, 0, -(far+near)/(far-near), -1,
    0, 0, -(2*far*near)/(far-near), 0,
  ]);

  return multiply4x4(proj, view);
}
function dot(a: [number, number, number], b: [number, number, number]) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
function multiply4x4(a: Float32Array, b: Float32Array) {
  // Column-major multiplication: out = a * b
  const out = new Float32Array(16);
  for (let col = 0; col < 4; col++) {
    for (let row = 0; row < 4; row++) {
      const ai0 = 0*4 + row, ai1 = 1*4 + row, ai2 = 2*4 + row, ai3 = 3*4 + row;
      const bi0 = col*4 + 0, bi1 = col*4 + 1, bi2 = col*4 + 2, bi3 = col*4 + 3;
      out[col*4 + row] = a[ai0]! * b[bi0]! + a[ai1]! * b[bi1]! + a[ai2]! * b[bi2]! + a[ai3]! * b[bi3]!;
    }
  }
  return out;
}

// Build view-projection using existing basis to ensure consistency with raytracer
function buildViewProjFromBasis(eye: [number, number, number], right: [number, number, number], up: [number, number, number], forward: [number, number, number], fovy: number, aspect: number, near: number, far: number) {
  // Camera looks along forward; for view matrix we need z = -forward in right-handed
  const z: [number, number, number] = [-forward[0], -forward[1], -forward[2]];
  const x = right;
  const y = up;
  const view = new Float32Array([
    x[0], x[1], x[2], 0,
    y[0], y[1], y[2], 0,
    z[0], z[1], z[2], 0,
    -dot(x, eye), -dot(y, eye), -dot(z, eye), 1,
  ]);
  const t = Math.tan(fovy / 2);
  const proj = new Float32Array([
    1/(aspect*t), 0, 0, 0,
    0, 1/t, 0, 0,
    0, 0, -(far+near)/(far-near), -1,
    0, 0, -(2*far*near)/(far-near), 0,
  ]);
  return multiply4x4(proj, view);
}


