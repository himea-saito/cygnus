(function(){const o=document.createElement("link").relList;if(o&&o.supports&&o.supports("modulepreload"))return;for(const t of document.querySelectorAll('link[rel="modulepreload"]'))s(t);new MutationObserver(t=>{for(const n of t)if(n.type==="childList")for(const c of n.addedNodes)c.tagName==="LINK"&&c.rel==="modulepreload"&&s(c)}).observe(document,{childList:!0,subtree:!0});function a(t){const n={};return t.integrity&&(n.integrity=t.integrity),t.referrerPolicy&&(n.referrerPolicy=t.referrerPolicy),t.crossOrigin==="use-credentials"?n.credentials="include":t.crossOrigin==="anonymous"?n.credentials="omit":n.credentials="same-origin",n}function s(t){if(t.ep)return;t.ep=!0;const n=a(t);fetch(t.href,n)}})();function N(e,o,a){const s=e.createShader(o);if(e.shaderSource(s,a),e.compileShader(s),!e.getShaderParameter(s,e.COMPILE_STATUS)){const t=e.getShaderInfoLog(s);throw e.deleteShader(s),new Error(`Shader compile error: ${t}
${a}`)}return s}function Q(e,o,a){const s=N(e,e.VERTEX_SHADER,o),t=N(e,e.FRAGMENT_SHADER,a),n=e.createProgram();if(e.attachShader(n,s),e.attachShader(n,t),e.bindAttribLocation(n,0,"aPos"),e.bindAttribLocation(n,1,"aTexCoord"),e.linkProgram(n),!e.getProgramParameter(n,e.LINK_STATUS)){const c=e.getProgramInfoLog(n);throw e.deleteProgram(n),new Error(`Program link error: ${c}`)}return e.deleteShader(s),e.deleteShader(t),n}const lr=`#version 300 es\r
precision highp float;\r
layout(location = 0) in vec2 aPos;\r
layout(location = 1) in vec2 aTexCoord;\r
out vec2 vUV;\r
void main() {\r
  gl_Position = vec4(aPos, 0.0, 1.0);\r
  vUV = aTexCoord;\r
}\r
\r
\r
`,fr=`#version 300 es\r
precision highp float;\r
\r
in vec2 vUV;\r
out vec4 fragColor;\r
\r
// Camera parameters (avoid struct to prevent driver packing issues)\r
uniform vec3 u_camPos;\r
uniform vec3 u_camRight;\r
uniform vec3 u_camUp;\r
uniform vec3 u_camForward;\r
uniform float u_tanHalfFov;\r
uniform float u_aspect;\r
uniform bool u_camMoving;\r
uniform int u_maxSteps;    // user-controlled max iteration count\r
\r
uniform vec2 u_resolution;  // internal pixel resolution\r
uniform float u_rs;         // Schwarzschild radius in scene units (we use 1.0)\r
uniform float u_disk_r1;    // inner disk radius (e.g. 2.2)\r
uniform float u_disk_r2;    // outer disk radius (e.g. 5.2)\r
uniform bool u_highSteps;   // request more steps while moving\r
uniform float u_time;       // seconds\r
uniform vec4 u_objPosRadius; // xyz = pos, w = radius\r
uniform vec3 u_objColor;    // object base color\r
uniform vec4 u_obj2PosRadius; // second object\r
uniform vec3 u_obj2Color;\r
uniform bool u_obj1Enabled;\r
uniform bool u_obj2Enabled;\r
// Background reference plane (distant dotted sheet) to show lensing\r
uniform float u_bgZ;            // z-plane (< 0) for background sheet\r
uniform float u_bgDotPeriod;    // spacing between dots\r
uniform float u_bgDotSize;      // radius of each dot in cell units\r
uniform vec3  u_bgColor;        // dot color\r
// Localized nebula cloud (static world object)\r
uniform vec4 u_cloudPosRadius; // xyz = position, w = radius\r
uniform vec3 u_cloudColor;\r
\r
// Integration constants\r
const float D_LAMBDA = 2e-3;  // step in affine parameter (scaled units)\r
uniform float u_escapeR;      // far radius to consider escaped\r
\r
// Globals for simple object hit (none for now)\r
\r
struct Ray {\r
  float x, y, z, r, theta, phi;\r
  float dr, dtheta, dphi;\r
  float E, L;\r
};\r
\r
Ray initRay(vec3 pos, vec3 dir, float rs) {\r
  Ray ray;\r
  ray.x = pos.x; ray.y = pos.y; ray.z = pos.z;\r
  ray.r = length(pos);\r
  // Match reference: polar axis along +Z, equatorial plane is Y=0\r
  ray.theta = acos(pos.z / max(ray.r, 1e-6));\r
  ray.phi = atan(pos.y, pos.x);\r
\r
  float dx = dir.x, dy = dir.y, dz = dir.z;\r
  float st = sin(ray.theta), ct = cos(ray.theta);\r
  float cp = cos(ray.phi), sp = sin(ray.phi);\r
  // From reference (axis = Z):\r
  ray.dr     = st*cp*dx + st*sp*dy + ct*dz;\r
  ray.dtheta = (ct*cp*dx + ct*sp*dy - st*dz) / max(ray.r, 1e-6);\r
  ray.dphi   = (-sp*dx + cp*dy) / max(ray.r * max(st, 1e-6), 1e-6);\r
\r
  ray.L = ray.r * ray.r * max(st, 1e-6) * ray.dphi;\r
  float f = 1.0 - rs / max(ray.r, 1e-6);\r
  float dt_dL = sqrt((ray.dr*ray.dr)/max(f,1e-6) + ray.r*ray.r*(ray.dtheta*ray.dtheta + st*st*ray.dphi*ray.dphi));\r
  ray.E = f * dt_dL;\r
  return ray;\r
}\r
\r
bool intercept(Ray ray, float rs) {\r
  return ray.r <= rs;\r
}\r
\r
bool crossesDisk(vec3 oldPos, vec3 newPos, float r1, float r2) {\r
  // Disk lies on X-Y plane (z = 0) in our world\r
  bool crossed = (oldPos.z * newPos.z < 0.0);\r
  float r = length(vec2(newPos.x, newPos.y));\r
  return crossed && (r >= r1 && r <= r2);\r
}\r
\r
void geodesicRHS(Ray ray, out vec3 d1, out vec3 d2, float rs) {\r
  float r = ray.r, theta = ray.theta;\r
  float dr = ray.dr, dtheta = ray.dtheta, dphi = ray.dphi;\r
  float f = 1.0 - rs / max(r, 1e-6);\r
  float dt_dL = ray.E / max(f, 1e-6);\r
\r
  d1 = vec3(dr, dtheta, dphi);\r
  d2.x = - (rs / (2.0 * r*r)) * f * dt_dL * dt_dL\r
       + (rs / (2.0 * r*r * max(f,1e-6))) * dr * dr\r
       + r * (dtheta*dtheta + sin(theta)*sin(theta)*dphi*dphi);\r
  d2.y = -2.0*dr*dtheta/max(r,1e-6) + sin(theta)*cos(theta)*dphi*dphi;\r
  d2.z = -2.0*dr*dphi/max(r,1e-6) - 2.0*cos(theta)/max(sin(theta),1e-6) * dtheta * dphi;\r
}\r
\r
void rk4Step(inout Ray ray, float dL, float rs) {\r
  vec3 k1a, k1b; geodesicRHS(ray, k1a, k1b, rs);\r
  ray.r      += dL * k1a.x;\r
  ray.theta  += dL * k1a.y;\r
  ray.phi    += dL * k1a.z;\r
  ray.dr     += dL * k1b.x;\r
  ray.dtheta += dL * k1b.y;\r
  ray.dphi   += dL * k1b.z;\r
  // Spherical to Cartesian with Z as polar axis\r
  ray.x = ray.r * sin(ray.theta) * cos(ray.phi);\r
  ray.y = ray.r * sin(ray.theta) * sin(ray.phi);\r
  ray.z = ray.r * cos(ray.theta);\r
}\r
\r
// --- Nebula background (cheap fBm over direction) ---\r
float hash31(vec3 p) {\r
  return fract(sin(dot(p, vec3(12.9898,78.233,37.719))) * 43758.5453);\r
}\r
float fbmDir(vec3 d) {\r
  // d is normalized ray direction; build a stable pseudo-3D fBm\r
  vec3 p = normalize(d) * 2.0;\r
  float v = 0.0;\r
  float a = 0.55;\r
  vec3 q = p;\r
  for (int i = 0; i < 5; ++i) {\r
    v += a * (0.5 + 0.5 * sin(dot(q, vec3(3.1, 5.2, 7.3))));\r
    q = q * 1.7 + 0.13;\r
    a *= 0.55;\r
  }\r
  return clamp(v, 0.0, 1.0);\r
}\r
\r
void main() {\r
  // Build primary ray\r
  float u = (2.0 * (gl_FragCoord.x + 0.5) / u_resolution.x - 1.0) * u_aspect * u_tanHalfFov;\r
  float v = (1.0 - 2.0 * (gl_FragCoord.y + 0.5) / u_resolution.y) * u_tanHalfFov;\r
  vec3 dir = normalize(u * u_camRight - v * u_camUp + u_camForward);\r
\r
  // Place camera in scaled units: inputs already in scaled scene units\r
  Ray ray = initRay(u_camPos, dir, u_rs);\r
\r
  vec4 color = vec4(0.0);\r
  vec3 prevPos = vec3(ray.x, ray.y, ray.z);\r
  bool hitBH = false;\r
  bool hitDisk = false;\r
  bool hitObj = false;\r
  bool hitBg = false;\r
  bool escaped = false;\r
  vec3 escDir = vec3(0.0);\r
  int iter = 0;\r
  bool hitCloud = false;\r
  vec3 cloudP = vec3(0.0);\r
  const float SKY_R = 1000.0; // fixed sphere radius for stable skybox sampling\r
  vec3 selCenter = vec3(0.0);\r
  vec3 selColor = vec3(1.0);\r
  bool selEmissive = false;\r
  float selRadius = 1.0;\r
\r
  // Adaptive steps and step size: fewer/larger when moving for performance\r
  int steps = u_maxSteps;\r
  if (u_highSteps) steps = int(float(steps) * 1.3);\r
\r
  for (int i = 0; i < 60000; ++i) {\r
    if (i >= steps) break;\r
    if (intercept(ray, u_rs)) { hitBH = true; break; }\r
    float dL = D_LAMBDA; // fixed step for stability\r
    rk4Step(ray, dL, u_rs);\r
    vec3 newPos = vec3(ray.x, ray.y, ray.z);\r
    // Object intersection (segment-sphere test) for two bodies\r
    if (!hitObj) {\r
      vec3 A = prevPos;\r
      vec3 B = newPos;\r
      vec3 AB = B - A;\r
      float a = dot(AB, AB);\r
      // Obj 1\r
      if (u_obj1Enabled) {\r
        vec3 C = u_objPosRadius.xyz;\r
        float R = u_objPosRadius.w;\r
        vec3 AC = A - C;\r
        float bq = 2.0 * dot(AB, AC);\r
        float cq = dot(AC, AC) - R*R;\r
        float disc2 = bq*bq - 4.0*a*cq;\r
        if (disc2 >= 0.0) {\r
          float s = sqrt(max(disc2, 0.0));\r
          float u1 = (-bq - s) / (2.0*a);\r
          float u2 = (-bq + s) / (2.0*a);\r
          if ((u1 >= 0.0 && u1 <= 1.0) || (u2 >= 0.0 && u2 <= 1.0)) { hitObj = true; selCenter = C; selColor = u_objColor; selEmissive = false; selRadius = R; break; }\r
        }\r
      }\r
      // Obj 2\r
      if (u_obj2Enabled && !hitObj) {\r
        vec3 C = u_obj2PosRadius.xyz;\r
        float R = u_obj2PosRadius.w;\r
        vec3 AC = A - C;\r
        float bq = 2.0 * dot(AB, AC);\r
        float cq = dot(AC, AC) - R*R;\r
        float disc2 = bq*bq - 4.0*a*cq;\r
        if (disc2 >= 0.0) {\r
          float s = sqrt(max(disc2, 0.0));\r
          float u1 = (-bq - s) / (2.0*a);\r
          float u2 = (-bq + s) / (2.0*a);\r
          if ((u1 >= 0.0 && u1 <= 1.0) || (u2 >= 0.0 && u2 <= 1.0)) { hitObj = true; selCenter = C; selColor = u_obj2Color; selEmissive = true; selRadius = R; break; }\r
        }\r
      }\r
    }\r
    if (crossesDisk(prevPos, newPos, u_disk_r1, u_disk_r2)) { hitDisk = true; break; }\r
    // Background dotted sheet at z = u_bgZ (behind BH)\r
    if (!hitBg && u_bgDotSize > 0.0) {\r
      float z0 = prevPos.z;\r
      float z1 = newPos.z;\r
      if ((z0 - u_bgZ) * (z1 - u_bgZ) <= 0.0) {\r
        float t = clamp((u_bgZ - z0) / max(z1 - z0, 1e-6), 0.0, 1.0);\r
        vec3 P = mix(prevPos, newPos, t);\r
        vec2 cell = fract(P.xy / u_bgDotPeriod) - 0.5;\r
        float d = length(cell);\r
        if (d < u_bgDotSize) {\r
          hitBg = true;\r
          ray.x = P.x; ray.y = P.y; ray.z = P.z; // store for shading below\r
          break;\r
        }\r
      }\r
    }\r
    // Localized cloud sphere\r
    if (!hitCloud) {\r
      float dC = distance(newPos, u_cloudPosRadius.xyz);\r
      if (dC <= u_cloudPosRadius.w) {\r
        hitCloud = true;\r
        cloudP = newPos;\r
        break;\r
      }\r
    }\r
    // Stable skybox capture when crossing a fixed radius\r
    if (!escaped) {\r
      float r0sq = dot(prevPos, prevPos);\r
      float r1sq = dot(newPos, newPos);\r
      float Rsq = SKY_R * SKY_R;\r
      if (r0sq < Rsq && r1sq >= Rsq) {\r
        vec3 A = prevPos;\r
        vec3 V = newPos - prevPos;\r
        float a = dot(V, V);\r
        float b = 2.0 * dot(A, V);\r
        float c = dot(A, A) - Rsq;\r
        float disc = b*b - 4.0*a*c;\r
        if (disc >= 0.0) {\r
          float s = sqrt(max(disc, 0.0));\r
          float t0 = (-b - s) / (2.0*a);\r
          float t1 = (-b + s) / (2.0*a);\r
          float t = (t0 >= 0.0 && t0 <= 1.0) ? t0 : ((t1 >= 0.0 && t1 <= 1.0) ? t1 : -1.0);\r
          if (t >= 0.0) {\r
            vec3 P = A + t * V;\r
            escaped = true;\r
            escDir = normalize(P);\r
            break;\r
          }\r
        }\r
      }\r
    }\r
\r
    prevPos = newPos;\r
    if (ray.r > u_escapeR) { escaped = true; escDir = normalize(newPos); break; }\r
    iter++;\r
  }\r
\r
  // If we ran out of iterations without a hit or escape, treat as skybox hit\r
  if (!hitBH && !hitDisk && !hitObj && !escaped) {\r
    escaped = true;\r
    escDir = normalize(vec3(ray.x, ray.y, ray.z));\r
  }\r
\r
\r
  if (hitDisk) {\r
    float r = length(vec3(ray.x, ray.y, ray.z)) / u_disk_r2;\r
    vec3 diskColor = vec3(1.0, r, 0.2);\r
    // simple limb brightening\r
    float a = clamp(1.5 - r, 0.0, 1.0);\r
    color = vec4(mix(diskColor * 0.6, diskColor, a), 1.0);\r
  } else if (hitBH) {\r
    // Solid alpha so BH occludes grid\r
    color = vec4(0.0, 0.0, 0.0, 1.0);\r
  } else if (hitObj) {\r
    vec3 P = vec3(ray.x, ray.y, ray.z);\r
    if (selEmissive) {\r
      // Emissive star with hotter core using screen-projected radial falloff\r
      vec3 toP = P - selCenter;\r
      vec3 Vdir = normalize(u_camPos - P);\r
      vec3 toPperp = toP - dot(toP, Vdir) * Vdir; // component perpendicular to view\r
      float rProj = length(toPperp) / max(selRadius, 1e-6); // 0=center, 1=edge\r
      // Wider orange core: weight stays high until ~80% of radius\r
      float w = smoothstep(0.8, 0.0, rProj);\r
      vec3 edgeCol = selColor;                // outer (yellowish)\r
      vec3 coreCol = vec3(1.0, 0.6, 0.2);     // inner (orange)\r
      vec3 grad = mix(edgeCol, coreCol, w);\r
      color = vec4(grad, 1.0);\r
    } else {\r
      // Lambert toward camera for non-emissive bodies\r
      vec3 N = normalize(P - selCenter);\r
      vec3 V = normalize(u_camPos - P);\r
      float ambient = 0.1;\r
      float diff = max(dot(N, V), 0.0);\r
      float intensity = ambient + (1.0 - ambient) * diff;\r
      vec3 shaded = selColor * intensity;\r
      color = vec4(shaded, 1.0);\r
    }\r
  } else if (hitBg) {\r
    color = vec4(u_bgColor, 0.9);\r
  } else if (hitCloud) {\r
    // Soft volumetric look via radial falloff in the cloud sphere\r
    float r = distance(cloudP, u_cloudPosRadius.xyz) / max(u_cloudPosRadius.w, 1e-6);\r
    float a = smoothstep(1.0, 0.2, r);\r
    vec3 c = u_cloudColor * (0.4 + 0.6 * a);\r
    color = vec4(c, a * 0.6);\r
  } else if (escaped) {\r
    // Minimal skybox: a few fixed world-space stars on the unit sphere\r
    vec3 d = normalize(escDir);\r
    const int NUM_STARS = 16;\r
    const vec3 STARS[16] = vec3[16](\r
      vec3( 0.923,  0.000,  0.384), vec3(-0.707,  0.000,  0.707),\r
      vec3( 0.000,  0.923,  0.384), vec3( 0.000, -0.707,  0.707),\r
      vec3( 0.577,  0.577,  0.577), vec3(-0.577,  0.577,  0.577),\r
      vec3( 0.577, -0.577,  0.577), vec3(-0.577, -0.577,  0.577),\r
      vec3( 0.965,  0.259,  0.000), vec3(-0.965,  0.259,  0.000),\r
      vec3( 0.259,  0.965,  0.000), vec3( 0.259, -0.965,  0.000),\r
      vec3( 0.866,  0.500,  0.000), vec3(-0.866,  0.500,  0.000),\r
      vec3( 0.500, -0.866,  0.000), vec3(-0.500, -0.866,  0.000)\r
    );\r
    float brightness = 0.0;\r
    for (int i = 0; i < NUM_STARS; ++i) {\r
      brightness = max(brightness, smoothstep(0.9998, 1.0, dot(d, normalize(STARS[i]))));\r
    }\r
    color = vec4(vec3(brightness), brightness * 0.9);\r
  } else {\r
    // Empty space\r
    color = vec4(0.0);\r
  }\r
\r
  fragColor = color;\r
}\r
\r
\r
`,ur=`#version 300 es\r
precision highp float;\r
\r
layout(location = 0) in vec3 aPos;\r
\r
// Use the exact same camera basis as the raytracer to avoid mismatch\r
uniform vec3 u_camPos;\r
uniform vec3 u_camRight;\r
uniform vec3 u_camUp;\r
uniform vec3 u_camForward;\r
uniform float u_aspect;\r
uniform float u_tanHalfFov;\r
\r
void main() {\r
  vec3 d = aPos - u_camPos;\r
  float cx = dot(d, u_camRight);\r
  float cy = dot(d, u_camUp);\r
  float cz = dot(d, u_camForward);\r
  // Convert to NDC using perspective with forward along -Z in view space\r
  float zView = -cz; // standard view space z (negative in front)\r
  // distance in front of camera (positive)\r
  float dist = max(-zView, 1e-3);\r
  float x_ndc = (cx / (dist * u_tanHalfFov * u_aspect));\r
  float y_ndc = (cy / (dist * u_tanHalfFov));\r
  gl_Position = vec4(x_ndc, y_ndc, 0.0, 1.0);\r
}\r
\r
\r
`,mr=`#version 300 es\r
precision highp float;\r
\r
out vec4 fragColor;\r
uniform vec4 u_color;\r
\r
void main() {\r
  fragColor = u_color;\r
}\r
\r
\r
`;class hr{radius;minRadius;maxRadius;azimuth=0;elevation=Math.PI/3;orbitSpeed;zoomSpeed;dragging=!1;moving=!1;lastX=0;lastY=0;constructor(o){this.radius=o.radius,this.minRadius=o.minRadius,this.maxRadius=o.maxRadius,this.orbitSpeed=o.orbitSpeed,this.zoomSpeed=o.zoomSpeed}attach(o){o.addEventListener("mousedown",a=>{a.button===0&&(this.dragging=!0,this.moving=!0,this.lastX=a.clientX,this.lastY=a.clientY)}),window.addEventListener("mouseup",()=>{this.dragging=!1}),window.addEventListener("mousemove",a=>{if(!this.dragging)return;const s=a.clientX-this.lastX,t=a.clientY-this.lastY;this.azimuth+=s*this.orbitSpeed,this.elevation=Y(this.elevation+t*this.orbitSpeed,.01,Math.PI-.01),this.lastX=a.clientX,this.lastY=a.clientY}),o.addEventListener("wheel",a=>{a.preventDefault(),a.deltaY<0?this.radius/=Math.pow(this.zoomSpeed,1):this.radius*=Math.pow(this.zoomSpeed,1),this.radius=Y(this.radius,this.minRadius,this.maxRadius),this.moving=!0},{passive:!1})}get position(){const o=this.radius,a=this.elevation,s=this.azimuth;return[o*Math.sin(a)*Math.cos(s),o*Math.sin(a)*Math.sin(s),o*Math.cos(a)]}endFrame(){this.dragging||(this.moving=!1)}}function Y(e,o,a){return Math.min(a,Math.max(o,e))}const P=document.getElementById("glcanvas"),Z=document.getElementById("err"),R=document.getElementById("overlay"),W=document.getElementById("startBtn"),X=document.getElementById("rememberStart");let G=!1;const rr=P.getContext("webgl2",{antialias:!1,preserveDrawingBuffer:!1});if(!rr)throw Z&&(Z.textContent="WebGL2 not supported or blocked. Try a different browser/device."),new Error("WebGL2 not supported");const r=rr,S=document.getElementById("maxSteps"),K=document.getElementById("maxStepsVal");let A=2e4;const b=document.getElementById("obj1Toggle"),y=document.getElementById("obj2Toggle");let z=!0,C=!0;if(S&&K){const e=()=>{A=parseInt(S.value,10);const o=Math.round(A/1e3);K.textContent=`${o}k`};S.addEventListener("input",e),e()}b&&(b.addEventListener("change",()=>{z=b.checked}),z=b.checked);y&&(y.addEventListener("change",()=>{C=y.checked}),C=y.checked);let p=200,g=150;function er(){const e=Math.max(100,Math.floor(window.innerWidth/8)),o=Math.max(75,Math.floor(window.innerHeight/8));p=Math.min(320,e),g=Math.min(240,o),P.width=p,P.height=g,r.viewport(0,0,p,g)}er();window.addEventListener("resize",er);const i=Q(r,lr,fr);r.useProgram(i);const m=Q(r,ur,mr),or=r.createVertexArray(),pr=r.createBuffer(),gr=r.createBuffer();let nr=0;function vr(){const s=[],t=[];for(let n=0;n<=25;n++)for(let c=0;c<=25;c++){const l=(c-12.5)*1,f=(n-25/2)*1,u=1,v=Math.hypot(l,f);let h=0;v>u?h=-10+2*Math.sqrt(Math.max(u*(v-u),0)):h=-10+2*Math.sqrt(u*u),s.push(l,f,h)}for(let n=0;n<25;n++)for(let c=0;c<25;c++){const l=n*26+c;t.push(l,l+1),t.push(l,l+25+1)}r.bindVertexArray(or),r.bindBuffer(r.ARRAY_BUFFER,pr),r.bufferData(r.ARRAY_BUFFER,new Float32Array(s),r.DYNAMIC_DRAW),r.enableVertexAttribArray(0),r.vertexAttribPointer(0,3,r.FLOAT,!1,12,0),r.bindBuffer(r.ELEMENT_ARRAY_BUFFER,gr),r.bufferData(r.ELEMENT_ARRAY_BUFFER,new Uint32Array(t),r.STATIC_DRAW),nr=t.length}vr();const tr=r.createVertexArray();r.bindVertexArray(tr);const br=new Float32Array([-1,1,0,1,-1,-1,0,0,1,-1,1,0,-1,1,0,1,1,-1,1,0,1,1,1,1]),yr=r.createBuffer();r.bindBuffer(r.ARRAY_BUFFER,yr);r.bufferData(r.ARRAY_BUFFER,br,r.STATIC_DRAW);r.enableVertexAttribArray(0);r.vertexAttribPointer(0,2,r.FLOAT,!1,16,0);r.enableVertexAttribArray(1);r.vertexAttribPointer(1,2,r.FLOAT,!1,16,8);const d={camPos:r.getUniformLocation(i,"u_camPos"),camRight:r.getUniformLocation(i,"u_camRight"),camUp:r.getUniformLocation(i,"u_camUp"),camForward:r.getUniformLocation(i,"u_camForward"),tanHalfFov:r.getUniformLocation(i,"u_tanHalfFov"),aspect:r.getUniformLocation(i,"u_aspect"),moving:r.getUniformLocation(i,"u_camMoving"),resolution:r.getUniformLocation(i,"u_resolution"),rs:r.getUniformLocation(i,"u_rs"),disk:{r1:r.getUniformLocation(i,"u_disk_r1"),r2:r.getUniformLocation(i,"u_disk_r2")},highSteps:r.getUniformLocation(i,"u_highSteps")},_=new hr({radius:6.34194,minRadius:1,maxRadius:60,orbitSpeed:.012,zoomSpeed:1.08});_.attach(P);let E=!1;window.addEventListener("keydown",e=>{e.key.toLowerCase()==="g"&&(E=!E)});function ar(){r.useProgram(i),r.clearColor(0,0,0,1),r.clear(r.COLOR_BUFFER_BIT);const e=_.position,a=w(Pr([0,0,0],e)),t=w(J(a,[0,0,1])),n=w(J(t,a));r.uniform3fv(d.camPos,e),r.uniform3fv(d.camRight,t),r.uniform3fv(d.camUp,n),r.uniform3fv(d.camForward,a),r.uniform1f(d.tanHalfFov,Math.tan(60*Math.PI/180*.5)),r.uniform1f(d.aspect,p/g),r.uniform1i(d.moving,_.moving?1:0),r.uniform2f(d.resolution,p,g),r.uniform1i(d.highSteps,E?1:0);const c=r.getUniformLocation(i,"u_escapeR");c&&r.uniform1f(c,5e4);const l=r.getUniformLocation(i,"u_maxSteps");l&&r.uniform1i(l,A);const f=performance.now()*.001,u=r.getUniformLocation(i,"u_time");u&&r.uniform1f(u,f);const v=20,h=[v*Math.cos(f*.2),v*Math.sin(f*.2),0],ir=1.5,k=r.getUniformLocation(i,"u_objPosRadius");k&&r.uniform4f(k,h[0],h[1],h[2],ir);const U=r.getUniformLocation(i,"u_objColor");U&&r.uniform3f(U,1,.9,.6);const B=r.getUniformLocation(i,"u_obj1Enabled");B&&r.uniform1i(B,z?1:0);const j=14,x=[j*Math.cos(f*.08),j*Math.sin(f*.08),0],sr=3,F=r.getUniformLocation(i,"u_obj2PosRadius");F&&r.uniform4f(F,x[0],x[1],x[2],sr);const D=r.getUniformLocation(i,"u_obj2Color");D&&r.uniform3f(D,1,.85,.4);const M=r.getUniformLocation(i,"u_obj2Enabled");M&&r.uniform1i(M,C?1:0);const T=r.getUniformLocation(i,"u_bgZ"),q=r.getUniformLocation(i,"u_bgDotPeriod"),I=r.getUniformLocation(i,"u_bgDotSize"),O=r.getUniformLocation(i,"u_bgColor");T&&r.uniform1f(T,-200),q&&r.uniform1f(q,8),I&&r.uniform1f(I,.17),O&&r.uniform3f(O,.9,.95,1);const L=[30,15,-50],cr=12,H=r.getUniformLocation(i,"u_cloudPosRadius"),V=r.getUniformLocation(i,"u_cloudColor");H&&r.uniform4f(H,L[0],L[1],L[2],cr),V&&r.uniform3f(V,.7,.4,.9),r.uniform1f(d.rs,1),r.uniform1f(d.disk.r1,2.2),r.uniform1f(d.disk.r2,5.2),r.useProgram(m);const dr=r.getUniformLocation(m,"u_color");r.uniform3fv(r.getUniformLocation(m,"u_camPos"),e),r.uniform3fv(r.getUniformLocation(m,"u_camRight"),t),r.uniform3fv(r.getUniformLocation(m,"u_camUp"),n),r.uniform3fv(r.getUniformLocation(m,"u_camForward"),a),r.uniform1f(r.getUniformLocation(m,"u_aspect"),p/g),r.uniform1f(r.getUniformLocation(m,"u_tanHalfFov"),Math.tan(60*Math.PI/180*.5)),r.uniform4f(dr,.6,.6,.6,.85),r.enable(r.BLEND),r.blendFunc(r.SRC_ALPHA,r.ONE_MINUS_SRC_ALPHA),r.disable(r.DEPTH_TEST),r.bindVertexArray(or),r.drawElements(r.LINES,nr,r.UNSIGNED_INT,0),r.enable(r.DEPTH_TEST),r.disable(r.BLEND),r.useProgram(i),r.enable(r.BLEND),r.blendFunc(r.SRC_ALPHA,r.ONE_MINUS_SRC_ALPHA),r.disable(r.DEPTH_TEST),r.bindVertexArray(tr),r.drawArrays(r.TRIANGLES,0,6),r.enable(r.DEPTH_TEST),r.disable(r.BLEND),_.endFrame(),requestAnimationFrame(ar)}function $(){G||(G=!0,requestAnimationFrame(ar))}function _r(){try{const e="cygnus_ack_gpu",o=localStorage.getItem(e)==="1";R&&(o&&(R.style.display="none",$()),W&&(W.onclick=()=>{R.style.display="none",X&&X.checked&&localStorage.setItem(e,"1"),$()}))}catch{}}_r();function Pr(e,o){return[e[0]-o[0],e[1]-o[1],e[2]-o[2]]}function J(e,o){return[e[1]*o[2]-e[2]*o[1],e[2]*o[0]-e[0]*o[2],e[0]*o[1]-e[1]*o[0]]}function w(e){const o=Math.hypot(e[0],e[1],e[2])||1;return[e[0]/o,e[1]/o,e[2]/o]}
