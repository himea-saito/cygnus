#version 300 es
precision highp float;

in vec2 vUV;
out vec4 fragColor;

// Camera parameters (avoid struct to prevent driver packing issues)
uniform vec3 u_camPos;
uniform vec3 u_camRight;
uniform vec3 u_camUp;
uniform vec3 u_camForward;
uniform float u_tanHalfFov;
uniform float u_aspect;
uniform bool u_camMoving;
uniform int u_maxSteps;    // user-controlled max iteration count

uniform vec2 u_resolution;  // internal pixel resolution
uniform float u_rs;         // Schwarzschild radius in scene units (we use 1.0)
uniform float u_disk_r1;    // inner disk radius (e.g. 2.2)
uniform float u_disk_r2;    // outer disk radius (e.g. 5.2)
uniform bool u_highSteps;   // request more steps while moving
uniform float u_time;       // seconds
uniform vec4 u_objPosRadius; // xyz = pos, w = radius
uniform vec3 u_objColor;    // object base color
uniform vec4 u_obj2PosRadius; // second object
uniform vec3 u_obj2Color;
uniform bool u_obj1Enabled;
uniform bool u_obj2Enabled;
// Background reference plane (distant dotted sheet) to show lensing
uniform float u_bgZ;            // z-plane (< 0) for background sheet
uniform float u_bgDotPeriod;    // spacing between dots
uniform float u_bgDotSize;      // radius of each dot in cell units
uniform vec3  u_bgColor;        // dot color
// Localized nebula cloud (static world object)
uniform vec4 u_cloudPosRadius; // xyz = position, w = radius
uniform vec3 u_cloudColor;

// Integration constants
const float D_LAMBDA = 2e-3;  // step in affine parameter (scaled units)
uniform float u_escapeR;      // far radius to consider escaped

// Globals for simple object hit (none for now)

struct Ray {
  float x, y, z, r, theta, phi;
  float dr, dtheta, dphi;
  float E, L;
};

Ray initRay(vec3 pos, vec3 dir, float rs) {
  Ray ray;
  ray.x = pos.x; ray.y = pos.y; ray.z = pos.z;
  ray.r = length(pos);
  // Match reference: polar axis along +Z, equatorial plane is Y=0
  ray.theta = acos(pos.z / max(ray.r, 1e-6));
  ray.phi = atan(pos.y, pos.x);

  float dx = dir.x, dy = dir.y, dz = dir.z;
  float st = sin(ray.theta), ct = cos(ray.theta);
  float cp = cos(ray.phi), sp = sin(ray.phi);
  // From reference (axis = Z):
  ray.dr     = st*cp*dx + st*sp*dy + ct*dz;
  ray.dtheta = (ct*cp*dx + ct*sp*dy - st*dz) / max(ray.r, 1e-6);
  ray.dphi   = (-sp*dx + cp*dy) / max(ray.r * max(st, 1e-6), 1e-6);

  ray.L = ray.r * ray.r * max(st, 1e-6) * ray.dphi;
  float f = 1.0 - rs / max(ray.r, 1e-6);
  float dt_dL = sqrt((ray.dr*ray.dr)/max(f,1e-6) + ray.r*ray.r*(ray.dtheta*ray.dtheta + st*st*ray.dphi*ray.dphi));
  ray.E = f * dt_dL;
  return ray;
}

bool intercept(Ray ray, float rs) {
  return ray.r <= rs;
}

bool crossesDisk(vec3 oldPos, vec3 newPos, float r1, float r2) {
  // Disk lies on X-Y plane (z = 0) in our world
  bool crossed = (oldPos.z * newPos.z < 0.0);
  float r = length(vec2(newPos.x, newPos.y));
  return crossed && (r >= r1 && r <= r2);
}

void geodesicRHS(Ray ray, out vec3 d1, out vec3 d2, float rs) {
  float r = ray.r, theta = ray.theta;
  float dr = ray.dr, dtheta = ray.dtheta, dphi = ray.dphi;
  float f = 1.0 - rs / max(r, 1e-6);
  float dt_dL = ray.E / max(f, 1e-6);

  d1 = vec3(dr, dtheta, dphi);
  d2.x = - (rs / (2.0 * r*r)) * f * dt_dL * dt_dL
       + (rs / (2.0 * r*r * max(f,1e-6))) * dr * dr
       + r * (dtheta*dtheta + sin(theta)*sin(theta)*dphi*dphi);
  d2.y = -2.0*dr*dtheta/max(r,1e-6) + sin(theta)*cos(theta)*dphi*dphi;
  d2.z = -2.0*dr*dphi/max(r,1e-6) - 2.0*cos(theta)/max(sin(theta),1e-6) * dtheta * dphi;
}

void rk4Step(inout Ray ray, float dL, float rs) {
  vec3 k1a, k1b; geodesicRHS(ray, k1a, k1b, rs);
  ray.r      += dL * k1a.x;
  ray.theta  += dL * k1a.y;
  ray.phi    += dL * k1a.z;
  ray.dr     += dL * k1b.x;
  ray.dtheta += dL * k1b.y;
  ray.dphi   += dL * k1b.z;
  // Spherical to Cartesian with Z as polar axis
  ray.x = ray.r * sin(ray.theta) * cos(ray.phi);
  ray.y = ray.r * sin(ray.theta) * sin(ray.phi);
  ray.z = ray.r * cos(ray.theta);
}

// --- Nebula background (cheap fBm over direction) ---
float hash31(vec3 p) {
  return fract(sin(dot(p, vec3(12.9898,78.233,37.719))) * 43758.5453);
}
float fbmDir(vec3 d) {
  // d is normalized ray direction; build a stable pseudo-3D fBm
  vec3 p = normalize(d) * 2.0;
  float v = 0.0;
  float a = 0.55;
  vec3 q = p;
  for (int i = 0; i < 5; ++i) {
    v += a * (0.5 + 0.5 * sin(dot(q, vec3(3.1, 5.2, 7.3))));
    q = q * 1.7 + 0.13;
    a *= 0.55;
  }
  return clamp(v, 0.0, 1.0);
}

void main() {
  // Build primary ray
  float u = (2.0 * (gl_FragCoord.x + 0.5) / u_resolution.x - 1.0) * u_aspect * u_tanHalfFov;
  float v = (1.0 - 2.0 * (gl_FragCoord.y + 0.5) / u_resolution.y) * u_tanHalfFov;
  vec3 dir = normalize(u * u_camRight - v * u_camUp + u_camForward);

  // Place camera in scaled units: inputs already in scaled scene units
  Ray ray = initRay(u_camPos, dir, u_rs);

  vec4 color = vec4(0.0);
  vec3 prevPos = vec3(ray.x, ray.y, ray.z);
  bool hitBH = false;
  bool hitDisk = false;
  bool hitObj = false;
  bool hitBg = false;
  bool escaped = false;
  vec3 escDir = vec3(0.0);
  int iter = 0;
  bool hitCloud = false;
  vec3 cloudP = vec3(0.0);
  const float SKY_R = 1000.0; // fixed sphere radius for stable skybox sampling
  vec3 selCenter = vec3(0.0);
  vec3 selColor = vec3(1.0);
  bool selEmissive = false;
  float selRadius = 1.0;

  // Adaptive steps and step size: fewer/larger when moving for performance
  int steps = u_maxSteps;
  if (u_highSteps) steps = int(float(steps) * 1.3);

  for (int i = 0; i < 60000; ++i) {
    if (i >= steps) break;
    if (intercept(ray, u_rs)) { hitBH = true; break; }
    float dL = D_LAMBDA; // fixed step for stability
    rk4Step(ray, dL, u_rs);
    vec3 newPos = vec3(ray.x, ray.y, ray.z);
    // Object intersection (segment-sphere test) for two bodies
    if (!hitObj) {
      vec3 A = prevPos;
      vec3 B = newPos;
      vec3 AB = B - A;
      float a = dot(AB, AB);
      // Obj 1
      if (u_obj1Enabled) {
        vec3 C = u_objPosRadius.xyz;
        float R = u_objPosRadius.w;
        vec3 AC = A - C;
        float bq = 2.0 * dot(AB, AC);
        float cq = dot(AC, AC) - R*R;
        float disc2 = bq*bq - 4.0*a*cq;
        if (disc2 >= 0.0) {
          float s = sqrt(max(disc2, 0.0));
          float u1 = (-bq - s) / (2.0*a);
          float u2 = (-bq + s) / (2.0*a);
          if ((u1 >= 0.0 && u1 <= 1.0) || (u2 >= 0.0 && u2 <= 1.0)) { hitObj = true; selCenter = C; selColor = u_objColor; selEmissive = false; selRadius = R; break; }
        }
      }
      // Obj 2
      if (u_obj2Enabled && !hitObj) {
        vec3 C = u_obj2PosRadius.xyz;
        float R = u_obj2PosRadius.w;
        vec3 AC = A - C;
        float bq = 2.0 * dot(AB, AC);
        float cq = dot(AC, AC) - R*R;
        float disc2 = bq*bq - 4.0*a*cq;
        if (disc2 >= 0.0) {
          float s = sqrt(max(disc2, 0.0));
          float u1 = (-bq - s) / (2.0*a);
          float u2 = (-bq + s) / (2.0*a);
          if ((u1 >= 0.0 && u1 <= 1.0) || (u2 >= 0.0 && u2 <= 1.0)) { hitObj = true; selCenter = C; selColor = u_obj2Color; selEmissive = true; selRadius = R; break; }
        }
      }
    }
    if (crossesDisk(prevPos, newPos, u_disk_r1, u_disk_r2)) { hitDisk = true; break; }
    // Background dotted sheet at z = u_bgZ (behind BH)
    if (!hitBg && u_bgDotSize > 0.0) {
      float z0 = prevPos.z;
      float z1 = newPos.z;
      if ((z0 - u_bgZ) * (z1 - u_bgZ) <= 0.0) {
        float t = clamp((u_bgZ - z0) / max(z1 - z0, 1e-6), 0.0, 1.0);
        vec3 P = mix(prevPos, newPos, t);
        vec2 cell = fract(P.xy / u_bgDotPeriod) - 0.5;
        float d = length(cell);
        if (d < u_bgDotSize) {
          hitBg = true;
          ray.x = P.x; ray.y = P.y; ray.z = P.z; // store for shading below
          break;
        }
      }
    }
    // Localized cloud sphere
    if (!hitCloud) {
      float dC = distance(newPos, u_cloudPosRadius.xyz);
      if (dC <= u_cloudPosRadius.w) {
        hitCloud = true;
        cloudP = newPos;
        break;
      }
    }
    // Stable skybox capture when crossing a fixed radius
    if (!escaped) {
      float r0sq = dot(prevPos, prevPos);
      float r1sq = dot(newPos, newPos);
      float Rsq = SKY_R * SKY_R;
      if (r0sq < Rsq && r1sq >= Rsq) {
        vec3 A = prevPos;
        vec3 V = newPos - prevPos;
        float a = dot(V, V);
        float b = 2.0 * dot(A, V);
        float c = dot(A, A) - Rsq;
        float disc = b*b - 4.0*a*c;
        if (disc >= 0.0) {
          float s = sqrt(max(disc, 0.0));
          float t0 = (-b - s) / (2.0*a);
          float t1 = (-b + s) / (2.0*a);
          float t = (t0 >= 0.0 && t0 <= 1.0) ? t0 : ((t1 >= 0.0 && t1 <= 1.0) ? t1 : -1.0);
          if (t >= 0.0) {
            vec3 P = A + t * V;
            escaped = true;
            escDir = normalize(P);
            break;
          }
        }
      }
    }

    prevPos = newPos;
    if (ray.r > u_escapeR) { escaped = true; escDir = normalize(newPos); break; }
    iter++;
  }

  // If we ran out of iterations without a hit or escape, treat as skybox hit
  if (!hitBH && !hitDisk && !hitObj && !escaped) {
    escaped = true;
    escDir = normalize(vec3(ray.x, ray.y, ray.z));
  }


  if (hitDisk) {
    float r = length(vec3(ray.x, ray.y, ray.z)) / u_disk_r2;
    vec3 diskColor = vec3(1.0, r, 0.2);
    // simple limb brightening
    float a = clamp(1.5 - r, 0.0, 1.0);
    color = vec4(mix(diskColor * 0.6, diskColor, a), 1.0);
  } else if (hitBH) {
    // Solid alpha so BH occludes grid
    color = vec4(0.0, 0.0, 0.0, 1.0);
  } else if (hitObj) {
    vec3 P = vec3(ray.x, ray.y, ray.z);
    if (selEmissive) {
      // Emissive star with hotter core using screen-projected radial falloff
      vec3 toP = P - selCenter;
      vec3 Vdir = normalize(u_camPos - P);
      vec3 toPperp = toP - dot(toP, Vdir) * Vdir; // component perpendicular to view
      float rProj = length(toPperp) / max(selRadius, 1e-6); // 0=center, 1=edge
      // Wider orange core: weight stays high until ~80% of radius
      float w = smoothstep(0.8, 0.0, rProj);
      vec3 edgeCol = selColor;                // outer (yellowish)
      vec3 coreCol = vec3(1.0, 0.6, 0.2);     // inner (orange)
      vec3 grad = mix(edgeCol, coreCol, w);
      color = vec4(grad, 1.0);
    } else {
      // Lambert toward camera for non-emissive bodies
      vec3 N = normalize(P - selCenter);
      vec3 V = normalize(u_camPos - P);
      float ambient = 0.1;
      float diff = max(dot(N, V), 0.0);
      float intensity = ambient + (1.0 - ambient) * diff;
      vec3 shaded = selColor * intensity;
      color = vec4(shaded, 1.0);
    }
  } else if (hitBg) {
    color = vec4(u_bgColor, 0.9);
  } else if (hitCloud) {
    // Soft volumetric look via radial falloff in the cloud sphere
    float r = distance(cloudP, u_cloudPosRadius.xyz) / max(u_cloudPosRadius.w, 1e-6);
    float a = smoothstep(1.0, 0.2, r);
    vec3 c = u_cloudColor * (0.4 + 0.6 * a);
    color = vec4(c, a * 0.6);
  } else if (escaped) {
    // Minimal skybox: a few fixed world-space stars on the unit sphere
    vec3 d = normalize(escDir);
    const int NUM_STARS = 16;
    const vec3 STARS[16] = vec3[16](
      vec3( 0.923,  0.000,  0.384), vec3(-0.707,  0.000,  0.707),
      vec3( 0.000,  0.923,  0.384), vec3( 0.000, -0.707,  0.707),
      vec3( 0.577,  0.577,  0.577), vec3(-0.577,  0.577,  0.577),
      vec3( 0.577, -0.577,  0.577), vec3(-0.577, -0.577,  0.577),
      vec3( 0.965,  0.259,  0.000), vec3(-0.965,  0.259,  0.000),
      vec3( 0.259,  0.965,  0.000), vec3( 0.259, -0.965,  0.000),
      vec3( 0.866,  0.500,  0.000), vec3(-0.866,  0.500,  0.000),
      vec3( 0.500, -0.866,  0.000), vec3(-0.500, -0.866,  0.000)
    );
    float brightness = 0.0;
    for (int i = 0; i < NUM_STARS; ++i) {
      brightness = max(brightness, smoothstep(0.9998, 1.0, dot(d, normalize(STARS[i]))));
    }
    color = vec4(vec3(brightness), brightness * 0.9);
  } else {
    // Empty space
    color = vec4(0.0);
  }

  fragColor = color;
}


