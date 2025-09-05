#version 300 es
precision highp float;

layout(location = 0) in vec3 aPos;

// Use the exact same camera basis as the raytracer to avoid mismatch
uniform vec3 u_camPos;
uniform vec3 u_camRight;
uniform vec3 u_camUp;
uniform vec3 u_camForward;
uniform float u_aspect;
uniform float u_tanHalfFov;

void main() {
  vec3 d = aPos - u_camPos;
  float cx = dot(d, u_camRight);
  float cy = dot(d, u_camUp);
  float cz = dot(d, u_camForward);
  // Convert to NDC using perspective with forward along -Z in view space
  float zView = -cz; // standard view space z (negative in front)
  // distance in front of camera (positive)
  float dist = max(-zView, 1e-3);
  float x_ndc = (cx / (dist * u_tanHalfFov * u_aspect));
  float y_ndc = (cy / (dist * u_tanHalfFov));
  gl_Position = vec4(x_ndc, y_ndc, 0.0, 1.0);
}


