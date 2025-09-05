type Vec3 = [number, number, number];

export class OrbitCamera {
  radius: number;
  minRadius: number;
  maxRadius: number;
  azimuth = 0;
  elevation = Math.PI / 3; // start above the grid (z > 0)
  orbitSpeed: number;
  zoomSpeed: number;
  dragging = false;
  moving = false;
  lastX = 0;
  lastY = 0;

  constructor(opts: { radius: number; minRadius: number; maxRadius: number; orbitSpeed: number; zoomSpeed: number; }) {
    this.radius = opts.radius;
    this.minRadius = opts.minRadius;
    this.maxRadius = opts.maxRadius;
    this.orbitSpeed = opts.orbitSpeed;
    this.zoomSpeed = opts.zoomSpeed;
  }

  attach(canvas: HTMLCanvasElement) {
    canvas.addEventListener('mousedown', (e) => {
      if (e.button === 0) {
        this.dragging = true;
        this.moving = true;
        this.lastX = e.clientX; this.lastY = e.clientY;
      }
    });
    window.addEventListener('mouseup', () => {
      this.dragging = false;
    });
    window.addEventListener('mousemove', (e) => {
      if (!this.dragging) return;
      const dx = e.clientX - this.lastX;
      const dy = e.clientY - this.lastY;
      // Typical orbit: azimuth increases dragging right; elevation increases dragging up
      this.azimuth += dx * this.orbitSpeed;
      this.elevation = clamp(this.elevation + dy * this.orbitSpeed, 0.01, Math.PI - 0.01);
      this.lastX = e.clientX; this.lastY = e.clientY;
    });
    canvas.addEventListener('wheel', (e) => {
      e.preventDefault();
      if (e.deltaY < 0) this.radius /= Math.pow(this.zoomSpeed, 1);
      else this.radius *= Math.pow(this.zoomSpeed, 1);
      this.radius = clamp(this.radius, this.minRadius, this.maxRadius);
      this.moving = true;
    }, { passive: false });
  }

  get position(): Vec3 {
    const r = this.radius;
    const el = this.elevation;
    const az = this.azimuth;
    // Z-up: orbit around Z axis
    return [
      r * Math.sin(el) * Math.cos(az),
      r * Math.sin(el) * Math.sin(az),
      r * Math.cos(el),
    ];
  }

  endFrame() {
    // simple inertia flag: if not dragging and no wheel this frame, mark still
    if (!this.dragging) this.moving = false;
  }
}

function clamp(x: number, a: number, b: number) {
  return Math.min(b, Math.max(a, x));
}


