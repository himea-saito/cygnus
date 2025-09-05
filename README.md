# Cygnus – WebGL Black Hole (Schwarzschild) Demo

Cygnus is a WebGL2 (Vite + TypeScript) demo that ray–traces null geodesics around a black hole using the Schwarzschild metric directly in the browser. It aims to replicate the look and spirit of kavan010’s excellent desktop project, but packaged as an easy‑to‑load web page.

- Reference project (CPU/GPU desktop): https://github.com/kavan010/black_hole
- Reference video (overview and results): https://www.youtube.com/watch?v=8-B6ryuBkCM

This project adapts those ideas to a fragment shader + grid overlay, with low‑res pixelated rendering, orbit camera, and a couple of orbiting bodies for frame of reference.

## Live site

- GitHub Pages: https://himea-saito.github.io/cygnus/

If your device is laggy or gets hot: lower the “Render distance” (top‑left), toggle off extra bodies, or zoom in.

## Controls

- Left‑drag: orbit camera
- Mouse wheel: zoom
- G: temporarily boost step count (render distance)
- HUD toggles: show/hide orbiting bodies
- HUD slider: adjust render distance (max iteration count)

## Requirements

- WebGL2 capable browser (Chrome/Edge/Firefox on desktop recommended)
- Node 18+ for local development

## Development

```bash
npm i
npm run dev
```

Open the local URL printed by Vite.

Build a production bundle:

```bash
npm run build
```

The output goes to `/dist/`.

## Deploying to GitHub Pages

This repo includes a GitHub Actions workflow that builds and deploys on pushes to `main`.

- Ensure the `base` path in `vite.config.ts` matches your repository name:

```ts
export default defineConfig({
  base: '/cygnus/',
})
```

- Push to `main`. In GitHub → Settings → Pages, set Build & deployment to “GitHub Actions”. The action will publish `/dist/` to Pages.

If hosting at the account root (username.github.io), set `base: '/'` instead.

## Notes on performance

- The simulation is GPU‑intensive by nature (ray marching per‑pixel). Lower the render distance if performance is poor.
- Mobile devices and iOS Safari may lack robust WebGL2 support; a desktop browser is recommended.

