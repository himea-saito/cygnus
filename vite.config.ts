import { defineConfig } from 'vite';

export default defineConfig({
  base: '/cygnus/',
  server: {
    host: true,
    port: 5173,
  },
  assetsInclude: ['**/*.glsl'],
});


