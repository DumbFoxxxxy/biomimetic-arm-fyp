import * as THREE from 'three';
import { OrbitControls } from './vendor/OrbitControls.js';

const host = document.querySelector('#canvas-host');
const viewport = document.querySelector('#viewport');
const loading = document.querySelector('#loading');
const layers = Object.fromEntries(['bones', 'muscles', 'wires', 'hardware'].map(key => [key, new THREE.Group()]));
const palette = { bones: '#d8cbb1', muscles: '#b84430', wires: '#387c91', hardware: '#78888d' };
let renderer, camera, controls, manifest, bounds;

function applyVisibility(values) {
  for (const [key, visible] of Object.entries(values)) {
    layers[key].visible = visible;
    document.querySelector(`[data-layer="${key}"]`).checked = visible;
  }
}

function fitView() {
  if (!bounds) return;
  const damping = controls.enableDamping;
  controls.enableDamping = false;
  controls.update();
  const size = bounds.getSize(new THREE.Vector3());
  const centre = bounds.getCenter(new THREE.Vector3());
  const vertical = THREE.MathUtils.degToRad(camera.fov / 2);
  const horizontal = Math.atan(Math.tan(vertical) * camera.aspect);
  const distance = size.length() / 2 / Math.sin(Math.min(vertical, horizontal)) * 1.1;
  camera.position.copy(centre).add(new THREE.Vector3(1.2, -2.8, .9).normalize().multiplyScalar(distance));
  camera.near = Math.max(.0001, distance / 1000);
  camera.far = distance * 50;
  camera.updateProjectionMatrix();
  controls.target.copy(centre);
  controls.minDistance = size.length() * .06;
  controls.maxDistance = distance * 8;
  controls.update();
  controls.enableDamping = damping;
}

function showError(error) {
  loading.hidden = false;
  loading.querySelector('.spinner')?.remove();
  loading.querySelector('strong').textContent = '暂时无法载入三维模型';
  loading.querySelector('span:last-child').textContent = location.protocol === 'file:'
    ? '请使用项目说明中的本地 HTTP 服务，或在启用 GitHub Pages 后打开网页。'
    : '请刷新重试，并确认已推送展示资源、使用支持 WebGL 的浏览器。';
  console.error(error);
  window.__armViewerError = String(error);
}

async function start() {
  renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false, powerPreference: 'high-performance' });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 1.75));
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 1.35;
  renderer.domElement.setAttribute('aria-label', '拖动旋转，滚轮缩放，右键平移当前肩部模型');
  renderer.domElement.tabIndex = 0;
  host.appendChild(renderer.domElement);
  const scene = new THREE.Scene();
  scene.background = new THREE.Color('#eeebe3');
  camera = new THREE.PerspectiveCamera(35, 1, .001, 20);
  camera.up.set(0, 0, 1);
  controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = .1;
  controls.autoRotateSpeed = .7;
  controls.screenSpacePanning = true;
  controls.listenToKeyEvents(renderer.domElement);
  scene.add(new THREE.HemisphereLight(0xffffff, 0xb7b3a3, 2.2));
  const key = new THREE.DirectionalLight(0xffffff, 3.4);
  key.position.set(1.5, -2, 3);
  scene.add(key);
  const fill = new THREE.DirectionalLight(0xe0ebf0, 1.7);
  fill.position.set(-2, 1, 1.2);
  scene.add(fill);
  Object.values(layers).forEach(group => scene.add(group));
  const resize = () => {
    const width = host.clientWidth, height = host.clientHeight;
    renderer.setSize(width, height, false);
    camera.aspect = width / Math.max(height, 1);
    camera.updateProjectionMatrix();
  };
  new ResizeObserver(resize).observe(host);
  resize();
  const metaResponse = await fetch(new URL('./scene.json', import.meta.url));
  if (!metaResponse.ok) throw new Error(`Model manifest: HTTP ${metaResponse.status}`);
  manifest = await metaResponse.json();
  const binaryResponse = await fetch(new URL(manifest.mesh_file, import.meta.url));
  if (!binaryResponse.ok) throw new Error(`Model geometry: HTTP ${binaryResponse.status}`);
  const buffer = await binaryResponse.arrayBuffer();
  if (buffer.byteLength !== manifest.mesh_bytes) throw new Error('Model byte count differs from exported manifest');
  for (const part of manifest.parts) {
    const geometry = new THREE.BufferGeometry();
    const vertices = new Float32Array(buffer, part.position_byte_offset, part.vertex_count * 3);
    const Indices = part.index_component_bytes === 2 ? Uint16Array : Uint32Array;
    const indices = new Indices(buffer, part.index_byte_offset, part.triangle_count * 3);
    geometry.setAttribute('position', new THREE.BufferAttribute(vertices, 3));
    geometry.setIndex(new THREE.BufferAttribute(indices, 1));
    geometry.computeVertexNormals();
    const material = new THREE.MeshStandardMaterial({
      color: palette[part.category], roughness: part.category === 'hardware' ? .5 : .76,
      metalness: part.category === 'hardware' ? .18 : 0, side: THREE.DoubleSide,
    });
    const mesh = new THREE.Mesh(geometry, material);
    mesh.name = part.key;
    layers[part.category].add(mesh);
  }
  bounds = new THREE.Box3(new THREE.Vector3(...manifest.bounds_m[0]), new THREE.Vector3(...manifest.bounds_m[1]));
  const size = bounds.getSize(new THREE.Vector3()).length();
  const grid = new THREE.GridHelper(size * 2.8, 22, 0xd1d2c5, 0xdcdcd0);
  grid.rotation.x = Math.PI / 2;
  grid.position.z = bounds.min.z - size * .05;
  grid.material.transparent = true;
  grid.material.opacity = .45;
  scene.add(grid);
  fitView();
  document.querySelector('#scene-count').textContent = `${manifest.part_count} 个对象 · 可自由旋转`;
  document.querySelectorAll('button:disabled, input:disabled').forEach(element => { element.disabled = false; });
  document.querySelectorAll('[data-layer]').forEach(input => input.addEventListener('change', () => {
    layers[input.dataset.layer].visible = input.checked;
    document.querySelectorAll('[data-preset]').forEach(button => button.setAttribute('aria-pressed', 'false'));
  }));
  const presets = {
    assembly: { bones: true, muscles: true, wires: true, hardware: true },
    skeleton: { bones: true, muscles: false, wires: false, hardware: false },
    muscles: { bones: true, muscles: true, wires: true, hardware: false },
  };
  document.querySelectorAll('[data-preset]').forEach(button => button.addEventListener('click', () => {
    applyVisibility(presets[button.dataset.preset]);
    document.querySelectorAll('[data-preset]').forEach(other => other.setAttribute('aria-pressed', String(other === button)));
  }));
  document.querySelector('#reset-view').addEventListener('click', fitView);
  document.querySelector('#auto-rotate').addEventListener('click', event => {
    controls.autoRotate = !controls.autoRotate;
    event.currentTarget.setAttribute('aria-pressed', String(controls.autoRotate));
  });
  document.querySelector('#fullscreen').addEventListener('click', async () => {
    try {
      if (document.fullscreenElement) await document.exitFullscreen();
      else if (viewport.requestFullscreen) await viewport.requestFullscreen();
    } catch { document.querySelector('#fullscreen').textContent = '当前浏览器不支持全屏'; }
  });
  if (!viewport.requestFullscreen) document.querySelector('#fullscreen').disabled = true;
  loading.hidden = true;
  let previous = 0;
  function animate(time) {
    requestAnimationFrame(animate);
    if (document.hidden || time - previous < 30) return;
    previous = time;
    controls.update();
    renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);
  window.__armViewer = { ready: true, renderer, scene, camera, controls, layers, manifest, fitView };
}

start().catch(showError);
