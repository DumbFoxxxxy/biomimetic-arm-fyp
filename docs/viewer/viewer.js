import * as THREE from 'three';
import { OrbitControls } from './vendor/OrbitControls.js';

const host = document.querySelector('#canvas-host');
const viewport = document.querySelector('#viewport');
const loading = document.querySelector('#loading');
const layers = Object.fromEntries(['bones', 'muscles', 'wires', 'hardware'].map(key => [key, new THREE.Group()]));
const palette = { bones: '#d8cbb1', muscles: '#b84430', wires: '#387c91', hardware: '#78888d' };
const parameters = new URLSearchParams(location.search);
const selected = parameters.get('model');
const combined = selected === 'combined';
const firstPass = selected === 'first-pass';
const motionMode = selected === 'motion';
const comparison = firstPass && parameters.get('compare') === 'opensim';
document.body.classList.toggle('opensim-mode', comparison);
document.body.classList.toggle('motion-mode', motionMode);
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
  const viewpoint = combined || firstPass || motionMode ? new THREE.Vector3(-1.8, -2.8, .8) : new THREE.Vector3(1.2, -2.8, .9);
  camera.position.copy(centre).add(viewpoint.normalize().multiplyScalar(distance));
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
  const shoulder = selected === 'shoulder';
  if (motionMode) {
    document.querySelector('#model-name').textContent = '整臂运动回放';
    document.querySelector('#scene-label').textContent = '整臂 · MuJoCo 轨迹';
    document.querySelector('#elbow-view').removeAttribute('aria-current');
    document.querySelector('#motion-view').setAttribute('aria-current', 'page');
    document.querySelector('#wire-label').textContent = '肩胛钢丝与导向';
    document.querySelector('#model-scope').textContent = 'K1 与 K7 显示几何插值生成的 TPU 形变。接触压力、材料应变和实际行程需另行验证。';
    document.querySelector('#design-description').textContent = '在播放中拖动模型可从任意方向查看动作。K1 可查看肘部拮抗晶格，K7 可查看八块肌腹。K6 与肩胛后缩采用钢丝牵引。内外旋的橙线为肱骨方向标记。';
    document.querySelector('#design-record').href = './showcase/#all-actions';
    document.querySelector('#design-record').textContent = '查看动作范围和验证记录 ↗';
    document.querySelector('[data-layer="muscles"]').closest('label').hidden = true;
    document.querySelector('.presets').hidden = true;
  }
  if (shoulder || combined || firstPass) {
    document.querySelector('#model-name').textContent = firstPass ? '桌面整臂首版' : combined ? '胸廓 · S1 肩部 · 肘部' : 'S1-Y2 肩部';
    document.querySelector('#scene-label').textContent = firstPass ? 'FIRST-PASS WHOLE ARM' : combined ? 'SHOULDER–ELBOW ASSEMBLY' : 'S1 SHOULDER';
    document.querySelector('#elbow-view').removeAttribute('aria-current');
    document.querySelector(firstPass ? '#first-pass-view' : combined ? '#combined-view' : '#shoulder-view').setAttribute('aria-current', 'page');
  }
  if (firstPass) {
    document.querySelector('#design-description').textContent = '从桌面底座、胸廓到肩胛、S1 肩部和肘部，观察八块 TPU 晶格、回位弹簧与八条固定端钢丝。肩袖活动件和中央肌腱见 K7 动画。';
    document.querySelector('#wire-label').textContent = '八条局部钢丝';
    document.querySelector('#design-record').href = 'https://github.com/DumbFoxxxxy/biomimetic-arm-fyp/tree/main/02_人体骨骼晶格肌肉方案/输出/whole_arm_first_pass';
    document.querySelector('#design-record').textContent = '查看首版装配与动作范围 ↗';
  } else if (combined) {
    document.querySelector('#design-description').textContent = '查看共用肱骨的 S1 肩部、前臂、串联滑块肘部和胸廓参考。肌肉按静置形态显示，部分肩部端部接口仍在整合。';
    document.querySelector('#wire-label').textContent = '肘部局部钢丝';
    document.querySelector('#design-record').href = 'https://github.com/DumbFoxxxxy/biomimetic-arm-fyp/tree/main/02_人体骨骼晶格肌肉方案/输出/integration_review/s1_series_elbow';
    document.querySelector('#design-record').textContent = '查看组合装配与验证记录 ↗';
  } else if (shoulder) {
    document.querySelector('#design-description').textContent = '查看 S1-Y2 的肩胛、锁骨、盂肱关节及 TPU 晶格肌肉。';
  }
  renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false, powerPreference: 'high-performance' });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 1.75));
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 1.35;
  renderer.domElement.setAttribute('aria-label', '拖动旋转，滚轮缩放，右键平移当前模型');
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
  const modelFile = motionMode ? 'whole-arm-actions.json' : firstPass ? 'first-pass.json'
    : combined ? 'combined.json' : shoulder ? 'shoulder.json' : 'scene.json';
  const metaUrl = new URL(`./${modelFile}`, import.meta.url);
  metaUrl.search = new URL(import.meta.url).search;
  const metaResponse = await fetch(metaUrl, { cache: 'no-cache' });
  if (!metaResponse.ok) throw new Error(`Model manifest: HTTP ${metaResponse.status}`);
  manifest = await metaResponse.json();
  const binaryUrl = new URL(manifest.mesh_file, import.meta.url);
  binaryUrl.searchParams.set('sha', manifest.mesh_sha256);
  const binaryResponse = await fetch(binaryUrl);
  if (!binaryResponse.ok) throw new Error(`Model geometry: HTTP ${binaryResponse.status}`);
  const buffer = await binaryResponse.arrayBuffer();
  if (buffer.byteLength !== manifest.mesh_bytes) throw new Error('Model byte count differs from exported manifest');
  const bodyGroups = motionMode ? Array.from({ length: Math.max(...manifest.part_body) + 1 }, () => new Map()) : [];
  for (const [partIndex, part] of manifest.parts.entries()) {
    const geometry = new THREE.BufferGeometry();
    const vertices = new Float32Array(buffer, part.position_byte_offset, part.vertex_count * 3);
    const Indices = part.index_component_bytes === 2 ? Uint16Array : Uint32Array;
    const indices = new Indices(buffer, part.index_byte_offset, part.triangle_count * 3);
    geometry.setAttribute('position', new THREE.BufferAttribute(vertices, 3));
    geometry.setIndex(new THREE.BufferAttribute(indices, 1));
    geometry.computeVertexNormals();
    const reference = combined && part.key.startsWith('thorax_reference:');
    const material = new THREE.MeshStandardMaterial({
      color: reference ? '#a2b6bd' : palette[part.category], roughness: part.category === 'hardware' ? .5 : .76,
      metalness: part.category === 'hardware' ? .18 : 0, side: THREE.DoubleSide,
      transparent: reference, opacity: reference ? .29 : 1, depthWrite: !reference,
    });
    const mesh = new THREE.Mesh(geometry, material);
    mesh.name = part.key;
    let parent = layers[part.category];
    if (motionMode) {
      const groups = bodyGroups[manifest.part_body[partIndex]];
      if (!groups.has(part.category)) {
        const group = new THREE.Group();
        parent.add(group);
        groups.set(part.category, group);
      }
      parent = groups.get(part.category);
    }
    parent.add(mesh);
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
  document.querySelectorAll('button:disabled, input:disabled:not([data-compare-layer])').forEach(element => { element.disabled = false; });
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
  let motion;
  if (motionMode) {
    const moduleUrl = new URL('./arm-motion.js', import.meta.url);
    moduleUrl.search = new URL(import.meta.url).search;
    try {
      const { startArmMotion } = await import(moduleUrl);
      motion = await startArmMotion(manifest, bodyGroups.map(groups => [...groups.values()]), layers, parameters.get('action'));
      document.querySelector('#arm-motion-controls').hidden = false;
    } catch (error) {
      document.querySelector('#scene-count').textContent = '动作数据加载失败 · 可旋转查看中立位';
      window.__armViewerMotionError = String(error);
    }
  }
  loading.hidden = true;
  let previous = 0;
  function animate(time) {
    requestAnimationFrame(animate);
    if (document.hidden || time - previous < 30) return;
    previous = time;
    motion?.update(time);
    controls.update();
    renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);
  window.__armViewer = { ready: true, renderer, scene, camera, controls, layers, manifest, fitView };
  if (comparison) {
    const panel = document.querySelector('#opensim-controls');
    panel.hidden = false;
    document.querySelector('#first-pass-view').removeAttribute('aria-current');
    document.querySelector('#opensim-view').setAttribute('aria-current', 'page');
    document.querySelector('#scene-label').textContent = 'OPENSIM · NEUTRAL';
    const moduleUrl = new URL('./opensim.js', import.meta.url);
    moduleUrl.search = new URL(import.meta.url).search;
    try {
      const { startOpenSim } = await import(moduleUrl);
      window.__armViewer.comparison = await startOpenSim(scene, manifest);
      applyVisibility(presets.skeleton);
      document.querySelectorAll('[data-preset]').forEach(button => {
        button.setAttribute('aria-pressed', String(button.dataset.preset === 'skeleton'));
      });
    } catch (error) {
      document.querySelector('#opensim-status').textContent = '暂时无法载入对照数据。可刷新重试，或使用下方装配图层。';
      window.__armViewer.comparisonError = String(error);
    }
  }
}

start().catch(showError);
