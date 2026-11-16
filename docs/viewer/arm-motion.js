import * as THREE from 'three';

const lineColors = ['#22846e', '#78972a', '#4575b7', '#8753a4'];
const up = new THREE.Vector3(0, 1, 0);

export async function startArmMotion(manifest, bodyGroups, wireLayer, initialAction) {
  const url = new URL(manifest.motion_file, import.meta.url);
  url.search = new URL(import.meta.url).search;
  const response = await fetch(url, { cache: 'no-cache' });
  if (!response.ok) throw new Error(`Arm motion: HTTP ${response.status}`);
  const data = await response.json();
  if (data.mesh_sha256 !== manifest.mesh_sha256) throw new Error('Motion belongs to different web geometry');

  const select = document.querySelector('#arm-action');
  const play = document.querySelector('#arm-play');
  const slider = document.querySelector('#arm-time');
  const readout = document.querySelector('#arm-measurement');
  const evidence = document.querySelector('#arm-evidence');
  const geometry = new THREE.CylinderGeometry(1, 1, 1, 8);
  const lines = lineColors.map((color, index) => {
    const material = new THREE.MeshBasicMaterial({ color });
    const length = Math.max(...Object.values(data.actions).map(entry => entry.samples[0][2][index]?.length ?? 1));
    return Array.from({ length: length - 1 }, () => {
      const mesh = new THREE.Mesh(geometry, material);
      wireLayer.add(mesh);
      return mesh;
    });
  });
  const a = new THREE.Vector3(), b = new THREE.Vector3(), direction = new THREE.Vector3();
  const endRotation = new THREE.Quaternion();
  let action, position, running = false, previous = 0;

  for (const [key, entry] of Object.entries(data.actions)) {
    const label = key.startsWith('K') ? `${key} ${entry.label}` : entry.label;
    select.add(new Option(label, key));
  }

  function setRunning(value) {
    running = value;
    previous = 0;
    play.textContent = running ? '暂停' : '播放';
  }

  function renderPose() {
    const frame = Math.min(position, action.samples.length - 1);
    const index = Math.floor(frame);
    const fraction = frame - index;
    const [extent, transforms, paths] = action.samples[index];
    const [nextExtent, nextTransforms, nextPaths] = action.samples[Math.min(index + 1, action.samples.length - 1)];
    bodyGroups.forEach((groups, id) => {
      const [x, y, z, qx, qy, qz, qw] = transforms[id];
      const [nx, ny, nz, nqx, nqy, nqz, nqw] = nextTransforms[id];
      endRotation.set(nqx, nqy, nqz, nqw).normalize();
      for (const group of groups) {
        group.position.set(x, y, z).lerp(a.set(nx, ny, nz), fraction);
        group.quaternion.set(qx, qy, qz, qw).normalize().slerp(endRotation, fraction);
      }
    });
    lines.forEach((segments, id) => {
      const path = paths[id];
      const next = nextPaths[id];
      const axial = id === 2 && ['K4', 'K5'].includes(select.value);
      segments.forEach((mesh, i) => {
        mesh.visible = Boolean(path && i + 1 < path.length);
        if (!mesh.visible) return;
        a.fromArray(path[i]).lerp(direction.fromArray(next[i]), fraction);
        b.fromArray(path[i + 1]).lerp(direction.fromArray(next[i + 1]), fraction);
        direction.subVectors(b, a);
        const distance = direction.length();
        mesh.visible = distance > 1e-7;
        if (!mesh.visible) return;
        const radius = id === 1 ? .0014 : .0012;
        mesh.position.copy(a).lerp(b, .5);
        mesh.quaternion.setFromUnitVectors(up, direction.normalize());
        mesh.scale.set(radius, distance, radius);
        if (id === 2) mesh.material.color.set(axial ? '#be701a' : lineColors[id]);
      });
    });
    slider.value = frame;
    const value = THREE.MathUtils.lerp(extent, nextExtent, fraction);
    readout.value = `${(Math.abs(value) < .005 ? 0 : value).toFixed(2)} ${action.units}`;
  }

  function choose() {
    action = data.actions[select.value];
    position = 0;
    slider.max = action.samples.length - 1;
    evidence.textContent = action.active ? '钢丝牵引与回位 · 理想重力预平衡' : '规定姿态 · 几何往返';
    renderPose();
  }

  select.value = Object.hasOwn(data.actions, initialAction) ? initialAction : 'K1';
  select.addEventListener('change', choose);
  play.addEventListener('click', () => {
    if (!running && position >= Number(slider.max)) position = 0;
    setRunning(!running);
  });
  slider.addEventListener('input', () => {
    setRunning(false);
    position = Number(slider.value);
    renderPose();
  });
  choose();
  return {
    update(time) {
      if (!running) return;
      if (previous) {
        position += Math.min((time - previous) / 1000, .1) / action.sample_interval;
        if (position >= action.samples.length - 1) {
          position = action.samples.length - 1;
          setRunning(false);
          renderPose();
          return;
        }
        renderPose();
      }
      previous = time;
    },
  };
}
