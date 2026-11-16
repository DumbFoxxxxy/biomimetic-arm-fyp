import * as THREE from 'three';

const palette = { reference: '#0085aa', engineering: '#cb6000', candidates: '#a13ace' };
const COLORS = Object.fromEntries(Object.entries(palette).map(([key, value]) =>
  [key, new THREE.MeshBasicMaterial({ color: value, toneMapped: false })]));
const segments = new THREE.CylinderGeometry(1, 1, 1, 8);
const markers = new THREE.SphereGeometry(1, 10, 8);
const up = new THREE.Vector3(0, 1, 0);

function tubes(paths, material, radius) {
  const spans = [];
  for (const points of paths) {
    for (let i = 1; i < points.length; i++) {
      const a = new THREE.Vector3(...points[i - 1]);
      const b = new THREE.Vector3(...points[i]);
      const delta = b.clone().sub(a);
      if (delta.lengthSq() > 1e-14) spans.push({ a, b, delta });
    }
  }
  const mesh = new THREE.InstancedMesh(segments, material, spans.length);
  const pose = new THREE.Object3D();
  spans.forEach(({ a, b, delta }, index) => {
    pose.position.copy(a).add(b).multiplyScalar(.5);
    pose.scale.set(radius, delta.length(), radius);
    pose.quaternion.setFromUnitVectors(up, delta.normalize());
    pose.updateMatrix();
    mesh.setMatrixAt(index, pose.matrix);
  });
  return mesh;
}

function dots(points, material, radius) {
  const mesh = new THREE.InstancedMesh(markers, material, points.length);
  const pose = new THREE.Object3D();
  pose.scale.setScalar(radius);
  points.forEach((point, index) => {
    pose.position.set(...point);
    pose.updateMatrix();
    mesh.setMatrixAt(index, pose.matrix);
  });
  return mesh;
}

export async function startOpenSim(scene, manifest) {
  const url = new URL('./opensim-overlay.json', import.meta.url);
  url.search = new URL(import.meta.url).search;
  const response = await fetch(url, { cache: 'no-cache' });
  if (!response.ok) throw new Error(`OpenSim overlay: HTTP ${response.status}`);
  const data = await response.json();
  if (data.schema_version !== 1 || !data.neutral || data.units !== 'm' ||
      data.assembly_sha256 !== manifest.assembly_sha256 || data.mesh_sha256 !== manifest.mesh_sha256) {
    throw new Error('OpenSim overlay and web mesh versions differ');
  }
  const groups = {
    elbow: ['FJ1486', 'FJ1477'],
    shoulder: ['FJ1468', 'FJ1467', 'FJ1513', 'FJ1507', 'FJ1500', 'FJ1504'],
    scapula: ['FJ1459_g3', 'FJ1536_g2'],
  };
  const notes = {
    elbow: 'K1 刚体与钢丝检查到 31.31°，实物行程待测。',
    shoulder: '冈下肌和肩胛下肌的橙线对应两条局部受拉跨段。',
    scapula: '现有前锯肌线在所查上回旋轨迹中放线。',
  };
  const selector = document.querySelector('#opensim-route');
  const status = document.querySelector('#opensim-status');
  const detail = document.querySelector('#opensim-detail');
  const visibility = Object.fromEntries(Object.keys(COLORS).map(key => [key, new THREE.Group()]));
  const geometry = Object.fromEntries(Object.entries(data.routes).map(([id, route]) => {
    const option = document.createElement('option');
    option.value = id;
    option.textContent = `${route.label} · ${id}`;
    selector.append(option);
    const content = {
      reference: [tubes(route.reference, COLORS.reference, .0010)],
      engineering: [tubes([route.engineering], COLORS.engineering, .0013),
        dots(route.engineering_ports, COLORS.engineering, .0020)],
      candidates: [dots(route.candidates.map(candidate => candidate.xyz), COLORS.candidates, .0030)],
    };
    for (const key of Object.keys(visibility)) visibility[key].add(...content[key]);
    return [id, content];
  }));
  scene.add(...Object.values(visibility));

  function showRoutes() {
    const selected = selector.value;
    const routes = groups[selected] ?? [selected];
    for (const [id, content] of Object.entries(geometry)) {
      for (const meshes of Object.values(content)) meshes.forEach(mesh => { mesh.visible = routes.includes(id); });
    }
    status.textContent = routes.map(id => data.routes[id].label).join('、');
    const route = data.routes[selected];
    detail.textContent = route
      ? `${route.source_model} · ${route.bundles.join('、')}。${route.zone}。`
      : notes[selected];
    detail.hidden = false;
  }

  selector.addEventListener('change', showRoutes);
  selector.disabled = false;
  document.querySelectorAll('[data-compare-layer]').forEach(input => {
    const group = visibility[input.dataset.compareLayer];
    group.visible = input.checked;
    input.addEventListener('change', () => { group.visible = input.checked; });
    input.disabled = false;
  });
  showRoutes();
  return { data, geometry, visibility };
}
