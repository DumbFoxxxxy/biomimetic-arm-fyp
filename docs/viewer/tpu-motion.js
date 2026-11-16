import * as THREE from 'three';

const supported = new Set(['K1', 'K7']);
const files = new Map();

async function getManifest(hash) {
  const url = new URL('./whole-arm-tpu.json', import.meta.url);
  url.search = new URL(import.meta.url).search;
  const response = await fetch(url, { cache: 'no-cache' });
  if (!response.ok) throw new Error(`TPU manifest HTTP ${response.status}`);
  const data = await response.json();
  if (data.mesh_sha256 !== hash) throw new Error('TPU asset belongs to a different rigid assembly');
  return data;
}

async function getBuffer(chunk) {
  const url = new URL(chunk.file, import.meta.url);
  url.searchParams.set('sha', chunk.sha256);
  const response = await fetch(url);
  if (!response.ok) throw new Error(`TPU geometry HTTP ${response.status}`);
  const buffer = await response.arrayBuffer();
  if (buffer.byteLength !== chunk.bytes) throw new Error('Incomplete TPU geometry');
  return buffer;
}

function buildGeometry(part, buffer) {
  const geometry = new THREE.BufferGeometry();
  const points = part.vertices * 3;
  geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(buffer, part.position_offset, points), 3));
  const Indices = part.index_bytes === 2 ? Uint16Array : Uint32Array;
  geometry.setIndex(new THREE.BufferAttribute(new Indices(buffer, part.index_offset, part.triangles * 3), 1));
  if (part.modes) {
    geometry.morphTargetsRelative = true;
    geometry.morphAttributes.position = part.modes.map(mode => {
      const encoded = new Int16Array(buffer, mode.offset, points);
      return new THREE.BufferAttribute(Float32Array.from(encoded, value => value * mode.scale_m), 3);
    });
  }
  return geometry;
}

function buildMuscle(part, buffer, key) {
  const color = key === 'tpu:elbow_extensor' ? '#a53b49' : '#c35835';
  const mesh = new THREE.Mesh(buildGeometry(part, buffer), new THREE.MeshStandardMaterial({
    color, roughness: .88, flatShading: true, side: THREE.DoubleSide,
  }));
  mesh.name = key;
  mesh.visible = false;
  return mesh;
}

function wireColor(key) {
  return key === 'FJ1486' ? '#117a7a' : key === 'FJ1477' ? '#165daf' : '#406d73';
}

export function createTPUMotion(hash, { muscles: muscleLayer, wires: wireLayer, hardware: hardwareLayer }, paint) {
  const checkbox = document.querySelector('#arm-tpu');
  const note = document.querySelector('#arm-tpu-status');
  const muscles = new Map();
  const wires = new Map();
  const cuffs = new Map();
  const end = new THREE.Vector3();
  const endRotation = new THREE.Quaternion();
  let activeKey, readyKey, data, generation = 0;

  async function choose(key) {
    const request = ++generation;
    activeKey = key;
    readyKey = undefined;
    checkbox.disabled = !supported.has(key);
    if (checkbox.disabled || !checkbox.checked) {
      note.textContent = checkbox.disabled ? '当前动作仅显示刚体与钢丝' : 'TPU 图层已隐藏';
      paint();
      return;
    }
    note.textContent = '正在加载 TPU 晶格…';
    paint();
    try {
      data ??= await getManifest(hash);
      const entries = Object.entries(data.parts).filter(([, part]) => part.coefficients[key]);
      const chunks = new Set(entries.map(([, part]) => part.chunk));
      await Promise.all(Array.from(chunks, async name => {
        if (!files.has(name)) {
          files.set(name, getBuffer(data.chunks[name]).catch(error => {
            files.delete(name);
            throw error;
          }));
        }
        await files.get(name);
      }));
      if (request !== generation) return;
      for (const [id, part] of entries) {
        if (!muscles.has(id)) {
          const mesh = buildMuscle(part, await files.get(part.chunk), id);
          muscleLayer.add(mesh);
          muscles.set(id, mesh);
        }
      }
      if (request !== generation) return;
      if (key === 'K7' && !cuffs.size) {
        const buffer = await files.get('shoulder');
        const material = new THREE.MeshStandardMaterial({ color: '#78888d', roughness: .52, flatShading: true, side: THREE.DoubleSide });
        for (const [id, part] of Object.entries(data.moving_cuffs)) {
          const mesh = new THREE.Mesh(buildGeometry(part, buffer), material);
          mesh.name = `K7 cuff ${id}`;
          mesh.visible = false;
          hardwareLayer.add(mesh);
          cuffs.set(id, mesh);
        }
      }
      if (request !== generation) return;
      readyKey = key;
      note.textContent = key === 'K7' ? '8 块 TPU 几何形变 · K7 既有走线' : '肘屈肌与伸肌 · TPU 几何形变';
      paint();
    } catch (error) {
      if (request === generation) {
        note.textContent = 'TPU 图层加载失败，可取消勾选后重试';
        console.error(error);
      }
    }
  }

  checkbox.addEventListener('change', () => choose(activeKey));

  return {
    choose,
    pose(key, frame) {
      const showing = readyKey === key;
      for (const mesh of muscles.values()) mesh.visible = false;
      for (const wire of wires.values()) wire.visible = false;
      for (const cuff of cuffs.values()) cuff.visible = false;
      if (!showing) return false;
      const { keyframes, indices, wires: samples } = data.actions[key];
      const webIndex = Math.floor(frame);
      const sourceFrame = THREE.MathUtils.lerp(indices[webIndex], indices[Math.min(webIndex + 1, indices.length - 1)], frame - webIndex);
      let i = 0;
      while (i + 1 < keyframes.length - 1 && keyframes[i + 1] < sourceFrame) i++;
      const fraction = (sourceFrame - keyframes[i]) / (keyframes[i + 1] - keyframes[i]);
      for (const [id, mesh] of muscles) {
        const coefficients = data.parts[id].coefficients[key];
        if (!coefficients) continue;
        mesh.visible = true;
        mesh.morphTargetInfluences?.forEach((_, j) => {
          mesh.morphTargetInfluences[j] = THREE.MathUtils.lerp(coefficients[i][j], coefficients[i + 1][j], fraction);
        });
      }
      if (key === 'K7') {
        for (const [id, mesh] of cuffs) {
          const { poses } = data.moving_cuffs[id];
          mesh.visible = true;
          mesh.position.fromArray(poses[i]).lerp(end.fromArray(poses[i + 1]), fraction);
          endRotation.fromArray(poses[i + 1], 3).normalize();
          mesh.quaternion.fromArray(poses[i], 3).normalize().slerp(endRotation, fraction);
        }
      }
      for (const [id, points] of Object.entries(samples[i])) {
        if (!wires.has(id)) {
          const geometry = new THREE.BufferGeometry();
          geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(points.length * 3), 3));
          const line = new THREE.Line(geometry, new THREE.LineBasicMaterial({ color: wireColor(id) }));
          line.frustumCulled = false;
          line.name = `TPU wire ${id}`;
          wireLayer.add(line);
          wires.set(id, line);
        }
        const line = wires.get(id);
        const position = line.geometry.attributes.position;
        points.forEach((point, p) => {
          point.forEach((coordinate, axis) => {
            position.array[p * 3 + axis] = THREE.MathUtils.lerp(coordinate, samples[i + 1][id][p][axis], fraction);
          });
        });
        position.needsUpdate = true;
        line.visible = true;
      }
      return true;
    },
  };
}
