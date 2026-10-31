"""Report basic geometry stats for the BodyParts3D OBJ meshes."""

import os
import struct


def read_obj(path):
    verts = []
    faces = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith("v "):
                parts = line.split()
                verts.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif line.startswith("f "):
                idx = []
                for tok in line.split()[1:]:
                    idx.append(int(tok.split("/")[0]) - 1)
                for k in range(1, len(idx) - 1):
                    faces.append((idx[0], idx[k], idx[k + 1]))
    return verts, faces


def edge_stats(faces):
    counts = {}
    for a, b, c in faces:
        for e in ((a, b), (b, c), (c, a)):
            key = (min(e), max(e))
            counts[key] = counts.get(key, 0) + 1
    boundary = sum(1 for v in counts.values() if v == 1)
    nonmanifold = sum(1 for v in counts.values() if v > 2)
    return len(counts), boundary, nonmanifold


def main():
    folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    print(f"{'file':<14}{'verts':>8}{'tris':>8}{'edges':>8}{'boundary':>10}{'nonman':>8}   dims (mm)")
    for name in sorted(os.listdir(folder)):
        if not name.lower().endswith(".obj"):
            continue
        verts, faces = read_obj(os.path.join(folder, name))
        xs = [v[0] for v in verts]
        ys = [v[1] for v in verts]
        zs = [v[2] for v in verts]
        dims = (max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))
        ne, bnd, nonman = edge_stats(faces)
        dims_txt = " x ".join(f"{d:8.2f}" for d in dims)
        print(f"{name:<14}{len(verts):>8}{len(faces):>8}{ne:>8}{bnd:>10}{nonman:>8}   {dims_txt}")


if __name__ == "__main__":
    main()
