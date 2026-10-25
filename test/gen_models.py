import math
import sys
import numpy as np


def write_obj(path, verts, faces, header=""):
    with open(path, "w") as f:
        if header:
            f.write("# " + header + "\n")
        for v in verts:
            f.write("v %.9g %.9g %.9g\n" % (v[0], v[1], v[2]))
        for t in faces:
            f.write("f %d %d %d\n" % (t[0] + 1, t[1] + 1, t[2] + 1))


def icosphere(subdiv=3, radius=1.0):
    t = (1.0 + math.sqrt(5.0)) / 2.0
    verts = [
        [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
        [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
        [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1],
    ]
    faces = [
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
    ]
    verts = [np.array(v, dtype=float) / np.linalg.norm(v) for v in verts]
    for _ in range(subdiv):
        cache = {}
        new_faces = []

        def midpoint(i, j):
            key = (min(i, j), max(i, j))
            if key in cache:
                return cache[key]
            m = verts[i] + verts[j]
            m = m / np.linalg.norm(m)
            verts.append(m)
            cache[key] = len(verts) - 1
            return cache[key]

        for a, b, c in faces:
            ab = midpoint(a, b)
            bc = midpoint(b, c)
            ca = midpoint(c, a)
            new_faces += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
        faces = new_faces
    verts = [v * radius for v in verts]
    return verts, faces


def torus(R=1.0, r=0.35, nu=60, nv=24):
    verts = []
    for i in range(nu):
        u = 2 * math.pi * i / nu
        for j in range(nv):
            v = 2 * math.pi * j / nv
            verts.append([(R + r * math.cos(v)) * math.cos(u),
                          (R + r * math.cos(v)) * math.sin(u),
                          r * math.sin(v)])
    faces = []
    for i in range(nu):
        for j in range(nv):
            a = i * nv + j
            b = ((i + 1) % nu) * nv + j
            c = ((i + 1) % nu) * nv + (j + 1) % nv
            d = i * nv + (j + 1) % nv
            faces += [[a, b, c], [a, c, d]]
    return verts, faces


def defect_soup(verts, faces, drop=40, flip=40, shift=0.0, seed=0):
    rng = np.random.default_rng(seed)
    faces = [f for f in faces]
    idx = rng.choice(len(faces), size=min(drop, len(faces)), replace=False)
    keep = [f for k, f in enumerate(faces) if k not in set(idx.tolist())]
    vs = [np.array(v, dtype=float) for v in verts]
    if shift > 0:
        for k in rng.choice(len(vs), size=len(vs) // 5, replace=False):
            vs[k] = vs[k] + rng.normal(scale=shift, size=3)
    flip_idx = rng.choice(len(keep), size=min(flip, len(keep)), replace=False)
    for k in flip_idx:
        keep[k] = [keep[k][0], keep[k][2], keep[k][1]]
    return vs, keep



def box(sx=2.0, sy=1.0, sz=0.6):
    x, y, z = sx / 2, sy / 2, sz / 2
    verts = [[-x,-y,-z],[x,-y,-z],[x,y,-z],[-x,y,-z],[-x,-y,z],[x,-y,z],[x,y,z],[-x,y,z]]
    faces = [[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[1,2,6],[1,6,5],
             [2,3,7],[2,7,6],[3,0,4],[3,4,7]]
    return verts, faces


def plane(n=20, size=2.0):
    verts, faces = [], []
    for i in range(n + 1):
        for j in range(n + 1):
            verts.append([size * i / n - size / 2, size * j / n - size / 2, 0.0])
    for i in range(n):
        for j in range(n):
            a = i * (n + 1) + j
            b = (i + 1) * (n + 1) + j
            c = (i + 1) * (n + 1) + j + 1
            d = i * (n + 1) + j + 1
            faces += [[a, b, c], [a, c, d]]
    return verts, faces

def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    v, f = icosphere(3, 1.0)
    write_obj(out + "/sphere.obj", v, f, "icosphere r=1")
    v, f = icosphere(4, 1.0)
    write_obj(out + "/sphere_fine.obj", v, f, "icosphere r=1 subdiv4")
    v, f = box()
    write_obj(out + "/box.obj", v, f, "box")
    v, f = plane()
    write_obj(out + "/plane.obj", v, f, "open plane")
    v, f = torus()
    write_obj(out + "/torus.obj", v, f, "torus")
    vs, fs = defect_soup(*icosphere(3, 1.0), drop=60, flip=40, shift=0.004, seed=7)
    write_obj(out + "/sphere_defect.obj", vs, fs, "defective sphere soup")
    print("written to", out)


if __name__ == "__main__":
    main()
