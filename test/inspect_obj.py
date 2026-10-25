import sys
import numpy as np

path = sys.argv[1]
verts = []
faces = 0
quads = 0
with open(path, "r", errors="ignore") as f:
    for line in f:
        if line.startswith("v "):
            p = line.split()
            verts.append((float(p[1]), float(p[2]), float(p[3])))
        elif line.startswith("f "):
            n = len(line.split()) - 1
            faces += 1
            if n > 3:
                quads += 1
V = np.asarray(verts)
print("verts", len(V), "faces", faces, "non-tri", quads)
print("bbox", V.min(axis=0), V.max(axis=0))
print("diag", np.linalg.norm(V.max(axis=0) - V.min(axis=0)))
print("center", (V.max(axis=0) + V.min(axis=0)) / 2)
