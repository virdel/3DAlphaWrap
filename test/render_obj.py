import sys
import numpy as np
from PIL import Image, ImageDraw

path = sys.argv[1]
out = sys.argv[2]
W = int(sys.argv[3]) if len(sys.argv) > 3 else 1400
H = int(sys.argv[4]) if len(sys.argv) > 4 else 1000

V = []
F = []
with open(path, "r", errors="ignore") as f:
    for line in f:
        if line.startswith("v "):
            p = line.split()
            V.append((float(p[1]), float(p[2]), float(p[3])))
        elif line.startswith("f "):
            idx = [int(t.split("/")[0]) - 1 for t in line.split()[1:]]
            for k in range(1, len(idx) - 1):
                F.append((idx[0], idx[k], idx[k + 1]))
V = np.asarray(V)
F = np.asarray(F)

center = (V.max(axis=0) + V.min(axis=0)) / 2
Vc = V - center
ext = np.abs(Vc).max()
scale = 0.92 * min(W, H) / (2 * ext)

th = np.deg2rad(35)
ph = np.deg2rad(22)
Ry = np.array([[np.cos(th), 0, np.sin(th)], [0, 1, 0], [-np.sin(th), 0, np.cos(th)]])
Rx = np.array([[1, 0, 0], [0, np.cos(ph), -np.sin(ph)], [0, np.sin(ph), np.cos(ph)]])
R = Rx @ Ry
P = Vc @ R.T
P[:, 2] = -P[:, 2]  # camera looks along -z; larger z = closer

sx = W / 2 + scale * P[:, 0]
sy = H / 2 - scale * P[:, 1]
depth = P[:, 2]

A = Vc[F[:, 0]]
B = Vc[F[:, 1]]
C = Vc[F[:, 2]]
N = np.cross(B - A, C - A)
ln = np.linalg.norm(N, axis=1)
ln[ln == 0] = 1
N = N / ln[:, None]
light = np.array([0.35, 0.65, 0.68])
light = light / np.linalg.norm(light)
shade = 0.30 + 0.70 * np.clip(np.abs(N @ light), 0, 1)

order = np.argsort(depth[F].mean(axis=1))
img = Image.new("RGB", (W, H), (245, 245, 242))
dr = ImageDraw.Draw(img)
base = np.array([176, 168, 152])
for ti in order:
    f = F[ti]
    s = shade[ti]
    col = tuple(int(c * (0.55 + 0.45 * s)) for c in base)
    dr.polygon([(sx[f[0]], sy[f[0]]), (sx[f[1]], sy[f[1]]), (sx[f[2]], sy[f[2]])], fill=col)
img.save(out)
print("saved", out, "faces", len(F))
