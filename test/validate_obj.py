"""Validate an OBJ wrap: manifoldness, orientation, volume, offset distance, enclosure."""
import sys
import numpy as np


def load_obj(path):
    verts, faces = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                verts.append([float(x) for x in line.split()[1:4]])
            elif line.startswith("f "):
                idx = [int(t.split("/")[0]) - 1 for t in line.split()[1:]]
                for k in range(1, len(idx) - 1):
                    faces.append([idx[0], idx[k], idx[k + 1]])
    return np.asarray(verts, dtype=np.float64), np.asarray(faces, dtype=np.int64)


def closest_point_to_tris(p, A, B, C):
    """Vectorized closest point on triangles (Ericson), p (3,), A/B/C (F,3)."""
    ab = B - A
    ac = C - A
    ap = p - A
    d1 = np.einsum("ij,ij->i", ab, ap)
    d2 = np.einsum("ij,ij->i", ac, ap)
    bp = p - B
    d3 = np.einsum("ij,ij->i", ab, bp)
    d4 = np.einsum("ij,ij->i", ac, bp)
    cp = p - C
    d5 = np.einsum("ij,ij->i", ab, cp)
    d6 = np.einsum("ij,ij->i", ac, cp)
    vc = d1 * d4 - d3 * d2
    vb = d5 * d2 - d1 * d6
    va = d3 * d6 - d5 * d4
    den = va + vb + vc
    den = np.where(np.abs(den) < 1e-30, 1.0, den)

    res = np.empty_like(A)
    v = np.clip(vb / den, 0, 1)
    w = np.clip(vc / den, 0, 1)
    res = A + ab * v[:, None] + ac * w[:, None]

    # vertex / edge regions
    m = (d1 <= 0) & (d2 <= 0)
    res[m] = A[m]
    m = (d3 >= 0) & (d4 <= d3)
    res[m] = B[m]
    m = (d6 >= 0) & (d5 <= d6)
    res[m] = C[m]
    m = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    t = (d1 / np.where(d1 - d3 == 0, 1, d1 - d3))[m]
    res[m] = A[m] + ab[m] * t[:, None]
    m = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    t = (d2 / np.where(d2 - d6 == 0, 1, d2 - d6))[m]
    res[m] = A[m] + ac[m] * t[:, None]
    m = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
    t = ((d4 - d3) / np.where((d4 - d3) + (d5 - d6) == 0, 1, (d4 - d3) + (d5 - d6)))[m]
    res[m] = B[m] + (C[m] - B[m]) * t[:, None]
    return res


def winding_numbers(points, A, B, C):
    out = np.zeros(len(points))
    for i, p in enumerate(points):
        a = A - p
        b = B - p
        c = C - p
        la = np.linalg.norm(a, axis=1)
        lb = np.linalg.norm(b, axis=1)
        lc = np.linalg.norm(c, axis=1)
        det = np.einsum("ij,ij->i", a, np.cross(b, c))
        den = la * lb * lc + np.einsum("ij,ij->i", a, b) * lc + np.einsum("ij,ij->i", b, c) * la + np.einsum("ij,ij->i", c, a) * lb
        ang = 2.0 * np.arctan2(det, den)
        out[i] = ang.sum() / (4.0 * np.pi)
    return out


def main():
    inp, out = sys.argv[1], sys.argv[2]
    V, F = load_obj(out)
    IV, IF = load_obj(inp)
    print(f"output: {len(V)} verts, {len(F)} faces")

    edges = {}
    for f in F:
        for k in range(3):
            a, b = int(f[k]), int(f[(k + 1) % 3])
            key = (min(a, b), max(a, b))
            edges.setdefault(key, []).append((a, b))
    cnt = {}
    for k, v in edges.items():
        cnt[len(v)] = cnt.get(len(v), 0) + 1
    odd = sum(v for k, v in cnt.items() if k % 2 == 1)
    print("edge incidence histogram:", dict(sorted(cnt.items())), " odd-incidence edges:", odd)

    # orientation: shared edges must be traversed in opposite directions
    bad_orient = 0
    for k, v in edges.items():
        if len(v) == 2 and v[0] == v[1]:
            bad_orient += 1
    print("orientation conflicts (2-face):", bad_orient)

    vol = np.einsum("ij,ij->i", V[F[:, 0]], np.cross(V[F[:, 1]], V[F[:, 2]])).sum() / 6.0
    print(f"signed volume: {vol:.6g}")

    # vertex distances to input (chunked)
    A = IV[IF[:, 0]]
    B = IV[IF[:, 1]]
    C = IV[IF[:, 2]]
    dmin, dmax, dsum = 1e30, 0.0, 0.0
    vsample = V if len(V) <= 1500 else V[:: max(1, len(V) // 1500)]
    for s in range(0, len(vsample), 32):
        chunk = vsample[s:s + 32]
        for p in chunk:
            q = closest_point_to_tris(p, A, B, C)
            d = np.linalg.norm(q - p, axis=1).min()
            dmin = min(dmin, d)
            dmax = max(dmax, d)
            dsum += d
    print(f"vertex distance ({len(vsample)} sampled): min {dmin:.6g}  mean {dsum/len(vsample):.6g}  max {dmax:.6g}")

    # enclosure: winding number of input vertices w.r.t. the output mesh
    OV = V
    OF = F
    sample = IV[:: max(1, len(IV) // 80)]
    wn = winding_numbers(sample, OV[OF[:, 0]], OV[OF[:, 1]], OV[OF[:, 2]])
    print(f"winding number of {len(sample)} input vertices: min {wn.min():.4f}  mean {wn.mean():.4f}  max {wn.max():.4f}")


if __name__ == "__main__":
    main()
