// alpha_wrap.cpp
// Standalone C++17 implementation of "Alpha Wrapping with an Offset"
// (Portaneri et al., ACM TOG 2022), no external libraries.
//
// Pipeline:
//   OBJ triangle soup -> AABB tree (distance / intersection oracle)
//   -> 3D Delaunay tetrahedralization (incremental Bowyer-Watson)
//   -> outside-in carving with alpha-empty-ball gates
//   -> Steiner points on the offset surface (rules R1/R2)
//   -> watertight triangle mesh -> OBJ
//
// Build:  g++ -O2 -std=c++17 alpha_wrap.cpp -o alpha_wrap
// Usage:  alpha_wrap in.obj out.obj [--alpha A] [--offset D] [--verbose]
//
// Limitations (documented in README): predicates use long double with a
// relative filter (no exact arithmetic); near-degenerate input can be
// perturbed by the insertion fallback. Single threaded.

#if defined(__MINGW32__) && !defined(__USE_MINGW_ANSI_STDIO)
#define __USE_MINGW_ANSI_STDIO 1
#endif

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// ----------------------------------------------------------------------------
// Basic geometry
// ----------------------------------------------------------------------------
struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};
static inline Vec3 operator+(const Vec3& a, const Vec3& b) { return Vec3(a.x + b.x, a.y + b.y, a.z + b.z); }
static inline Vec3 operator-(const Vec3& a, const Vec3& b) { return Vec3(a.x - b.x, a.y - b.y, a.z - b.z); }
static inline Vec3 operator*(double s, const Vec3& a) { return Vec3(s * a.x, s * a.y, s * a.z); }
static inline Vec3 operator/(const Vec3& a, double s) { return Vec3(a.x / s, a.y / s, a.z / s); }
static inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
static inline double norm2(const Vec3& a) { return dot(a, a); }
static inline double norm(const Vec3& a) { return std::sqrt(norm2(a)); }

// ----------------------------------------------------------------------------
// Robust-ish predicates: long double + static relative filter.
// ----------------------------------------------------------------------------
#ifdef _MSC_VER
static const long double REL_EPS = 1e-12L;
#else
static const long double REL_EPS = 1e-15L;
#endif

static inline int signFilter(long double det, long double mag) {
    if (det > REL_EPS * mag) return 1;
    if (det < -REL_EPS * mag) return -1;
    return 0;
}

// det[b-a, c-a, d-a]
static inline long double orient3dVal(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    long double ax = (long double)b.x - a.x, ay = (long double)b.y - a.y, az = (long double)b.z - a.z;
    long double bx = (long double)c.x - a.x, by = (long double)c.y - a.y, bz = (long double)c.z - a.z;
    long double cx = (long double)d.x - a.x, cy = (long double)d.y - a.y, cz = (long double)d.z - a.z;
    long double t1 = ax * (by * cz - bz * cy);
    long double t2 = ay * (bx * cz - bz * cx);
    long double t3 = az * (bx * cy - by * cx);
    return t1 - t2 + t3;
}

static inline int orientSign(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    long double ax = (long double)b.x - a.x, ay = (long double)b.y - a.y, az = (long double)b.z - a.z;
    long double bx = (long double)c.x - a.x, by = (long double)c.y - a.y, bz = (long double)c.z - a.z;
    long double cx = (long double)d.x - a.x, cy = (long double)d.y - a.y, cz = (long double)d.z - a.z;
    long double t1 = ax * (by * cz - bz * cy);
    long double t2 = ay * (bx * cz - bz * cx);
    long double t3 = az * (bx * cy - by * cx);
    long double det = t1 - t2 + t3;
    long double mag = fabsl(t1) + fabsl(t2) + fabsl(t3);
    return signFilter(det, mag);
}

static inline long double det3(long double m[3][3], int c0, int c1, int c2) {
    return m[0][c0] * (m[1][c1] * m[2][c2] - m[1][c2] * m[2][c1])
         - m[0][c1] * (m[1][c0] * m[2][c2] - m[1][c2] * m[2][c0])
         + m[0][c2] * (m[1][c0] * m[2][c1] - m[1][c1] * m[2][c0]);
}

// det of 4x4 matrix (rows given). Also accumulates |terms| for filtering.
static inline long double det4(long double m[4][4], long double& mag) {
    long double m3[3][3];
    long double total = 0;
    long double acc = 0;
    for (int c = 0; c < 4; ++c) {
        for (int r = 1; r < 4; ++r)
            for (int k = 0, cc = 0; k < 4; ++k) {
                if (k == c) continue;
                m3[r - 1][cc++] = m[r][k];
            }
        long double minor = det3(m3, 0, 1, 2);
        long double term = m[0][c] * minor;
        long double s = (c % 2 == 0) ? 1.0L : -1.0L;
        total += s * term;
        // crude magnitude: |m0c| * sum of |3x3 terms|
        long double sub = 0;
        sub += fabsl(m3[0][0] * (m3[1][1] * m3[2][2] - m3[1][2] * m3[2][1]));
        sub += fabsl(m3[0][1] * (m3[1][0] * m3[2][2] - m3[1][2] * m3[2][0]));
        sub += fabsl(m3[0][2] * (m3[1][0] * m3[2][1] - m3[1][1] * m3[2][0]));
        acc += fabsl(m[0][c]) * sub;
    }
    mag = acc;
    return total;
}

// Returns +1 if p is inside the circumsphere of positively oriented (a,b,c,d).
// Convention verified by selfTest(): inside  =>  det4 < 0.
static inline int inSphereSign(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d, const Vec3& p) {
    long double m[4][4];
    const Vec3* q[4] = {&a, &b, &c, &d};
    for (int i = 0; i < 4; ++i) {
        m[i][0] = (long double)q[i]->x - p.x;
        m[i][1] = (long double)q[i]->y - p.y;
        m[i][2] = (long double)q[i]->z - p.z;
        m[i][3] = m[i][0] * m[i][0] + m[i][1] * m[i][1] + m[i][2] * m[i][2];
    }
    long double mag = 0;
    long double det = det4(m, mag);
    int s = signFilter(det, mag);
    int o = orientSign(a, b, c, d);
    if (o == 0) return 0;
    return (o > 0) ? -s : s;
}

// ----------------------------------------------------------------------------
// Triangle soup + AABB tree
// ----------------------------------------------------------------------------
struct Tri {
    Vec3 a, b, c;
    Vec3 lo, hi;
};

static std::vector<Tri> g_tris;

struct BVHNode {
    Vec3 lo, hi;
    int left = -1, right = -1;
    int start = 0, count = 0;
};
static std::vector<BVHNode> g_nodes;
static std::vector<int> g_order;

static inline void triBBox(const Tri& t, Vec3& lo, Vec3& hi) {
    lo.x = std::min({t.a.x, t.b.x, t.c.x});
    lo.y = std::min({t.a.y, t.b.y, t.c.y});
    lo.z = std::min({t.a.z, t.b.z, t.c.z});
    hi.x = std::max({t.a.x, t.b.x, t.c.x});
    hi.y = std::max({t.a.y, t.b.y, t.c.y});
    hi.z = std::max({t.a.z, t.b.z, t.c.z});
}

static inline double boxDist2(const Vec3& lo, const Vec3& hi, const Vec3& p) {
    double dx = std::max({lo.x - p.x, 0.0, p.x - hi.x});
    double dy = std::max({lo.y - p.y, 0.0, p.y - hi.y});
    double dz = std::max({lo.z - p.z, 0.0, p.z - hi.z});
    return dx * dx + dy * dy + dz * dz;
}

static inline bool boxOverlap(const Vec3& a0, const Vec3& a1, const Vec3& b0, const Vec3& b1) {
    return a0.x <= b1.x && b0.x <= a1.x && a0.y <= b1.y && b0.y <= a1.y && a0.z <= b1.z && b0.z <= a1.z;
}

static Vec3 closestOnTri(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    Vec3 ab = b - a, ac = c - a, ap = p - a;
    double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) return a;
    Vec3 bp = p - b;
    double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return b;
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        double v = d1 / (d1 - d3);
        return a + v * ab;
    }
    Vec3 cp = p - c;
    double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return c;
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        double w = d2 / (d2 - d6);
        return a + w * ac;
    }
    double va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + w * (c - b);
    }
    double den = 1.0 / (va + vb + vc);
    double v = vb * den, w = vc * den;
    return a + v * ab + w * ac;
}

static int buildBVH(int start, int end) {
    int id = (int)g_nodes.size();
    g_nodes.push_back(BVHNode());
    Vec3 lo(1e300, 1e300, 1e300), hi(-1e300, -1e300, -1e300);
    for (int i = start; i < end; ++i) {
        const Tri& t = g_tris[g_order[i]];
        lo.x = std::min(lo.x, t.lo.x); lo.y = std::min(lo.y, t.lo.y); lo.z = std::min(lo.z, t.lo.z);
        hi.x = std::max(hi.x, t.hi.x); hi.y = std::max(hi.y, t.hi.y); hi.z = std::max(hi.z, t.hi.z);
    }
    g_nodes[id].lo = lo;
    g_nodes[id].hi = hi;
    int n = end - start;
    if (n <= 4) {
        g_nodes[id].start = start;
        g_nodes[id].count = n;
        return id;
    }
    Vec3 ext = hi - lo;
    int axis = (ext.x > ext.y && ext.x > ext.z) ? 0 : (ext.y > ext.z ? 1 : 2);
    int mid = start + n / 2;
    std::nth_element(g_order.begin() + start, g_order.begin() + mid, g_order.begin() + end,
                     [axis](int i, int j) {
                         const Tri& a = g_tris[i];
                         const Tri& b = g_tris[j];
                         double ca = axis == 0 ? (a.lo.x + a.hi.x) : (axis == 1 ? (a.lo.y + a.hi.y) : (a.lo.z + a.hi.z));
                         double cb = axis == 0 ? (b.lo.x + b.hi.x) : (axis == 1 ? (b.lo.y + b.hi.y) : (b.lo.z + b.hi.z));
                         return ca < cb;
                     });
    int l = buildBVH(start, mid);
    int r = buildBVH(mid, end);
    g_nodes[id].left = l;
    g_nodes[id].right = r;
    g_nodes[id].count = 0;
    return id;
}

static void closestRec(int node, const Vec3& p, double& best2, Vec3& bestP) {
    const BVHNode& nd = g_nodes[node];
    if (boxDist2(nd.lo, nd.hi, p) >= best2) return;
    if (nd.count > 0) {
        for (int i = nd.start; i < nd.start + nd.count; ++i) {
            const Tri& t = g_tris[g_order[i]];
            Vec3 q = closestOnTri(p, t.a, t.b, t.c);
            double d2 = norm2(q - p);
            if (d2 < best2) { best2 = d2; bestP = q; }
        }
        return;
    }
    closestRec(nd.left, p, best2, bestP);
    closestRec(nd.right, p, best2, bestP);
}

static double distToInput(const Vec3& p) {
    double best2 = 1e300;
    Vec3 bestP;
    closestRec(0, p, best2, bestP);
    return std::sqrt(best2);
}

static Vec3 closestToInput(const Vec3& p) {
    double best2 = 1e300;
    Vec3 bestP = p;
    closestRec(0, p, best2, bestP);
    return bestP;
}

// ----------------------------------------------------------------------------
// Triangle / tetrahedron intersection (conservative)
// ----------------------------------------------------------------------------
static bool segTri(const Vec3& p0, const Vec3& p1, const Vec3& a, const Vec3& b, const Vec3& c, double eps) {
    Vec3 dir = p1 - p0;
    Vec3 e1 = b - a, e2 = c - a;
    Vec3 h = cross(dir, e2);
    double det = dot(e1, h);
    if (std::fabs(det) < 1e-18) return false;  // parallel: ignored (conservative enough in practice)
    double inv = 1.0 / det;
    Vec3 s = p0 - a;
    double u = inv * dot(s, h);
    if (u < -eps || u > 1.0 + eps) return false;
    Vec3 q = cross(s, e1);
    double v = inv * dot(dir, q);
    if (v < -eps || u + v > 1.0 + eps) return false;
    double t = inv * dot(e2, q);
    return t >= -eps && t <= 1.0 + eps;
}

// ----------------------------------------------------------------------------
// Tetrahedral mesh
// ----------------------------------------------------------------------------
struct Tet {
    int v[4] = {0, 0, 0, 0};
    int nb[4] = {-1, -1, -1, -1};
    int nbi[4] = {-1, -1, -1, -1};
    unsigned char state = 0;  // 0 inside, 1 outside
    bool alive = true;
};

static std::vector<Vec3> g_V;
static std::vector<Tet> g_T;
static double g_diag = 1.0;
static double g_epsLen = 1e-9;

static inline Vec3 P(int i) { return g_V[i]; }

// Outward-oriented triangle of face i (opposite vertex i) of a positively
// oriented tet, validated in selfTest().
static const int FACE[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};

static inline void faceVerts(const Tet& t, int f, Vec3& a, Vec3& b, Vec3& c) {
    a = g_V[t.v[FACE[f][0]]];
    b = g_V[t.v[FACE[f][1]]];
    c = g_V[t.v[FACE[f][2]]];
}

static int makeTet(int a, int b, int c, int d) {
    if (orientSign(g_V[a], g_V[b], g_V[c], g_V[d]) < 0) std::swap(b, c);
    Tet t;
    t.v[0] = a; t.v[1] = b; t.v[2] = c; t.v[3] = d;
    g_T.push_back(t);
    return (int)g_T.size() - 1;
}

static inline uint64_t fkey(int a, int b, int c) {
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
    return ((uint64_t)a << 42) | ((uint64_t)b << 21) | (uint64_t)c;
}

static Vec3 triCircumcenter(const Vec3& a, const Vec3& b, const Vec3& c, double* R) {
    Vec3 A = b - a, B = c - a;
    double aa = norm2(A), bb = norm2(B), ab = dot(A, B);
    double det = aa * bb - ab * ab;
    if (std::fabs(det) < 1e-300) {
        if (R) *R = 1e300;
        return (a + b + c) / 3.0;
    }
    double rhs1 = aa * 0.5, rhs2 = bb * 0.5;
    double al = (rhs1 * bb - ab * rhs2) / det;
    double be = (aa * rhs2 - ab * rhs1) / det;
    Vec3 ctr = a + al * A + be * B;
    if (R) *R = norm(ctr - a);
    return ctr;
}

static inline double det3v(const Vec3& u, const Vec3& v, const Vec3& w) {
    return dot(u, cross(v, w));
}

static Vec3 tetCircumcenter(int ti, double* R) {
    const Tet& t = g_T[ti];
    const Vec3& a = g_V[t.v[0]];
    const Vec3& b = g_V[t.v[1]];
    const Vec3& c = g_V[t.v[2]];
    const Vec3& d = g_V[t.v[3]];
    Vec3 A = b - a, B = c - a, C = d - a;
    double det = det3v(A, B, C);
    if (std::fabs(det) < 1e-300) {
        if (R) *R = 1e300;
        return (a + b + c + d) / 4.0;
    }
    Vec3 rhs(0.5 * norm2(A), 0.5 * norm2(B), 0.5 * norm2(C));
    Vec3 x = (rhs.x * cross(B, C) + rhs.y * cross(C, A) + rhs.z * cross(A, B)) / det;
    Vec3 ctr = a + x;
    if (R) *R = norm(x);
    return ctr;
}

// ----------------------------------------------------------------------------
// Point location
// ----------------------------------------------------------------------------
static int locateCell(const Vec3& p, int hint) {
    int c = hint;
    if (c < 0 || c >= (int)g_T.size() || !g_T[c].alive) {
        c = -1;
        for (int i = 0; i < (int)g_T.size(); ++i) {
            if (g_T[i].alive) { c = i; break; }
        }
        if (c < 0) return -1;
    }
    for (int it = 0; it < 200000; ++it) {
        const Tet& t = g_T[c];
        bool moved = false;
        for (int i = 0; i < 4; ++i) {
            Vec3 a, b, cc;
            faceVerts(t, i, a, b, cc);
            if (orientSign(a, b, cc, p) > 0) {
                int nb = t.nb[i];
                if (nb < 0) return -1;
                c = nb;
                moved = true;
                break;
            }
        }
        if (!moved) return c;
    }
    return -2;
}

// ----------------------------------------------------------------------------
// Incremental Delaunay insertion (Bowyer-Watson)
// ----------------------------------------------------------------------------
struct BoundaryFace {
    int cell, face, nb, nbf;
};

static std::vector<int> g_stamp;
static int g_curStamp = 0;
static double totalAliveVolume();
static long long g_failLoc = 0, g_failDup = 0, g_failStar = 0, g_failLink = 0, g_failCav = 0;
static bool g_verbose = false;
static bool g_debug = false;
static long long g_insCalls = 0;

static bool insertPoint(const Vec3& p0, int& hint, std::vector<int>& newCells) {
    newCells.clear();
    ++g_insCalls;
    if (g_debug && g_insCalls == 3) {
        std::printf("[enter %d] size=%d total=%.9g\n", (int)g_insCalls, (int)g_T.size(), totalAliveVolume());
        for (int ti = 0; ti < (int)g_T.size(); ++ti) {
            if (!g_T[ti].alive) continue;
            const Tet& tt = g_T[ti];
            double vv = std::fabs((double)orient3dVal(g_V[tt.v[0]], g_V[tt.v[1]], g_V[tt.v[2]], g_V[tt.v[3]])) / 6.0;
            std::printf("  E %d v=(%d %d %d %d) vol=%.9g\n", ti, tt.v[0], tt.v[1], tt.v[2], tt.v[3], vv);
        }
    }
    Vec3 p = p0;
    for (int attempt = 0; attempt < 4; ++attempt) {
        int c0 = locateCell(p, hint);
        if (c0 < 0) { ++g_failLoc; return false; }

        bool dup = false;
        for (int i = 0; i < 4; ++i) {
            if (norm(p - g_V[g_T[c0].v[i]]) < 1e-9 * g_diag) { dup = true; break; }
        }
        if (dup) { ++g_failDup; return false; }

        if ((int)g_stamp.size() < (int)g_T.size()) g_stamp.resize(g_T.size(), -1);
        ++g_curStamp;
        std::vector<int> cavity;
        cavity.reserve(64);
        g_stamp[c0] = g_curStamp;
        cavity.push_back(c0);
        bool overflow = false;
        for (size_t k = 0; k < cavity.size(); ++k) {
            const Tet& t = g_T[cavity[k]];
            for (int i = 0; i < 4; ++i) {
                int nb = t.nb[i];
                if (nb < 0 || g_stamp[nb] == g_curStamp) continue;
                const Tet& u = g_T[nb];
                if (inSphereSign(g_V[u.v[0]], g_V[u.v[1]], g_V[u.v[2]], g_V[u.v[3]], p) >= 0) {
                    g_stamp[nb] = g_curStamp;
                    cavity.push_back(nb);
                    if (cavity.size() > 200000) { overflow = true; break; }
                }
            }
            if (overflow) break;
        }
        if (overflow) { ++g_failCav; return false; }

        std::vector<BoundaryFace> bf;
        bf.reserve(cavity.size() * 2);
        for (int c : cavity) {
            const Tet& t = g_T[c];
            for (int i = 0; i < 4; ++i) {
                int nb = t.nb[i];
                if (nb < 0 || g_stamp[nb] != g_curStamp) bf.push_back({c, i, nb, t.nbi[i]});
            }
        }
        if (bf.empty()) { ++g_failCav; return false; }

        bool star = true;
        for (const BoundaryFace& f : bf) {
            Vec3 a, b, cc;
            faceVerts(g_T[f.cell], f.face, a, b, cc);
            if (orientSign(a, b, cc, p) > 0) { star = false; break; }
        }
        if (!star) {
            ++g_failStar;
            double s = (attempt + 1) * 1e-8 * g_diag;
            p = p0 + Vec3(s, -s, s);
            continue;
        }

        double volOldDbg = 0;
        if (g_debug && g_curStamp <= 6) {
            double vOld = 0, vNew = 0;
            double minOr = 1e300;
            for (int c : cavity) {
                const Tet& t = g_T[c];
                vOld += std::fabs((double)orient3dVal(g_V[t.v[0]], g_V[t.v[1]], g_V[t.v[2]], g_V[t.v[3]])) / 6.0;
            }
            for (const BoundaryFace& f : bf) {
                Vec3 a, b, cc;
                faceVerts(g_T[f.cell], f.face, a, b, cc);
                double o = (double)orient3dVal(a, b, cc, p);
                vNew += std::fabs(o) / 6.0;
                minOr = std::min(minOr, o);
            }
            volOldDbg = vOld;
            std::printf("[insert %d] cav=%d bf=%d volOld=%.9g volNew=%.9g minOrient=%.3g\n",
                        g_curStamp, (int)cavity.size(), (int)bf.size(), vOld, vNew, minOr);
        }
        // create new tets; p is vertex 3, boundary face is face 3
        int base = (int)g_T.size();
        for (const BoundaryFace& f : bf) {
            int j0 = FACE[f.face][0], j1 = FACE[f.face][1], j2 = FACE[f.face][2];
            int ia = g_T[f.cell].v[j0], ib = g_T[f.cell].v[j1], ic = g_T[f.cell].v[j2];
            if (orientSign(g_V[ia], g_V[ib], g_V[ic], p) < 0) std::swap(ib, ic);
            Tet t;
            t.v[0] = ia; t.v[1] = ib; t.v[2] = ic; t.v[3] = -1;
            g_T.push_back(t);
        }
        int pIdx = (int)g_V.size();
        g_V.push_back(p);

        // link new tets among themselves (faces 0,1,2 contain p)
        std::unordered_map<uint64_t, std::pair<int, int>> fmap;
        fmap.reserve(bf.size() * 3);
        for (size_t k = 0; k < bf.size(); ++k) {
            int nt = base + (int)k;
            g_T[nt].v[3] = pIdx;
            for (int f = 0; f < 3; ++f) {
                int j0 = FACE[f][0], j1 = FACE[f][1], j2 = FACE[f][2];
                uint64_t key = fkey(g_T[nt].v[j0], g_T[nt].v[j1], g_T[nt].v[j2]);
                auto it = fmap.find(key);
                if (it == fmap.end()) {
                    fmap[key] = {nt, f};
                } else {
                    int ot = it->second.first, of = it->second.second;
                    g_T[nt].nb[f] = ot; g_T[nt].nbi[f] = of;
                    g_T[ot].nb[of] = nt; g_T[ot].nbi[of] = f;
                    fmap.erase(it);
                }
            }
        }
        if (!fmap.empty()) {
            // inconsistent cavity: undo new tets, keep old mesh intact
            for (size_t k = 0; k < bf.size(); ++k) g_T[base + k].alive = false;
            ++g_failLink;
            return false;
        }

        // boundary links and old-neighbour updates
        for (size_t k = 0; k < bf.size(); ++k) {
            int nt = base + (int)k;
            int nb = bf[k].nb;
            g_T[nt].nb[3] = nb;
            g_T[nt].nbi[3] = bf[k].nbf;
            if (nb >= 0) {
                g_T[nb].nb[bf[k].nbf] = nt;
                g_T[nb].nbi[bf[k].nbf] = 3;
            }
        }
        if (g_debug && g_curStamp <= 8) {
            double vv = 0;
            for (size_t k = 0; k < bf.size(); ++k) {
                const Tet& nt = g_T[base + (int)k];
                vv += std::fabs((double)orient3dVal(g_V[nt.v[0]], g_V[nt.v[1]], g_V[nt.v[2]], g_V[nt.v[3]])) / 6.0;
            }
            std::printf("[created %d] n=%d volNewTets=%.9g\n", g_curStamp, (int)bf.size(), vv);
        }
        if (g_debug && g_curStamp <= 8) {
            double vAll = 0;
            for (int ti = 0; ti < (int)g_T.size(); ++ti) {
                if (!g_T[ti].alive) continue;
                const Tet& tt = g_T[ti];
                vAll += std::fabs((double)orient3dVal(g_V[tt.v[0]], g_V[tt.v[1]], g_V[tt.v[2]], g_V[tt.v[3]])) / 6.0;
            }
            std::printf("[before destroy %d] aliveVol=%.9g cavityVol=%.9g\n", g_curStamp, vAll, volOldDbg);
        }
        if (g_debug && g_curStamp <= 8) {
            double sumCav = 0;
            int badIdx = 0, dup = 0;
            std::vector<int> cs = cavity;
            std::sort(cs.begin(), cs.end());
            for (size_t i = 0; i < cs.size(); ++i) {
                if (i && cs[i] == cs[i-1]) ++dup;
                if (cs[i] >= base) ++badIdx;
                const Tet& tt = g_T[cs[i]];
                sumCav += std::fabs((double)orient3dVal(g_V[tt.v[0]], g_V[tt.v[1]], g_V[tt.v[2]], g_V[tt.v[3]])) / 6.0;
            }
            int oldAlive = 0;
            double oldSum = 0;
            for (int ti = 0; ti < base; ++ti) {
                if (!g_T[ti].alive) continue;
                ++oldAlive;
                const Tet& tt = g_T[ti];
                oldSum += std::fabs((double)orient3dVal(g_V[tt.v[0]], g_V[tt.v[1]], g_V[tt.v[2]], g_V[tt.v[3]])) / 6.0;
            }
            std::printf("[pre-destroy %d] base=%d oldAlive=%d oldSum=%.9g sumCav=%.9g dup=%d badIdx=%d\n",
                        g_curStamp, base, oldAlive, oldSum, sumCav, dup, badIdx);
            if (g_curStamp == 3) {
                for (int ti = 0; ti < base; ++ti) {
                    if (!g_T[ti].alive) continue;
                    const Tet& tt = g_T[ti];
                    double vv = std::fabs((double)orient3dVal(g_V[tt.v[0]], g_V[tt.v[1]], g_V[tt.v[2]], g_V[tt.v[3]])) / 6.0;
                    std::printf("  P %d v=(%d %d %d %d) vol=%.9g\n", ti, tt.v[0], tt.v[1], tt.v[2], tt.v[3], vv);
                }
            }
        }
        for (int c : cavity) g_T[c].alive = false;

        if (g_debug && g_curStamp <= 8) {
            std::printf("[post-destroy %d] total=%.9g\n", g_curStamp, totalAliveVolume());
        }
        newCells.resize(bf.size());
        for (size_t k = 0; k < bf.size(); ++k) newCells[k] = base + (int)k;
        hint = newCells[0];
        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------
// Wrapping
// ----------------------------------------------------------------------------
static double g_alpha = 0.0;
static double g_delta = 0.0;

static bool tetIntersectsInput(int ti) {
    const Tet& t = g_T[ti];
    Vec3 lo(1e300, 1e300, 1e300), hi(-1e300, -1e300, -1e300);
    for (int i = 0; i < 4; ++i) {
        const Vec3& p = g_V[t.v[i]];
        lo.x = std::min(lo.x, p.x); lo.y = std::min(lo.y, p.y); lo.z = std::min(lo.z, p.z);
        hi.x = std::max(hi.x, p.x); hi.y = std::max(hi.y, p.y); hi.z = std::max(hi.z, p.z);
    }
    double eps = 1e-9;
    // recursive query
    std::vector<int> stack;
    stack.push_back(0);
    while (!stack.empty()) {
        int n = stack.back();
        stack.pop_back();
        const BVHNode& nd = g_nodes[n];
        if (!boxOverlap(nd.lo, nd.hi, lo, hi)) continue;
        if (nd.count > 0) {
            for (int i = nd.start; i < nd.start + nd.count; ++i) {
                const Tri& tr = g_tris[g_order[i]];
                bool hit = false;
                // 1. any triangle vertex inside the tet
                for (int vi = 0; vi < 3 && !hit; ++vi) {
                    const Vec3& q = (vi == 0 ? tr.a : (vi == 1 ? tr.b : tr.c));
                    bool in = true;
                    for (int k = 0; k < 4; ++k) {
                        Vec3 a, b, c;
                        faceVerts(t, k, a, b, c);
                        if (orientSign(a, b, c, q) > 0) { in = false; break; }
                    }
                    if (in) hit = true;
                }
                // 2. any tet vertex on the triangle
                for (int vi = 0; vi < 4 && !hit; ++vi) {
                    const Vec3& q = g_V[t.v[vi]];
                    if (norm2(closestOnTri(q, tr.a, tr.b, tr.c) - q) < eps * eps) hit = true;
                }
                // 3. tet edges vs triangle
                if (!hit) {
                    static const int EDGES[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
                    for (int e = 0; e < 6 && !hit; ++e) {
                        const Vec3& q0 = g_V[t.v[EDGES[e][0]]];
                        const Vec3& q1 = g_V[t.v[EDGES[e][1]]];
                        if (segTri(q0, q1, tr.a, tr.b, tr.c, eps)) hit = true;
                    }
                }
                // 4. triangle edges vs tet faces
                if (!hit) {
                    Vec3 tv[4];
                    for (int k = 0; k < 4; ++k) tv[k] = g_V[t.v[k]];
                    static const int TFACES[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
                    Vec3 te[3] = {tr.a, tr.b, tr.c};
                    for (int e = 0; e < 3 && !hit; ++e) {
                        const Vec3& q0 = te[e];
                        const Vec3& q1 = te[(e + 1) % 3];
                        for (int f = 0; f < 4 && !hit; ++f) {
                            if (segTri(q0, q1, tv[TFACES[f][0]], tv[TFACES[f][1]], tv[TFACES[f][2]], eps))
                                hit = true;
                        }
                    }
                }
                if (hit) return true;
            }
        } else {
            stack.push_back(nd.left);
            stack.push_back(nd.right);
        }
    }
    return false;
}

static double totalAliveVolume() {
    double v = 0;
    for (int ti = 0; ti < (int)g_T.size(); ++ti) {
        if (!g_T[ti].alive) continue;
        const Tet& t = g_T[ti];
        v += std::fabs((double)orient3dVal(g_V[t.v[0]], g_V[t.v[1]], g_V[t.v[2]], g_V[t.v[3]])) / 6.0;
    }
    return v;
}

static void validateMesh(const char* tag) {
    double vol = 0;
    long long badOrient = 0, badNb = 0, alive = 0;
    for (int ti = 0; ti < (int)g_T.size(); ++ti) {
        const Tet& t = g_T[ti];
        if (!t.alive) continue;
        ++alive;
        double o = (double)orient3dVal(g_V[t.v[0]], g_V[t.v[1]], g_V[t.v[2]], g_V[t.v[3]]);
        if (!(o > 0)) ++badOrient;
        vol += std::fabs(o) / 6.0;
        for (int f = 0; f < 4; ++f) {
            int nb = t.nb[f];
            if (nb < 0) continue;
            if (!g_T[nb].alive) { ++badNb; continue; }
            bool ok = false;
            for (int k = 0; k < 4; ++k) if (g_T[nb].nb[k] == ti) { ok = true; break; }
            if (!ok) ++badNb;
        }
    }
    std::printf("[validate %s] alive=%d vol=%.9g badOrient=%d badNb=%d\n", tag, (int)alive, vol,
                (int)badOrient, (int)badNb);
}
static double gatePhi(int ti, int fi, Vec3* outKo, Vec3* outKn) {
    const Tet& t = g_T[ti];
    Vec3 a, b, c;
    faceVerts(t, fi, a, b, c);
    Vec3 kn = tetCircumcenter(ti, nullptr);
    Vec3 ko;
    int nb = t.nb[fi];
    if (nb < 0) {
        double rf = 0;
        ko = triCircumcenter(a, b, c, &rf);
    } else {
        ko = tetCircumcenter(nb, nullptr);
    }
    if (outKo) *outKo = ko;
    if (outKn) *outKn = kn;
    Vec3 d = kn - ko;
    double dd = norm2(d);
    double tt = 0.0;
    if (dd > 1e-300) tt = dot(a - ko, d) / dd;
    if (tt < 0) tt = 0;
    if (tt > 1) tt = 1;
    Vec3 m = ko + tt * d;
    return norm(m - a);
}

// first crossing of segment [A,B] with the offset level set, from A to B,
// where d(A) > delta and d(B) < delta (first + -> - transition)
static bool segmentCrossesOffset(const Vec3& A, const Vec3& B, Vec3& out) {
    double ga = distToInput(A) - g_delta;
    if (ga <= 0) return false;
    double gb = distToInput(B) - g_delta;
    if (gb < 0) {
        double lo = 0.0, hi = 1.0;
        for (int k = 0; k < 40; ++k) {
            double mid = 0.5 * (lo + hi);
            Vec3 xm = A + mid * (B - A);
            if (distToInput(xm) - g_delta > 0) lo = mid; else hi = mid;
        }
        out = A + lo * (B - A);
        return true;
    }
    const int N = 16;
    double gPrev = ga, tPrev = 0.0;
    for (int i = 1; i <= N; ++i) {
        double t = (double)i / N;
        double g = distToInput(A + t * (B - A)) - g_delta;
        if (gPrev > 0 && g < 0) {
            double lo = tPrev, hi = t;
            for (int k = 0; k < 40; ++k) {
                double mid = 0.5 * (lo + hi);
                if (distToInput(A + mid * (B - A)) - g_delta > 0) lo = mid; else hi = mid;
            }
            out = A + lo * (B - A);
            return true;
        }
        gPrev = g;
        tPrev = t;
    }
    return false;
}

static bool offsetPointTowardCenter(const Vec3& kn, Vec3& out) {
    Vec3 p = closestToInput(kn);
    double dk = norm(kn - p);
    if (dk <= g_delta * 1.0000001) return false;
    double lo = 0.0, hi = 1.0;
    Vec3 xm;
    for (int k = 0; k < 60; ++k) {
        double mid = 0.5 * (lo + hi);
        xm = kn + mid * (p - kn);
        double g = distToInput(xm) - g_delta;
        if (g > 0) lo = mid; else hi = mid;
    }
    out = kn + lo * (p - kn);
    return true;
}

struct Gate {
    double phi;
    int cell;
    int face;
    bool operator<(const Gate& o) const { return phi < o.phi; }
};

// ----------------------------------------------------------------------------
// OBJ IO
// ----------------------------------------------------------------------------
static bool loadOBJ(const char* path, std::vector<Vec3>& verts, std::vector<std::array<int, 3>>& faces) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string tag;
        ss >> tag;
        if (tag == "v") {
            double x, y, z;
            if (!(ss >> x >> y >> z)) continue;
            verts.emplace_back(x, y, z);
        } else if (tag == "f") {
            std::vector<int> idx;
            std::string tok;
            while (ss >> tok) {
                size_t slash = tok.find('/');
                std::string s = (slash == std::string::npos) ? tok : tok.substr(0, slash);
                if (s.empty()) continue;
                int i = std::atoi(s.c_str());
                if (i > 0) i -= 1;
                else if (i < 0) i = (int)verts.size() + i;
                else continue;
                idx.push_back(i);
            }
            if (idx.size() < 3) continue;
            for (size_t k = 1; k + 1 < idx.size(); ++k) faces.push_back({idx[0], idx[k], idx[k + 1]});
        }
    }
    return !verts.empty() && !faces.empty();
}

static void writeOBJ(const char* path, const std::vector<Vec3>& verts, const std::vector<std::array<int, 3>>& faces) {
    std::ofstream out(path);
    out << "# alpha wrap (standalone)\n";
    out << "# vertices " << verts.size() << "  faces " << faces.size() << "\n";
    for (const Vec3& v : verts) {
        char buf[128];
        std::snprintf(buf, sizeof(buf), "v %.9g %.9g %.9g\n", v.x, v.y, v.z);
        out << buf;
    }
    for (const auto& f : faces) out << "f " << f[0] + 1 << " " << f[1] + 1 << " " << f[2] + 1 << "\n";
}

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------
static bool selfTest() {
    // regular tet
    Vec3 a(0, 0, 0), b(1, 0, 0), c(0, 1, 0), d(0, 0, 1);
    if (orientSign(a, b, c, d) <= 0) { std::printf("selftest: orient sign failed\n"); return false; }
    Vec3 in(0.25, 0.25, 0.25), out(2, 2, 2);
    if (inSphereSign(a, b, c, d, in) <= 0) { std::printf("selftest: insphere inside failed\n"); return false; }
    if (inSphereSign(a, b, c, d, out) >= 0) { std::printf("selftest: insphere outside failed\n"); return false; }
    // outward faces
    for (int f = 0; f < 4; ++f) {
        Vec3 q;
        if (f == 0) q = a; else if (f == 1) q = b; else if (f == 2) q = c; else q = d;
        int i0 = FACE[f][0], i1 = FACE[f][1], i2 = FACE[f][2];
        Vec3 vv[4] = {a, b, c, d};
        if (orientSign(vv[i0], vv[i1], vv[i2], q) >= 0) {
            std::printf("selftest: face %d orientation failed\n", f);
            return false;
        }
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::printf("usage: %s input.obj output.obj [--alpha A] [--offset D] [--verbose]\n", argv[0]);
        return 1;
    }
    const char* inPath = argv[1];
    const char* outPath = argv[2];
    double alphaArg = -1, deltaArg = -1;
    bool verbose = false;
    for (int i = 3; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--alpha") && i + 1 < argc) alphaArg = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--offset") && i + 1 < argc) deltaArg = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--verbose")) verbose = true;
        else if (!std::strcmp(argv[i], "--debug")) g_debug = true;
    }
    g_verbose = verbose;
    if (!selfTest()) return 2;

    auto t0 = std::chrono::steady_clock::now();

    std::vector<Vec3> srcV;
    std::vector<std::array<int, 3>> srcF;
    if (!loadOBJ(inPath, srcV, srcF)) {
        std::printf("error: failed to load %s\n", inPath);
        return 1;
    }
    std::printf("input: %d vertices, %d triangles\n", (int)srcV.size(), (int)srcF.size());

    // build triangle soup
    g_tris.reserve(srcF.size());
    Vec3 bbLo(1e300, 1e300, 1e300), bbHi(-1e300, -1e300, -1e300);
    for (const auto& f : srcF) {
        Tri t;
        t.a = srcV[f[0]]; t.b = srcV[f[1]]; t.c = srcV[f[2]];
        triBBox(t, t.lo, t.hi);
        g_tris.push_back(t);
        bbLo.x = std::min(bbLo.x, t.lo.x); bbLo.y = std::min(bbLo.y, t.lo.y); bbLo.z = std::min(bbLo.z, t.lo.z);
        bbHi.x = std::max(bbHi.x, t.hi.x); bbHi.y = std::max(bbHi.y, t.hi.y); bbHi.z = std::max(bbHi.z, t.hi.z);
    }
    g_diag = norm(bbHi - bbLo);
    if (!(g_diag > 0)) { std::printf("error: degenerate input\n"); return 1; }
    g_epsLen = 1e-9 * g_diag;

    g_alpha = (alphaArg > 0) ? alphaArg : 0.02 * g_diag;
    g_delta = (deltaArg > 0) ? deltaArg : 0.002 * g_diag;
    if (g_delta >= g_alpha) std::printf("warning: offset >= alpha, output will be coarse\n");
    std::printf("alpha = %.6g  offset = %.6g  diag = %.6g\n", g_alpha, g_delta, g_diag);

    g_order.resize(g_tris.size());
    for (size_t i = 0; i < g_order.size(); ++i) g_order[i] = (int)i;
    buildBVH(0, (int)g_order.size());
    if (g_debug) {
        std::printf("[dist] origin=%.6g  (0,0,1.1)=%.6g  (0,0,2)=%.6g\n",
                    distToInput(Vec3(0,0,0)), distToInput(Vec3(0,0,1.1)), distToInput(Vec3(0,0,2)));
    }

    // initial loose box
    double pad = 0.03 * g_diag + 2.0 * g_delta;
    Vec3 lo = bbLo - Vec3(pad, pad, pad);
    Vec3 hi = bbHi + Vec3(pad, pad, pad);
    for (int i = 0; i < 8; ++i) {
        g_V.emplace_back((i & 1) ? hi.x : lo.x, (i & 2) ? hi.y : lo.y, (i & 4) ? hi.z : lo.z);
    }
    static const int BOXT[6][4] = {{0,1,3,7},{0,1,5,7},{0,4,5,7},{0,4,6,7},{0,2,6,7},{0,2,3,7}};
    for (int i = 0; i < 6; ++i) makeTet(BOXT[i][0], BOXT[i][1], BOXT[i][2], BOXT[i][3]);

    // link initial tets
    {
        std::unordered_map<uint64_t, std::pair<int, int>> m;
        m.reserve(g_T.size() * 4);
        for (int ti = 0; ti < (int)g_T.size(); ++ti) {
            for (int f = 0; f < 4; ++f) {
                int i0 = FACE[f][0], i1 = FACE[f][1], i2 = FACE[f][2];
                uint64_t key = fkey(g_T[ti].v[i0], g_T[ti].v[i1], g_T[ti].v[i2]);
                auto it = m.find(key);
                if (it == m.end()) {
                    m[key] = {ti, f};
                } else {
                    int ot = it->second.first, of = it->second.second;
                    g_T[ti].nb[f] = ot; g_T[ti].nbi[f] = of;
                    g_T[ot].nb[of] = ti; g_T[ot].nbi[of] = f;
                    m.erase(it);
                }
            }
        }
    }

    if (g_debug) validateMesh("init");
    std::priority_queue<Gate> pq;
    for (int ti = 0; ti < (int)g_T.size(); ++ti) {
        for (int f = 0; f < 4; ++f) {
            if (g_T[ti].nb[f] < 0) pq.push({gatePhi(ti, f, nullptr, nullptr), ti, f});
        }
    }

    if (g_debug) {
        for (int ti = 0; ti < (int)g_T.size(); ++ti) {
            const Tet& t = g_T[ti];
            std::printf("tet %d v=(%d %d %d %d) orient=%d intersectsInput=%d phi0=%.6g\n", ti,
                        t.v[0], t.v[1], t.v[2], t.v[3],
                        orientSign(g_V[t.v[0]], g_V[t.v[1]], g_V[t.v[2]], g_V[t.v[3]]),
                        (int)tetIntersectsInput(ti), gatePhi(ti, 0, nullptr, nullptr));
        }
        for (int i = 0; i < 3; ++i) {
            std::printf("input v%d = (%.3g %.3g %.3g) locate=%d\n", i, srcV[i].x, srcV[i].y, srcV[i].z,
                        locateCell(srcV[i], 0));
        }
    }

    int hint = 0;
    long long r1 = 0, r2 = 0, carved = 0, stalls = 0, inserted = 0;
    std::vector<int> newCells;
    while (!pq.empty()) {
        Gate g = pq.top();
        pq.pop();
        if (!g_T[g.cell].alive || g_T[g.cell].state != 0) continue;
        int nb = g_T[g.cell].nb[g.face];
        if (nb >= 0 && g_T[nb].state != 1) continue;
        Vec3 ko, kn;
        double phi = gatePhi(g.cell, g.face, &ko, &kn);
        if (!(phi > g_alpha)) continue;
        Vec3 ps;
        bool doInsert = false;
        bool r1hit = segmentCrossesOffset(ko, kn, ps);
        double rn = norm(kn - g_V[g_T[g.cell].v[0]]);
        if (g_debug && carved < 3) {
            std::printf("[gate] cell=%d face=%d phi=%.6g ko=(%.3g %.3g %.3g) kn=(%.3g %.3g %.3g) d(ko)=%.6g d(kn)=%.6g r1=%d psIn=%.6g/%.6g\n",
                        g.cell, g.face, phi, ko.x, ko.y, ko.z, kn.x, kn.y, kn.z,
                        distToInput(ko), distToInput(kn), (int)r1hit,
                        r1hit ? norm(ps - kn) : -1.0, rn);
        }
        if (r1hit) {
            if (norm(ps - kn) < rn * (1.0 - 1e-9)) {
                doInsert = true;
                ++r1;
            }
        }
        if (!doInsert && tetIntersectsInput(g.cell)) {
            if (offsetPointTowardCenter(kn, ps)) {
                double rn = norm(kn - g_V[g_T[g.cell].v[0]]);
                if (norm(ps - kn) < rn * (1.0 - 1e-9)) {
                    doInsert = true;
                    ++r2;
                }
            }
        }
        if (doInsert) {
            bool okIns = insertPoint(ps, hint, newCells);
            if (okIns) {
                ++inserted;
                static int dbgIns = 0;
                if (g_debug && dbgIns < 3) { ++dbgIns; validateMesh("after-insert"); }
                for (int nc : newCells) {
                    for (int f = 0; f < 4; ++f) {
                        int nb2 = g_T[nc].nb[f];
                        if (nb2 < 0 || g_T[nb2].state == 1) pq.push({gatePhi(nc, f, nullptr, nullptr), nc, f});
                    }
                }
                continue;
            }
            ++stalls;
            if (tetIntersectsInput(g.cell)) continue;  // cannot carve safely
        }
        // carve
        g_T[g.cell].state = 1;
        ++carved;
        for (int f = 0; f < 4; ++f) {
            int nb2 = g_T[g.cell].nb[f];
            if (nb2 >= 0 && g_T[nb2].state == 0) {
                int back = g_T[g.cell].nbi[f];
                pq.push({gatePhi(nb2, back, nullptr, nullptr), nb2, back});
            }
        }
    }

    // output
    std::vector<std::array<int, 3>> outF;
    std::vector<int> remap(g_V.size(), -1);
    std::vector<Vec3> outV;
    long long insideCells = 0, outsideCells = 0;
    for (int ti = 0; ti < (int)g_T.size(); ++ti) {
        if (!g_T[ti].alive) continue;
        if (g_T[ti].state == 0) ++insideCells; else ++outsideCells;
    }
    for (int ti = 0; ti < (int)g_T.size(); ++ti) {
        if (!g_T[ti].alive || g_T[ti].state != 0) continue;
        const Tet& t = g_T[ti];
        for (int f = 0; f < 4; ++f) {
            int nb = t.nb[f];
            if (nb >= 0 && g_T[nb].state != 1) continue;
            int ids[3] = {t.v[FACE[f][0]], t.v[FACE[f][1]], t.v[FACE[f][2]]};
            std::array<int, 3> tri;
            for (int k = 0; k < 3; ++k) {
                if (remap[ids[k]] < 0) {
                    remap[ids[k]] = (int)outV.size();
                    outV.push_back(g_V[ids[k]]);
                }
                tri[k] = remap[ids[k]];
            }
            outF.push_back(tri);
        }
    }

    // basic watertight check
    std::unordered_map<uint64_t, int> edgeCount;
    auto ekey = [](int a, int b) {
        if (a > b) std::swap(a, b);
        return ((uint64_t)a << 32) | (uint32_t)b;
    };
    for (const auto& f : outF) {
        for (int k = 0; k < 3; ++k) edgeCount[ekey(f[k], f[(k + 1) % 3])]++;
    }
    long long nonManifold = 0;
    for (auto& kv : edgeCount) if (kv.second != 2) ++nonManifold;

    writeOBJ(outPath, outV, outF);
    auto t1 = std::chrono::steady_clock::now();
    double sec = std::chrono::duration<double>(t1 - t0).count();
    std::printf("output: %d vertices, %d triangles\n", (int)outV.size(), (int)outF.size());
    std::printf("tets: %d (alive inside %d, outside %d), inserted %d, R1 %d, R2 %d, carved %d, stalls %d\n",
                (int)g_T.size(), (int)insideCells, (int)outsideCells, (int)inserted, (int)r1, (int)r2, (int)carved, (int)stalls);
    std::printf("non-manifold edges: %d\n", (int)nonManifold);
    std::printf("insert failures: locate %d dup %d star %d link %d cavity %d\n",
                (int)g_failLoc, (int)g_failDup, (int)g_failStar, (int)g_failLink, (int)g_failCav);
    std::printf("time: %.2f s\n", sec);
    if (verbose) {
        double dmin = 1e300, dmax = 0, dsum = 0;
        for (const Vec3& v : outV) {
            double d = distToInput(v);
            dmin = std::min(dmin, d);
            dmax = std::max(dmax, d);
            dsum += d;
        }
        std::printf("output vertex distance to input: min %.6g  mean %.6g  max %.6g (offset %.6g)\n",
                    dmin, dsum / std::max<size_t>(1, outV.size()), dmax, g_delta);
    }
    return 0;
}
