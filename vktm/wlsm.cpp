// wlsm: weighted local search for (v,k,t,m) lottery coverings.
//
// A (v,k,t,m) covering is a set of k-blocks such that every m-subset M of {1..v} meets some block in >= t points.
// State: cnt[M] = number of blocks meeting M in >= t points; weights w[M] on m-sets (breakout weighting).
//
// A move swaps x_out -> x_in in block B. Only m-sets M = S + {x} + R change status, where S is a (t-1)-subset of
// B - x_out and R an (m-t)-subset of the points outside B + x_in:
//   x = x_out: |M cap B| = t before, t-1 after  -> uncovered if cnt[M] == 1
//   x = x_in : |M cap B| = t-1 before, t after  -> covered   if cnt[M] == 0
// so the weighted gain of a move is one pass over the pairs (S, R).
//
// Each step: pick a heavy uncovered m-set M, sample candidate moves that make some block meet M in t points
// (block with |B cap M| = t-1, x_in in M - B, x_out in B - M), apply the best non-tabu one, bump the weights of
// uncovered m-sets. Shrink mode as in wls2: on a full cover, save it, delete the block covering the fewest m-sets
// uniquely, repair.
//
// Usage: wlsm v k t m b seconds seed [init|-] [cands=8] [tenure=5] [targetB]
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
using namespace std;
typedef uint64_t u64;
typedef unsigned __int128 M;
static inline int ctz128(M x) { u64 lo = (u64)x; return lo ? __builtin_ctzll(lo) : 64 + __builtin_ctzll((u64)(x >> 64)); }
static inline int pc128(M x) { return __builtin_popcountll((u64)x) + __builtin_popcountll((u64)(x >> 64)); }
#define BIT(p) (((M)1) << (p))

int v, k, t, m, b;
u64 Cb[130][24];
inline u64 binom(int n, int r) { return (r < 0 || n < r || n < 0) ? 0 : Cb[n][r]; }

vector<uint16_t> cnt;
vector<uint32_t> w;
vector<uint32_t> unc; vector<int32_t> upos;
vector<M> bm;
vector<int> tabuUntil;  // [block*128 + point]: re-adding point forbidden until step
long step = 0;

inline void add_unc(uint32_t r) { upos[r] = unc.size(); unc.push_back(r); }
inline void rem_unc(uint32_t r) { int32_t p = upos[r]; uint32_t l = unc.back(); unc[p] = l; upos[l] = p; unc.pop_back(); upos[r] = -1; }

inline int pts_of(M x, int* p) { int n = 0; for (; x; x &= x - 1) p[n++] = ctz128(x); return n; }
inline uint32_t rank_sorted(const int* p, int n) { u64 r = 0; for (int i = 0; i < n; ++i) r += binom(p[i], i + 1); return (uint32_t)r; }
// colex rank by 16-bit chunks: RT[(chunk*65536 + bits)*(m+1) + (#points below the chunk)]
vector<uint32_t> RT; int NCH;
void build_rank_tables() {
    NCH = (v + 15) / 16; RT.assign((size_t)NCH * 65536 * (m + 1), 0);
    for (int c = 0; c < NCH; ++c)
        for (int bits = 1; bits < 65536; ++bits)
            for (int off = 0; off <= m; ++off) {
                u64 r = 0; int idx = off; bool ok = true;
                for (int p = 0; p < 16; ++p) if ((bits >> p) & 1) { if (++idx > m) { ok = false; break; } r += binom(16 * c + p, idx); }
                if (ok) RT[((size_t)c * 65536 + bits) * (m + 1) + off] = (uint32_t)r;
            }
}
inline uint32_t rank_mask(M x) {
    uint32_t r = 0; int off = 0;
    for (int c = 0; c < NCH; ++c) { unsigned bits = (unsigned)(x >> (16 * c)) & 0xffff; if (bits) { r += RT[((size_t)c * 65536 + bits) * (m + 1) + off]; off += __builtin_popcount(bits); } }
    return r;
}
inline M unrank(u64 r) {
    M x = 0;
    for (int i = m; i >= 1; --i) { int c = i - 1; while (binom(c + 1, i) <= r) ++c; x |= BIT(c); r -= binom(c, i); }
    return x;
}

template <class F> inline void subsets(int n, int r, F f) {
    if (r > n) return;
    int idx[24]; for (int i = 0; i < r; ++i) idx[i] = i;
    while (true) {
        f(idx);
        int i = r - 1; while (i >= 0 && idx[i] == n - r + i) --i;
        if (i < 0) return;
        ++idx[i]; for (int l = i + 1; l < r; ++l) idx[l] = idx[l - 1] + 1;
    }
}

// call f(mask_lose, mask_gain) for every pair (S, R): lose = S + xout + R, gain = S + xin + R
template <class F> inline void pairs(M B, int xout, int xin, F f) {
    int Sp[24]; int ns = pts_of(B & ~BIT(xout), Sp);
    M outside = ~B & ~BIT(xin);
    if (v < 128) outside &= (BIT(v) - 1);
    int Rp[128]; int nr = pts_of(outside, Rp);
    // precompute the R masks once
    static vector<M> Rm; Rm.clear();
    subsets(nr, m - t, [&](int* idx) { M x = 0; for (int i = 0; i < m - t; ++i) x |= BIT(Rp[idx[i]]); Rm.push_back(x); });
    subsets(ns, t - 1, [&](int* idx) {
        M s = 0; for (int i = 0; i < t - 1; ++i) s |= BIT(Sp[idx[i]]);
        M so = s | BIT(xout), si = s | BIT(xin);
        for (M r : Rm) f(so | r, si | r);
    });
}

// all m-sets met by block B in >= t points
template <class F> inline void covered_by(M B, F f) {
    int Bp[24]; int nb = pts_of(B, Bp);
    M outside = ~B; if (v < 128) outside &= (BIT(v) - 1);
    int Op[128]; int no = pts_of(outside, Op);
    for (int j = t; j <= min(k, m); ++j) {
        subsets(nb, j, [&](int* ia) {
            M s = 0; for (int i = 0; i < j; ++i) s |= BIT(Bp[ia[i]]);
            subsets(no, m - j, [&](int* ib) { M x = s; for (int i = 0; i < m - j; ++i) x |= BIT(Op[ib[i]]); f(x); });
        });
    }
}

void save(const char* fn) {
    FILE* f = fopen(fn, "w");
    for (M x : bm) { int c = 0; for (int p = 0; p < v; ++p) if ((x >> p) & 1) fprintf(f, "%d%c", p + 1, ++c < k ? ' ' : '\n'); }
    fclose(f);
}

int main(int argc, char** argv) {
    if (argc < 8) { fprintf(stderr, "usage: wlsm v k t m b seconds seed [init|-] [cands=8] [tenure=5] [targetB]\n"); return 1; }
    v = atoi(argv[1]); k = atoi(argv[2]); t = atoi(argv[3]); m = atoi(argv[4]); b = atoi(argv[5]);
    double secs = atof(argv[6]); unsigned seed = atoi(argv[7]);
    const char* init = (argc > 8 && argv[8][0] != '-') ? argv[8] : nullptr;
    int ncand = argc > 9 ? atoi(argv[9]) : 8;
    int tenure = argc > 10 ? atoi(argv[10]) : 5;
    int targetB = argc > 11 ? atoi(argv[11]) : b;
    for (int n = 0; n < 130; ++n) { Cb[n][0] = 1; for (int r = 1; r < 24; ++r) Cb[n][r] = n ? Cb[n - 1][r - 1] + Cb[n - 1][r] : 0; }
    u64 NM = binom(v, m);
    build_rank_tables();
    if (NM > 400000000ULL) { fprintf(stderr, "C(v,m) = %llu too large\n", (unsigned long long)NM); return 1; }
    cnt.assign(NM, 0); w.assign(NM, 1); upos.assign(NM, -1);
    mt19937_64 rng(seed);
    bm.assign(b, 0);
    int loaded = 0;
    if (init) {
        FILE* f = fopen(init, "r");
        for (int i = 0; i < b && f; ++i) {
            M x = 0; int c = 0;
            for (int j = 0; j < k; ++j) { int p; if (fscanf(f, "%d", &p) == 1) { x |= BIT(p - 1); ++c; } }
            if (c == k) { bm[i] = x; ++loaded; } else break;
        }
        if (f) fclose(f);
    }
    for (int i = loaded; i < b; ++i) {
        vector<int> pts(v); for (int p = 0; p < v; ++p) pts[p] = p; shuffle(pts.begin(), pts.end(), rng);
        M x = 0; for (int j = 0; j < k; ++j) x |= BIT(pts[j]); bm[i] = x;
    }
    for (int i = 0; i < b; ++i) covered_by(bm[i], [&](M x) { uint32_t r = rank_mask(x); cnt[r]++; });
    for (u64 r = 0; r < NM; ++r) if (!cnt[r]) add_unc((uint32_t)r);
    tabuUntil.assign((size_t)b * 128, 0);
    long best = unc.size();
    fprintf(stderr, "init: loaded %d blocks, C(v,m)=%llu, uncovered %ld\n", loaded, (unsigned long long)NM, best);
    auto start = chrono::steady_clock::now();
    int bestFull = -1;
    vector<int> cand;
    while (true) {
        if (unc.empty()) {
            bestFull = b;
            char fn[256]; sprintf(fn, "mcover_%d_%d_%d_%d_%d_%u.txt", v, k, t, m, b, seed); save(fn);
            double el = chrono::duration<double>(chrono::steady_clock::now() - start).count();
            printf("FULL %d blocks at %.1fs (%ld steps)\n", b, el, step); fflush(stdout);
            if (b <= targetB) break;
            int bi = -1; long bl = 1L << 60;
            for (int i = 0; i < b; ++i) { long l = 0; covered_by(bm[i], [&](M x) { if (cnt[rank_mask(x)] == 1) ++l; }); if (l < bl) { bl = l; bi = i; } }
            covered_by(bm[bi], [&](M x) { uint32_t r = rank_mask(x); if (--cnt[r] == 0) add_unc(r); });
            bm[bi] = bm[b - 1]; bm.pop_back();
            memcpy(&tabuUntil[(size_t)bi * 128], &tabuUntil[(size_t)(b - 1) * 128], 128 * sizeof(int));
            --b;
            for (uint32_t r : unc) w[r] = 1;
            best = unc.size();
            continue;
        }
        if ((step & 63) == 0) {
            double el = chrono::duration<double>(chrono::steady_clock::now() - start).count();
            if (el > secs) break;
        }
        ++step;
        // target: heaviest of a small sample of uncovered m-sets
        uint32_t rT = unc[rng() % unc.size()];
        for (int s = 0; s < 3; ++s) { uint32_t r2 = unc[rng() % unc.size()]; if (w[r2] > w[rT]) rT = r2; }
        M T = unrank(rT);
        cand.clear();
        for (int i = 0; i < b; ++i) if (pc128(bm[i] & T) == t - 1) cand.push_back(i);
        if (cand.empty()) {  // every block meets T in < t-1 points: pull a random block towards T
            int i = rng() % b; M B = bm[i];
            int xin = -1, xout = -1;
            M inT = T & ~B; M outB = B & ~T;
            int a[128]; int na = pts_of(inT, a); int c[128]; int nc = pts_of(outB, c);
            xin = a[rng() % na]; xout = c[rng() % nc];
            pairs(B, xout, xin, [&](M lo, M gi) {
                uint32_t rl = rank_mask(lo); if (--cnt[rl] == 0) add_unc(rl);
                uint32_t rg = rank_mask(gi); if (cnt[rg]++ == 0) rem_unc(rg);
            });
            bm[i] = (B & ~BIT(xout)) | BIT(xin);
            continue;
        }
        // exact evaluation of every move (i, xin in T-B, xout in B-T) over up to ncand candidate blocks:
        //   gain(xin, xout) = G[xin] - g1[xin][xout]   (m-sets S+xin+R, S in B, uncovered; minus those with xout in S)
        //   loss(xout, xin) = L[xout] - l1[xout][xin]  (m-sets S+xout+R, S in B-xout, cnt==1; minus those with xin in R)
        if ((int)cand.size() > ncand) { shuffle(cand.begin(), cand.end(), rng); cand.resize(ncand); }
        long bestGain = -(1L << 60); int bi = -1, bout = -1, bin_ = -1; int nties = 0;
        static int64_t g1[128][128], l1[128][128]; int64_t G[128], L[128];
        for (int i : cand) {
            M B = bm[i];
            int Bp[24]; int nb = pts_of(B, Bp);
            M outB = ~B; if (v < 128) outB &= (BIT(v) - 1);
            int a[24]; int na = pts_of(T & ~B, a); int o[24]; int no = pts_of(B & ~T, o);
            for (int ia = 0; ia < na; ++ia) {          // gain side, one pass per xin
                int xin = a[ia]; G[xin] = 0; for (int j = 0; j < nb; ++j) g1[xin][Bp[j]] = 0;
                int Op[128]; int nO = pts_of(outB & ~BIT(xin), Op);
                subsets(nb, t - 1, [&](int* is) {
                    M s = BIT(xin); for (int q = 0; q < t - 1; ++q) s |= BIT(Bp[is[q]]);
                    subsets(nO, m - t, [&](int* ir) {
                        M x = s; for (int q = 0; q < m - t; ++q) x |= BIT(Op[ir[q]]);
                        uint32_t r = rank_mask(x);
                        if (cnt[r] == 0) { int64_t ww = w[r]; G[xin] += ww; for (int q = 0; q < t - 1; ++q) g1[xin][Bp[is[q]]] += ww; }
                    });
                });
            }
            int Op[128]; int nO = pts_of(outB, Op);
            for (int io = 0; io < no; ++io) {          // loss side, one pass per xout
                int xout = o[io]; L[xout] = 0; for (int j = 0; j < na; ++j) l1[xout][a[j]] = 0;
                int Sp[24]; int ns = pts_of(B & ~BIT(xout), Sp);
                subsets(ns, t - 1, [&](int* is) {
                    M s = BIT(xout); for (int q = 0; q < t - 1; ++q) s |= BIT(Sp[is[q]]);
                    subsets(nO, m - t, [&](int* ir) {
                        M x = s; for (int q = 0; q < m - t; ++q) x |= BIT(Op[ir[q]]);
                        uint32_t r = rank_mask(x);
                        if (cnt[r] == 1) { int64_t ww = w[r]; L[xout] += ww; for (int q = 0; q < m - t; ++q) { int p = Op[ir[q]]; if ((T >> p) & 1) l1[xout][p] += ww; } }
                    });
                });
            }
            for (int ia = 0; ia < na; ++ia) {
                int xin = a[ia];
                if (tabuUntil[(size_t)i * 128 + xin] > step) continue;
                for (int io = 0; io < no; ++io) {
                    int xout = o[io];
                    long g = (G[xin] - g1[xin][xout]) - (L[xout] - l1[xout][xin]);
                    if (g > bestGain) { bestGain = g; bi = i; bout = xout; bin_ = xin; nties = 1; }
                    else if (g == bestGain && rng() % (++nties) == 0) { bi = i; bout = xout; bin_ = xin; }
                }
            }
        }
        if (bi < 0) continue;
        M B = bm[bi];
        // apply: lose first (cnt can only drop for lose-sets, rise for gain-sets; the two families are disjoint)
        pairs(B, bout, bin_, [&](M lo, M gi) {
            uint32_t rl = rank_mask(lo); if (--cnt[rl] == 0) add_unc(rl);
            uint32_t rg = rank_mask(gi); if (cnt[rg]++ == 0) rem_unc(rg);
        });
        bm[bi] = (B & ~BIT(bout)) | BIT(bin_);
        tabuUntil[(size_t)bi * 128 + bout] = step + tenure + (int)(rng() % 3);
        for (uint32_t r : unc) w[r] += 1;
        if ((long)unc.size() < best) best = unc.size();
    }
    if (bestFull > 0) printf("SHRINK result: full cover with %d blocks\n", bestFull);
    double el = chrono::duration<double>(chrono::steady_clock::now() - start).count();
    printf("WLSM C(%d,%d,%d,%d) b=%d: best %ld uncovered, %ld steps, %.1fs (%.0f steps/s)\n", v, k, t, m, b, best, step, el, step / max(el, 1e-9));
}
