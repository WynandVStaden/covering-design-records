// Weighted local search (row-weighting + tabu/configuration checking) for covering designs.
// Ideas imported from SAT/MaxSAT local search (clause weighting, "breakout") and from
// RWLS for unicost set covering (Gao et al. 2015).
//
// State: b blocks (k-subsets of {0..v-1}). Each t-set S has a weight w_S (init 1).
// Step: pick an uncovered t-set T (mostly the heaviest of a small sample); among all blocks B with
// |B cap T| = t-1 and all xout in B \ T, choose the swap (xout -> xin = T \ B) maximising
//    gain = sum_{S newly covered} w_S - sum_{S newly uncovered} w_S,
// excluding tabu moves (re-adding a point recently removed from that block). Apply it (even if gain<0),
// then add +1 to the weight of every still-uncovered t-set. Track the true number uncovered.
// cov[S] = number of blocks containing S (a t-set with cov == 1 is "uniquely" covered).
//
// Usage: wls128 v k t b seconds seed [init_file|-] [tenure=7] [noise=0.05] [targetB]
// Shrink mode: whenever the cover is complete it is saved as shrink_cover_v_k_t_b_seed.txt, the block covering
// the fewest t-sets uniquely is deleted, and the search repairs the hole; this repeats until b <= targetB or time
// runs out. Compile: g++ -O3 -march=native -o wls128 wls128.cpp (GCC or Clang; needs unsigned __int128).
#include <cstdio>
#include <cstdlib>
#include <cstdint>
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

int v, k, t, b;
u64 Cb[130][20];
inline u64 binom(int n, int r) { return (r < 0 || n < r || n < 0) ? 0 : Cb[n][r]; }

vector<uint16_t> cov;
vector<uint32_t> w;
vector<u64> unc; vector<int64_t> upos;
vector<M> bm;                   // block bitmasks (128-bit)
vector<vector<int>> tabuUntil;  // tabuUntil[block][point] = step until which adding point is tabu
long step = 0;

inline void add_unc(u64 r) { upos[r] = unc.size(); unc.push_back(r); }
inline void rem_unc(u64 r) { int64_t p = upos[r]; u64 l = unc.back(); unc[p] = l; upos[l] = p; unc.pop_back(); upos[r] = -1; }

// rank of t-set given as bitmask (bits ascending)
inline u64 rank_mask(M m) {
    u64 r = 0; int i = 1;
    while (m) { int p = ctz128(m); r += binom(p, i); ++i; m &= m - 1; }
    return r;
}
inline M unrank(u64 r) {
    M m = 0;
    for (int i = t; i >= 1; --i) { int c = i - 1; while (binom(c + 1, i) <= r) ++c; m |= BIT(c); r -= binom(c, i); }
    return m;
}

// enumerate t-subsets of mask B that contain point x, call f(rank)
template <class F> inline void each_with(M B, int x, F f) {
    M rest = B & ~BIT(x);
    int pts[128]; int m = 0; for (M q = rest; q; q &= q - 1) pts[m++] = ctz128(q);
    int tt = t - 1; int idx[20];
    for (int i = 0; i < tt; ++i) idx[i] = i;
    while (true) {
        M S = BIT(x); for (int i = 0; i < tt; ++i) S |= BIT(pts[idx[i]]);
        f(rank_mask(S));
        if (tt == 0) return;
        int i = tt - 1; while (i >= 0 && idx[i] == m - tt + i) --i;
        if (i < 0) return;
        ++idx[i]; for (int l = i + 1; l < tt; ++l) idx[l] = idx[l - 1] + 1;
    }
}
template <class F> inline void each_all(M B, F f) {
    int pts[128]; int m = 0; for (M q = B; q; q &= q - 1) pts[m++] = ctz128(q);
    int idx[20]; for (int i = 0; i < t; ++i) idx[i] = i;
    while (true) {
        M S = 0; for (int i = 0; i < t; ++i) S |= BIT(pts[idx[i]]);
        f(rank_mask(S));
        int i = t - 1; while (i >= 0 && idx[i] == m - t + i) --i;
        if (i < 0) return;
        ++idx[i]; for (int l = i + 1; l < t; ++l) idx[l] = idx[l - 1] + 1;
    }
}

int main(int argc, char** argv) {
    v = atoi(argv[1]); k = atoi(argv[2]); t = atoi(argv[3]); b = atoi(argv[4]);
    double secs = atof(argv[5]); unsigned seed = atoi(argv[6]);
    const char* init = (argc > 7 && argv[7][0] != '-') ? argv[7] : nullptr;
    int tenure = argc > 8 ? atoi(argv[8]) : 7;
    double noise = argc > 9 ? atof(argv[9]) : 0.05;
    for (int n = 0; n < 130; ++n) { Cb[n][0] = 1; for (int r = 1; r < 20; ++r) Cb[n][r] = n ? Cb[n - 1][r - 1] + Cb[n - 1][r] : 0; }
    u64 NT = binom(v, t);
    cov.assign(NT, 0); w.assign(NT, 1); upos.assign(NT, -1);
    mt19937_64 rng(seed);
    bm.assign(b, 0);
    int loaded = 0;
    if (init) {
        FILE* f = fopen(init, "r");
        for (int i = 0; i < b && f; ++i) {
            M m = 0; int cnt = 0;
            for (int j = 0; j < k; ++j) { int p; if (fscanf(f, "%d", &p) == 1) { m |= BIT(p - 1); ++cnt; } }
            if (cnt == k) { bm[i] = m; ++loaded; } else break;
        }
        if (f) fclose(f);
    }
    for (int i = loaded; i < b; ++i) {
        vector<int> pts(v); for (int p = 0; p < v; ++p) pts[p] = p; shuffle(pts.begin(), pts.end(), rng);
        M m = 0; for (int j = 0; j < k; ++j) m |= BIT(pts[j]); bm[i] = m;
    }
    for (int i = 0; i < b; ++i) each_all(bm[i], [&](u64 r) { cov[r]++; });
    for (u64 r = 0; r < NT; ++r) if (cov[r] == 0) add_unc(r);
    tabuUntil.assign(b, vector<int>(v, 0));
    long best = unc.size(); vector<M> bestB = bm;
    fprintf(stderr, "init: loaded %d blocks, uncovered %ld\n", loaded, best);
    auto start = chrono::steady_clock::now();
    uniform_real_distribution<double> U01(0, 1);
    
    int targetB = argc > 10 ? atoi(argv[10]) : b;
    int bestFull = -1; vector<M> bestFullB;
    double levelStart = 0;
    while (true) {
        if (unc.empty()) {
            bestFull = b; bestFullB = bm;
            {   // save every full cover reached
                char fn[256]; sprintf(fn, "shrink_cover_%d_%d_%d_%d_%u.txt", v, k, t, b, seed);
                FILE* f = fopen(fn, "w");
                for (M m : bm) { int c = 0; for (int p = 0; p < v; ++p) if ((m >> p) & 1) fprintf(f, "%d%c", p + 1, ++c < k ? ' ' : '\n'); }
                fclose(f);
                double el = chrono::duration<double>(chrono::steady_clock::now() - start).count();
                printf("FULL %d blocks at %.1fs\n", b, el); fflush(stdout);
            }
            if (b <= targetB) break;
            // delete the block whose removal uncovers fewest t-sets
            int bi = -1; long bl = 1L << 60;
            for (int i = 0; i < b; ++i) { long l = 0; each_all(bm[i], [&](u64 r) { if (cov[r] == 1) ++l; }); if (l < bl) { bl = l; bi = i; } }
            each_all(bm[bi], [&](u64 r) { if (--cov[r] == 0) add_unc(r); });
            bm[bi] = bm[b - 1]; bm.pop_back(); tabuUntil[bi] = tabuUntil[b - 1]; tabuUntil.pop_back(); --b;
            for (u64 r : unc) w[r] = 1;
            best = unc.size(); bestB = bm;
            continue;
        }
        if ((step & 1023) == 0) {
            double el = chrono::duration<double>(chrono::steady_clock::now() - start).count();
            if (el > secs) break;
        }
        ++step;
        // pick target uncovered t-set: heaviest of a small sample (or random with prob noise)
        u64 rT;
        if (U01(rng) < noise) rT = unc[rng() % unc.size()];
        else {
            rT = unc[rng() % unc.size()];
            for (int s = 0; s < 3; ++s) { u64 r2 = unc[rng() % unc.size()]; if (w[r2] > w[rT]) rT = r2; }
        }
        M T = unrank(rT);
        // evaluate moves
        long bestGain = -(1L << 60); int bi = -1, bout = -1, bin_ = -1; int nties = 0;
        for (int i = 0; i < b; ++i) {
            M I = bm[i] & T;
            if (pc128(I) != t - 1) continue;
            int xin = ctz128(T & ~bm[i]);
            if (tabuUntil[i][xin] > step) continue;
            // make: t-subsets of B' containing xin that are uncovered (independent of xout except they must not contain xout)
            for (M q = bm[i] & ~T; q; q &= q - 1) {
                int xout = ctz128(q);
                M B2 = (bm[i] & ~BIT(xout)) | BIT(xin);
                long brk = 0, mk = 0;
                each_with(bm[i], xout, [&](u64 r) { if (cov[r] == 1) brk += w[r]; });
                each_with(B2, xin, [&](u64 r) { if (cov[r] == 0) mk += w[r]; });
                long g = mk - brk;
                if (g > bestGain) { bestGain = g; bi = i; bout = xout; bin_ = xin; nties = 1; }
                else if (g == bestGain && (rng() % (++nties)) == 0) { bi = i; bout = xout; bin_ = xin; }
            }
        }
        if (bi < 0) {  // all candidate moves tabu: random move on target
            vector<int> cands; for (int i = 0; i < b; ++i) if (pc128(bm[i] & T) == t - 1) cands.push_back(i);
            if (cands.empty()) continue;
            bi = cands[rng() % cands.size()]; bin_ = ctz128(T & ~bm[bi]);
            M rest = bm[bi] & ~T; int nr = pc128(rest); int pick = rng() % nr;
            for (M q = rest; q; q &= q - 1, --pick) if (pick == 0) { bout = ctz128(q); break; }
        }
        // apply
        M B = bm[bi];
        each_with(B, bout, [&](u64 r) { if (--cov[r] == 0) add_unc(r); });
        M B2 = (B & ~BIT(bout)) | BIT(bin_);
        each_with(B2, bin_, [&](u64 r) { if (cov[r]++ == 0) rem_unc(r); });
        bm[bi] = B2;
        tabuUntil[bi][bout] = step + tenure + (int)(rng() % 4);
        // weighting: bump uncovered t-sets
        for (u64 r : unc) w[r] += 1;
        long c = unc.size();
        if (c < best) { best = c; bestB = bm; }
    }
    if (bestFull > 0) printf("SHRINK result: full cover with %d blocks\n", bestFull);
    double el = chrono::duration<double>(chrono::steady_clock::now() - start).count();
    printf("WLS C(%d,%d,%d) b=%d seed=%u tenure=%d noise=%.2f: best %ld uncovered, %ld steps, %.1fs (%.0f steps/s)\n",
           v, k, t, b, seed, tenure, noise, best, step, el, step / max(el, 1e-9));
    if (best <= 5) {
        char fn[256];
        if (best == 0) sprintf(fn, "wls_cover_%d_%d_%d_%d_%u.txt", v, k, t, b, seed);
        else sprintf(fn, "wls_near_%d_%d_%d_%d_%u_c%ld.txt", v, k, t, b, seed, best);
        FILE* f = fopen(fn, "w");
        for (M m : bestB) { int c = 0; for (int p = 0; p < v; ++p) if ((m >> p) & 1) fprintf(f, "%d%c", p + 1, ++c < k ? ' ' : '\n'); }
        fclose(f); printf("saved %s\n", fn);
    }
}
