// sa_cover: simulated annealing for covering designs, after Nurmela and Ostergard (1993), run in the same
// decremental (shrink) mode as wls2 so the two can be compared from the same start covers and time budget.
//
// State: b blocks of size k; cost = number of uncovered t-subsets. Move: replace a random point x of a random block B
// by a random point y outside B. The change in cost is exact: the t-subsets of B that contain x and are covered only
// by B become uncovered, the uncovered t-subsets of B-x+y that contain y become covered (the two groups are disjoint).
// Metropolis acceptance at temperature T; T is multiplied by ALPHA every 1000 moves and reset to T0 when it falls
// below TMIN (reheating). Whenever the cost is 0 the cover is saved and the block covering the fewest t-subsets alone
// is deleted, as in wls2.
//
// Usage: sa_cover v k t b seconds seed init_file [targetB]      env: T0 (0.15), TMIN (0.02), ALPHA (0.995)
// T0 was tuned on C(57,12,3) over {0.05, 0.15, 0.5} x ALPHA {0.995, 0.9995}; 0.15 / 0.995 left the fewest uncovered.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
using namespace std;
typedef uint64_t u64;

int v, k, t, b;
u64 Cb[130][24];
vector<uint16_t> cov;
vector<vector<int>> blk;   // sorted point lists
long unc = 0;

inline u64 rank_sorted(const int* s) { u64 r = 0; for (int i = 0; i < t; ++i) r += Cb[s[i]][i + 1]; return r; }

// call f(rank) for every t-subset of the sorted point list P (size n) that contains point x (x in P)
template <class F> void tsets_with(const vector<int>& P, int x, F f) {
    int others[128], n = 0; for (int p : P) if (p != x) others[n++] = p;
    int idx[24]; for (int i = 0; i < t - 1; ++i) idx[i] = i;
    int s[24];
    while (true) {
        // merge x into the chosen t-1 others (sorted)
        int o = 0; bool placed = false;
        for (int i = 0; i < t - 1; ++i) { int q = others[idx[i]]; if (!placed && x < q) { s[o++] = x; placed = true; } s[o++] = q; }
        if (!placed) s[o++] = x;
        f(rank_sorted(s));
        if (t - 1 == 0) break;
        int i = t - 2; while (i >= 0 && idx[i] == n - (t - 1) + i) --i;
        if (i < 0) break;
        ++idx[i]; for (int j = i + 1; j < t - 1; ++j) idx[j] = idx[j - 1] + 1;
    }
}
template <class F> void tsets_all(const vector<int>& P, F f) {
    int n = P.size(); int idx[24]; for (int i = 0; i < t; ++i) idx[i] = i; int s[24];
    while (true) {
        for (int i = 0; i < t; ++i) s[i] = P[idx[i]];
        f(rank_sorted(s));
        int i = t - 1; while (i >= 0 && idx[i] == n - t + i) --i;
        if (i < 0) break;
        ++idx[i]; for (int j = i + 1; j < t; ++j) idx[j] = idx[j - 1] + 1;
    }
}

void save(const char* fn) {
    FILE* f = fopen(fn, "w");
    for (auto& B : blk) { for (int i = 0; i < k; ++i) fprintf(f, i ? " %d" : "%d", B[i] + 1); fputc('\n', f); }
    fclose(f);
}

int main(int argc, char** argv) {
    if (argc < 8) { fprintf(stderr, "usage: sa_cover v k t b seconds seed init_file [targetB]\n"); return 1; }
    v = atoi(argv[1]); k = atoi(argv[2]); t = atoi(argv[3]); b = atoi(argv[4]);
    double secs = atof(argv[5]); unsigned seed = atoi(argv[6]); int targetB = argc > 8 ? atoi(argv[8]) : 1;
    double T0 = getenv("T0") ? atof(getenv("T0")) : 0.15, TMIN = getenv("TMIN") ? atof(getenv("TMIN")) : 0.02;
    double ALPHA = getenv("ALPHA") ? atof(getenv("ALPHA")) : 0.995;
    for (int n = 0; n < 130; ++n) { Cb[n][0] = 1; for (int r = 1; r < 24; ++r) Cb[n][r] = n ? Cb[n - 1][r - 1] + Cb[n - 1][r] : 0; }
    cov.assign(Cb[v][t], 0);
    FILE* fi = fopen(argv[7], "r"); if (!fi) { perror(argv[7]); return 1; }
    int x; vector<int> cur; char line[4096];
    while (fgets(line, sizeof line, fi)) {
        cur.clear(); char* p = line; int nread;
        while (sscanf(p, "%d%n", &x, &nread) == 1) { cur.push_back(x - 1); p += nread; }
        if ((int)cur.size() == k) { sort(cur.begin(), cur.end()); blk.push_back(cur); }
    }
    fclose(fi);
    b = blk.size();
    for (auto& B : blk) tsets_all(B, [&](u64 r) { ++cov[r]; });
    for (u64 r = 0; r < cov.size(); ++r) unc += cov[r] == 0;
    printf("init: %d blocks, %ld uncovered\n", b, unc);
    mt19937_64 rng(seed);
    uniform_real_distribution<double> U01(0, 1);
    auto start = chrono::steady_clock::now();
    double T = T0; long moves = 0; int bestFull = 0;
    vector<char> inB(v);
    while (true) {
        if ((moves & 1023) == 0 && chrono::duration<double>(chrono::steady_clock::now() - start).count() > secs) break;
        if (unc == 0) {
            bestFull = b;
            char fn[256]; sprintf(fn, "sa_cover_%d_%d_%d_%d_%u.txt", v, k, t, b, seed); save(fn);
            printf("FULL %d blocks at %.1fs\n", b, chrono::duration<double>(chrono::steady_clock::now() - start).count()); fflush(stdout);
            if (b <= targetB) break;
            int bi = -1; long bl = 1L << 60;
            for (int i = 0; i < b; ++i) { long l = 0; tsets_all(blk[i], [&](u64 r) { l += cov[r] == 1; }); if (l < bl) { bl = l; bi = i; } }
            tsets_all(blk[bi], [&](u64 r) { if (--cov[r] == 0) ++unc; });
            blk[bi] = blk.back(); blk.pop_back(); --b;
            continue;
        }
        ++moves;
        if (moves % 1000 == 0) { T *= ALPHA; if (T < TMIN) T = T0; }
        int bi = rng() % b; vector<int>& B = blk[bi];
        int xo = B[rng() % k];
        fill(inB.begin(), inB.end(), 0); for (int p : B) inB[p] = 1;
        int yi; do { yi = rng() % v; } while (inB[yi]);
        long lost = 0, gained = 0;
        tsets_with(B, xo, [&](u64 r) { lost += cov[r] == 1; });
        vector<int> B2; B2.reserve(k); for (int p : B) if (p != xo) B2.push_back(p);
        B2.insert(upper_bound(B2.begin(), B2.end(), yi), yi);
        tsets_with(B2, yi, [&](u64 r) { gained += cov[r] == 0; });
        long d = lost - gained;
        if (d <= 0 || U01(rng) < exp(-d / T)) {
            tsets_with(B, xo, [&](u64 r) { if (--cov[r] == 0) ++unc; });
            tsets_with(B2, yi, [&](u64 r) { if (cov[r]++ == 0) --unc; });
            B = B2;
        }
    }
    if (bestFull > 0) printf("SHRINK result: full cover with %d blocks\n", bestFull);
    double el = chrono::duration<double>(chrono::steady_clock::now() - start).count();
    printf("SA C(%d,%d,%d) b=%d seed=%u: %ld uncovered, %ld moves, %.1fs (%.0f moves/s)\n", v, k, t, b, seed, unc, moves, el, moves / max(el, 1e-9));
}
