// wlsm2big: wlsm2 (weighted local search for (v,k,t,m) coverings with t-level move scores) engineered for large
// covers: C(v,m) up to ~1.6e9 and hundreds of thousands of blocks. Same search, same shrink mode; additions:
//   - witness byte sw[M] (sum of local indices of covered t-subsets), so Bc lives only on the unique covered t-subset
//     and a 1<->2 transition costs one update instead of C(m,t); O(m-t) colex ranking of t-subsets
//   - parallel (OpenMP) setup, score rebuild, superset enumeration per flip, candidate scan and candidate scoring
//   - POOL / POOL_REFRESH: deletion candidates from a pool of the lowest-scoring blocks, early exit at score 0
//   - NORESET + BUMP_GROWTH (per mille per level): keep weights between levels, geometric bump = decay of old weights
//   - LEVEL_BUDGET / GROWTH / REVERT_SECS: revert to the last full cover (undo log) and delete another block when a
//     deletion cannot be repaired; STALL stops after that many seconds without a new best cover
//   - SAVE_SECS throttles saving; the best full cover is always saved at exit, also on SIGTERM/SIGINT
//   - cost guard: refuses (exit 3) when one t-set flip would cost more than MAXFLIP (default 1e8) updates, or when
//     C(m,t) > 256 (the witness byte could not name the covered t-subset)
// Build with -DCHECK for a brute-force recount of all bookkeeping every 2000 steps, -DPROF for phase timings.
//
// Settings used for the published re-shrinks: POOL=2000 POOL_REFRESH=500 NORESET=1 BUMP_GROWTH=30 LEVEL_BUDGET=3000
// GROWTH=50 REVERT_SECS=120 SAVE_SECS=300 (STALL=60 for batches), cands=1000, tenure=5.
//
// An m-set M is covered iff one of its t-subsets lies inside a block. So coverage is tracked at two levels:
//   ct[T]  = number of blocks containing t-set T           (as in the (v,k,t) solver)
//   cnt[M] = number of covered t-subsets of m-set M        (M covered iff cnt[M] > 0)
// Breakout weights live on m-sets and are pushed down to t-sets:
//   A[T]  = sum of w[M] over uncovered M containing T      (gain from covering T)
//   Bc[T] = sum of w[M] over M with cnt[M] == 1 containing T (for a covered T: loss from uncovering it)
// A swap x_out -> x_in in block B is scored with C(k-1,t-1)-size sums over t-sets of B, like wls2:
//   loss(x_out) = sum Bc[T] over T in B containing x_out with ct[T] == 1
//   gain(x_in)  = sum A[T]  over T in B - x_out + x_in containing x_in with ct[T] == 0
// (an m-set containing two flipping t-sets is counted twice; the applied move is always exact).
// m-sets are only touched when a t-set flips between covered and uncovered.
//
// Usage: wlsm2big v k t m b seconds seed [init|-] [cands=12] [tenure=5] [targetB]
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <omp.h>
#include <unordered_map>
#include <csignal>
#include <mutex>
static std::mutex uncMu;
static volatile sig_atomic_t stopFlag = 0;
static void on_signal(int) { stopFlag = 1; }
using namespace std;
typedef uint64_t u64;
typedef unsigned __int128 M;
static inline int ctz128(M x) { u64 lo = (u64)x; return lo ? __builtin_ctzll(lo) : 64 + __builtin_ctzll((u64)(x >> 64)); }
static inline int pc128(M x) { return __builtin_popcountll((u64)x) + __builtin_popcountll((u64)(x >> 64)); }
#define BIT(p) (((M)1) << (p))

int v, k, t, m, b;
u64 Cb[130][24];
inline u64 binom(int n, int r) { return (r < 0 || n < r || n < 0) ? 0 : Cb[n][r]; }

vector<uint16_t> ct;            // per t-set
vector<int64_t> A, Bc;          // per t-set; A holds sum of base weights of uncovered supersets
vector<int32_t> U;              // per t-set: number of uncovered m-sets containing it
int64_t E = 0;                  // global bump counter: an uncovered m-set weighs w[M] + E
int64_t delta = 1;
#ifdef PROF
double tScan = 0, tScore = 0, tApply = 0; unsigned long long nVisit = 0, n01 = 0, n10 = 0, n12 = 0, n21 = 0;
#define PC(x) x
#else
#define PC(x)
#endif
static inline double nowS() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }              // bump per step; grows geometrically per level (= decaying older weights)
vector<uint16_t> cnt;           // per m-set
vector<int64_t> w;              // per m-set: weight if covered, base (weight - E) if uncovered
vector<uint32_t> unc; vector<M> uncM; vector<int32_t> upos;
vector<M> bm;
vector<int> tabuUntil;
long step = 0;

inline void add_unc(uint32_t r, M x) { upos[r] = unc.size(); unc.push_back(r); uncM.push_back(x); }
inline void rem_unc(uint32_t r) { int32_t p = upos[r]; uint32_t l = unc.back(); unc[p] = l; uncM[p] = uncM.back(); upos[l] = p; unc.pop_back(); uncM.pop_back(); upos[r] = -1; }
inline int pts_of(M x, int* p) { int n = 0; for (; x; x &= x - 1) p[n++] = ctz128(x); return n; }

// chunked colex rank for m-sets
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
inline uint32_t rankM(M x) {
    uint32_t r = 0; int off = 0;
    for (int c = 0; c < NCH; ++c) { unsigned bits = (unsigned)(x >> (16 * c)) & 0xffff; if (bits) { r += RT[((size_t)c * 65536 + bits) * (m + 1) + off]; off += __builtin_popcount(bits); } }
    return r;
}
inline u64 rankT(M x) { u64 r = 0; int i = 1; while (x) { int p = ctz128(x); r += binom(p, i); ++i; x &= x - 1; } return r; }
inline M unrankM(u64 r) {
    M x = 0;
    for (int i = m; i >= 1; --i) { int c = i - 1; while (binom(c + 1, i) <= r) ++c; x |= BIT(c); r -= binom(c, i); }
    return x;
}

// next m-set in colex order (= next larger mask with the same popcount), Gosper's hack on 128 bits
static inline M next_comb(M x) { M c = x & (~x + 1); M r = x + c; return r | (((r ^ x) >> 2) >> ctz128(c)); }
template <class F> inline void subsets(int n, int r, F f) {
    if (r > n || r < 0) return;
    int idx[24]; for (int i = 0; i < r; ++i) idx[i] = i;
    while (true) {
        f(idx);
        if (r == 0) return;
        int i = r - 1; while (i >= 0 && idx[i] == n - r + i) --i;
        if (i < 0) return;
        ++idx[i]; for (int l = i + 1; l < r; ++l) idx[l] = idx[l - 1] + 1;
    }
}
// t-subsets of m-set mask: f(rankT). Each t-subset is M minus R = m-t positions; its colex rank is a sum of R+1
// segment sums of S_d[j] = sum_{i<j} binom(P[i], i+1-d), so each subset costs O(R) instead of O(t) bit scans.
int Rm;                           // m - t
vector<uint32_t> LU;              // local index (colex rank of removed positions) -> bitmask of removed positions
template <class F> inline void tsubs(M x, F f) {
    int P[24]; pts_of(x, P);
    u64 S[14][25];
    for (int d = 0; d <= Rm; ++d) { S[d][0] = 0; for (int j = 0; j < m; ++j) S[d][j + 1] = S[d][j] + binom(P[j], j + 1 - d); }
    subsets(m, Rm, [&](int* rm) {
        u64 r = 0; int prev = 0;
        for (int d = 0; d < Rm; ++d) { r += S[d][rm[d]] - S[d][prev]; prev = rm[d] + 1; }
        r += S[Rm][m] - S[Rm][prev];
        f(r);
    });
}
// witness byte: sw[M] = sum (mod 256) of the local indices of the covered t-subsets of M; when cnt[M] == 1 it is
// exactly the local index of the single covered t-subset, so Bc needs one update instead of C(m,t)
vector<uint8_t> sw;
inline int local_idx(M x, M T) {  // colex rank of the positions of M \ T inside M
    int idx = 0, i = 0;
    for (M c = x & ~T; c; c &= c - 1, ++i) { int p = ctz128(c); int pos = pc128(x & (BIT(p) - 1)); idx += (int)binom(pos, i + 1); }
    return idx;
}
inline u64 unique_rank(M x, int li) {  // rankT of the t-subset of M with local index li
    uint32_t rmask = LU[li]; M T = 0; int pos = 0;
    for (M c = x; c; c &= c - 1, ++pos) if (!((rmask >> pos) & 1)) T |= BIT(ctz128(c));
    return rankT(T);
}
// m-sets containing t-set mask T: f(rankM, mask)
M FULL;
template <class F> inline void msupers(M T, F f) {
    int O[128]; int n = pts_of(FULL & ~T, O);
    subsets(n, m - t, [&](int* id) { M x = T; for (int i = 0; i < m - t; ++i) x |= BIT(O[id[i]]); f(rankM(x), x); });
}

// m-sets containing t-set T, enumerated in parallel chunks (each m-set exactly once, so per-m-set state needs no
// atomics; shared t-set scores use atomic adds and the uncovered list a critical section)
template <class F> inline void msupers_par(M T, F f) {
    int O[128]; int n = pts_of(FULL & ~T, O); const int R = m - t;
    u64 N = binom(n, R); const u64 CH = 128;
    #pragma omp parallel for schedule(dynamic, 1) if (N >= 1024)
    for (u64 c0 = 0; c0 < N; c0 += CH) {
        int id[24]; u64 rr = c0;
        for (int i = R; i >= 1; --i) { int c = i - 1; while (binom(c + 1, i) <= rr) ++c; id[i - 1] = c; rr -= binom(c, i); }
        u64 c1 = min(N, c0 + CH);
        for (u64 q = c0; q < c1; ++q) {
            M x = T; for (int i = 0; i < R; ++i) x |= BIT(O[id[i]]);
            f(rankM(x), x);
            if (R > 0) { int i = 0; while (i < R - 1 && id[i] + 1 == id[i + 1]) { id[i] = i; ++i; } ++id[i]; }  // next colex
        }
    }
}
#define AADD(a, v) __atomic_fetch_add(&(a), (v), __ATOMIC_RELAXED)
inline void flip_on(M T, u64 rT) {   // t-set T (rank rT) becomes covered
    msupers_par(T, [&](uint32_t r, M x) {
        uint16_t old = cnt[r]++; int li = local_idx(x, T); PC(++nVisit; if (old == 0) ++n01; else if (old == 1) ++n12;)
        if (old == 0) {
            { std::lock_guard<std::mutex> g(uncMu); rem_unc(r); }
            int64_t base = w[r]; int64_t ww = base + E; w[r] = ww;
            tsubs(x, [&](u64 u) { AADD(A[u], -base); AADD(U[u], -1); }); AADD(Bc[rT], ww);
        }
        else if (old == 1) { AADD(Bc[unique_rank(x, sw[r])], -w[r]); }
        sw[r] = (uint8_t)(sw[r] + li);
    });
}
inline void flip_off(M T, u64 rT) {  // t-set T (rank rT) becomes uncovered
    msupers_par(T, [&](uint32_t r, M x) {
        uint16_t old = cnt[r]--; int li = local_idx(x, T); PC(++nVisit; if (old == 1) ++n10; else if (old == 2) ++n21;)
        sw[r] = (uint8_t)(sw[r] - li);
        if (old == 1) {
            { std::lock_guard<std::mutex> g(uncMu); add_unc(r, x); }
            int64_t ww = w[r]; int64_t base = ww - E; w[r] = base;
            tsubs(x, [&](u64 u) { AADD(A[u], base); AADD(U[u], 1); }); AADD(Bc[rT], -ww);
        }
        else if (old == 2) { AADD(Bc[unique_rank(x, sw[r])], w[r]); }
    });
}
// t-sets of block B containing point x: f(mask, rank)
template <class F> inline void tsets_with(M B, int x, F f) {
    int P[128]; int n = pts_of(B & ~BIT(x), P);
    subsets(n, t - 1, [&](int* id) { M s = BIT(x); for (int i = 0; i < t - 1; ++i) s |= BIT(P[id[i]]); f(s, rankT(s)); });
}
template <class F> inline void tsets_all(M B, F f) {
    int P[128]; int n = pts_of(B, P);
    subsets(n, t, [&](int* id) { M s = 0; for (int i = 0; i < t; ++i) s |= BIT(P[id[i]]); f(s, rankT(s)); });
}
inline void block_remove_point(M B, int x) { tsets_with(B, x, [&](M s, u64 r) { if (--ct[r] == 0) flip_off(s, r); }); }
inline void block_add_point(M B2, int x) { tsets_with(B2, x, [&](M s, u64 r) { if (ct[r]++ == 0) flip_on(s, r); }); }


// exact search for one k-block meeting every m-set in Us in >= t points (each m-set needs one of its
// t-subsets inside the block); returns true and the block points in `out`
long dfsNodes;
bool dfs_cover(const vector<M>& Us, M Bk, M& out) {
    if (++dfsNodes > 200000) return false;
    const M* first = nullptr;
    for (const M& x : Us) if (pc128(x & Bk) < t) { first = &x; break; }
    if (!first) { out = Bk; return true; }
    int P[24]; int n = pts_of(*first, P);
    bool found = false;
    subsets(n, t, [&](int* id) {
        if (found) return;
        M s = Bk; for (int i = 0; i < t; ++i) s |= BIT(P[id[i]]);
        if (pc128(s) <= k && s != Bk && dfs_cover(Us, s, out)) found = true;
    });
    return found;
}
inline void remove_block(M B) { tsets_all(B, [&](M s, u64 r) { if (--ct[r] == 0) flip_off(s, r); }); }
inline void add_block(M B) { tsets_all(B, [&](M s, u64 r) { if (ct[r]++ == 0) flip_on(s, r); }); }

void rebuild_scores() {  // recompute A and Bc from cnt and w, in parallel over chunks of m-sets
    fill(A.begin(), A.end(), 0); fill(Bc.begin(), Bc.end(), 0); fill(U.begin(), U.end(), 0);
    u64 NM = cnt.size(); const u64 CH = 1 << 20;
    #pragma omp parallel for schedule(dynamic, 1)
    for (u64 c0 = 0; c0 < NM; c0 += CH) {
        u64 c1 = min(NM, c0 + CH); M x = unrankM(c0);
        for (u64 r = c0; r < c1; ++r, x = next_comb(x)) {
            if (cnt[r] > 1) continue;
            int64_t ww = w[r];
            if (cnt[r] == 0) tsubs(x, [&](u64 u) { __atomic_fetch_add(&A[u], ww, __ATOMIC_RELAXED); __atomic_fetch_add(&U[u], 1, __ATOMIC_RELAXED); });
            else __atomic_fetch_add(&Bc[unique_rank(x, sw[r])], ww, __ATOMIC_RELAXED);
        }
    }
}

void save(const char* fn) {
    FILE* f = fopen(fn, "w");
    for (M x : bm) { int c = 0; for (int p = 0; p < v; ++p) if ((x >> p) & 1) fprintf(f, "%d%c", p + 1, ++c < k ? ' ' : '\n'); }
    fclose(f);
}

int main(int argc, char** argv) {
    if (argc < 8) { fprintf(stderr, "usage: wlsm2 v k t m b seconds seed [init|-] [cands=12] [tenure=5] [targetB]\n"); return 1; }
    v = atoi(argv[1]); k = atoi(argv[2]); t = atoi(argv[3]); m = atoi(argv[4]); b = atoi(argv[5]);
    double secs = atof(argv[6]); unsigned seed = atoi(argv[7]);
    const char* init = (argc > 8 && argv[8][0] != '-') ? argv[8] : nullptr;
    int ncand = argc > 9 ? atoi(argv[9]) : 12;
    int tenure = argc > 10 ? atoi(argv[10]) : 5;
    int targetB = argc > 11 ? atoi(argv[11]) : b;
    for (int n = 0; n < 130; ++n) { Cb[n][0] = 1; for (int r = 1; r < 24; ++r) Cb[n][r] = n ? Cb[n - 1][r - 1] + Cb[n - 1][r] : 0; }
    FULL = (v < 128) ? (BIT(v) - 1) : ~(M)0;
    u64 NM = binom(v, m), NT = binom(v, t);
    if (NM > 2000000000ULL) { fprintf(stderr, "C(v,m) too large\n"); return 1; }
    {   // cost guard: one t-set flip touches C(v-t,m-t) m-sets and each 0<->1 transition updates C(m,t) t-set scores
        double flipCost = (double)binom(v - t, m - t) * (double)binom(m, t);
        double maxFlip = getenv("MAXFLIP") ? atof(getenv("MAXFLIP")) : 1e8;
        if (flipCost > maxFlip) { fprintf(stderr, "TOO_EXPENSIVE: flip cost %.3g > %.3g, use zerocov instead\n", flipCost, maxFlip); return 3; }
    }
    // the witness byte holds a local index mod 256, so it names the unique covered t-subset only when C(m,t) <= 256
    if (binom(m, t) > 256) { fprintf(stderr, "TOO_EXPENSIVE: C(m,t) = %llu > 256 overflows the witness byte, use zerocov instead\n", (unsigned long long)binom(m, t)); return 3; }
    build_rank_tables();
    Rm = m - t; LU.assign(binom(m, Rm), 0);
    subsets(m, Rm, [&](int* id) { int li = 0; uint32_t mk = 0; for (int i = 0; i < Rm; ++i) { li += (int)binom(id[i], i + 1); mk |= 1u << id[i]; } LU[li] = mk; });
    ct.assign(NT, 0); A.assign(NT, 0); Bc.assign(NT, 0); U.assign(NT, 0);
    cnt.assign(NM, 0); w.assign(NM, 1); upos.assign(NM, -1); sw.assign(NM, 0);  // E = 0 at start, so base = weight
    mt19937_64 rng(seed);
    signal(SIGTERM, on_signal); signal(SIGINT, on_signal);
    bm.assign(b, 0); int loaded = 0;
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
    // (re)build the whole state from the blocks in bm, with all weights reset
    auto init_state = [&]() {
        fill(ct.begin(), ct.end(), 0); fill(cnt.begin(), cnt.end(), 0); fill(sw.begin(), sw.end(), 0); E = 0;
        unc.clear(); uncM.clear();
        // count coverage first, then list only the m-sets that are actually uncovered (listing all of them first
        // costs ~20 bytes per m-set, too much for C(v,m) near 1e9)
        for (u64 r = 0; r < NM; ++r) { w[r] = 1; upos[r] = -1; }
        vector<M> newT;  // distinct covered t-sets, in the order first seen
        for (int i = 0; i < b; ++i) tsets_all(bm[i], [&](M s, u64 r) { if (ct[r]++ == 0) newT.push_back(s); });
        #pragma omp parallel for schedule(dynamic, 256)
        for (size_t j = 0; j < newT.size(); ++j) msupers(newT[j], [&](uint32_t q, M x) { __atomic_fetch_add(&cnt[q], 1, __ATOMIC_RELAXED); __atomic_fetch_add(&sw[q], (uint8_t)local_idx(x, newT[j]), __ATOMIC_RELAXED); });
        { M x = 0; for (u64 r = 0; r < NM; ++r, x = r < NM ? next_comb(x) : x) { if (r == 0) x = unrankM(0); if (cnt[r] == 0) add_unc((uint32_t)r, x); } }
        rebuild_scores();
        tabuUntil.assign((size_t)b * 128, 0);
    };
    init_state();
    long levelSteps = getenv("LEVEL_STEPS") ? atol(getenv("LEVEL_STEPS")) : 0;  // 0 = never restart a level
    int delPool = getenv("DEL_POOL") ? atoi(getenv("DEL_POOL")) : 20;        // deletion choices to cycle through
    // undo log of the moves since the last full cover: 'd' = block i (mask a) deleted, 's' = block i changed a -> c
    struct Op { char ty; int i; M a, c; }; vector<Op> undo; int fullB = -1;
    auto fullCoverCopy = [&]() { vector<M> cp = bm;   // the last full cover, rebuilt from bm by undoing the log
        for (size_t q = undo.size(); q-- > 0;) { const Op& o = undo[q];
            if (o.ty == 's') cp[o.i] = o.a; else { if (o.i == (int)cp.size()) cp.push_back(o.a); else { cp.push_back(cp[o.i]); cp[o.i] = o.a; } } }
        return cp; };
    vector<M> fullCover; long levelStart = 0; int delTry = 0; vector<int> delOrder;
    double stallSecs = getenv("STALL") ? atof(getenv("STALL")) : 0, lastFullTime = 0;
    int pairPool = getenv("PAIRX") ? atoi(getenv("PAIRX")) : 0;   // blocks considered for two-for-one exchange
    if ((getenv("POOL") || getenv("LAZY")) && (pairPool > 0 || levelSteps > 0)) {  // inconsistent with the pool/heap indices
        fprintf(stderr, "PAIRX/LEVEL_STEPS ignored with POOL/LAZY\n"); pairPool = 0; levelSteps = 0;
    }
    int maxU = getenv("MAXU") ? atoi(getenv("MAXU")) : 5000;
    bool noReset = getenv("NORESET") != nullptr;  // keep breakout weights between shrink levels
    double saveSecs = getenv("SAVE_SECS") ? atof(getenv("SAVE_SECS")) : 0, lastSave = -1e9, lastStatus = 0; int lastFullB = -1, lastSavedB = INT32_MAX;
    long best = unc.size();
    fprintf(stderr, "init: loaded %d blocks, C(v,m)=%llu, uncovered %ld\n", loaded, (unsigned long long)NM, best);
    auto start = chrono::steady_clock::now();
    int bestFull = -1;
    vector<int> cand;
    // LAZY: deletion candidates in a min-heap of (score, index, mask); stale entries are skipped or re-scored on pop
    bool lazy = getenv("LAZY") != nullptr;
    struct HE { int64_t s; int i; M x; bool operator<(const HE& o) const { return s > o.s; } };
    vector<HE> heap; auto hpush = [&](int64_t sc, int i) { heap.push_back({sc, i, bm[i]}); push_heap(heap.begin(), heap.end()); };
    auto bscore = [&](int i) { int64_t l = 0; tsets_all(bm[i], [&](M, u64 r) { if (ct[r] == 1) l += Bc[r]; }); return l; };
    int64_t bumpGrowth = getenv("BUMP_GROWTH") ? atoll(getenv("BUMP_GROWTH")) : 0;  // per mille per level
    double lastProgress = 0; int progB = INT32_MAX;  // time of the last new best full cover (reverts do not count)
    long fullCount = 0; int resetEvery = getenv("RESET_EVERY") ? atoi(getenv("RESET_EVERY")) : 0;
    double revertSecs = getenv("REVERT_SECS") ? atof(getenv("REVERT_SECS")) : 0; int reverts = 0;
    long levelBudget = getenv("LEVEL_BUDGET") ? atol(getenv("LEVEL_BUDGET")) : 0;   // max steps to repair one deletion
    double growth = getenv("GROWTH") ? atof(getenv("GROWTH")) : 0;                    // revert if uncovered > growth x post-deletion
    size_t postDelUnc = 0;
    auto reset_weights = [&]() {  // all covered m-sets weigh 1, uncovered ones 1 (base 1 - E); parallel rebuild
        #pragma omp parallel for schedule(static)
        for (u64 r = 0; r < NM; ++r) w[r] = cnt[r] ? 1 : 1 - E;
        rebuild_scores();
    };
    int poolSize = getenv("POOL") ? atoi(getenv("POOL")) : 0, poolRefresh = getenv("POOL_REFRESH") ? atoi(getenv("POOL_REFRESH")) : 500, sinceRefresh = 0;
    struct MH { size_t operator()(const M& x) const { return std::hash<u64>()((u64)x ^ (u64)(x >> 64) * 0x9E3779B97F4A7C15ULL); } };
    vector<M> pool; vector<int64_t> poolS; std::unordered_map<M, int, MH> posOf;
    while (true) {
        if (unc.empty()) {
            bestFull = b;
            double el = chrono::duration<double>(chrono::steady_clock::now() - start).count();
            if (el - lastSave >= saveSecs || b <= targetB) {  // writing a huge cover after every removal is slow
                char fn[256]; sprintf(fn, "mcover_%d_%d_%d_%d_%d_%u.txt", v, k, t, m, b, seed); save(fn); lastSave = el; lastSavedB = b; }
            lastFullB = b;
            if (fullCount % 100 == 0 || b <= targetB) { printf("FULL %d blocks at %.1fs (%ld steps)\n", b, el, step); fflush(stdout); }
            lastFullTime = el; if (b < progB) { progB = b; lastProgress = el; }
            if (b <= targetB) break;
            // reset weights, then delete the block whose unique t-sets carry the least Bc (= m-sets it alone covers)
            ++fullCount;
            if (bumpGrowth > 0) delta += max<int64_t>(1, delta * bumpGrowth / 1000);
            if (!noReset || (resetEvery > 0 && fullCount % resetEvery == 0) || delta > (1LL << 24)) { reset_weights(); delta = 1; }
            undo.clear(); fullB = b; fullCover.assign(1, 0); delTry = 0;  // fullCover is only a 'have one' marker now
            if (poolSize > 0) {
              // POOL: keep the poolSize lowest-scoring blocks (by mask); re-score only those after each removal and
              // refresh the pool from all blocks (in parallel) every poolRefresh removals
              if (pool.empty() || ++sinceRefresh >= poolRefresh) {
                  vector<pair<int64_t,int>> sc(b);
                  #pragma omp parallel for schedule(static)
                  for (int i = 0; i < b; ++i) sc[i] = {bscore(i), i};
                  int P = min(b, poolSize); nth_element(sc.begin(), sc.begin() + (P - 1), sc.end()); sort(sc.begin(), sc.begin() + P);
                  pool.clear(); poolS.clear(); for (int q = 0; q < P; ++q) { pool.push_back(bm[sc[q].second]); poolS.push_back(sc[q].first); }  // ascending stale score
                  posOf.clear(); for (int i = 0; i < b; ++i) posOf[bm[i]] = i;
                  sinceRefresh = 0;
              }
              int pick = -1; int64_t bestS = INT64_MAX; size_t pq = 0;
              for (size_t q = 0; q < pool.size(); ++q) {
                  if (poolS[q] >= bestS) break;   // pool is sorted by stale score: the rest are (approximately) no better
                  auto it = posOf.find(pool[q]);
                  if (it == posOf.end()) { pool.erase(pool.begin() + q); poolS.erase(poolS.begin() + q); --q; continue; }
                  int64_t sc = bscore(it->second); poolS[q] = sc;
                  if (sc < bestS) { bestS = sc; pick = it->second; pq = q; if (sc == 0) break; }  // 0 is the minimum: stop
              }
              if (pick >= 0) {   // keep the pool sorted: move the picked entry's updated neighbours into place lazily
                  for (size_t q = 1; q < pool.size() && q < 64; ++q) if (poolS[q] < poolS[q - 1]) { swap(poolS[q], poolS[q - 1]); swap(pool[q], pool[q - 1]); }
              }
              (void)pq;
              if (pick < 0) { sinceRefresh = poolRefresh; pick = 0; }
              delOrder.assign(1, pick);
            }
            else if (!lazy) { vector<pair<int64_t,int>> sc; for (int i = 0; i < b; ++i) { int64_t l = 0; tsets_all(bm[i], [&](M, u64 r) { if (ct[r] == 1) l += Bc[r]; }); sc.push_back({l, i}); }
              sort(sc.begin(), sc.end()); delOrder.clear(); for (auto& p : sc) delOrder.push_back(p.second); }
            else {
              if (heap.empty()) { heap.reserve(2 * (size_t)b); for (int i = 0; i < b; ++i) heap.push_back({bscore(i), i, bm[i]}); make_heap(heap.begin(), heap.end()); }
              int pick = -1;
              while (!heap.empty()) {
                  HE h = heap.front(); pop_heap(heap.begin(), heap.end()); heap.pop_back();
                  if (h.i >= b || bm[h.i] != h.x) continue;            // block changed or moved: stale entry
                  int64_t sc = bscore(h.i);
                  if (heap.empty() || sc <= heap.front().s) { pick = h.i; break; }
                  hpush(sc, h.i);
              }
              if (pick < 0) { for (int i = 0; i < b; ++i) hpush(bscore(i), i); pick = 0; }
              delOrder.assign(1, pick);
            }
            if (pairPool > 0) {  // block exchange: remove two blocks, look for one block covering what they leave
                bool improved = false; int P = min(b, pairPool); long pairsTried = 0;
                for (int a = 0; a < P && !improved; ++a) for (int c = a + 1; c < P && !improved; ++c) {
                    int i = delOrder[a], j = delOrder[c];
                    remove_block(bm[i]); remove_block(bm[j]); ++pairsTried;
                    M nb = 0; bool ok = false;
                    if ((int)unc.size() <= maxU) { vector<M> Us(uncM.begin(), uncM.end()); dfsNodes = 0; ok = dfs_cover(Us, 0, nb); }
                    if (ok) {
                        vector<int> free_; for (int p = 0; p < v; ++p) if (!((nb >> p) & 1)) free_.push_back(p);
                        shuffle(free_.begin(), free_.end(), rng);
                        for (int q = 0; pc128(nb) < k; ++q) nb |= BIT(free_[q]);
                        add_block(nb);
                        vector<M> nbm; for (int q = 0; q < b; ++q) if (q != i && q != j) nbm.push_back(bm[q]);
                        nbm.push_back(nb); bm = nbm; b = bm.size(); improved = true;
                    } else { add_block(bm[i]); add_block(bm[j]); }
                }
                fprintf(stderr, "pair exchange at b=%d: %s after %ld pairs\n", fullB, improved ? "SUCCESS" : "none", pairsTried);
                if (improved) continue;
            }
            int bi = delOrder[0];
            levelStart = step;
            M delM = bm[bi]; undo.push_back({'d', bi, delM, 0});
            tsets_all(bm[bi], [&](M s, u64 r) { if (--ct[r] == 0) flip_off(s, r); });
            bm[bi] = bm[b - 1]; bm.pop_back();
            memcpy(&tabuUntil[(size_t)bi * 128], &tabuUntil[(size_t)(b - 1) * 128], 128 * sizeof(int));
            --b; best = unc.size(); postDelUnc = unc.size();
            if (lazy && bi < b) hpush(INT64_MIN / 2, bi);  // re-scored when popped
            if (poolSize > 0) { posOf.erase(delM); if (bi < b) posOf[bm[bi]] = bi; }
            continue;
        }
        if (levelSteps > 0 && !fullCover.empty() && step - levelStart > levelSteps) {
            ++delTry;
            int pick = delTry < (int)delOrder.size() && delTry < delPool ? delOrder[delTry] : delOrder[rng() % delOrder.size()];
            bm = fullCoverCopy(); undo.clear(); b = bm.size(); init_state();
            tsets_all(bm[pick], [&](M s, u64 r) { if (--ct[r] == 0) flip_off(s, r); });
            bm[pick] = bm[b - 1]; bm.pop_back(); --b;
            levelStart = step;
        }
        if (true) {   // checked every step: one step can take long on expensive instances
            double el = chrono::duration<double>(chrono::steady_clock::now() - start).count();
            if (el > secs || stopFlag) break;
            if (el - lastStatus >= 30) { lastStatus = el; fprintf(stderr, "status %.0fs: b=%d uncovered=%zu steps=%ld\n", el, b, unc.size(), step);
#ifdef PROF
                { double S = (double)max(step, 1L); fprintf(stderr, "  prof: scan %.1fs score %.1fs apply %.1fs | per step: visits %.0f 0to1 %.0f 1to0 %.0f 1to2 %.0f 2to1 %.0f\n", tScan, tScore, tApply, nVisit / S, n01 / S, n10 / S, n12 / S, n21 / S); }
#endif
            }
            if (stallSecs > 0 && el - lastProgress > stallSecs) break;  // no shrink for stallSecs: give up
            bool stalled = (revertSecs > 0 && el - lastFullTime > revertSecs) || (levelBudget > 0 && step - levelStart > levelBudget)
                        || (growth > 0 && (double)unc.size() > max(1000.0, growth * (double)postDelUnc));
            if (stalled && !fullCover.empty()) {
                // stalled: go back to the last full cover (apply the difference block by block), reset the weights,
                // and delete a different low-score block instead
                for (size_t q = undo.size(); q-- > 0;) {   // undo every move since the last full cover
                    const Op& o = undo[q];
                    if (o.ty == 's') {
                        int xin = ctz128(o.c & ~o.a), xout = ctz128(o.a & ~o.c);
                        block_remove_point(o.c, xin); bm[o.i] = o.a; block_add_point(o.a, xout);
                        if (poolSize > 0) { posOf.erase(o.c); posOf[o.a] = o.i; }
                    } else {
                        if (o.i == b) bm.push_back(o.a);
                        else { bm.push_back(bm[o.i]); memcpy(&tabuUntil[(size_t)b * 128], &tabuUntil[(size_t)o.i * 128], 128 * sizeof(int)); bm[o.i] = o.a;
                               if (poolSize > 0) posOf[bm[b]] = b; }
                        memset(&tabuUntil[(size_t)o.i * 128], 0, 128 * sizeof(int));
                        ++b; add_block(o.a); if (poolSize > 0) posOf[o.a] = o.i;
                    }
                }
                undo.clear();
                if (getenv("REVERT_RESET")) reset_weights();
                vector<pair<int64_t,int>> sc(b);
                #pragma omp parallel for schedule(static)
                for (int i = 0; i < b; ++i) sc[i] = {bscore(i), i};
                ++reverts; int P = min(b, 64); nth_element(sc.begin(), sc.begin() + (P - 1), sc.end()); sort(sc.begin(), sc.begin() + P);
                int bi = sc[(reverts) % P].second;   // skip the block whose deletion stalled
                M delM = bm[bi]; undo.push_back({'d', bi, delM, 0});
                remove_block(bm[bi]); bm[bi] = bm[b - 1]; memcpy(&tabuUntil[(size_t)bi * 128], &tabuUntil[(size_t)(b - 1) * 128], 128 * sizeof(int)); bm.pop_back(); --b;
                if (poolSize > 0) { posOf.erase(delM); if (bi < b) posOf[bm[bi]] = bi; }
                if (lazy) heap.clear();
                fprintf(stderr, "revert %d at %.0fs: back to %d blocks, deleted rank-%d block, uncovered=%zu\n", reverts, el, fullB, reverts % P, unc.size());
                lastFullTime = el; levelStart = step; postDelUnc = unc.size();
                continue;  // the deletion may leave nothing uncovered: let the loop top handle a full cover
            }
        }
        ++step;
#ifdef CHECK
        if (step % 2000 == 0) {
            vector<uint16_t> ct2(NT, 0); for (int i = 0; i < b; ++i) tsets_all(bm[i], [&](M, u64 r) { ct2[r]++; });
            for (u64 r = 0; r < NT; ++r) if (ct2[r] != ct[r]) { fprintf(stderr, "CHECK ct mismatch\n"); exit(2); }
            vector<int64_t> A2(NT, 0), B2(NT, 0);
            for (u64 r = 0; r < NM; ++r) { M x = unrankM(r); int c = 0; tsubs(x, [&](u64 u) { if (ct[u]) ++c; });
                if (c != cnt[r]) { fprintf(stderr, "CHECK cnt mismatch\n"); exit(2); }
                if (c == 0) tsubs(x, [&](u64 u) { A2[u] += w[r] + E; }); else if (c == 1) tsubs(x, [&](u64 u) { if (ct[u]) B2[u] += w[r]; }); }
            for (u64 r = 0; r < NT; ++r) if (A2[r] != A[r] + E * U[r] || B2[r] != Bc[r]) { fprintf(stderr, "CHECK A/Bc mismatch at step %ld\n", step); exit(2); }
            if (step % 20000 == 0) fprintf(stderr, "CHECK ok step %ld\n", step);
        }
#endif
        PC(double t0 = nowS();)
        uint32_t rT = unc[rng() % unc.size()];
        for (int s = 0; s < 3; ++s) { uint32_t r2 = unc[rng() % unc.size()]; if (w[r2] > w[rT]) rT = r2; }
        M T = uncM[upos[rT]];
        cand.clear();
        if (b >= 20000) {  // parallel scan, per-thread lists concatenated in index order
            static vector<vector<int>> tl; int nt = omp_get_max_threads(); if ((int)tl.size() < nt) tl.resize(nt);
            #pragma omp parallel
            {
                vector<int>& my = tl[omp_get_thread_num()]; my.clear();
                #pragma omp for schedule(static) nowait
                for (int i = 0; i < b; ++i) if (pc128(bm[i] & T) == t - 1) my.push_back(i);
            }
            for (int q = 0; q < nt; ++q) cand.insert(cand.end(), tl[q].begin(), tl[q].end());
        } else for (int i = 0; i < b; ++i) if (pc128(bm[i] & T) == t - 1) cand.push_back(i);
        if (cand.empty()) {  // no block is one point away: pull a random block towards T
            int i = rng() % b; M B = bm[i];
            int a[24]; int na = pts_of(T & ~B, a); int o[24]; int no = pts_of(B & ~T, o);
            if (!na || !no) continue;
            int xin = a[rng() % na], xout = o[rng() % no];
            block_remove_point(B, xout); M B2 = (B & ~BIT(xout)) | BIT(xin); bm[i] = B2; block_add_point(B2, xin);
            undo.push_back({'s', i, B, B2});
            if (lazy) hpush(INT64_MIN / 2, i);
            if (poolSize > 0) { posOf.erase(B); posOf[B2] = i; }
            continue;
        }
        if ((int)cand.size() > ncand) { shuffle(cand.begin(), cand.end(), rng); cand.resize(ncand); }
        PC(double t1 = nowS(); tScan += t1 - t0;)
        int64_t bestGain = INT64_MIN; int bi = -1, bout = -1, bin_ = -1; long nties = 0;
        u64 stepSeed = rng();
        #pragma omp parallel if (cand.size() >= 48)
        {
            static thread_local int64_t g1[128][128]; int64_t G[128], L[128];
            mt19937_64 trng(stepSeed + 0x9E3779B97F4A7C15ULL * (omp_get_thread_num() + 1));
            int64_t myBest = INT64_MIN; int mi = -1, mo = -1, mn = -1; long mt = 0;
            #pragma omp for schedule(dynamic, 4)
            for (size_t ci = 0; ci < cand.size(); ++ci) {
                int i = cand[ci]; M B = bm[i];
                int a[24]; int na = pts_of(T & ~B, a); int o[24]; int no = pts_of(B & ~T, o);
                for (int io = 0; io < no; ++io) {  // loss of removing xout
                    int xo = o[io]; int64_t l = 0;
                    tsets_with(B, xo, [&](M, u64 r) { if (ct[r] == 1) l += Bc[r]; });
                    L[xo] = l;
                }
                int Bp[24]; int nb = pts_of(B, Bp);
                for (int ia = 0; ia < na; ++ia) {  // gain of adding xin, minus the part using xout
                    int xi = a[ia]; G[xi] = 0; for (int j = 0; j < nb; ++j) g1[xi][Bp[j]] = 0;
                    subsets(nb, t - 1, [&](int* id) {
                        M s = BIT(xi); for (int q = 0; q < t - 1; ++q) s |= BIT(Bp[id[q]]);
                        u64 r = rankT(s);
                        if (ct[r] == 0) { int64_t aa = A[r] + E * U[r]; G[xi] += aa; for (int q = 0; q < t - 1; ++q) g1[xi][Bp[id[q]]] += aa; }
                    });
                }
                for (int ia = 0; ia < na; ++ia) {
                    int xi = a[ia];
                    if (tabuUntil[(size_t)i * 128 + xi] > step) continue;
                    for (int io = 0; io < no; ++io) {
                        int xo = o[io];
                        int64_t g = (G[xi] - g1[xi][xo]) - L[xo];
                        if (g > myBest) { myBest = g; mi = i; mo = xo; mn = xi; mt = 1; }
                        else if (g == myBest && trng() % (++mt) == 0) { mi = i; mo = xo; mn = xi; }
                    }
                }
            }
            #pragma omp critical
            {
                if (mi >= 0) {
                    if (myBest > bestGain) { bestGain = myBest; bi = mi; bout = mo; bin_ = mn; nties = mt; }
                    else if (myBest == bestGain) { nties += mt; if ((long)(trng() % nties) < mt) { bi = mi; bout = mo; bin_ = mn; } }
                }
            }
        }
        PC(double t2 = nowS(); tScore += t2 - t1;)
        if (bi < 0) continue;
        M B = bm[bi];
        block_remove_point(B, bout);
        M B2 = (B & ~BIT(bout)) | BIT(bin_); bm[bi] = B2;
        block_add_point(B2, bin_);
        undo.push_back({'s', bi, B, B2});
        if (lazy) hpush(INT64_MIN / 2, bi);
        if (poolSize > 0) { posOf.erase(B); posOf[B2] = bi; }
        tabuUntil[(size_t)bi * 128 + bout] = step + tenure + (int)(rng() % 3);
        PC(tApply += nowS() - t2;)
        E += delta;  // bump every uncovered m-set's weight
        if ((long)unc.size() < best) best = unc.size();
    }
    if (bestFull > 0) printf("SHRINK result: full cover with %d blocks\n", bestFull);
    if (fullB > 0 && fullB < lastSavedB) {  // best full cover not yet on disk
        bm = fullCoverCopy(); char fn[256]; sprintf(fn, "mcover_%d_%d_%d_%d_%d_%u.txt", v, k, t, m, (int)bm.size(), seed); save(fn);
        fprintf(stderr, "saved %s\n", fn);
    }
    double el = chrono::duration<double>(chrono::steady_clock::now() - start).count();
    printf("WLSM2 C(%d,%d,%d,%d) b=%d: best %ld uncovered, %ld steps, %.1fs (%.0f steps/s)\n", v, k, t, m, b, best, step, el, step / max(el, 1e-9));
}
