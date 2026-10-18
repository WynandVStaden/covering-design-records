// wlsm2: weighted local search for (v,k,t,m) coverings with t-level move scores.
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
// Usage: wlsm2 v k t m b seconds seed [init|-] [cands=12] [tenure=5] [targetB]
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

vector<uint16_t> ct;            // per t-set
vector<int64_t> A, Bc;          // per t-set; A holds sum of base weights of uncovered supersets
vector<int32_t> U;              // per t-set: number of uncovered m-sets containing it
int64_t E = 0;                  // global bump counter: an uncovered m-set weighs w[M] + E
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
// t-subsets of m-set mask: f(rankT)
template <class F> inline void tsubs(M x, F f) {
    int P[24]; int n = pts_of(x, P);
    subsets(n, t, [&](int* id) { M s = 0; for (int i = 0; i < t; ++i) s |= BIT(P[id[i]]); f(rankT(s)); });
}
// m-sets containing t-set mask T: f(rankM, mask)
M FULL;
template <class F> inline void msupers(M T, F f) {
    int O[128]; int n = pts_of(FULL & ~T, O);
    subsets(n, m - t, [&](int* id) { M x = T; for (int i = 0; i < m - t; ++i) x |= BIT(O[id[i]]); f(rankM(x), x); });
}

inline void flip_on(M T) {   // t-set T becomes covered
    msupers(T, [&](uint32_t r, M x) {
        uint16_t old = cnt[r]++;
        if (old == 0) { rem_unc(r); int64_t base = w[r]; int64_t ww = base + E; w[r] = ww; tsubs(x, [&](u64 u) { A[u] -= base; U[u]--; Bc[u] += ww; }); }
        else if (old == 1) { int64_t ww = w[r]; tsubs(x, [&](u64 u) { Bc[u] -= ww; }); }
    });
}
inline void flip_off(M T) {  // t-set T becomes uncovered
    msupers(T, [&](uint32_t r, M x) {
        uint16_t old = cnt[r]--;
        if (old == 1) { add_unc(r, x); int64_t ww = w[r]; int64_t base = ww - E; w[r] = base; tsubs(x, [&](u64 u) { Bc[u] -= ww; A[u] += base; U[u]++; }); }
        else if (old == 2) { int64_t ww = w[r]; tsubs(x, [&](u64 u) { Bc[u] += ww; }); }
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
inline void block_remove_point(M B, int x) { tsets_with(B, x, [&](M s, u64 r) { if (--ct[r] == 0) flip_off(s); }); }
inline void block_add_point(M B2, int x) { tsets_with(B2, x, [&](M s, u64 r) { if (ct[r]++ == 0) flip_on(s); }); }


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
inline void remove_block(M B) { tsets_all(B, [&](M s, u64 r) { if (--ct[r] == 0) flip_off(s); }); }
inline void add_block(M B) { tsets_all(B, [&](M s, u64 r) { if (ct[r]++ == 0) flip_on(s); }); }

void rebuild_scores() {  // recompute A and Bc from cnt and w
    fill(A.begin(), A.end(), 0); fill(Bc.begin(), Bc.end(), 0); fill(U.begin(), U.end(), 0);
    u64 NM = cnt.size();
    for (u64 r = 0; r < NM; ++r) {
        if (cnt[r] > 1) continue;
        M x = unrankM(r); int64_t ww = w[r];
        if (cnt[r] == 0) tsubs(x, [&](u64 u) { A[u] += ww; U[u]++; }); else tsubs(x, [&](u64 u) { Bc[u] += ww; });
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
    if (NM > 300000000ULL) { fprintf(stderr, "C(v,m) too large\n"); return 1; }
    build_rank_tables();
    ct.assign(NT, 0); A.assign(NT, 0); Bc.assign(NT, 0); U.assign(NT, 0);
    cnt.assign(NM, 0); w.assign(NM, 1); upos.assign(NM, -1);  // E = 0 at start, so base = weight
    mt19937_64 rng(seed);
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
        fill(ct.begin(), ct.end(), 0); fill(cnt.begin(), cnt.end(), 0); E = 0;
        unc.clear(); uncM.clear();
        for (u64 r = 0; r < NM; ++r) { w[r] = 1; add_unc((uint32_t)r, unrankM(r)); }
        for (int i = 0; i < b; ++i) tsets_all(bm[i], [&](M s, u64 r) { if (ct[r]++ == 0) { msupers(s, [&](uint32_t q, M) { if (cnt[q]++ == 0) rem_unc(q); }); } });
        rebuild_scores();
        tabuUntil.assign((size_t)b * 128, 0);
    };
    init_state();
    long levelSteps = getenv("LEVEL_STEPS") ? atol(getenv("LEVEL_STEPS")) : 0;  // 0 = never restart a level
    int delPool = getenv("DEL_POOL") ? atoi(getenv("DEL_POOL")) : 20;        // deletion choices to cycle through
    vector<M> fullCover; long levelStart = 0; int delTry = 0; vector<int> delOrder;
    double stallSecs = getenv("STALL") ? atof(getenv("STALL")) : 0, lastFullTime = 0;
    int pairPool = getenv("PAIRX") ? atoi(getenv("PAIRX")) : 0;   // blocks considered for two-for-one exchange
    int maxU = getenv("MAXU") ? atoi(getenv("MAXU")) : 5000;
    long best = unc.size();
    fprintf(stderr, "init: loaded %d blocks, C(v,m)=%llu, uncovered %ld\n", loaded, (unsigned long long)NM, best);
    auto start = chrono::steady_clock::now();
    int bestFull = -1;
    vector<int> cand;
    static int64_t g1[128][128]; int64_t G[128], L[128];
    while (true) {
        if (unc.empty()) {
            bestFull = b;
            char fn[256]; sprintf(fn, "mcover_%d_%d_%d_%d_%d_%u.txt", v, k, t, m, b, seed); save(fn);
            double el = chrono::duration<double>(chrono::steady_clock::now() - start).count();
            printf("FULL %d blocks at %.1fs (%ld steps)\n", b, el, step); fflush(stdout);
            lastFullTime = el;
            if (b <= targetB) break;
            // reset weights, then delete the block whose unique t-sets carry the least Bc (= m-sets it alone covers)
            for (u64 r = 0; r < NM; ++r) w[r] = cnt[r] ? 1 : 1 - E; rebuild_scores();
            fullCover = bm; delTry = 0;
            { vector<pair<int64_t,int>> sc; for (int i = 0; i < b; ++i) { int64_t l = 0; tsets_all(bm[i], [&](M, u64 r) { if (ct[r] == 1) l += Bc[r]; }); sc.push_back({l, i}); }
              sort(sc.begin(), sc.end()); delOrder.clear(); for (auto& p : sc) delOrder.push_back(p.second); }
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
                fprintf(stderr, "pair exchange at b=%d: %s after %ld pairs\n", (int)fullCover.size(), improved ? "SUCCESS" : "none", pairsTried);
                if (improved) continue;
            }
            int bi = delOrder[0];
            levelStart = step;
            tsets_all(bm[bi], [&](M s, u64 r) { if (--ct[r] == 0) flip_off(s); });
            bm[bi] = bm[b - 1]; bm.pop_back();
            memcpy(&tabuUntil[(size_t)bi * 128], &tabuUntil[(size_t)(b - 1) * 128], 128 * sizeof(int));
            --b; best = unc.size();
            continue;
        }
        if (levelSteps > 0 && !fullCover.empty() && step - levelStart > levelSteps) {
            ++delTry;
            int pick = delTry < (int)delOrder.size() && delTry < delPool ? delOrder[delTry] : delOrder[rng() % delOrder.size()];
            bm = fullCover; b = bm.size(); init_state();
            tsets_all(bm[pick], [&](M s, u64 r) { if (--ct[r] == 0) flip_off(s); });
            bm[pick] = bm[b - 1]; bm.pop_back(); --b;
            levelStart = step;
        }
        if ((step & 255) == 0) {
            double el = chrono::duration<double>(chrono::steady_clock::now() - start).count();
            if (el > secs) break;
            if (stallSecs > 0 && el - lastFullTime > stallSecs) break;  // no shrink for stallSecs: give up
        }
        ++step;
#ifdef CHECK
        if (step % 2000 == 0) {
            vector<uint16_t> ct2(NT, 0); for (int i = 0; i < b; ++i) tsets_all(bm[i], [&](M, u64 r) { ct2[r]++; });
            for (u64 r = 0; r < NT; ++r) if (ct2[r] != ct[r]) { fprintf(stderr, "CHECK ct mismatch\n"); exit(2); }
            vector<int64_t> A2(NT, 0), B2(NT, 0);
            for (u64 r = 0; r < NM; ++r) { M x = unrankM(r); int c = 0; tsubs(x, [&](u64 u) { if (ct[u]) ++c; });
                if (c != cnt[r]) { fprintf(stderr, "CHECK cnt mismatch\n"); exit(2); }
                if (c == 0) tsubs(x, [&](u64 u) { A2[u] += w[r] + E; }); else if (c == 1) tsubs(x, [&](u64 u) { B2[u] += w[r]; }); }
            for (u64 r = 0; r < NT; ++r) if (A2[r] != A[r] + E * U[r] || B2[r] != Bc[r]) { fprintf(stderr, "CHECK A/Bc mismatch at step %ld\n", step); exit(2); }
            if (step % 20000 == 0) fprintf(stderr, "CHECK ok step %ld\n", step);
        }
#endif
        uint32_t rT = unc[rng() % unc.size()];
        for (int s = 0; s < 3; ++s) { uint32_t r2 = unc[rng() % unc.size()]; if (w[r2] > w[rT]) rT = r2; }
        M T = uncM[upos[rT]];
        cand.clear();
        for (int i = 0; i < b; ++i) if (pc128(bm[i] & T) == t - 1) cand.push_back(i);
        if (cand.empty()) {  // no block is one point away: pull a random block towards T
            int i = rng() % b; M B = bm[i];
            int a[24]; int na = pts_of(T & ~B, a); int o[24]; int no = pts_of(B & ~T, o);
            if (!na || !no) continue;
            int xin = a[rng() % na], xout = o[rng() % no];
            block_remove_point(B, xout); M B2 = (B & ~BIT(xout)) | BIT(xin); bm[i] = B2; block_add_point(B2, xin);
            continue;
        }
        if ((int)cand.size() > ncand) { shuffle(cand.begin(), cand.end(), rng); cand.resize(ncand); }
        int64_t bestGain = INT64_MIN; int bi = -1, bout = -1, bin_ = -1; int nties = 0;
        for (int i : cand) {
            M B = bm[i];
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
                    if (g > bestGain) { bestGain = g; bi = i; bout = xo; bin_ = xi; nties = 1; }
                    else if (g == bestGain && rng() % (++nties) == 0) { bi = i; bout = xo; bin_ = xi; }
                }
            }
        }
        if (bi < 0) continue;
        M B = bm[bi];
        block_remove_point(B, bout);
        M B2 = (B & ~BIT(bout)) | BIT(bin_); bm[bi] = B2;
        block_add_point(B2, bin_);
        tabuUntil[(size_t)bi * 128 + bout] = step + tenure + (int)(rng() % 3);
        ++E;  // bump every uncovered m-set's weight by 1
        if ((long)unc.size() < best) best = unc.size();
    }
    if (bestFull > 0) printf("SHRINK result: full cover with %d blocks\n", bestFull);
    double el = chrono::duration<double>(chrono::steady_clock::now() - start).count();
    printf("WLSM2 C(%d,%d,%d,%d) b=%d: best %ld uncovered, %ld steps, %.1fs (%.0f steps/s)\n", v, k, t, m, b, best, step, el, step / max(el, 1e-9));
}
