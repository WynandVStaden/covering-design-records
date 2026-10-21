// Coverage analysis and zero-coverage block removal for C(v,k,t,m) coverings.
//
// An m-set M is covered by block B when |B & M| >= t. The exclusive coverage of a block is the number of m-sets
// covered by that block and by no other (WSC "Coverage analysis"). A block with exclusive coverage 0 can be deleted
// without uncovering anything. Deleting one block can raise others' exclusive coverage, so blocks are removed one at
// a time with exact counts.
//
// usage: zerocov v k t m in.txt [out.txt]
//   prints the coverage analysis of the input, removes duplicates and then zero-coverage blocks (scanning in file
//   order, repeating until none is left), checks the result is still 100% covered and writes it to out.txt.
// Memory: 2 bytes per m-subset of the v points (C(v,m) <= ~4e9 on a 16 GB machine).
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <set>
#include <algorithm>
#include <omp.h>
using namespace std;
typedef uint64_t u64;

int v, k, t, m;
u64 C[130][40];
vector<vector<int>> blocks;
vector<uint16_t> cnt;

// colex rank of a sorted m-subset of 0..v-1
inline u64 rankset(const int* s) { u64 r = 0; for (int i = 0; i < m; i++) r += C[s[i]][i + 1]; return r; }

// call f(rank) for every m-set M with |B & M| >= t
template <class F> void forCovered(const vector<int>& B, F f) {
    vector<int> in = B, out; sort(in.begin(), in.end());
    out.clear(); vector<char> isB(v, 0); for (int x : in) isB[x] = 1; for (int x = 0; x < v; x++) if (!isB[x]) out.push_back(x);
    int K = in.size(), R = out.size();
    vector<int> a, b, s(m);
    for (int j = max(t, m - R); j <= min(K, m); j++) {
        a.assign(j, 0); for (int i = 0; i < j; i++) a[i] = i;
        while (true) {
            int r = m - j; b.assign(r, 0); for (int i = 0; i < r; i++) b[i] = i;
            while (true) {
                // merge in[a] and out[b] into sorted s
                int p = 0, q = 0, o = 0;
                while (p < j || q < r) {
                    if (q >= r || (p < j && in[a[p]] < out[b[q]])) s[o++] = in[a[p++]]; else s[o++] = out[b[q++]];
                }
                f(rankset(s.data()));
                int i = r - 1; while (i >= 0 && b[i] == R - r + i) i--; if (i < 0) break;
                b[i]++; for (int x = i + 1; x < r; x++) b[x] = b[x - 1] + 1;
            }
            int i = j - 1; while (i >= 0 && a[i] == K - j + i) i--; if (i < 0) break;
            a[i]++; for (int x = i + 1; x < j; x++) a[x] = a[x - 1] + 1;
        }
    }
}

u64 exclusive(const vector<int>& B) { u64 e = 0; forCovered(B, [&](u64 r) { e += cnt[r] == 1; }); return e; }

int main(int argc, char** argv) {
    if (argc < 6) { fprintf(stderr, "usage: zerocov v k t m in.txt [out.txt]\n"); return 2; }
    v = atoi(argv[1]); k = atoi(argv[2]); t = atoi(argv[3]); m = atoi(argv[4]);
    for (int n = 0; n < 130; n++) { C[n][0] = 1; for (int r = 1; r < 40; r++) C[n][r] = n ? C[n - 1][r - 1] + C[n - 1][r] : 0; }
    FILE* fin = fopen(argv[5], "r"); if (!fin) { perror(argv[5]); return 2; }
    set<vector<int>> seen; int dup = 0, x; vector<int> cur;
    char line[4096];
    while (fgets(line, sizeof line, fin)) {
        cur.clear(); char* p = line; int n;
        while (sscanf(p, "%d%n", &x, &n) == 1) { cur.push_back(x - 1); p += n; }
        if (cur.empty()) continue;
        sort(cur.begin(), cur.end());
        if ((int)cur.size() != k || cur.front() < 0 || cur.back() >= v || unique(cur.begin(), cur.end()) != cur.end()) {
            fprintf(stderr, "bad block: %s", line); return 2;
        }
        if (seen.insert(cur).second) blocks.push_back(cur); else dup++;
    }
    u64 N = C[v][m]; cnt.assign(N, 0);
    #pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < blocks.size(); i++) forCovered(blocks[i], [&](u64 r) { __atomic_fetch_add(&cnt[r], 1, __ATOMIC_RELAXED); });
    u64 unc = 0;
#pragma omp parallel for reduction(+:unc)
    for (u64 i = 0; i < N; i++) unc += cnt[i] == 0;
    printf("C(%d,%d,%d,%d): %zu blocks (+%d duplicates removed), %llu of %llu m-sets uncovered\n", v, k, t, m,
           blocks.size(), dup, (unsigned long long)unc, (unsigned long long)N);
    if (unc) { printf("not a covering, nothing removed\n"); return 1; }
    // coverage analysis of the (deduplicated) input
    u64 zero0 = 0; vector<u64> ex(blocks.size());
    #pragma omp parallel for schedule(dynamic, 64) reduction(+:zero0)
    for (size_t i = 0; i < blocks.size(); i++) { ex[i] = exclusive(blocks[i]); zero0 += ex[i] == 0; }
    printf("coverage analysis: %llu blocks with exclusive coverage 0\n", (unsigned long long)zero0);
    if (getenv("ANALYSIS")) { FILE* fa = fopen(getenv("ANALYSIS"), "w");
        for (size_t i = 0; i < blocks.size(); i++) { for (int y : blocks[i]) fprintf(fa, "%d ", y + 1); fprintf(fa, ": %llu\n", (unsigned long long)ex[i]); }
        fclose(fa); }
    // remove zero-coverage blocks one at a time, in file order. Deleting a block only lowers counts, so other blocks'
    // exclusive coverage can only grow: a block above zero now can never become removable. Only the blocks at zero in
    // the analysis above are candidates; each is re-checked against the current counts just before it is removed.
    vector<char> alive(blocks.size(), 1); size_t removed = 0;
    for (size_t i = 0; i < blocks.size(); i++) if (ex[i] == 0 && exclusive(blocks[i]) == 0) {
        forCovered(blocks[i], [&](u64 r) { cnt[r]--; }); alive[i] = 0; removed++;
    }
    unc = 0;
#pragma omp parallel for reduction(+:unc)
    for (u64 i = 0; i < N; i++) unc += cnt[i] == 0;
    printf("removed %zu zero-coverage blocks: %zu -> %zu blocks, uncovered after: %llu\n", removed, blocks.size(),
           blocks.size() - removed, (unsigned long long)unc);
    if (unc) { printf("BUG: result not covered\n"); return 3; }
    if (argc > 6) { FILE* fo = fopen(argv[6], "w");
        for (size_t i = 0; i < blocks.size(); i++) if (alive[i]) { for (int j = 0; j < k; j++) fprintf(fo, j ? " %d" : "%d", blocks[i][j] + 1); fputc('\n', fo); }
        fclose(fo); }
    return 0;
}
