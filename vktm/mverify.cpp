// Independent brute-force checker for (v,k,t,m) coverings: every m-subset must meet some block in >= t points.
// usage: mverify v k t m file
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
typedef unsigned __int128 M;
static inline int pc(M x) { return __builtin_popcountll((uint64_t)x) + __builtin_popcountll((uint64_t)(x >> 64)); }
int main(int argc, char** argv) {
    int v = atoi(argv[1]), k = atoi(argv[2]), t = atoi(argv[3]), m = atoi(argv[4]);
    FILE* f = fopen(argv[5], "r"); std::vector<M> B; int p; M x = 0; int c = 0; bool bad = false;
    while (fscanf(f, "%d", &p) == 1) { if (p < 1 || p > v || ((x >> (p - 1)) & 1)) bad = true; x |= ((M)1) << (p - 1); if (++c == k) { B.push_back(x); x = 0; c = 0; } }
    if (c || bad) { printf("INVALID format\n"); return 1; }
    for (size_t i = 0; i < B.size(); ++i) for (size_t j = i + 1; j < B.size(); ++j) if (B[i] == B[j]) { printf("INVALID duplicate block\n"); return 1; }
    int idx[24]; for (int i = 0; i < m; ++i) idx[i] = i;
    unsigned long long n = 0, unc = 0;
    while (true) {
        M s = 0; for (int i = 0; i < m; ++i) s |= ((M)1) << idx[i];
        bool ok = false; for (M b : B) if (pc(b & s) >= t) { ok = true; break; }
        if (!ok) ++unc; ++n;
        int i = m - 1; while (i >= 0 && idx[i] == v - m + i) --i;
        if (i < 0) break;
        ++idx[i]; for (int l = i + 1; l < m; ++l) idx[l] = idx[l - 1] + 1;
    }
    printf("C(%d,%d,%d,%d) with %zu blocks: %s (%llu of %llu m-subsets uncovered)\n", v, k, t, m, B.size(), unc ? "INVALID" : "VALID", unc, n);
    return unc ? 1 : 0;
}
