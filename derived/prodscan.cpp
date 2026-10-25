// Product/interval construction scan for lottery designs C(v,k,t,m).
// Term: D1 x D2 where D1 = C(v1,k1,t1,m1) on V1, D2 = C(v2,k2,t2,m2) on V2, k1+k2=k, t1+t2=t.
// It covers every m-set M with |M&V1| in [m1, m-m2]. A set of terms covering all feasible a=|M&V1| gives a cover.
// usage: prodscan ctable_out table1 [table2 ...] > out   (tables: lines "v k t m b")
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const int V=90,K=25,T=20,M=20; const ll INF=(ll)4e18, CAP=(ll)1e15;
ll c[V+1][K+1][T+1][M+1]; ll tab[V+1][K+1][T+1][M+1];
int rule[V+1][K+1][T+1][M+1]; // 0 base/table, 1 product (v1 stored), 2 m-1, 3 k-1, 4 t+1, 5 v+1, 6 v-1,m-1, 7 dual
int rv1[V+1][K+1][T+1][M+1];
static inline ll mul(ll a,ll b){ if(a>=INF||b>=INF) return INF; if(a==0||b==0) return 0; if(a>CAP/b) return CAP; return a*b;}
bool upd(int v,int k,int t,int m,ll x,int r,int v1=0){ if(x<c[v][k][t][m]){c[v][k][t][m]=x;rule[v][k][t][m]=r;rv1[v][k][t][m]=v1;return true;} return false;}
int main(int argc,char**argv){
  for(int v=0;v<=V;v++)for(int k=0;k<=K;k++)for(int t=0;t<=T;t++)for(int m=0;m<=M;m++){
    ll x=INF; tab[v][k][t][m]=INF;
    if(m<t) {c[v][k][t][m]=INF;continue;}
    if(t==0) x=(m>v?0:1);
    else if(m>v) x=0;
    else if(t>k) x=INF;
    else if(v<=k) x=1;
    else if(t==1) x=(v-m+1+k-1)/k;
    c[v][k][t][m]=x; rule[v][k][t][m]=0;
  }
  for(int i=2;i<argc;i++){ FILE*f=fopen(argv[i],"r"); char line[512];
    while(fgets(line,512,f)){ int v,k,t,m; ll b; if(sscanf(line,"%d %d %d %d %lld",&v,&k,&t,&m,&b)!=5) continue;
      if(v>V||k>K||t>T||m>M||m<t) continue; tab[v][k][t][m]=min(tab[v][k][t][m],b); if(b<c[v][k][t][m]){c[v][k][t][m]=b;rule[v][k][t][m]=0;} }
    fclose(f);}
  for(int pass=0;pass<6;pass++){
    long changes=0;
    for(int v=1;v<=V;v++){
      // closure at this v first
      for(int it=0;it<2;it++)
      for(int k=0;k<=K;k++)for(int t=0;t<=T;t++)for(int m=t;m<=M;m++){
        if(m>t) changes+=upd(v,k,t,m,c[v][k][t][m-1],2);
        if(k>0) changes+=upd(v,k,t,m,c[v][k-1][t][m],3);
        if(t<T&&t+1<=m) changes+=upd(v,k,t,m,c[v][k][t+1][m],4);
        if(v<V) changes+=upd(v,k,t,m,c[v+1][k][t][m],5);
        if(m>0&&m-1>=t) changes+=upd(v,k,t,m,c[v-1][k][t][m-1],6);
      }
      vector<pair<int,int>> kts; for(int k=2;k<=K;k++)for(int t=2;t<=k;t++) if(v>k) kts.push_back({k,t});
      long ch=0;
      #pragma omp parallel for schedule(dynamic) reduction(+:ch)
      for(size_t q=0;q<kts.size();q++){ int k=kts[q].first,t=kts[q].second;
        ll best_m[M+1]; int bv1[M+1]; for(int m=0;m<=M;m++){best_m[m]=INF;bv1[m]=0;}
        static thread_local ll W[M+1][M+1];
        for(int v1=1;v1<=v/2;v1++){ int v2=v-v1; int A=min(M,v1),B=min(M,v2);
          for(int a=0;a<=A;a++)for(int b=0;b<=B;b++) W[a][b]=INF;
          for(int k1=0;k1<=k;k1++){int k2=k-k1;
            for(int t1=max(0,t-k2);t1<=min(t,k1);t1++){int t2=t-t1;
              for(int m1=t1;m1<=A;m1++){ ll x1=c[v1][k1][t1][m1]; if(x1>=INF) continue;
                for(int m2=t2;m2<=B;m2++){ ll x=mul(x1,c[v2][k2][t2][m2]); if(x<W[m1][m2]) W[m1][m2]=x; }}}}
          // monotone fix: W decreasing in both args
          for(int a=0;a<=A;a++)for(int b=1;b<=B;b++) W[a][b]=min(W[a][b],W[a][b-1]);
          for(int a=1;a<=A;a++)for(int b=0;b<=B;b++) W[a][b]=min(W[a][b],W[a-1][b]);
          for(int m=t;m<=min(M,v);m++){ int lo=max(0,m-v2),hi=min(m,v1); ll bst[M+2];
            for(int x=0;x<=M+1;x++) bst[x]=INF; bst[lo]=0;
            for(int r=lo;r<=hi;r++){ ll y=INF; for(int l=lo;l<=r;l++){ if(bst[l]>=INF) continue; ll w=W[l][m-r]; if(w>=INF) continue; ll z=bst[l]+w; if(z<y) y=z;} bst[r+1]=min(y,CAP);}
            if(bst[hi+1]<best_m[m]){best_m[m]=bst[hi+1];bv1[m]=v1;}
          }
        }
        for(int m=t;m<=min(M,v);m++) if(best_m[m]<c[v][k][t][m]){c[v][k][t][m]=best_m[m];rule[v][k][t][m]=1;rv1[v][k][t][m]=bv1[m];ch++;}
      }
      changes+=ch;
      // complement duality: C(v,k,t,m) = C(v, v-k, v-k-m+t, v-m)
      for(int it=0;it<2;it++)
      for(int k=0;k<=K;k++)for(int t=0;t<=T;t++)for(int m=t;m<=M;m++){
        int k2=v-k,t2=v-k-m+t,m2=v-m; if(k2<0||k2>K||m2<0||m2>M) continue; if(t2<0) t2=0; if(t2>m2||t2>T) continue;
        if(c[v][k2][t2][m2]<c[v][k][t][m]){c[v][k][t][m]=c[v][k2][t2][m2];rule[v][k][t][m]=7;rv1[v][k][t][m]=0;changes++;}
      }
    }
    fprintf(stderr,"pass %d changes %ld\n",pass,changes);
    if(!changes) break;
  }
  // report: cells where c < table value (table entry present)
  for(int v=1;v<=V;v++)for(int k=0;k<=K;k++)for(int t=2;t<=T;t++)for(int m=t;m<=M;m++)
    if(tab[v][k][t][m]<INF && c[v][k][t][m]<tab[v][k][t][m])
      printf("%d %d %d %d %lld %lld %d %d\n",v,k,t,m,tab[v][k][t][m],c[v][k][t][m],rule[v][k][t][m],rv1[v][k][t][m]);
  // dump full c table
  FILE*g=fopen(argv[1],"w");
  for(int v=0;v<=V;v++)for(int k=0;k<=K;k++)for(int t=0;t<=T;t++)for(int m=t;m<=M;m++) if(c[v][k][t][m]<INF)
    fprintf(g,"%d %d %d %d %lld %d %d\n",v,k,t,m,c[v][k][t][m],rule[v][k][t][m],rv1[v][k][t][m]);
  fclose(g);
}
