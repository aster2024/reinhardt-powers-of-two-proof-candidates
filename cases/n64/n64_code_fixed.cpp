#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <vector>
#include <chrono>
#include <string>
#include <fstream>
using namespace std;
struct Item { int64_t sum; uint32_t support; uint32_t neg; };
struct Rec { uint32_t support; uint32_t neg; int64_t sum; };
static void gen_rec(const int64_t* w,int idx,int64_t sum,uint32_t supp,uint32_t neg,vector<Item>&out,int offset){
 if(idx==16){out.push_back({sum,supp,neg});return;}
 uint32_t bit=1u<<(offset+idx); int64_t v=w[idx];
 gen_rec(w,idx+1,sum,supp,neg,out,offset);
 gen_rec(w,idx+1,sum+v,supp|bit,neg,out,offset);
 gen_rec(w,idx+1,sum-v,supp|bit,neg|bit,out,offset);
}
static vector<Rec> near_ternary(const int64_t w[32],int64_t Q,const char*lab){
 size_t sz=1;for(int i=0;i<16;i++)sz*=3;
 vector<Item>A;A.reserve(sz);vector<Item>B;B.reserve(sz);
 auto t0=chrono::steady_clock::now();
 gen_rec(w,0,0,0,0,A,0);gen_rec(w+16,0,0,0,0,B,16);
 fprintf(stderr,"%s generated\n",lab);
 auto cmp=[](const Item&a,const Item&b){return a.sum<b.sum;};
 sort(A.begin(),A.end(),cmp);sort(B.begin(),B.end(),cmp);
 vector<Rec>out;out.reserve(4000000);size_t lo=0,hi=0;
 for(size_t ii=A.size();ii-->0;){auto&a=A[ii];int64_t L=-a.sum-Q,U=-a.sum+Q;
  while(lo<B.size()&&B[lo].sum<L)++lo;if(hi<lo)hi=lo;while(hi<B.size()&&B[hi].sum<=U)++hi;
  for(size_t k=lo;k<hi;k++)out.push_back({a.support|B[k].support,a.neg|B[k].neg,a.sum+B[k].sum});
 }
 fprintf(stderr,"%s outputs %zu in %.2fs\n",lab,out.size(),chrono::duration<double>(chrono::steady_clock::now()-t0).count());
 sort(out.begin(),out.end(),[](const Rec&a,const Rec&b){if(a.support!=b.support)return a.support<b.support;return a.sum<b.sum;});
 return out;
}
static pair<size_t,size_t> range_support(const vector<Rec>&v,uint32_t s){
 auto lo=lower_bound(v.begin(),v.end(),s,[](const Rec&a,uint32_t b){return a.support<b;});
 auto hi=upper_bound(lo,v.end(),s,[](uint32_t b,const Rec&a){return b<a.support;});
 return {(size_t)(lo-v.begin()),(size_t)(hi-v.begin())};
}
int main(){
 const int64_t X[32]={490824570458246,1471291271993348,2448213503984324,3419237775206025,4382024803137396,5334255149497968,6273634807977830,7197900730699763,8104826280099797,8992226593092132,9857963844595681,10699952397741944,11516163828356906,12304631811612537,13063456859075535,13790810894741338,14484941659029338,15144176930129691,15766928552532125,16351696263031674,16897071304994141,17401739822174228,17864486023910306,18284195114070613,18659855976694778,18990563611860733,19275521315908797,19514042600770571,19705552847778825,19849590691974200,19945809133573804,19993976373924084};
 const int64_t Y[32]={19993976373924084,19945809133573804,19849590691974200,19705552847778825,19514042600770571,19275521315908797,18990563611860733,18659855976694778,18284195114070613,17864486023910306,17401739822174228,16897071304994141,16351696263031674,15766928552532125,15144176930129691,14484941659029338,13790810894741338,13063456859075535,12304631811612537,11516163828356906,10699952397741944,9857963844595681,8992226593092132,8104826280099797,7197900730699763,6273634807977830,5334255149497968,4382024803137396,3419237775206025,2448213503984324,1471291271993348,490824570458246};
 const int64_t Q=170000064; // ceil(1e16*1.7e-8)+64
 auto A=near_ternary(X,Q,"X"),B=near_ternary(Y,Q,"Y");
 vector<string>codes; const uint32_t ALL=0xffffffffu; uint32_t last=0;pair<size_t,size_t>rg;bool first=true;
 __int128 Q2=(__int128)Q*Q;
 for(auto&x:A){uint32_t v=ALL^x.support;if(first||v!=last){rg=range_support(B,v);last=v;first=false;}
  for(size_t k=rg.first;k<rg.second;k++){__int128 r2=(__int128)x.sum*x.sum+(__int128)B[k].sum*B[k].sum;if(r2<=Q2){
    string c(64,'?');for(int j=0;j<32;j++){uint32_t bit=1u<<j;int kk=63-j;if(x.support&bit){int d=(x.neg&bit)?-1:1;c[j]=d>0?'+':'-';c[kk]=c[j];}else{int f=(B[k].neg&bit)?-1:1;int e=-f;c[j]=e>0?'+':'-';c[kk]=e>0?'-':'+';}}codes.push_back(c);
  }}
 }
 sort(codes.begin(),codes.end());codes.erase(unique(codes.begin(),codes.end()),codes.end());
 ofstream f("n64_fixed_survivors.txt");for(auto&s:codes)f<<s<<'\n';
 printf("fixed survivors %zu\n",codes.size());
}
