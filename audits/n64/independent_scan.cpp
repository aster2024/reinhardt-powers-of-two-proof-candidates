#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <vector>
#include <chrono>
#include <fstream>
using namespace std;
struct Half { int64_t sum; uint32_t support; uint32_t minus; };
struct Near { uint32_t support; uint32_t minus; int64_t sum; };

static void enumerate_half(const int64_t w[32], const int idx[16], int pos,
                           int64_t sum, uint32_t support, uint32_t minus,
                           vector<Half>& out) {
    if (pos==16) { out.push_back({sum,support,minus}); return; }
    int j=idx[pos]; uint32_t bit=uint32_t(1)<<j;
    enumerate_half(w,idx,pos+1,sum,support,minus,out);
    enumerate_half(w,idx,pos+1,sum+w[j],support|bit,minus,out);
    enumerate_half(w,idx,pos+1,sum-w[j],support|bit,minus|bit,out);
}

static vector<Near> near_zero(const int64_t w[32], int64_t radius, const char* tag) {
    int even[16], odd[16];
    for(int k=0;k<16;k++){even[k]=2*k; odd[k]=2*k+1;}
    size_t total=1; for(int k=0;k<16;k++) total*=3;
    vector<Half> E,O; E.reserve(total); O.reserve(total);
    auto t0=chrono::steady_clock::now();
    enumerate_half(w,even,0,0,0,0,E);
    enumerate_half(w,odd,0,0,0,0,O);
    auto bysum=[](const Half&a,const Half&b){return a.sum<b.sum;};
    sort(E.begin(),E.end(),bysum); sort(O.begin(),O.end(),bysum);
    vector<Near> out; out.reserve(4000000);
    size_t lo=0, hi=0;
    // E is traversed from largest to smallest, so the target interval in O moves right.
    for(size_t ii=E.size(); ii-->0;){
        const auto &e=E[ii];
        int64_t low=-e.sum-radius, high=-e.sum+radius;
        while(lo<O.size() && O[lo].sum<low) ++lo;
        if(hi<lo) hi=lo;
        while(hi<O.size() && O[hi].sum<=high) ++hi;
        for(size_t k=lo;k<hi;k++) out.push_back({e.support|O[k].support,e.minus|O[k].minus,e.sum+O[k].sum});
    }
    sort(out.begin(),out.end(),[](const Near&a,const Near&b){
        if(a.support!=b.support) return a.support<b.support;
        if(a.sum!=b.sum) return a.sum<b.sum;
        return a.minus<b.minus;
    });
    fprintf(stderr,"%s interleaved outputs %zu in %.2fs\n",tag,out.size(),
            chrono::duration<double>(chrono::steady_clock::now()-t0).count());
    return out;
}

static pair<size_t,size_t> support_range(const vector<Near>& v,uint32_t support){
    auto lo=lower_bound(v.begin(),v.end(),support,[](const Near&a,uint32_t b){return a.support<b;});
    auto hi=upper_bound(lo,v.end(),support,[](uint32_t b,const Near&a){return b<a.support;});
    return {size_t(lo-v.begin()),size_t(hi-v.begin())};
}

int main(){
    const int64_t X[32]={490824570458246,1471291271993348,2448213503984324,3419237775206025,4382024803137396,5334255149497968,6273634807977830,7197900730699763,8104826280099797,8992226593092132,9857963844595681,10699952397741944,11516163828356906,12304631811612537,13063456859075535,13790810894741338,14484941659029338,15144176930129691,15766928552532125,16351696263031674,16897071304994141,17401739822174228,17864486023910306,18284195114070613,18659855976694778,18990563611860733,19275521315908797,19514042600770571,19705552847778825,19849590691974200,19945809133573804,19993976373924084};
    const int64_t Y[32]={19993976373924084,19945809133573804,19849590691974200,19705552847778825,19514042600770571,19275521315908797,18990563611860733,18659855976694778,18284195114070613,17864486023910306,17401739822174228,16897071304994141,16351696263031674,15766928552532125,15144176930129691,14484941659029338,13790810894741338,13063456859075535,12304631811612537,11516163828356906,10699952397741944,9857963844595681,8992226593092132,8104826280099797,7197900730699763,6273634807977830,5334255149497968,4382024803137396,3419237775206025,2448213503984324,1471291271993348,490824570458246};
    const int64_t Q=170000064;
    auto HX=near_zero(X,Q,"X"), HY=near_zero(Y,Q,"Y");
    const uint32_t ALL=0xffffffffu;
    __int128 q2=(__int128)Q*Q;
    vector<uint64_t> codes; codes.reserve(1024);
    uint32_t previous=0; pair<size_t,size_t> rg{0,0}; bool first=true;
    for(const auto &x:HX){
        uint32_t ysupp=ALL^x.support;
        if(first || ysupp!=previous){rg=support_range(HY,ysupp);previous=ysupp;first=false;}
        for(size_t k=rg.first;k<rg.second;k++){
            const auto &y=HY[k];
            __int128 norm2=(__int128)x.sum*x.sum+(__int128)y.sum*y.sum;
            if(norm2>q2) continue;
            uint64_t bits=0;
            for(int j=0;j<32;j++){
                uint32_t bit=uint32_t(1)<<j; int jj=63-j;
                int cj, ck;
                if(x.support&bit){
                    int p=(x.minus&bit)?-1:1;
                    cj=p; ck=p;
                }else{
                    int f=(y.minus&bit)?-1:1; // f multiplies -2 sin(beta), hence q=-f
                    cj=-f; ck=f;
                }
                if(cj>0) bits|=uint64_t(1)<<j;
                if(ck>0) bits|=uint64_t(1)<<jj;
            }
            codes.push_back(bits);
        }
    }
    sort(codes.begin(),codes.end()); codes.erase(unique(codes.begin(),codes.end()),codes.end());
    ofstream out("/tmp/n64_independent_masks.txt");
    char buf[32];
    for(uint64_t x:codes){snprintf(buf,sizeof(buf),"%016llx",(unsigned long long)x);out<<buf<<'\n';}
    printf("independent fixed survivors %zu\n",codes.size());
}
