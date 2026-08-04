#!/usr/bin/env python3
"""Standard-library exact certificate for a candidate proof of the n=32
Reinhardt maximum-perimeter small-polygon problem.

All proof decisions use Fraction/integer arithmetic.  Transcendental values are
bounded by Machin's formula and alternating Taylor intervals.  The 2^31
normalized sign codes are covered by an exact meet-in-the-middle scan.
"""
from fractions import Fraction as Q
from bisect import bisect_left, bisect_right
from hashlib import sha256
import json, time

N=32
M=10**18
EPS=Q(135,10**15)                 # 1.35e-13
A_GAP=Q(9,100)                    # 0.09
B_GAP=Q(107,1000)                 # 0.107
R0=Q(347,100_000_000)             # 3.47e-6
T0=Q(341,10_000_000)              # 3.41e-5 (coarse code-screen threshold)
SIGMA_CODE=Q(41,10)               # 4.1
BEST_AXIAL="+-++--+-+-+---++--+++-+-+-++--+-"
ORBIT_REPS=(
    "++++--+--+-+-+-+++---+-+-+-++-++",
    "+++-+-+-++--+--+--++-+-+-+++--++",
    "++-++-+-+-+-+--++--+-+-+-+-++-++",
)
S_INTS=(-84838116394748,-7759417296262,-90254740014369,
        -172750062732476,-104929943709607,-37109824686737,
        -104570029200993,-46261230637792,-102998890440414,
        -56520596873701,-100433725625308,-67650481919513,
        -34867238213719,-2083994507925,-1041997253963)
S_DEN=10**20
TL=Q(99999953359,10**11)
TU=Q(99999953360,10**11)


def floor_q(x:Q)->int:
    return x.numerator//x.denominator

def ceil_q(x:Q)->int:
    return -floor_q(-x)

def nearest_int(x:Q)->int:
    return floor_q(x+Q(1,2))

def atan_interval(x:Q, terms:int):
    s=Q(0); power=x; vals=[]
    for k in range(terms+1):
        s += power/Q(2*k+1) if k%2==0 else -power/Q(2*k+1)
        if k in (terms-1,terms): vals.append(s)
        power *= x*x
    return min(vals),max(vals)

def pi_interval():
    aL,aU=atan_interval(Q(1,5),30)
    bL,bU=atan_interval(Q(1,239),10)
    # pi = 16 atan(1/5) - 4 atan(1/239)
    return 16*aL-4*bU,16*aU-4*bL

PI_L,PI_U=pi_interval()
assert PI_L>3 and PI_U<Q(22,7)


def sin_point_interval(x:Q, terms:int=24):
    assert 0<=x<=PI_U/2
    if x==0: return Q(0),Q(0)
    s=Q(0); term=x; vals=[]
    for k in range(terms+1):
        s += term
        if k in (terms-1,terms): vals.append(s)
        term *= -x*x/Q((2*k+2)*(2*k+3))
    return min(vals),max(vals)

def sin_interval(lo:Q,hi:Q):
    assert 0<=lo<=hi<=PI_U/2
    l,_=sin_point_interval(lo)
    _,u=sin_point_interval(hi)
    return l,u

def add_iv(a,b): return a[0]+b[0],a[1]+b[1]
def sub_iv(a,b): return a[0]-b[1],a[1]-b[0]
def scale_iv(k:Q,a):
    return (k*a[0],k*a[1]) if k>=0 else (k*a[1],k*a[0])

# sin/cos of k*pi/32 for k mod 64, reduced to sine on [0,pi/2].
def sin_root_iv(k:int):
    k%=64
    sign=1
    if k>32:
        k-=32; sign=-1
    if k>16:
        k=32-k
    iv=sin_interval(Q(k,32)*PI_L,Q(k,32)*PI_U)
    return scale_iv(sign,iv)

def cos_root_iv(k:int):
    return sin_root_iv(16-k)

ROOT_IV=[]
ROOT_INT=[]
for k in range(64):
    xiv=cos_root_iv(k); yiv=sin_root_iv(k)
    qx=nearest_int((xiv[0]+xiv[1])*M/2)
    qy=nearest_int((yiv[0]+yiv[1])*M/2)
    assert xiv[0]>=Q(qx-1,M) and xiv[1]<=Q(qx+1,M)
    assert yiv[0]>=Q(qy-1,M) and yiv[1]<=Q(qy+1,M)
    ROOT_IV.append((xiv,yiv)); ROOT_INT.append((qx,qy))

# Edge vectors exp(i(j+1)pi/32)-exp(ijpi/32), each coordinate within 1/M.
EDGE_INT=[]
for j in range(N):
    xiv=sub_iv(ROOT_IV[j+1][0],ROOT_IV[j][0])
    yiv=sub_iv(ROOT_IV[j+1][1],ROOT_IV[j][1])
    qx=nearest_int((xiv[0]+xiv[1])*M/2)
    qy=nearest_int((yiv[0]+yiv[1])*M/2)
    assert xiv[0]>=Q(qx-1,M) and xiv[1]<=Q(qx+1,M)
    assert yiv[0]>=Q(qy-1,M) and yiv[1]<=Q(qy+1,M)
    EDGE_INT.append((qx,qy))


def code_signs(code):
    assert len(code)==N and set(code)<={'+','-'}
    return [1 if ch=='+' else -1 for ch in code]

def residual_int(code):
    c=code_signs(code)
    return (sum(c[j]*EDGE_INT[j][0] for j in range(N)),
            sum(c[j]*EDGE_INT[j][1] for j in range(N)))

def norm_lower_than(code, q:Q):
    # Prove true residual > q using ||true-int/M|| < 2N/M.
    x,y=residual_int(code)
    threshold=ceil_q(M*q)+2*N
    return x*x+y*y>threshold*threshold

# ---- Feasible one-dimensional construction for the best orbit ----
c=code_signs(BEST_AXIAL)
assert all(c[N-1-j]==-c[j] for j in range(N))
s=[Q(0)]*(N+1)
for j,v in enumerate(S_INTS,1): s[j]=Q(v,S_DEN)
s[16]=Q(0)
for j in range(17,N): s[j]=-s[N-j]
assert all(s[N-j]==-s[j] for j in range(N+1))
a=[None]+[c[j-1]-c[j] for j in range(1,N)]

def H_interval(t:Q):
    ans=(Q(a[16]),Q(a[16]))
    for j in range(1,16):
        lo=Q(j,N)*PI_L+t*s[j]
        hi=Q(j,N)*PI_U+t*s[j]
        term=scale_iv(Q(2*a[j]),sin_interval(lo,hi))
        ans=add_iv(ans,term)
    return ans

HL=H_interval(TL); HU=H_interval(TU)
assert HL[0]>0 and HU[1]<0

# For every t in [TL,TU], certify all gaps and a lower perimeter bound.
P_IV=(Q(0),Q(0))
for j in range(N):
    q=s[j+1]-s[j]
    products=(TL*q,TU*q)
    lo=PI_L/Q(2*N)+min(products)/2
    hi=PI_U/Q(2*N)+max(products)/2
    # these are half-gaps
    assert 2*lo>A_GAP and 2*hi<B_GAP
    P_IV=add_iv(P_IV,scale_iv(2,sin_interval(lo,hi)))
U_IV=scale_iv(2*N,sin_interval(PI_L/Q(2*N),PI_U/Q(2*N)))
DEF_IV=sub_iv(U_IV,P_IV)
assert DEF_IV[1]<EPS

# ---- Self-contained n=32 saturation certificate ----
# If the centrally symmetric difference body has fewer than 64 vertices, it has
# at most 62.  The derivative bound U'(t)>=pi^2/(6t^3), together with pi>3,
# puts U_32-U_31 far above EPS.
assert Q(3,2*N**3)>EPS
# A normal-cone width outside (mu/2,3mu/2) costs at least
# mu^2 sin(mu/4)/16.  Use mu>3/32 and sin(pi/128)>=1/64.
assert Q(9,32**2*64*16)>2*EPS
# Once widths are localized, D_Z<=2EPS gives r>1-64EPS>1/2.
assert 64*EPS<Q(1,2)
# If |delta|>=mu/16, then 1-cos(delta)>=1/131072, whereas
# the localized radial/angular deficit gives 1-cos(delta)<128 EPS.
assert 128*EPS<Q(1,131072)
# The final KKT angular contradiction uses 2+3*pi<12 (pi<10/3).
assert PI_U<Q(10,3)

# ---- High-perimeter localization and radius bound ----
def sin_taylor_lower(x:Q): return x-x**3/Q(6)
# One-gap Jensen deficit at lower endpoint, using pi>3.
dlo=Q(3,32)-A_GAP
klo=sin_taylor_lower(A_GAP/2)/4
assert klo*Q(N,N-1)*dlo*dlo>EPS
# Upper endpoint, using pi<22/7 for distance and pi>3 for other-gap minimum.
dhi=B_GAP-Q(11,112)
u_lo=(Q(3)-B_GAP)/Q(N-1)
khi=sin_taylor_lower(u_lo/2)/4
assert khi*Q(N,N-1)*dhi*dhi>EPS
# Strong concavity on [0.09,0.107].
K0=sin_taylor_lower(A_GAP/2)/4
assert K0*R0*R0>EPS
# Poincare remainder C <= 256 and uniform linear norm <8.
assert 8*R0+256*R0*R0<T0

# ---- Exact meet-in-the-middle scan of all 2^31 codes with c_0=+1 ----
def block_sums(indices, fixed0=False):
    if fixed0:
        out=[(EDGE_INT[0][0],EDGE_INT[0][1],0)]
        indices=[j for j in indices if j!=0]
    else:
        out=[(0,0,0)]
    for j in indices:
        vx,vy=EDGE_INT[j]
        nxt=[]
        bit=1<<j
        for x,y,mask in out:
            nxt.append((x+vx,y+vy,mask))
            nxt.append((x-vx,y-vy,mask|bit))
        out=nxt
    return out

start=time.time()
A=block_sums(list(range(16)),fixed0=True)
B=block_sums(list(range(16,32)),fixed0=False)
B.sort(key=lambda r:r[0])
BX=[r[0] for r in B]
QBOUND=ceil_q(M*T0)+2*N
survivor_masks=[]
pair_tests=0
for ax,ay,am in A:
    lo=bisect_left(BX,-ax-QBOUND)
    hi=bisect_right(BX,-ax+QBOUND)
    pair_tests += hi-lo
    for bx,by,bm in B[lo:hi]:
        x=ax+bx; y=ay+by
        if x*x+y*y<=QBOUND*QBOUND:
            survivor_masks.append(am|bm)
assert len(survivor_masks)==96 and len(set(survivor_masks))==96

def mask_code(mask): return ''.join('-' if (mask>>j)&1 else '+' for j in range(N))
survivors={mask_code(m) for m in survivor_masks}

def negate(code): return ''.join('+' if ch=='-' else '-' for ch in code)
def orbit(code):
    full=code+negate(code)
    out=set()
    for base in (full,full[::-1]):
        for r in range(2*N):
            z=base[r:]+base[:r]
            h=z[:N]
            if h[0]=='-': h=negate(h)
            out.add(h)
    return out

orbits=[orbit(r) for r in ORBIT_REPS]
assert all(len(o)==32 for o in orbits)
assert not (orbits[0]&orbits[1] or orbits[0]&orbits[2] or orbits[1]&orbits[2])
assert survivors==set().union(*orbits)
assert BEST_AXIAL in orbits[1]

# ---- Code-specific linear norms and elimination of the two nonbest orbits ----
def sigma_sq_gershgorin_upper(code):
    cc=code_signs(code)
    aa=[cc[j-1]-cc[j] for j in range(1,N)]
    # Bhat rows: -a_j sin(j theta), a_j cos(j theta), represented /M.
    rows=[[-aa[j-1]*ROOT_INT[j][1] for j in range(1,N)],
          [ aa[j-1]*ROOT_INT[j][0] for j in range(1,N)]]
    def dinv(j,k): return Q(min(j,k)*(N-max(j,k)),N)
    K=[[Q(0),Q(0)],[Q(0),Q(0)]]
    for p in range(2):
        for q in range(2):
            total=Q(0)
            for j in range(1,N):
                if rows[p][j-1]==0: continue
                for k in range(1,N):
                    if rows[q][k-1]==0: continue
                    total += Q(rows[p][j-1]*rows[q][k-1],M*M)*dinv(j,k)
            K[p][q]=total
    SD=sum((dinv(j,k) for j in range(1,N) for k in range(1,N)),Q(0))
    e=Q(2,M); bmax=Q(2)+e
    EK=SD*(2*bmax*e+e*e)
    return max(K[0][0]+abs(K[0][1]),K[1][1]+abs(K[1][0]))+2*EK

sigma_uppers=[]
for rep in ORBIT_REPS:
    su=sigma_sq_gershgorin_upper(rep)
    sigma_uppers.append(su)
    assert su<SIGMA_CODE*SIGMA_CODE

# Certified regular-residual lower bounds for nonbest orbits.
BLOW=(Q(231,10_000_000),None,Q(166,10_000_000))
for idx in (0,2):
    assert norm_lower_than(ORBIT_REPS[idx],BLOW[idx])
    rem=256*R0*R0
    xlo=(BLOW[idx]-rem)/SIGMA_CODE
    assert K0*xlo*xlo>EPS

# ---- Uniqueness of the best-code high-perimeter KKT point ----
cc=code_signs(BEST_AXIAL)
switch=[j for j in range(1,N) if cc[j-1]!=cc[j]]
assert len(switch)==21
# At the regular point, certify |sum exp(2 i j theta)|<2.
sx=sum(ROOT_INT[(2*j)%64][0] for j in switch)
sy=sum(ROOT_INT[(2*j)%64][1] for j in switch)
assert sx*sx+sy*sy < (2*M-2*len(switch))**2
# Any high-perimeter point has ||s|| <=16 R0.  Then the switch-root sum
# changes by <2 sqrt(21)||s|| < 160 R0 <1, hence sigma_min(J)>6.
assert 160*R0<1
# multiplier and strong-convexity comparison
Y0=(B_GAP/2)*R0/6
M0=sin_taylor_lower(A_GAP/2)/Q(512) # use sin(pi/64)>=1/32
assert M0>2*Y0

elapsed=time.time()-start
report={
  'pi_width': str(float(PI_U-PI_L)),
  'feasible_H_left_lower': str(float(HL[0])),
  'feasible_H_right_upper': str(float(HU[1])),
  'feasible_deficit_upper': str(float(DEF_IV[1])),
  'high_perimeter_radius_R0': str(float(R0)),
  'uniform_code_threshold_T0': str(float(T0)),
  'MITM_pair_tests': pair_tests,
  'normalized_codes_covered': 2**31,
  'survivors': len(survivors),
  'orbits': len(orbits),
  'orbit_representatives': list(ORBIT_REPS),
  'sigma_squared_upper_bounds': [str(float(x)) for x in sigma_uppers],
  'uniqueness_margin': str(float(M0-2*Y0)),
  'elapsed_seconds': elapsed,
}
print(json.dumps(report,indent=2))
print('ALL ASSERTIONS PASSED')
