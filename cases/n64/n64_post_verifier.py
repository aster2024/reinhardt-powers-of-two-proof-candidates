#!/usr/bin/env python3
from fractions import Fraction as Q
from collections import defaultdict
import json
# Reuse rigorously generated intervals/constants. Import prints its report by design.
import n64_analytic_verifier as A
N=64; M=10**16; R0=A.R0; EPS=A.EPS; A_GAP=A.A_GAP; B_GAP=A.B_GAP
REPS=(
'++++++--++++++-----+--+-+++--+---+--+++---++-+++-++---+-++-+++++',
'++++--++-++--+---+--+-+-+-+-+-+-++-+++-++--+--++----++-++--++-++',
'++++-+++-+--+-+---++--+-+----+--++--++-++++-+-++--+++-+-++-+---+',
'++++-+++-+-+----++--+---+-++-+----++++-+--+-+++-++--++++-+-+---+',
'++++-+--+++-+---++-+----++-+---+-++-++-+--++++--+-++-++-+---+-++',
'+++++---+++--+++--++-+-+--+--++--++-+--+-+--+-+-+-+++---+++-++++')
BEST='-++++++-----+--+-+++--+---+--+++---++-+++-++---+-++-+++++------+'
BLOW=(None,Q(486,10**12),Q(508,10**12),Q(628,10**12),Q(764,10**12),Q(823,10**12))

def floor_q(x):return x.numerator//x.denominator
def ceil_q(x):return -floor_q(-x)
def nearest(x):return floor_q(x+Q(1,2))
def neg(s):return ''.join('+' if x=='-' else '-' for x in s)
def canon(code):
    full=code+neg(code); out=[]
    for base in (full,full[::-1]):
        for r in range(2*N):
            z=base[r:]+base[:r]; h=z[:N]
            if h[0]=='-':h=neg(h)
            out.append(h)
    return min(out)

survivors=[z.strip() for z in open('n64_fixed_survivors.txt') if z.strip()]
assert len(survivors)==896 and len(set(survivors))==896
g=defaultdict(list)
for z in survivors:g[canon(z)].append(z)
assert set(g)==set(REPS)
assert sorted(len(v) for v in g.values())==[128,128,128,128,128,256]
assert canon(BEST)==REPS[0]

# Roots exp(i*k*pi/64), each coordinate within 1/M of ROOT_INT.
def sin_root_iv(k):
    k%=128; sign=1
    if k>64:k-=64;sign=-1
    if k>32:k=64-k
    iv=A.sin_interval(Q(k,64)*A.PI_L,Q(k,64)*A.PI_U)
    return A.scale_iv(sign,iv)
def cos_root_iv(k):return sin_root_iv(32-k)
ROOT_INT=[]
for k in range(128):
    xi=cos_root_iv(k); yi=sin_root_iv(k)
    qx=nearest((xi[0]+xi[1])*M/2); qy=nearest((yi[0]+yi[1])*M/2)
    assert xi[0]>=Q(qx-1,M) and xi[1]<=Q(qx+1,M)
    assert yi[0]>=Q(qy-1,M) and yi[1]<=Q(qy+1,M)
    ROOT_INT.append((qx,qy))
EDGE_INT=[(ROOT_INT[j+1][0]-ROOT_INT[j][0],ROOT_INT[j+1][1]-ROOT_INT[j][1]) for j in range(N)]

def signs(code):return [1 if z=='+' else -1 for z in code]
def residual_int(code):
    c=signs(code)
    return sum(c[j]*EDGE_INT[j][0] for j in range(N)),sum(c[j]*EDGE_INT[j][1] for j in range(N))
def prove_residual_gt(code,q):
    x,y=residual_int(code); t=ceil_q(M*q)+2*N
    return x*x+y*y>t*t

def sigma_sq_upper(code):
    c=signs(code); aa=[c[j-1]-c[j] for j in range(1,N)]
    rows=[[-aa[j-1]*ROOT_INT[j][1] for j in range(1,N)],
          [ aa[j-1]*ROOT_INT[j][0] for j in range(1,N)]]
    def dinv(j,k):return Q(min(j,k)*(N-max(j,k)),N)
    K=[[Q(0),Q(0)],[Q(0),Q(0)]]
    for p in range(2):
      for q in range(2):
        total=Q(0)
        for j in range(1,N):
          if rows[p][j-1]==0:continue
          for k in range(1,N):
            if rows[q][k-1]==0:continue
            total+=Q(rows[p][j-1]*rows[q][k-1],M*M)*dinv(j,k)
        K[p][q]=total
    SD=sum((dinv(j,k) for j in range(1,N) for k in range(1,N)),Q(0))
    e=Q(2,M); bmax=Q(2)+e; EK=SD*(2*bmax*e+e*e)
    return max(K[0][0]+abs(K[0][1]),K[1][1]+abs(K[1][0]))+2*EK

sigmas=[]
K0=(A_GAP/2-(A_GAP/2)**3/Q(6))/4
rem=512*R0*R0
for i,rep in enumerate(REPS):
    su=sigma_sq_upper(rep);sigmas.append(su);assert su<36
    if i:
        assert prove_residual_gt(rep,BLOW[i])
        xlo=(BLOW[i]-rem)/6
        assert K0*xlo*xlo>EPS

# Continuous uniqueness for best orbit.
c=signs(BEST); switch=[j for j in range(1,N) if c[j-1]!=c[j]]
assert len(switch)==27
sx=sum(ROOT_INT[(2*j)%128][0] for j in switch);sy=sum(ROOT_INT[(2*j)%128][1] for j in switch)
assert sx*sx+sy*sy < (5*M-2*len(switch))**2
# Poincare gives ||s||<23 R0. Root-sum change <2 sqrt(27)||s||<276 R0<1.
assert 276*R0<1
# Hence switch-root norm<6 and sigma_min(Dg)>6.
Y0=(B_GAP/2)*R0/6
# sin(pi/128)>1/45, so lambda_min(path)>4/45^2.
root_l,_=A.sin_interval(A.PI_L/Q(2*N),A.PI_U/Q(2*N));assert root_l>Q(1,45)
M0=2*(A_GAP/2-(A_GAP/2)**3/Q(6))/Q(45**2)
assert M0>2*Y0

print(json.dumps({
 'survivors':len(survivors),'orbits':len(g),'orbit_sizes':sorted(len(v) for v in g.values()),
 'sigma_sq_uppers':[float(x) for x in sigmas],
 'nonbest_deficit_margins':[float(K0*((BLOW[i]-rem)/6)**2-EPS) for i in range(1,6)],
 'uniqueness_margin':float(M0-2*Y0)
},indent=2))
print('ALL POST-SCREEN ASSERTIONS PASSED')
