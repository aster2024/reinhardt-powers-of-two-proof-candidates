#!/usr/bin/env python3
"""Exact/interval audit for the proposed n=16 Reinhardt proof.

Standard library only.  All finite-code decisions are made with integer
interval arithmetic.  Trigonometric values at multiples of pi/16 are
constructed from nested square roots, so the code scan uses no floating point.
The scalar inequalities are checked with exact Fractions and alternating
Taylor/Machin bounds.
"""
from fractions import Fraction as Q
from itertools import product
from math import isqrt
import hashlib

# ---------------------------------------------------------------------------
# Exact rational elementary-function bounds
# ---------------------------------------------------------------------------

def atan_alt_bounds(x: Q, terms: int):
    """Bounds atan(x) for 0<x<1 by an alternating series."""
    assert 0 < x < 1 and terms >= 2
    s = Q(0)
    lo = hi = None
    for k in range(terms):
        term = x ** (2*k + 1) / (2*k + 1)
        s += term if k % 2 == 0 else -term
        # even last index -> upper; odd last index -> lower
        if k % 2 == 0:
            hi = s
        else:
            lo = s
    assert lo is not None and hi is not None and lo <= hi
    return lo, hi


def add_iv(x, y): return x[0] + y[0], x[1] + y[1]
def sub_iv(x, y): return x[0] - y[1], x[1] - y[0]
def mul_iv_q(x, y):
    vals = (x[0]*y[0], x[0]*y[1], x[1]*y[0], x[1]*y[1])
    return min(vals), max(vals)
def scale_iv(x, a: Q):
    return (x[0]*a, x[1]*a) if a >= 0 else (x[1]*a, x[0]*a)
def div_iv_pos_q(x, a: Q):
    assert a > 0
    return x[0]/a, x[1]/a


def sin_point_bounds(x: Q, terms: int = 10):
    """Exact alternating-series bounds for sin(x), 0<=x<=0.3."""
    assert 0 <= x <= Q(3,10)
    if x == 0: return Q(0), Q(0)
    s = Q(0); lo = hi = None
    for k in range(terms):
        # recurrence-free exact term
        import math
        term = x ** (2*k+1) / math.factorial(2*k+1)
        s += term if k % 2 == 0 else -term
        if k % 2 == 0: hi = s
        else: lo = s
    return lo, hi


def cos_point_bounds(x: Q, terms: int = 10):
    """Exact alternating-series bounds for cos(x), 0<=x<=0.3."""
    assert 0 <= x <= Q(3,10)
    s = Q(0); lo = hi = None
    for k in range(terms):
        import math
        term = x ** (2*k) / math.factorial(2*k)
        s += term if k % 2 == 0 else -term
        if k % 2 == 0: hi = s
        else: lo = s
    return lo, hi


def sin_iv_small(x):
    # sin is increasing on all intervals used here.
    return sin_point_bounds(x[0])[0], sin_point_bounds(x[1])[1]

def cos_iv_small(x):
    # cos is decreasing on all intervals used here.
    return cos_point_bounds(x[1])[0], cos_point_bounds(x[0])[1]

# Machin formula, exact interval.
a5 = atan_alt_bounds(Q(1,5), 34)
a239 = atan_alt_bounds(Q(1,239), 12)
PI = sub_iv(scale_iv(a5, Q(16)), scale_iv(a239, Q(4)))
assert PI[0] > Q(314159265358979323846, 10**20)
assert PI[1] < Q(314159265358979323847, 10**20)

# ---------------------------------------------------------------------------
# Fixed-point interval arithmetic for nested radicals and the code scan
# ---------------------------------------------------------------------------
SCALE = 10**75

def fdiv(a, b): return a // b
def cdiv(a, b): return -((-a) // b)
def Iadd(x,y): return x[0]+y[0], x[1]+y[1]
def Ineg(x): return -x[1], -x[0]
def Isub(x,y): return Iadd(x,Ineg(y))
def Imul_int(x,k): return (x[0]*k,x[1]*k) if k>=0 else (x[1]*k,x[0]*k)
def Idiv_int(x,k): return fdiv(x[0],k), cdiv(x[1],k)
def Imul(x,y):
    vals=(x[0]*y[0],x[0]*y[1],x[1]*y[0],x[1]*y[1])
    return fdiv(min(vals),SCALE), cdiv(max(vals),SCALE)
def Isqrt(x):
    assert x[0] >= 0
    lo=isqrt(x[0]*SCALE)
    hi=isqrt(x[1]*SCALE)
    if hi*hi < x[1]*SCALE: hi += 1
    return lo,hi
def Irecip(x):
    assert x[0] > 0
    return fdiv(SCALE*SCALE,x[1]), cdiv(SCALE*SCALE,x[0])
def Idiv(x,y): return Imul(x,Irecip(y))
def Cmul(z,w):
    a,b=z; c,d=w
    return Isub(Imul(a,c),Imul(b,d)), Iadd(Imul(a,d),Imul(b,c))
def Iminabs(x):
    if x[0] <= 0 <= x[1]: return 0
    return min(abs(x[0]),abs(x[1]))
def Imaxabs(x): return max(abs(x[0]),abs(x[1]))
def rawprod(x,y):
    vals=(x[0]*y[0],x[0]*y[1],x[1]*y[0],x[1]*y[1])
    return min(vals),max(vals)
def ceil_sqrt(n):
    r=isqrt(n)
    return r if r*r==n else r+1

ONE=(SCALE,SCALE); ZERO=(0,0); HALF=(SCALE//2,SCALE//2)
COS_PI_4=Isqrt(HALF)
COS_PI_8=Isqrt(Idiv_int(Iadd(ONE,COS_PI_4),2))
COS_PI_16=Isqrt(Idiv_int(Iadd(ONE,COS_PI_8),2))
SIN_PI_32=Isqrt(Idiv_int(Isub(ONE,COS_PI_16),2))
SIN_PI_16=Isqrt(Idiv_int(Isub(ONE,COS_PI_8),2))

# Exact algebraic evaluation of Bingane's C_16 lower-bound construction.
# If A=atan(sec(pi/8)-1), B=asin(sin(A)cos(pi/16)), delta=A-B,
# then all quantities below determine cos(delta/2) using only radicals.
t=Isub(Irecip(COS_PI_8),ONE)
den=Isqrt(Iadd(ONE,Imul(t,t)))
sinA=Idiv(t,den); cosA=Irecip(den)
y=Imul(sinA,COS_PI_16)
cosB=Isqrt(Isub(ONE,Imul(y,y)))
cosDelta=Iadd(Imul(cosA,cosB),Imul(sinA,y))
cosHalfDelta=Isqrt(Idiv_int(Iadd(ONE,cosDelta),2))
L_C16=Imul_int(Imul(SIN_PI_32,cosHalfDelta),32)
L0=Q(6273095,2000000)  # 3.1365475
assert Q(L_C16[0],SCALE) > L0

# Roots z_j=exp(i*j*pi/16), rigorously enclosed by repeated multiplication.
z=[(ONE,ZERO)]
base=(COS_PI_16,SIN_PI_16)
for _ in range(16): z.append(Cmul(z[-1],base))
assert z[16][0][0] <= -SCALE <= z[16][0][1]
assert z[16][1][0] <= 0 <= z[16][1][1]

# ---------------------------------------------------------------------------
# Scalar inequalities used in the analytic proof
# ---------------------------------------------------------------------------

def qiv(x: Q): return (x,x)
def lt_iv(x,y): return x[1] < y[0]
def gt_iv(x,y): return x[0] > y[1]

# 30 sin(pi/30) < L0.
pi30=div_iv_pos_q(PI,Q(30))
assert scale_iv(sin_iv_small(pi30),Q(30))[1] < L0

# Deficit bounds from the algebraic sin(pi/32) interval.
UP=(Q(32*SIN_PI_32[0],SCALE), Q(32*SIN_PI_32[1],SCALE))
UZ=scale_iv(UP,Q(2))
EPS=sub_iv(UP,qiv(L0))
DZ=sub_iv(UZ,qiv(2*L0))
assert EPS[1] < Q(991,10**9)
assert DZ[1] < Q(1982,10**9)

# F_omega(t)=2sin(t/2)+62sin((2pi-t)/62), at t=.18,.21.
def Fomega(tq):
    t=qiv(tq)
    first=scale_iv(sin_iv_small(div_iv_pos_q(t,Q(2))),Q(2))
    second_arg=div_iv_pos_q(sub_iv(scale_iv(PI,Q(2)),t),Q(62))
    second=scale_iv(sin_iv_small(second_arg),Q(62))
    return add_iv(first,second)
assert Fomega(Q(18,100))[1] < 2*L0
assert Fomega(Q(21,100))[1] < 2*L0

# r>0.99998 and |delta|<.005.
sin009=sin_point_bounds(Q(9,100))
cos0005=cos_point_bounds(Q(5,1000))
assert Q(2)*Q(2,100000)*sin009[0] > DZ[1]
r0=Q(99998,100000)
delta_loss=Q(2)*r0*sin009[0]*(Q(1)-cos0005[1])
assert delta_loss > DZ[1]

# KKT angular contradiction: ratio < .00585 and total projective distance < .0142.
ratio_hi=sin_point_bounds(Q(105,1000))[1] / sin009[0] * sin_point_bounds(Q(5,1000))[1]
assert ratio_hi < Q(585,100000)
assert PI[1]/2*ratio_hi + Q(5,1000) < Q(142,10000)
assert Q(142,10000) < Q(17,100)

# F_alpha(t)=2sin(t/2)+30sin((pi-t)/30), t=.189,.204.
def Falpha(tq):
    t=qiv(tq)
    first=scale_iv(sin_iv_small(div_iv_pos_q(t,Q(2))),Q(2))
    second_arg=div_iv_pos_q(sub_iv(PI,t),Q(30))
    second=scale_iv(sin_iv_small(second_arg),Q(30))
    return add_iv(first,second)
assert Falpha(Q(189,1000))[1] < L0
assert Falpha(Q(204,1000))[1] < L0

# ||x||_2 < .0065 and Poincare constant < 26.1.
sin00945=sin_point_bounds(Q(189,2000))
assert Q(4)*EPS[1] / sin00945[0] < Q(65,10000)**2
sinpi32_lo=Q(SIN_PI_32[0],SCALE)
assert Q(1)/(Q(4)*sinpi32_lo*sinpi32_lo) < Q(261,10)

# Uniform constants for the KKT uniqueness proof.
m_lo=Q(2)*sin00945[0]*sinpi32_lo*sinpi32_lo
assert m_lo > Q(18,10000)
grad_hi=sin_point_bounds(Q(102,1000))[1]*Q(65,10000)
assert grad_hi < Q(663,10**6)
assert Q(65,10000)/(Q(2)*sinpi32_lo) < Q(332,10000)
assert Q(3317,1000)**2 > 11
assert Q(2)*Q(3317,1000)*Q(332,10000) < Q(2203,10000)
assert Q(2)*(Q(11)-Q(12203,10000)) > Q(442,100)**2
assert Q(663,10**6)/Q(442,100) < Q(151,10**6)
assert Q(18,10000) > Q(2)*Q(151,10**6)

# ---------------------------------------------------------------------------
# Exact finite scan of all normalized half-codes c_0=+1
# ---------------------------------------------------------------------------
N=16
K=[[min(j,k)*(N-max(j,k)) for k in range(1,N)] for j in range(1,N)]

def quad_num(v,w):
    # Encloses v^T D^{-1} w as N/(16*SCALE^2).
    lo=hi=0
    for j in range(15):
        for k in range(15):
            p=rawprod(v[j],w[k]); kk=K[j][k]
            lo += kk*p[0]; hi += kk*p[1]
    return lo,hi

R_num,R_den=13,2000       # R=.0065
C_num,C_den=261,10        # C=26.1
survivors=[]
worst=None
common_den=4*SCALE*R_den*R_den*C_den

for bits in range(1<<15):
    c=[1]+[1 if ((bits>>(j-1))&1)==0 else -1 for j in range(1,16)]

    # b_c=sum c_j(z_{j+1}-z_j), interval scaled by SCALE.
    br=(0,0); bi=(0,0)
    for j in range(16):
        br=Iadd(br,Imul_int(Isub(z[j+1][0],z[j][0]),c[j]))
        bi=Iadd(bi,Imul_int(Isub(z[j+1][1],z[j][1]),c[j]))
    brmin=Iminabs(br); bimin=Iminabs(bi)
    bnorm_num=isqrt(brmin*brmin+bimin*bimin)  # lower bound / SCALE

    # B columns q_j=i(c_{j-1}-c_j)z_j, j=1,...,15.
    vr=[]; vi=[]
    for j in range(1,16):
        a=c[j-1]-c[j]
        vr.append(Imul_int(Ineg(z[j][1]),a))
        vi.append(Imul_int(z[j][0],a))

    M00=quad_num(vr,vr); M11=quad_num(vi,vi); M01=quad_num(vr,vi)
    trU=M00[1]+M11[1]
    diff=Isub(M00,M11)
    rad=ceil_sqrt(Imaxabs(diff)**2+4*Imaxabs(M01)**2)
    lam_num=cdiv(trU+rad,2)   # lambda_max <= lam_num/(16*SCALE^2)
    op_num=ceil_sqrt(lam_num) # sigma_c <= op_num/(4*SCALE)

    # Exact lower numerator for |b|-R*sigma-C*R^2.
    margin_num=(bnorm_num*4*R_den*R_den*C_den
                -R_num*op_num*R_den*C_den
                -C_num*R_num*R_num*4*SCALE)
    code=tuple(c)
    if margin_num <= 0:
        survivors.append(code)
    elif worst is None or margin_num < worst[0]:
        worst=(margin_num,code)

assert len(survivors)==16
REP=tuple(1 if ch=='+' else -1 for ch in '+--+-++-+--+-++-')
FULL=REP+tuple(-x for x in REP)
orbit=set()
for k in range(32):
    rot=FULL[k:]+FULL[:k]
    reflected=tuple(rot[(-j)%32] for j in range(32))
    for arr in (rot,reflected):
        half=arr[:16]
        if half[0]<0: half=tuple(-x for x in half)
        orbit.add(tuple(half))
assert set(survivors)==orbit and len(orbit)==16
assert worst is not None
worst_margin=Q(worst[0],common_den)
assert worst_margin > Q(2111,10**6)  # >0.002111

# Human-readable output.
def code_string(c): return ''.join('+' if x==1 else '-' for x in c)
print('Bingane C16 lower interval:')
print('  [', Q(L_C16[0],SCALE), ',', Q(L_C16[1],SCALE), ']')
print('  hence L(C16) > 3.1365475')
print('Scalar analytic inequalities: PASS')
print('Normalized half-codes scanned:',1<<15)
print('Survivors:',len(survivors))
for c in sorted(survivors): print(' ',code_string(c))
print('All survivors form the dihedral orbit of',code_string(REP))
print('Worst excluded code:',code_string(worst[1]))
print('Certified exclusion margin >',float(worst_margin))
print('ALL ASSERTIONS PASSED')
