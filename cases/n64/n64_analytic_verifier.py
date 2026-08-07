#!/usr/bin/env python3
from fractions import Fraction as Q
import json
from pathlib import Path
import re
N=64
EPS=Q(284,10**25) # 2.84e-23
CODE='-++++++-----+--+-+++--+---+--+++---++-+++-++---+-++-+++++------+'
S_DEN=10**40
S_INTS=[113553801169565658142384194265, -9880196998990015085044055410, -133314195167545688312472305084, -256748193336101361539900554759, -380182191504657034767328804433, -503616189673212707994757054108, -627050187841768381222185303782, -527080630656387284719934645129, -427111073471006188217683986476, -327141516285625091715433327823, -227171959100243995213182669170, -127202401914862898710932010517, -224518637743365449610080769943, -131254926381658612514512915081, -37991215019951775418945060220, -120535648883632559553294480236, -35302307098282933183053776754, -109412180028470470567041786873, -183522052958658007951029796992, -257631925888845545335017807111, -199919646362187568673185756696, -142207366835529592011353706281, -196345090812107605367512622967, -149035276097607335865970026399, -101725461383107066364427429832, -54415646668606796862884833264, -75982724031176814146623311785, -39897017973167708772317662954, -3811311915158603398012014123, -2540874610105735598674676082, -1270437305052867799337338041]
TL=Q(999999999987873350,10**18)
TU=Q(999999999987873353,10**18)
A_GAP=Q(45,1000)
B_GAP=Q(54,1000)
R0=Q(711,10**13) # 7.11e-11
T_S=Q(17,10**9)  # coarse regular-root sum threshold 1.7e-8

def atan_interval(x,terms):
    s=Q(0); p=x; vals=[]
    for k in range(terms+1):
        s += p/Q(2*k+1) if k%2==0 else -p/Q(2*k+1)
        if k in (terms-1,terms): vals.append(s)
        p*=x*x
    return min(vals),max(vals)
def pi_interval():
    aL,aU=atan_interval(Q(1,5),30); bL,bU=atan_interval(Q(1,239),10)
    return 16*aL-4*bU,16*aU-4*bL
PI_L,PI_U=pi_interval()
assert PI_L>3 and PI_U<Q(22,7)

def sin_point_interval(x,terms=24):
    assert Q(0)<=x<=PI_U/2
    s=Q(0); term=x; vals=[]
    for k in range(terms+1):
        s+=term
        if k in (terms-1,terms): vals.append(s)
        term *= -x*x/Q((2*k+2)*(2*k+3))
    return min(vals),max(vals)
def sin_interval(lo,hi):
    assert Q(0)<=lo<=hi<=PI_U/2
    l,_=sin_point_interval(lo); _,u=sin_point_interval(hi)
    return l,u
def add_iv(a,b):return a[0]+b[0],a[1]+b[1]
def sub_iv(a,b):return a[0]-b[1],a[1]-b[0]
def scale_iv(k,a):return (k*a[0],k*a[1]) if k>=0 else (k*a[1],k*a[0])

def signs(s): return [1 if ch=='+' else -1 for ch in s]
c=signs(CODE)
assert len(c)==N and all(c[N-1-j]==-c[j] for j in range(N))
s=[Q(0)]*(N+1)
for j,v in enumerate(S_INTS,1):s[j]=Q(v,S_DEN)
s[32]=Q(0)
for j in range(33,N):s[j]=-s[N-j]
assert all(s[N-j]==-s[j] for j in range(N+1))
a=[None]+[c[j-1]-c[j] for j in range(1,N)]

def H_iv(t):
    ans=(Q(a[32]),Q(a[32]))
    for j in range(1,32):
        lo=Q(j,N)*PI_L+t*s[j]; hi=Q(j,N)*PI_U+t*s[j]
        ans=add_iv(ans,scale_iv(Q(2*a[j]),sin_interval(lo,hi)))
    return ans
HL=H_iv(TL); HU=H_iv(TU)
assert HL[0]>0 and HU[1]<0

# Every t in bracket has ordered gaps and a certified lower perimeter.
P=(Q(0),Q(0))
for j in range(N):
    q=s[j+1]-s[j]
    prods=(TL*q,TU*q)
    lo=PI_L/Q(2*N)+min(prods)/2
    hi=PI_U/Q(2*N)+max(prods)/2
    assert 2*lo>A_GAP and 2*hi<B_GAP
    P=add_iv(P,scale_iv(Q(2),sin_interval(lo,hi)))
U=scale_iv(Q(2*N),sin_interval(PI_L/Q(2*N),PI_U/Q(2*N)))
DEF=sub_iv(U,P)
assert DEF[1]<EPS

# Exact endpoint checks used by the geometric localization argument.
U63=scale_iv(Q(126),sin_interval(PI_L/Q(126),PI_U/Q(126)))
U64_MINUS_U63=sub_iv(U,U63)
assert U64_MINUS_U63[0]>Q(3,2*N**3)

def width_deficit(which):
    t_half=(Q(which,4*N)*PI_L,Q(which,4*N)*PI_U)
    coef=(Q(1)-Q(which,4*N))/Q(127)
    avg_half=(coef*PI_L,coef*PI_U)
    ans=scale_iv(Q(2),U)
    ans=sub_iv(ans,scale_iv(Q(2),sin_interval(*t_half)))
    return sub_iv(ans,scale_iv(Q(2*127),sin_interval(*avg_half)))
WIDTH_DEFICITS=(width_deficit(1),width_deficit(3))
WIDTH_COARSE=Q(9,N**2*(2*N)*16)
assert all(iv[0]>WIDTH_COARSE for iv in WIDTH_DEFICITS)

def gap_deficit(alpha):
    first=scale_iv(Q(2),sin_interval(alpha/2,alpha/2))
    lo=(PI_L-alpha)/(2*(N-1)); hi=(PI_U-alpha)/(2*(N-1))
    rest=scale_iv(Q(2*(N-1)),sin_interval(lo,hi))
    return sub_iv(sub_iv(U,first),rest)
GAP_DEFICITS=(gap_deficit(A_GAP),gap_deficit(B_GAP))
assert all(iv[0]>EPS for iv in GAP_DEFICITS)

# Certify every fixed-point weight embedded in the exhaustive C++ scan.
cpp=(Path(__file__).resolve().parent/'n64_code_fixed.cpp').read_text()
fixed={}
for name in ('X','Y'):
    match=re.search(rf'const int64_t {name}\[32\]=\{{([^}}]+)\}};',cpp)
    assert match is not None
    fixed[name]=[int(x) for x in match.group(1).split(',')]
    assert len(fixed[name])==32
for j in range(32):
    x=scale_iv(Q(2),sin_interval(Q(2*j+1,128)*PI_L,Q(2*j+1,128)*PI_U))
    y=scale_iv(Q(2),sin_interval(Q(63-2*j,128)*PI_L,Q(63-2*j,128)*PI_U))
    xi,yi=fixed['X'][j],fixed['Y'][j]
    assert x[0]>=Q(xi-1,10**16) and x[1]<=Q(xi+1,10**16)
    assert y[0]>=Q(yi-1,10**16) and y[1]<=Q(yi+1,10**16)

# Analytic constants for localization and uniform code screen.
def sin_lower(x):return x-x**3/Q(6)
K0=sin_lower(A_GAP/2)/4
assert K0*R0*R0>EPS
# Dirichlet Poincare constant <512 follows from 2*sin(pi/128)>0.049 below.
# Twisted linear coefficient sqrt(128)cos(pi/128)<12.
assert Q(23,2)*R0+512*R0*R0 < Q(82,10**11) # b threshold 8.2e-10
# Divide by |xi-1|=2sin(pi/128)>0.049.
root_gap_l,_=sin_interval(PI_L/Q(2*N),PI_U/Q(2*N))
assert 2*root_gap_l>Q(49,1000)
assert Q(82,10**11)/Q(49,1000) < T_S

# Saturation constants become extremely loose at this deficit.
# Missing two difference-body vertices costs > 3/(2N^3).
assert Q(3,2*N**3)>EPS
# Normal width outside (mu/2,3mu/2) costs > 9/(N^2*(2N)*16).
assert Q(9,N**2*(2*N)*16)>2*EPS
# Radial and angular errors are far below the same general contradiction constants.
assert 2*N*EPS<Q(1,2)
assert 4*N*EPS<Q(1,2**20)

report={
 'pi_width': float(PI_U-PI_L),
 'H_left_lower': float(HL[0]),
 'H_right_upper': float(HU[1]),
 'deficit_upper': float(DEF[1]),
 'U64_minus_U63_lower':float(U64_MINUS_U63[0]),
 'width_deficit_boundary_lowers':[float(iv[0]) for iv in WIDTH_DEFICITS],
 'gap_deficit_boundary_lowers':[float(iv[0]) for iv in GAP_DEFICITS],
 'fixed_weights_certified':64,
 'epsilon':float(EPS),
 'R0':float(R0),
 'S_screen':float(T_S),
}
print(json.dumps(report,indent=2))
print('ALL ANALYTIC ASSERTIONS PASSED')
