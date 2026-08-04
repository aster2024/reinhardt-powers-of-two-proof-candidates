from fractions import Fraction as Q
from pathlib import Path
import os, sys
import importlib.util, math, re

ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT / 'cases' / 'n64'
os.chdir(CASE)
sys.path.insert(0, str(CASE))

spec=importlib.util.spec_from_file_location('A', CASE / 'n64_analytic_verifier.py')
A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
N=64; M=10**16

def add(a,b): return a[0]+b[0],a[1]+b[1]
def sub(a,b): return a[0]-b[1],a[1]-b[0]
def scale(k,a): return A.scale_iv(k,a)
def sin_rat_pi(num,den): return A.sin_interval(Q(num,den)*A.PI_L,Q(num,den)*A.PI_U)
def sin_general(lo,hi): return A.sin_interval(lo,hi)

# Missing-two-vertices comparison: U_64-U_63.
U64=scale(128,sin_rat_pi(1,128))
U63=scale(126,sin_rat_pi(1,126))
diff=sub(U64,U63)
assert diff[0] > Q(3,2*N**3)

# Width-Jensen deficit with m=128, regular width mu=pi/64.
# Analytically D'(t) has the sign of t-mu, so outside the box minima are at endpoints.
def width_deficit(which):
    # t=which*mu/2, which=1 or 3
    # h(t)=2sin(t/2); other average=(2pi-t)/(127)
    t_half=(Q(which,4*N)*A.PI_L,Q(which,4*N)*A.PI_U)
    # avg/2=(2pi-t)/(2*127)=(1-which/(4N))*pi/127
    coef=(Q(1)-Q(which,4*N))/Q(127)
    avg_half=(coef*A.PI_L,coef*A.PI_U)
    val=scale(2,U64)
    val=sub(val,scale(2,sin_general(*t_half)))
    val=sub(val,scale(2*127,sin_general(*avg_half)))
    return val
wd1=width_deficit(1); wd3=width_deficit(3)
coarse=Q(9,N**2*(2*N)*16)
assert wd1[0]>coarse and wd3[0]>coarse

# One-gap Jensen deficits at alpha=.045 and .054. D'(alpha) changes sign at theta.
def gap_deficit(alpha):
    h1=scale(2,A.sin_interval(alpha/2,alpha/2))
    # other gap average=(pi-alpha)/(N-1), half of it below
    lo=(A.PI_L-alpha)/(2*(N-1)); hi=(A.PI_U-alpha)/(2*(N-1))
    h2=scale(2*(N-1),A.sin_interval(lo,hi))
    return sub(sub(U64,h1),h2)
g45=gap_deficit(Q(45,1000)); g54=gap_deficit(Q(54,1000))
assert g45[0]>A.EPS and g54[0]>A.EPS

# Certify all C++ fixed-point weights, omitted from the bundled verifier.
text=(CASE / 'n64_code_fixed.cpp').read_text(); arrays={}
for name in ('X','Y'):
    m=re.search(rf'const int64_t {name}\[32\]=\{{([^}}]+)\}};',text)
    arrays[name]=[int(x) for x in m.group(1).split(',')]
for j in range(32):
    x=A.scale_iv(2,sin_rat_pi(2*j+1,128))
    y=A.scale_iv(2,sin_rat_pi(63-2*j,128))
    xi,yi=arrays['X'][j],arrays['Y'][j]
    assert x[0]>=Q(xi-1,M) and x[1]<=Q(xi+1,M)
    assert y[0]>=Q(yi-1,M) and y[1]<=Q(yi+1,M)

# Corrected fully generic residual padding: coordinate errors <=2N, Euclidean <=2N*sqrt(2)<2N*3/2=192.
# Use 256 integer units and recheck all five nonwinning lower bounds.
spec2=importlib.util.spec_from_file_location('P', CASE / 'n64_post_verifier.py')
P=importlib.util.module_from_spec(spec2);spec2.loader.exec_module(P)
for i in range(1,6):
    x,y=P.residual_int(P.REPS[i])
    t=P.ceil_q(M*P.BLOW[i])+256
    assert x*x+y*y>t*t

print({
 'U64_minus_U63_lower':float(diff[0]),
 'width_deficit_boundary_lowers':[float(wd1[0]),float(wd3[0])],
 'gap_deficit_boundary_lowers':[float(g45[0]),float(g54[0])],
 'fixed_weights_certified':64,
 'nonbest_residuals_pass_padding':256,
})
print('ALL SUPPLEMENTAL ASSERTIONS PASSED')
