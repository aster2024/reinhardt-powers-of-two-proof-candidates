import numpy as np
from mpmath import mp, mpf, sin, cos, pi, findroot

# ---------------------------------------------------------------- sigma_c
def sigma_sq(code, N):
    c = np.array([1 if ch == '+' else -1 for ch in code], dtype=float)
    a = c[:-1] - c[1:]
    th = np.pi / N
    j = np.arange(1, N)
    B = np.vstack([-a * np.sin(j * th), a * np.cos(j * th)])
    jj, kk = np.meshgrid(j, j, indexing='ij')
    Dinv = np.minimum(jj, kk) * (N - np.maximum(jj, kk)) / N
    K = B @ Dinv @ B.T
    return max(np.linalg.eigvalsh(K))

def negate(c):
    return ''.join('+' if ch == '-' else '-' for ch in c)

def orbit(code, N):
    full = code + negate(code)
    out = set()
    for base in (full, full[::-1]):
        for r in range(2 * N):
            z = base[r:] + base[:r]
            h = z[:N]
            if h[0] == '-': h = negate(h)
            out.add(h)
    return out

REPS32 = ("++++--+--+-+-+-+++---+-+-+-++-++",
          "+++-+-+-++--+--+--++-+-+-+++--++",
          "++-++-+-+-+-+--++--+-+-+-+-++-++")
print("=== sigma_c^2, independent eigenvalue computation (n=32) ===")
print("verifier's certified Gershgorin upper bounds: 16.364240, 16.453073, 16.254762")
for nm, r in zip("ABC", REPS32):
    o = orbit(r, 32)
    vals = [sigma_sq(x, 32) for x in o]
    print(f"  orbit {nm}: sigma^2 = {max(vals):.12f}"
          f"   (spread over orbit: {np.ptp(vals):.2e})")

# ------------------------------------------------- section 3 feasible point
mp.dps = 80
N = 32
AX = "+-++--+-+-+---++--+++-+-+-++--+-"
c = [1 if ch == '+' else -1 for ch in AX]
S_INTS = (-84838116394748, -7759417296262, -90254740014369,
          -172750062732476, -104929943709607, -37109824686737,
          -104570029200993, -46261230637792, -102998890440414,
          -56520596873701, -100433725625308, -67650481919513,
          -34867238213719, -2083994507925, -1041997253963)
s = [mpf(0)] * (N + 1)
for i, v in enumerate(S_INTS, 1): s[i] = mpf(v) / mpf(10) ** 20
s[16] = mpf(0)
for j in range(17, N): s[j] = -s[N-j]
s[32] = mpf(0)
a = [None] + [c[j-1] - c[j] for j in range(1,N)]
th = pi/N

def H(t):
    return 2*sum(a[j]*sin(j*th+t*s[j]) for j in range(1,16)) + a[16]

tstar = findroot(H, (mpf('0.99999953359'), mpf('0.99999953360')), solver='bisect')
print("\n=== Section 3 feasible construction (n=32) ===")
print("root t* =", mp.nstr(tstar, 30))
phi = [j*th+tstar*s[j] for j in range(N+1)]
gaps = [phi[j+1]-phi[j] for j in range(N)]
perim = sum(2*sin(g/2) for g in gaps)
U = 2*N*sin(pi/(2*N))
G = sum(c[j]*(mp.e**(1j*phi[j+1])-mp.e**(1j*phi[j])) for j in range(N))
print("min gap =", mp.nstr(min(gaps),20), "max gap =", mp.nstr(max(gaps),20))
print("closure residual =", mp.nstr(abs(G),20))
print("U =", mp.nstr(U,50))
print("perimeter =",mp.nstr(perim,50))
print("deficit =",mp.nstr(U-perim,30))

# --------------------------------------------- true optimum for fixed code
def solve_code(code, n, init_phi=None):
    cs = [1 if ch == '+' else -1 for ch in code]
    thn = pi/n
    def eqs(*v):
        ph=[mpf(0)]+list(v[:n-1])+[pi]
        lam=v[n-1:n+1]
        out=[]
        for k in range(1,n):
            grad = cos((ph[k]-ph[k-1])/2)-cos((ph[k+1]-ph[k])/2)
            ak=cs[k-1]-cs[k]
            out.append(grad + lam[0]*(-ak*sin(ph[k])) + lam[1]*(ak*cos(ph[k])))
        gr=sum(cs[j]*(cos(ph[j+1])-cos(ph[j])) for j in range(n))
        gi=sum(cs[j]*(sin(ph[j+1])-sin(ph[j])) for j in range(n))
        out += [gr,gi]
        return tuple(out)
    if init_phi is None:
        start=[mpf(k)*thn for k in range(1,n)]+[mpf('0'),mpf('0')]
    else:
        start=list(init_phi[1:-1])+[mpf('0'),mpf('0')]
    sol=findroot(eqs, tuple(start), tol=mpf('1e-55'), maxsteps=100)
    ph=[mpf(0)]+[sol[i] for i in range(n-1)]+[pi]
    vals=eqs(*sol)
    return sum(2*sin((ph[j+1]-ph[j])/2) for j in range(n)), ph, max(abs(v) for v in vals)

print("\n=== fixed-code KKT solves ===")
p32,ph32,res32=solve_code(AX,32,phi)
print("n32 p =",mp.nstr(p32,60))
print("n32 deficit =",mp.nstr(U-p32,30),"max eq res",mp.nstr(res32,10))
REP16="+--+-++-+--+-++-"
p16,ph16,res16=solve_code(REP16,16)
U16=2*16*sin(pi/32)
print("n16 p =",mp.nstr(p16,50))
print("n16 deficit =",mp.nstr(U16-p16,30),"max eq res",mp.nstr(res16,10))
print("n16 first gaps", [mp.nstr(ph16[j+1]-ph16[j],25) for j in range(4)])
