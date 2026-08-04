from pathlib import Path
import numpy as np, subprocess, sys

HERE = Path(__file__).resolve().parent
N = 32
th = np.pi / N
xi = np.exp(1j * np.arange(N + 1) * th)
E = xi[1:] - xi[:-1]

REPS = ("++++--+--+-+-+-+++---+-+-+-++-++",
        "+++-+-+-++--+--+--++-+-+-+++--++",
        "++-++-+-+-+-+--++--+-+-+-+-++-++")
AXIAL = "+-++--+-+-+---++--+++-+-+-++--+-"


def signs(code):
    return np.array([1 if ch == '+' else -1 for ch in code])


def b(code):
    return complex(np.sum(signs(code) * E))


def negate(c):
    return ''.join('+' if ch == '-' else '-' for ch in c)


def orbit(code):
    full = code + negate(code)
    out = set()
    for base in (full, full[::-1]):
        for r in range(2 * N):
            z = base[r:] + base[:r]
            h = z[:N]
            if h[0] == '-':
                h = negate(h)
            out.add(h)
    return out


mine = set()
vals = {}
for line in (HERE / 'independent_survivors.txt').open():
    c, v = line.split()
    mine.add(c)
    vals[c] = float(v)

orbs = [orbit(r) for r in REPS]
print("independent survivors:", len(mine))
print("union of the 3 claimed orbits equals my survivor set:",
      mine == set().union(*orbs))
print("orbit sizes:", [len(o) for o in orbs],
      "| pairwise disjoint:",
      not (orbs[0] & orbs[1] or orbs[0] & orbs[2] or orbs[1] & orbs[2]))
print("axial code of section 3 lies in orbit B:", AXIAL in orbs[1])
print()
for name, o, r in zip("ABC", orbs, REPS):
    bs = [abs(b(c)) for c in o]
    print(f"orbit {name}: rep={r}  |b_c| = {min(bs):.6e} .. {max(bs):.6e}"
          f"  (constant over orbit: {np.ptp(bs) < 1e-14})")
print()
print("Section 8 claims |b_A| > 2.31e-5 and |b_C| > 1.66e-5:")
print("   |b_A| =", f"{abs(b(REPS[0])):.8e}", " |b_B| =", f"{abs(b(REPS[1])):.8e}",
      " |b_C| =", f"{abs(b(REPS[2])):.8e}")
print()
print("necessary bound from the proof: |b_c| <= 8*R0 + 256*R0^2 = %.6e"
      % (8 * 3.47e-6 + 256 * 3.47e-6 ** 2))
print("threshold actually used in the scan T0 = 3.41e-5")
