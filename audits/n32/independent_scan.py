"""Independent re-check of the n=32 code screening.

Different implementation from the bundled verifier:
  * float64 (error ~1e-15, threshold ~3.4e-5, so the margin is ~1e10)
  * DIFFERENT meet-in-the-middle split: A = indices 0..7 and 16..23,
                                        B = indices 8..15 and 24..31
  * full |A| x |B| product (2^31 pairs) with no range pruning
"""
from pathlib import Path
import numpy as np, time, sys

HERE = Path(__file__).resolve().parent
N = 32
th = np.pi / N
xi = np.exp(1j * np.arange(N + 1) * th)
E = xi[1:] - xi[:-1]              # e_j = xi^{j+1} - xi^j, j=0..31

T0 = 3.41e-5

IA = [j for j in range(N) if (j // 8) % 2 == 0]   # 0..7, 16..23
IB = [j for j in range(N) if (j // 8) % 2 == 1]   # 8..15, 24..31
assert sorted(IA + IB) == list(range(N))
assert 0 in IA


def block(indices, fix_first_plus):
    """All signed sums over `indices`, plus the mask of which got '-'."""
    sums = np.zeros(1, dtype=np.complex128)
    masks = np.zeros(1, dtype=np.int64)
    idx = list(indices)
    if fix_first_plus:
        sums = sums + E[idx[0]]
        idx = idx[1:]
    for j in idx:
        sums = np.concatenate([sums + E[j], sums - E[j]])
        masks = np.concatenate([masks, masks | (1 << j)])
    return sums, masks


t0 = time.time()
SA, MA = block(IA, True)    # 2^15 entries (c_0 pinned to +1)
SB, MB = block(IB, False)   # 2^16 entries
print("block sizes:", SA.size, SB.size, "-> pairs:", SA.size * SB.size)

ax = SA.real.copy(); ay = SA.imag.copy()
bx = SB.real.copy(); by = SB.imag.copy()

thr2 = T0 * T0
hits = []
CH = 64
for s in range(0, SA.size, CH):
    e = min(s + CH, SA.size)
    dx = ax[s:e, None] + bx[None, :]
    dy = ay[s:e, None] + by[None, :]
    d2 = dx * dx + dy * dy
    ii, jj = np.nonzero(d2 <= thr2)
    for i, j in zip(ii, jj):
        hits.append((MA[s + i] | MB[j], float(np.sqrt(d2[i, j]))))
    if (s // CH) % 128 == 0:
        print(f"  {e}/{SA.size} A-rows, {len(hits)} hits, {time.time()-t0:.0f}s",
              flush=True)

print("elapsed %.1fs" % (time.time() - t0))
print("survivors found:", len(hits))


def mask_code(m):
    return ''.join('-' if (m >> j) & 1 else '+' for j in range(N))


codes = sorted({mask_code(m) for m, _ in hits})
print("distinct codes:", len(codes))
np.save(HERE / 'independent_survivor_norms.npy',
        np.array([v for _, v in hits]))
with (HERE / 'independent_survivors.txt').open('w') as f:
    for m, v in sorted(hits, key=lambda r: r[1]):
        f.write(f"{mask_code(m)} {v:.6e}\n")
for m, v in sorted(hits, key=lambda r: r[1])[:8]:
    print(f"  {mask_code(m)}  |b|={v:.6e}")
