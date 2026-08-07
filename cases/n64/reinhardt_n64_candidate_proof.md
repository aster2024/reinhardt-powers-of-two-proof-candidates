# A computer-assisted candidate proof for the maximum-perimeter small 64-gon

**Status:** complete candidate proof, 5 August 2026; independently unreviewed.

This document gives a candidate proof that the maximum perimeter among convex
64-gons of diameter at most one is attained by a unique congruence class.  The
proof uses exact rational interval arithmetic and an exhaustive integer code
certificate.  No binary floating-point comparison is used in a proof decision.

The winning saturated half-code is the dihedral orbit of

```text
-++++++-----+--+-+++--+---+--+++---++-+++-++---+-++-+++++------+
```

which is the axially symmetric extension of the published quarter code

```text
-++++++-----+--+-+++--+---+--+++
```

A non-rigorous high-precision evaluation of the unique maximizer is

\[
3.14127725093277286806199141550246829795626209630809641\ldots .
\]

The proof does not rely on this decimal expansion.

Put

\[
 n=64,\qquad \theta=\frac\pi{64},\qquad
 U=128\sin\frac\pi{128},\qquad
 \varepsilon=2.84\cdot10^{-23}.
\]

## 1. Difference-body reduction and reconstruction

For a convex polygon `P`, let `Z=P-P`.  Then

\[
 Z\subseteq \overline B(0,1),\qquad
 \operatorname{per}(Z)=2\operatorname{per}(P),
\]

and `Z` is centrally symmetric with at most 128 edges.

### Lemma 1.1 (reconstruction)

Let `Z` be a centrally symmetric convex 128-gon with all listed vertices
genuine.  Label one half cyclically by

\[
 z_0,\ldots,z_{63},\qquad z_{64}=-z_0,
\]

and put `e_j=z_{j+1}-z_j`.  If a code `c_j in {+-1}` satisfies

\[
 \sum_{j=0}^{63}c_je_j=0,                                      \tag{1.1}
\]

then the vectors `c_j e_j`, sorted by direction, are the edges of a closed
convex 64-gon `P_c`, and

\[
 P_c-P_c=Z.
\]

**Proof.**  The selected vectors contain exactly one vector from every
antipodal pair `{e_j,-e_j}`.  They are nonzero and have distinct directions.
Equation (1.1) gives closure.  Sorting by direction therefore produces a
strictly convex polygon.  Merging its edge list with the edge list of its
negative gives the full cyclic edge list of `Z`; centered polygons with the
same cyclic edge list coincide.  ∎

Summation by parts rewrites (1.1) as

\[
 -(c_0+c_{63})z_0+\sum_{j=1}^{63}(c_{j-1}-c_j)z_j=0,             \tag{1.2}
\]

whose coefficients lie in `{0,+-2}`.

### Lemma 1.2 (at most one interior half-vertex)

At a local maximum with a 128-vertex difference body, at most one of
`z_0,...,z_63` lies strictly inside the unit disk, and such a vertex has a
nonzero coefficient in (1.2).

**Proof.**  If two interior vertices have nonzero coefficients `a_r,a_s`, use

\[
 \delta z_r=a_sh,\qquad \delta z_s=-a_rh.
\]

If an interior vertex has zero coefficient, move it alone.  For sufficiently
small perturbations of either sign, the disk constraints, cyclic order and
genuine-vertex conditions remain valid.  Along the perturbation the perimeter
is a sum of norms of affine functions, hence convex.  A generic `h` makes one
affected norm strictly convex, contradicting a two-sided local maximum.  ∎

## 2. A rigorously feasible near-regular 64-gon

Use the winning axial code `c`, satisfying

\[
 c_{63-j}=-c_j.
\]

Let `s_0=s_32=s_64=0`, `s_{64-j}=-s_j`, and set

\[
 s_j=m_j/10^{40}\qquad(1\le j\le31),
\]

where the integer list `m_j` is contained in `n64_analytic_verifier.py`.  Put

\[
 \phi_j(t)=j\theta+t s_j.
\]

Writing `a_j=c_{j-1}-c_j`, axial symmetry makes the closure residual purely
imaginary:

\[
 G(t)=iH(t),\qquad
 H(t)=2\sum_{j=1}^{31}a_j\sin(j\theta+t s_j)+a_{32}.             \tag{2.1}
\]

Exact rational interval arithmetic proves

\[
 H(0.999999999987873350)>5.67\cdot10^{-28},
\]

\[
 H(0.999999999987873353)<-6.00\cdot10^{-28}.                    \tag{2.2}
\]

Hence an exact root exists in that bracket.  At every point in the bracket all
angle gaps lie in `(0.045,0.054)`.  The same interval calculation gives

\[
 U-\operatorname{per}(P(t_*))
 <2.835602\cdot10^{-23}<\varepsilon.                            \tag{2.3}
\]

Lemma 1.1 therefore reconstructs a genuinely feasible small 64-gon.

The verifier obtains `pi` from Machin's formula and bounds sine by alternating
Taylor intervals.  Its certified `pi` interval has width below
`6.05*10^-44`.

## 3. Saturation of every global maximizer

### 3.1 The difference body has 128 genuine vertices

If `Z` had fewer than 128 vertices, central symmetry would give at most 126,
so

\[
 \operatorname{per}(P)\le126\sin\frac\pi{126}=U_{63}.
\]

For `V(t)=2t sin(pi/(2t))`, the standard derivative estimate gives

\[
 U-U_{63}>\frac{3}{2\cdot64^3}>\varepsilon,
\]

contradicting (2.3).  Thus `Z` has 128 genuine vertices.

### 3.2 Near-regular normal cones

For a difference-body vertex write

\[
 z_j=r_je^{i\phi_j},
\]

and let `omega_j` be its normal-cone width, `eta_j` its normal-cone midpoint,
and `delta_j=eta_j-phi_j`.  Cauchy's formula is

\[
 \operatorname{per}(Z)=
 \sum_{j=0}^{127}2r_j\sin\frac{\omega_j}{2}\cos\delta_j,
 \qquad \sum_j\omega_j=2\pi.                                  \tag{3.1}
\]

Put `mu=pi/64`.  Fixing one width and averaging the other 127 widths with
Jensen's inequality, together with strong concavity, shows that a width outside

\[
 \frac\mu2<\omega_j<\frac{3\mu}{2}                              \tag{3.2}
\]

would cost more than `2 epsilon`.  The verifier checks a coarse lower bound

\[
 \frac{9}{64^2\cdot128\cdot16}>2\varepsilon.
\]

The nonnegative deficit decomposition in (3.1) then yields

\[
 1-r_j<128\varepsilon,
 \qquad 1-\cos\delta_j<256\varepsilon.                          \tag{3.3}
\]

Consequently

\[
 |\delta_j|<\frac\mu{16}.                                      \tag{3.4}
\]

Successive radial directions, and hence any two different half-vertices in
projective angle, are separated by more than

\[
 \frac\mu2-2\frac\mu{16}=\frac{3\mu}{8}.                       \tag{3.5}
\]

### 3.3 KKT excludes the last interior radius

By Lemma 1.2 suppose only `z_r` is interior.  Its coefficient is `a_r=+-2`.
MFCQ holds: move active circle vertices radially inward and use the free
2-dimensional displacement of `z_r` to repair (1.2).

For the objective `per(Z)/2`, the half-vertex gradient is

\[
 g_j=2\sin\frac{\omega_j}{2}e^{i\eta_j}.                        \tag{3.6}
\]

At `r`, stationarity gives `g_r=a_r lambda`, so the projective direction of
`lambda` is `eta_r`.  There is another `j != r` with `a_j !=0`; otherwise
(1.2) would force `z_r=0`.  This second vertex is on the unit circle.  Project
its KKT equation onto its circle tangent to obtain

\[
 |\sin(\eta_r-\phi_j)|
 =\frac{\sin(\omega_j/2)}{\sin(\omega_r/2)}|\sin\delta_j|.      \tag{3.7}
\]

By (3.2) the ratio is less than three.  Equations (3.4), (3.7) and
`sin t>=2t/pi` give

\[
 d_{\mathbb {RP}^1}(\phi_r,\phi_j)
 <\frac{(2+3\pi)\mu}{32}<\frac{3\mu}{8},
\]

contradicting (3.5).  Thus all 128 difference-body vertices lie on the unit
circle.  The original saturation conjecture is therefore proved for `n=64`.

## 4. Saturated angle model and localization

Rotate so that `phi_0=0`, `phi_64=pi`, and put

\[
 \alpha_j=\phi_{j+1}-\phi_j,
 \qquad
 \phi_j=j\theta+s_j,
 \qquad x_j=s_{j+1}-s_j.
\]

The objective is

\[
 F=\sum_{j=0}^{63}2\sin\frac{\alpha_j}{2},                     \tag{4.1}
\]

and closure is

\[
 G_c=\sum_{j=0}^{63}c_j(e^{i\phi_{j+1}}-e^{i\phi_j})=0.         \tag{4.2}
\]

A one-gap Jensen bound proves that every feasible point with deficit at most
`epsilon` satisfies

\[
 0.045<\alpha_j<0.054.                                          \tag{4.3}
\]

On this box, strong concavity gives

\[
 U-F\ge \kappa\|x\|_2^2,
 \qquad
 \kappa>\frac14\left(0.0225-\frac{0.0225^3}{6}\right).
\]

The exact verifier checks

\[
 \boxed{\|x\|_2<R_0:=7.11\cdot10^{-11}.}                       \tag{4.4}
\]

## 5. Uniform regular-residual screen

Let `xi=e^{i theta}` and

\[
 b_c=(\xi-1)\sum_{j=0}^{63}c_j\xi^j.                            \tag{5.1}
\]

Taylor expansion gives

\[
 G_c=b_c+L_c(s)+R_c(s).
\]

The Dirichlet Poincare inequality and `sin(pi/128)>1/45` give

\[
 |R_c(s)|\le512\|x\|_2^2.                                      \tag{5.2}
\]

The twisted-Dirichlet operator estimate, valid for general `n`, is

\[
 |L_c(s)|\le\sqrt{2n}\cos\frac\pi{2n}\,\|x\|_2.               \tag{5.3}
\]

For `n=64`, the coefficient is less than `23/2`.  Thus every competitive code
satisfies

\[
 |b_c|<8.2\cdot10^{-10}.                                       \tag{5.4}
\]

Since `|xi-1|=2sin(pi/128)>0.049`, it must satisfy

\[
 \left|\sum_{j=0}^{63}c_j\xi^j\right|<1.7\cdot10^{-8}.          \tag{5.5}
\]

## 6. Orthogonal reflection-pair decomposition

This is the new device that makes a full `n=64` code certificate feasible.
Pair indices `j` and `63-j`, for `0<=j<=31`, and put

\[
 \beta_j=\left(j-\frac{63}{2}\right)\theta,
\]

\[
 p_j=\frac{c_j+c_{63-j}}2,
 \qquad
 q_j=\frac{c_j-c_{63-j}}2.
\]

For each `j`, exactly one of `p_j,q_j` is zero and the other belongs to
`{+-1}`.  Direct calculation gives

\[
 \sum_{j=0}^{63}c_j\xi^j
 =2e^{63i\theta/2}
 \left(\sum_{j=0}^{31}p_j\cos\beta_j
       +i\sum_{j=0}^{31}q_j\sin\beta_j\right).                 \tag{6.1}
\]

Thus the real and imaginary coordinates after a fixed rotation are two
orthogonal one-dimensional ternary subset sums, with complementary supports.
Conversely every pair of ternary coefficient vectors with complementary
supports determines exactly one 64-bit code.  Hence (6.1) is a bijective
reparameterization of all

\[
 4^{32}=2^{64}
\]

half-codes, not a symmetry restriction.

## 7. Exact exhaustive code certificate

At scale `M=10^16`, the core analytic verifier proves that each of the 32 horizontal and 32
vertical weights in (6.1) lies within `1/M` of the integer printed in
`n64_code_fixed.cpp`.  Split each ternary sum into two blocks of 16 terms.
Each block has only

\[
 3^{16}=43,046,721
\]

states.  Sorting the two block lists and using a moving interval window
enumerates all coordinate sums within the conservative fixed-point radius.
Horizontal and vertical records are then matched exactly by complementary
support masks.

If the true complex sum satisfies (5.5), its fixed-point coordinate vector is
inside the integer circle used by the program; the accumulated rounding error
is less than 64 integer units.  Therefore the scan cannot omit a competitive
code.

The exact integer result is

\[
 \boxed{896\text{ half-codes survive}.}                         \tag{7.1}
\]

They form exactly six full-dihedral orbits, with orbit sizes

\[
 128,128,128,128,128,256.                                      \tag{7.2}
\]

Five orbits contain reflection-symmetric representatives; the last orbit is
generic.  The winning published axial code is the first orbit.

The recorded packaging run took 20.43 seconds and 1.47 GB; runtime is
environment dependent and is not a proof assertion.  The decomposition
replaces an infeasible direct scan of `2^63` normalized codes.

## 8. Exact elimination of the five nonwinning orbits

For a fixed code define

\[
 \sigma_c=\sup_{s\ne0}\frac{|L_c(s)|}{\|x\|_2}.
\]

As in the `n=32` proof,

\[
 \sigma_c^2=\lambda_{\max}(B_cD^{-1}B_c^T),                    \tag{8.1}
\]

where `D` is the 63-by-63 Dirichlet path matrix and

\[
 (D^{-1})_{jk}=\frac{\min(j,k)(64-\max(j,k))}{64}.              \tag{8.2}
\]

Exact fixed-point Gershgorin bounds prove for all six representatives

\[
 \sigma_c^2<36.                                                 \tag{8.3}
\]

For the five nonwinning representatives, exact regular-residual lower bounds
are respectively

\[
 |b_c|>
 4.86,\ 5.08,\ 6.28,\ 7.64,\ 8.23
 \quad\text{times }10^{-10}.                                   \tag{8.4}
\]

For these residual comparisons the core post-verifier uses a 256-unit
fixed-point padding, larger than the generic Euclidean rounding bound
`128*sqrt(2)`.  If one of these codes had deficit at most `epsilon`, equations (4.4), (5.2),
(8.3) would imply

\[
 \|x\|_2>
 \frac{|b_c|-512R_0^2}{6}.                                     \tag{8.5}
\]

Combining (8.5) with strong concavity gives deficits exceeding `epsilon`; the
smallest certified excess is

\[
 8.50\cdot10^{-24}.                                             \tag{8.6}
\]

Therefore only the winning orbit can contain a global maximizer.

## 9. Continuous global uniqueness inside the winning orbit

Use the axial representative from Section 2.  It has 27 switch indices

\[
 J=\{j:c_{j-1}\ne c_j\}.
\]

Let `f=-F` and `g=(Re G_c,Im G_c)`.

### 9.1 Strong convexity

The Hessian of `f` is the weighted Dirichlet path matrix with weights

\[
 w_j=\frac12\sin\frac{\alpha_j}{2}.
\]

Using (4.3) and `sin(pi/128)>1/45`, throughout the high-perimeter ball

\[
 \nabla^2f\succeq m_0I,
 \qquad
 m_0>
 \frac{2}{45^2}
 \left(0.0225-\frac{0.0225^3}{6}\right).                       \tag{9.1}
\]

### 9.2 Constraint Jacobian and multiplier bound

The nonzero columns of `Dg` are `a_j i e^{i phi_j}`.  Therefore

\[
 \sigma_{\min}(Dg)^2
 =2\left(27-\left|\sum_{j\in J}e^{2i\phi_j}\right|\right).     \tag{9.2}
\]

At the regular point the exact certificate proves

\[
 \left|\sum_{j\in J}e^{2ij\theta}\right|<5.                    \tag{9.3}
\]

Poincare gives `||s||<23R_0`, so the root sum changes by less than
`276R_0<1`.  Hence its norm remains below six and

\[
 \sigma_{\min}(Dg)>6.                                          \tag{9.4}
\]

The objective gradient obeys

\[
 \|\nabla f\|_2
 \le\sin(0.054/2)\|x\|_2
 <(0.054/2)R_0.
\]

At a KKT point,

\[
 \|y\|<Y_0:=\frac{(0.054/2)R_0}{6}.                            \tag{9.5}
\]

### 9.3 Two KKT points cannot coexist

If `u,v` are feasible high-perimeter KKT points and `d=u-v`, strong convexity
gives

\[
 d^T(\nabla f(u)-\nabla f(v))\ge m_0\|d\|^2.                   \tag{9.6}
\]

The scalar exponential remainder and `|a_j|<=2` give

\[
 \|g(u)-g(v)-Dg(v)d\|\le\|d\|^2,                               \tag{9.7}
\]

and the analogous reverse estimate.  Subtracting the two stationarity
equations and using feasibility yields

\[
 m_0\|d\|^2\le2Y_0\|d\|^2.                                   \tag{9.8}
\]

The exact verifier proves

\[
 m_0-2Y_0>2.22\cdot10^{-5}>0.                                 \tag{9.9}
\]

Thus `d=0`.  The winning orbit has at most one high-perimeter KKT point after
fixing rotation.

## 10. Completion

The original problem has a global maximizer by compactness.  Sections 2--3
show every global maximizer has a saturated 128-vertex difference body.
Sections 4--8 prove that its code belongs to the winning orbit.  Its gaps are
strictly inside (4.3), and (9.4) gives LICQ, so it is a KKT point.  Section 9
proves that this point is unique after normalization.  Dihedral code changes,
rotation, reflection and translation give congruent polygons.

Therefore, subject to independent checking of the geometric arguments and the
finite certificate, there is exactly one maximizing congruence class of small
64-gons.

## 11. Machine-assisted components

Analytic:

1. difference-body reconstruction;
2. at-most-one-interior-vertex perturbation;
3. saturation via normal-cone localization and KKT;
4. high-perimeter angle localization;
5. uniform regular-residual screen;
6. reflection-pair orthogonal decomposition;
7. nonbest-orbit lower-deficit inequalities;
8. best-orbit KKT uniqueness.

Finite certificate:

1. rigorous one-dimensional feasible root and deficit;
2. root-of-unity fixed-point enclosures;
3. exhaustive coverage of all `2^64` codes via two ternary MITM scans;
4. partition of 896 survivors into six dihedral orbits;
5. residual lower bounds and operator-norm upper bounds;
6. switch-root bound for continuous uniqueness.

## 12. External-validation requirements

This is not yet a literature-accepted solution.  Before public announcement it
should receive:

1. an independent implementation, preferably with Arb/FLINT or a proof
   assistant;
2. line-by-line review of Lemmas 1.1, 1.2 and the saturation KKT normalization;
3. a separate implementation of the orthogonal-pair exhaustive scan;
4. interval Newton only if certified decimal coordinates/perimeter are desired.
