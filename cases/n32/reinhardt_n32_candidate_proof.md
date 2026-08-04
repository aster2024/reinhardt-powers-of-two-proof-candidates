# A computer-assisted candidate proof for the maximum-perimeter small 32-gon

**Status:** complete candidate proof, 5 August 2026; not peer reviewed.

The proof below is self-contained apart from standard elementary facts about
convex polygons (Cauchy's perimeter formula, the perimeter maximum of an
`m`-gon in the unit disk, and the KKT theorem under MFCQ).  Every numerical
inequality and every finite-code assertion is checked by the accompanying
standard-library verifier.  No floating-point comparison is used in a proof
decision.

## 1. Statement

Let `P` range over convex 32-gons of diameter at most one.  Then the global
maximum of `per(P)` is attained by a unique congruence class.  Its fully
saturated difference body has sign-code orbit represented by

```text
+++-+-+-++--+--+--++-+-+-+++--++
```

or, in an axially symmetric representative,

```text
+-++--+-+-+---++--+++-+-+-++--+-.
```

The unique maximizing point is the unique high-perimeter KKT point of the
fixed-code problem for this orbit.  A non-rigorous high-precision Newton
evaluation gives

\[
 3.1403311569546193658254013805774586723120530983\ldots,
\]

in agreement with Mulansky--Potschka.  The proof itself does not depend on
this decimal expansion.

Throughout, put

\[
 n=32,\qquad \theta=\frac{\pi}{32},\qquad
 U=64\sin\frac{\pi}{64},\qquad
 \varepsilon=1.35\cdot10^{-13}.
\]

Here `U` is Reinhardt's perimeter upper bound for a small 32-gon.

## 2. Difference bodies and reconstruction

For a convex polygon `P`, let

\[
 Z=P-P.
\]

Then

\[
 Z\subseteq \overline B(0,1),\qquad
 \operatorname{per}(Z)=2\operatorname{per}(P),
\]

and `Z` is centrally symmetric with at most 64 edges.

### Lemma 2.1 (reconstruction)

Suppose `Z` is a nondegenerate centrally symmetric convex 64-gon, meaning that all 64 listed points are genuine vertices.  Label one half
of its vertices cyclically by

\[
 z_0,z_1,\ldots,z_{31},\qquad z_{32}=-z_0,
\]

and put `e_j=z_{j+1}-z_j`.  If a code `c_j in {+1,-1}` satisfies

\[
 \sum_{j=0}^{31}c_je_j=0,                                      \tag{2.1}
\]

then the vectors `c_j e_j`, sorted by direction, are the edge vectors of a
closed convex 32-gon `P_c`, and

\[
 P_c-P_c=Z.
\]

**Proof.**  The selected set contains exactly one vector from each antipodal
pair `{e_j,-e_j}`.  Equation (2.1) gives closure.  Sorting nonzero vectors of distinct directions gives a convex polygon with positive exterior angles.  The cyclic merge of its edge
list with the edge list of its negative is precisely the complete cyclic edge
list of `Z`.  Convex polygons with the same cyclic edge list differ only by a
translation; since both difference bodies are centered at the origin, they
coincide.  ∎

Consequently, sufficiently small perturbations of the `z_j` that preserve
(2.1), the disk constraints, genuine-vertex inequalities and cyclic order remain feasible for the
original small-polygon problem.

Summation by parts rewrites (2.1) as

\[
 -(c_0+c_{31})z_0+\sum_{j=1}^{31}(c_{j-1}-c_j)z_j=0.             \tag{2.2}
\]

All variable coefficients therefore belong to `{−2,0,2}`.

### Lemma 2.2 (at most one interior half-vertex)

At a local perimeter maximum whose difference body has 64 genuine vertices,
at most one of `z_0,...,z_31` lies strictly inside the unit disk.  If it exists,
its coefficient in (2.2) is nonzero.

**Proof.**  Write (2.2) as `sum a_j z_j=0`.  If two interior vertices `z_r,z_s`
have nonzero coefficients, use the two-sided perturbation

\[
 z_r(t)=z_r+t a_s h,\qquad z_s(t)=z_s-t a_r h.                  \tag{2.3}
\]

If an interior vertex has zero coefficient, move that vertex alone.  In all
cases the equality is preserved.  Since the moved vertices are interior and
cyclic order and the genuine-vertex inequalities are open conditions, both signs of
small `t` are feasible.  The perimeter is a sum of functions of the form
`||u+tv||`, hence is convex in `t`.  Choosing `h` outside the finite collection
of lines parallel to affected edges makes at least one summand strictly
convex.  A two-sided local maximum is impossible.  The same one-variable
argument shows that a unique interior vertex cannot have zero coefficient. ∎

This proof covers adjacent moved vertices and the wrap-around edge; possible
cancellation on their common edge does not affect the other incident edges.

## 3. A rigorously feasible near-regular polygon

Use the axially symmetric code

```text
c = +-++--+-+-+---++--+++-+-+-++--+-.
```

It satisfies `c_{31-j}=-c_j`.  Define rational perturbations `s_0=s_16=s_32=0`,
`s_{32-j}=-s_j`, and for `j=1,...,15`, set `s_j=m_j/10^20`, where

```text
m = (-84838116394748, -7759417296262, -90254740014369,
     -172750062732476, -104929943709607, -37109824686737,
     -104570029200993, -46261230637792, -102998890440414,
     -56520596873701, -100433725625308, -67650481919513,
     -34867238213719, -2083994507925, -1041997253963).
```

For

\[
 \phi_j(t)=j\theta+t s_j,
\]

symmetry makes the closure residual purely imaginary.  If
`a_j=c_{j-1}-c_j`, its imaginary component is

\[
 H(t)=2\sum_{j=1}^{15}a_j\sin(j\theta+t s_j)+a_{16}.             \tag{3.1}
\]

The exact interval verifier proves

\[
 H(0.99999953359)>5.53\cdot10^{-17},
\]

\[
 H(0.99999953360)<-7.84\cdot10^{-17}.                            \tag{3.2}
\]

Hence the intermediate value theorem gives a root `t_*` in this interval.
All angle gaps at every point of the bracket lie in `(0.09,0.107)`.  At the
root, Lemma 2.1 reconstructs a feasible small 32-gon, and exact interval
evaluation gives

\[
 U-\operatorname{per}(P(t_*))
 <1.336299789852041\cdot10^{-13}<\varepsilon.                   \tag{3.3}
\]

Thus every global maximizer has deficit at most `epsilon`.

The interval implementation obtains a rational enclosure of `pi` from
Machin's formula

\[
 \pi=16\arctan\frac15-4\arctan\frac1{239}
\]

and uses alternating Taylor bounds for sine.  Its certified `pi` interval has
width below `6.05*10^{-44}`.

## 4. Saturation for the global maximizer

This section is independent of any sign-code enumeration.

### 4.1 The difference body has 64 genuine vertices

If `Z` had fewer than 64 vertices, central symmetry would give at most 62.
The maximum perimeter of an `m`-gon in the unit disk is
`2m sin(pi/m)`, so

\[
 \operatorname{per}(P)\le62\sin\frac{\pi}{62}=U_{31}.
\]

For `U(t)=2t sin(pi/(2t))`,

\[
 U'(t)=2(\sin a-a\cos a)
      =2\int_0^a u\sin u\,du
      \ge \frac{\pi^2}{6t^3},\qquad a=\frac{\pi}{2t}.
\]

Therefore

\[
 U-U_{31}\ge\frac{\pi^2}{6\cdot32^3}
 >\frac{3}{2\cdot32^3}>\varepsilon.                             \tag{4.1}
\]

This contradicts (3.3).  Hence `Z` has exactly 64 genuine vertices.

### 4.2 Uniform localization of its normal cones

For a vertex `z_j`, write

- `r_j=||z_j||`;
- `omega_j` for its outward normal-cone width;
- `eta_j` for the cone midpoint direction;
- `phi_j` for the polar direction of `z_j`;
- `delta_j=eta_j-phi_j`.

Cauchy's formula is

\[
 \operatorname{per}(Z)
 =\sum_{j=0}^{63}2r_j\sin\frac{\omega_j}{2}\cos\delta_j,
 \qquad \sum_j\omega_j=2\pi.                                   \tag{4.2}
\]

Put `mu=pi/32`.  If one normal width is `t`, Jensen's inequality bounds the
right side by

\[
 F(t)=2\sin\frac t2
 +126\sin\frac{2\pi-t}{126}.
\]

Strong concavity of `2 sin(t/2)` shows that at either `t=mu/2` or
`t=3mu/2`,

\[
 2U-F(t)\ge\frac{\mu^2\sin(\mu/4)}{16}
 >\frac{9}{32^2\cdot64\cdot16}>2\varepsilon.                   \tag{4.3}
\]

Since `per(Z)>=2U-2epsilon`, monotonicity of `F` on each side of `mu` gives

\[
 \frac\mu2<\omega_j<\frac{3\mu}{2}.                             \tag{4.4}
\]

The nonnegative deficit decomposition from (4.2) now gives

\[
 1-r_j<64\varepsilon<\frac12,                                   \tag{4.5}
\]

and

\[
 1-\cos\delta_j<128\varepsilon.                                 \tag{4.6}
\]

If `|delta_j|>=mu/16=pi/512`, then

\[
 1-\cos\delta_j\ge\frac{2\delta_j^2}{\pi^2}
 \ge\frac1{131072}>128\varepsilon,
\]

contradicting (4.6).  Thus

\[
 |\delta_j|<\frac\mu{16}.                                       \tag{4.7}
\]

Consecutive normal midpoints differ by
`(omega_j+omega_{j+1})/2>mu/2`; hence consecutive radial directions differ by
more than

\[
 \frac\mu2-2\frac\mu{16}=\frac{3\mu}{8}.                        \tag{4.8}
\]

The same lower bound holds for the projective angular distance between any two
distinct half-vertices.

### 4.3 KKT excludes the last possible interior radius

By Lemma 2.2, suppose there is exactly one interior half-vertex `z_r`, with
`a_r=+-2`.  MFCQ holds explicitly: move every active circle vertex radially
inward and use the free displacement of `z_r` to repair the two components of
(2.2).  The equality Jacobian has rank two because the `z_r` block is
`a_r I_2`.

For the objective `per(P)=per(Z)/2`, the gradient with respect to a half-vertex
is

\[
 g_j=2\sin\frac{\omega_j}{2}e^{i\eta_j}.                         \tag{4.9}
\]

At the interior vertex, stationarity gives `g_r=a_r lambda`, so

\[
 ||\lambda||=\sin\frac{\omega_r}{2}
\]

and the projective direction of `lambda` is `eta_r`.  There is another index
`j != r` with `a_j !=0`; otherwise (2.2) would force `z_r=0`, impossible for a
vertex of the full-dimensional difference body.  This second vertex is on the
unit circle.  Projecting its KKT equation onto its tangent gives

\[
 |\sin(\eta_r-\phi_j)|
 =\frac{\sin(\omega_j/2)}{\sin(\omega_r/2)}|\sin\delta_j|.       \tag{4.10}
\]

By (4.4), the ratio is less than

\[
 \frac{\sin(3\mu/4)}{\sin(\mu/4)}<3.
\]

If `vartheta` is the projective distance from `eta_r` to `phi_j`, then
`sin(vartheta)>=2vartheta/pi` and (4.7)--(4.10) imply

\[
 \vartheta<\frac{3\pi\mu}{32}.
\]

Consequently

\[
 d_{\mathbb{RP}^1}(\phi_r,\phi_j)
 <\frac\mu{16}+\frac{3\pi\mu}{32}
 =\frac{(2+3\pi)\mu}{32}
 <\frac{3\mu}{8},                                                \tag{4.11}
\]

using `pi<10/3`.  This contradicts (4.8).  Hence every one of the 64 difference-
body vertices is on the unit circle.

## 5. The saturated angle/code model

Rotate so that

\[
 \phi_0=0,\qquad \phi_{32}=\pi,
\]

and put

\[
 \alpha_j=\phi_{j+1}-\phi_j>0.
\]

For a code `c in {+-1}^{32}`, the problem is

\[
 \max \sum_{j=0}^{31}2\sin\frac{\alpha_j}{2},                   \tag{5.1}
\]

subject to

\[
 G_c(\phi)=\sum_{j=0}^{31}c_j
 (e^{i\phi_{j+1}}-e^{i\phi_j})=0.                               \tag{5.2}
\]

Write

\[
 \phi_j=j\theta+s_j,\quad s_0=s_{32}=0,\quad
 x_j=s_{j+1}-s_j.
\]

Then `sum x_j=0`.

### Lemma 5.1 (all high-perimeter gaps are localized)

Every feasible point with deficit at most `epsilon` satisfies

\[
 0.09<\alpha_j<0.107.                                            \tag{5.3}
\]

**Proof.**  If one gap equals `t`, Jensen bounds the other 31 gaps by making
them equal.  The resulting one-variable upper envelope is increasing for
`t<theta` and decreasing for `t>theta`.  At `t=0.09`, strong concavity and
`pi>3` give a deficit larger than

\[
 \frac14\left(0.045-\frac{0.045^3}{6}\right)
 \frac{32}{31}\left(\frac3{32}-0.09\right)^2>\varepsilon.
\]

At `t=0.107`, use `theta<11/112` and the lower bound
`(3-0.107)/31` for the common remaining gap; the analogous rational
inequality is again larger than `epsilon`.  ∎

On this box, strong concavity yields

\[
 U-\operatorname{per}(P)
 \ge \kappa ||x||_2^2,\qquad
 \kappa=\frac14\sin(0.045)
 >\frac14\left(0.045-\frac{0.045^3}{6}\right).                  \tag{5.4}
\]

The verifier checks that the right-hand lower bound times
`(3.47*10^-6)^2` exceeds `epsilon`.  Therefore

\[
 ||x||_2<R_0:=3.47\cdot10^{-6}.                                 \tag{5.5}
\]

## 6. A uniform code-screening inequality

Let `xi=e^{i theta}`.  At the regular point,

\[
 b_c=(\xi-1)\sum_{j=0}^{31}c_j\xi^j.                             \tag{6.1}
\]

Taylor expansion gives

\[
 G_c(s)=b_c+L_c(s)+R_c(s),                                       \tag{6.2}
\]

where

\[
 L_c(s)=i\sum_{j=0}^{31}c_j
 (\xi^{j+1}s_{j+1}-\xi^j s_j).                                  \tag{6.3}
\]

The scalar exponential remainder and the Dirichlet Poincare inequality give

\[
 |R_c(s)|\le||s||_2^2
 \le\frac{||x||_2^2}{4\sin^2(\pi/64)}
 \le256||x||_2^2.                                                \tag{6.4}
\]

The following uniform linear estimate is useful beyond `n=32`.

### Lemma 6.1 (twisted-Dirichlet operator bound)

For general `n`,

\[
 |L_c(s)|\le\sqrt{2n}\cos\frac{\pi}{2n}\,||x||_2.               \tag{6.5}
\]

**Proof.**  Cauchy--Schwarz gives

\[
 |L_c(s)|^2\le n\sum_{j=0}^{n-1}|\xi s_{j+1}-s_j|^2.
\]

The twisted Dirichlet quadratic on `s_1,...,s_{n-1}` has eigenvalues

\[
 \mu_k=2-2\cos\theta\cos(k\theta),
\]

whereas the ordinary Dirichlet energy `sum x_j^2` has eigenvalues

\[
 \lambda_k=2-2\cos(k\theta).
\]

For `k>=1`,

\[
 \mu_k\le2\cos^2(\theta/2)\lambda_k,
\]

because half the difference between the right and left sides is
`cos(theta)-cos(k theta)>=0`.  This proves (6.5). ∎

For `n=32`, (6.5) is less than `8||x||`.  Thus every code capable of attaining
deficit at most `epsilon` must satisfy

\[
 |b_c|<8R_0+256R_0^2<3.41\cdot10^{-5}.                           \tag{6.6}
\]

## 7. Exact meet-in-the-middle code certificate

Overall negation of a code is immaterial, so normalize `c_0=+1`.  There are
`2^31=2,147,483,648` normalized half-codes.

The verifier encloses every coordinate of every edge

\[
 e_j=\xi^{j+1}-\xi^j
\]

by a rational interval and chooses an integer vector `E_j` at scale
`M=10^18` such that each coordinate differs from `E_j/M` by less than `1/M`.
It splits the code into blocks `0,...,15` and `16,...,31`, forms all `2^15` and
`2^16` signed integer sums, sorts the second block by its first coordinate,
and performs an exact integer range search.  Since the total Euclidean
rounding error is less than `64/M`, every true residual satisfying (6.6) is
included in the integer search circle.

The exact result is:

\[
 \boxed{96\text{ normalized codes survive}.}                    \tag{7.1}
\]

They are exactly three disjoint dihedral orbits, each of size 32, represented
by

```text
A = ++++--+--+-+-+-+++---+-+-+-++-++
B = +++-+-+-++--+--+--++-+-+-+++--++
C = ++-++-+-+-+-+--++--+-+-+-+-++-++.
```

The axially symmetric code of Section 3 belongs to orbit `B`.

This certificate scans every normalized code, not only axially symmetric
ones.  It therefore does not assume Mossinghoff's symmetry conjecture.

## 8. Exact exclusion of orbits A and C

For a fixed code let

\[
 \sigma_c=\sup_{s\ne0}\frac{|L_c(s)|}{||x||_2}.
\]

Writing `a_j=c_{j-1}-c_j`, the squared norm is the largest eigenvalue of the
2-by-2 matrix

\[
 K_c=B_cD^{-1}B_c^T,                                             \tag{8.1}
\]

where `D` is the Dirichlet path matrix and the rows of `B_c` are

\[
 (-a_j\sin(j\theta))_{j=1}^{31},\qquad
 ( a_j\cos(j\theta))_{j=1}^{31}.                                \tag{8.2}
\]

The verifier uses

\[
 (D^{-1})_{jk}=\frac{\min(j,k)(32-\max(j,k))}{32}                \tag{8.3}
\]

and exact fixed-point error bounds.  Gershgorin gives the certified squared-
norm upper bounds

\[
 \sigma_A^2<16.36425,\qquad
 \sigma_B^2<16.45308,\qquad
 \sigma_C^2<16.25477,                                           \tag{8.4}
\]

so in particular every survivor has `sigma_c<4.1`.

The quantities `|b_c|` and `sigma_c` are invariant under the full dihedral action and overall sign change: rotations multiply the complex residual and linear map by a unit complex number, while reflections conjugate them and reverse the Dirichlet path.  It is therefore enough to check one representative of each orbit.

The same exact residual intervals prove

\[
 |b_A|>2.31\cdot10^{-5},\qquad
 |b_C|>1.66\cdot10^{-5}.                                        \tag{8.5}
\]

For any high-perimeter feasible point of either code, (6.2), (6.4), (5.5),
and (8.4) imply

\[
 ||x||_2>
 \frac{|b_c|-256R_0^2}{4.1}.                                    \tag{8.6}
\]

Combining (8.6) with (5.4), the exact rational lower deficits are

\[
 U-\operatorname{per}(P)>3.5689976\cdot10^{-13}
 \quad\text{for orbit A},                                       \tag{8.7}
\]

and

\[
 U-\operatorname{per}(P)>1.8428631\cdot10^{-13}
 \quad\text{for orbit C}.                                       \tag{8.8}
\]

Both exceed `epsilon`, so neither orbit can contain a global maximizer.  Orbit
`B` is the unique surviving code class.

## 9. Continuous global uniqueness inside orbit B

Use the axial representative of Section 3.  It has 21 switch indices

\[
 J=\{j: c_{j-1}\ne c_j\}.
\]

Let `f=-per(P)` and let `g=(Re G_c,Im G_c)`.

### 9.1 Strong convexity

The Hessian of `f` is the weighted Dirichlet path matrix with edge weights

\[
 w_j=\frac12\sin\frac{\alpha_j}{2}.
\]

Using (5.3), `sin(pi/64)>=1/32`, and the alternating lower bound for
`sin(0.045)`,

\[
 \nabla^2 f\succeq m_0 I,
 \qquad
 m_0:=\frac{0.045-0.045^3/6}{512}
 >8.78609\cdot10^{-5}.                                          \tag{9.1}
\]

### 9.2 Constraint qualification and multiplier bound

At a point `phi`, the nonzero columns of `Dg` are `a_j i e^{i phi_j}`.  Hence

\[
 \sigma_{\min}(Dg)^2
 =2\left(21-\left|\sum_{j\in J}e^{2i\phi_j}\right|\right).       \tag{9.2}
\]

The exact fixed-point certificate shows that at the regular point

\[
 \left|\sum_{j\in J}e^{2ij\theta}\right|<2.                     \tag{9.3}
\]

Furthermore, (6.4) gives `||s||<16R_0`.  Since `sqrt(21)<5`,

\[
 \left|\sum_{j\in J}
 (e^{2i\phi_j}-e^{2ij\theta})\right|
 \le2\sqrt{21}\,||s||<160R_0<1.                                 \tag{9.4}
\]

Equations (9.2)--(9.4) imply

\[
 \sigma_{\min}(Dg)>6.                                           \tag{9.5}
\]

Thus LICQ holds throughout the high-perimeter ball.

At the regular point the objective gradient vanishes.  If
`r_j=cos(alpha_j/2)-cos(theta/2)`, then the path difference operator and
(5.3) give

\[
 ||\nabla f||_2\le\sin(0.107/2)||x||_2
 <\frac{0.107}{2}R_0.                                           \tag{9.6}
\]

At a KKT point, `nabla f+(Dg)^T y=0`, so (9.5)--(9.6) give

\[
 ||y||<Y_0:=\frac{(0.107/2)R_0}{6}
 <3.095\cdot10^{-8}.                                             \tag{9.7}
\]

### 9.3 Two KKT points are impossible

Suppose `u,v` are two high-perimeter feasible KKT points and put `d=u-v`.
For the complex closure constraint, `|a_j|<=2` and the scalar exponential
remainder imply

\[
 ||g(u)-g(v)-Dg(v)d||_2\le||d||_2^2,                             \tag{9.8}
\]

and the analogous bound holds with `u` and `v` interchanged.
Subtract the two stationarity equations and take the inner product with `d`.
Strong convexity, feasibility, (9.8), and (9.7) give

\[
 m_0||d||_2^2\le2Y_0||d||_2^2.                                  \tag{9.9}
\]

But the exact verifier proves

\[
 m_0-2Y_0>8.7799\cdot10^{-5}>0.                                 \tag{9.10}
\]

Therefore `d=0`.  Orbit `B` has at most one high-perimeter KKT point after
fixing rotation.

## 10. Completion of the global argument

The original polygon problem has a global maximizer by compactness.  Sections
3--4 show that every global maximizer has a saturated 64-vertex difference
body.  Sections 5--8 show that its code must lie in orbit `B`.  Its angle gaps
are strictly between `0.09` and `0.107`, so no order inequality is active, and
Section 9 gives LICQ; hence it is a KKT point.  Section 9 proves that such a
point is unique for a normalized representative.  Dihedral changes of code,
rotation, reflection and translation produce congruent polygons.

Thus there is exactly one maximizing congruence class, and its code is orbit
`B`.  This proves the theorem, subject only to independent checking of the
written geometric lemmas and execution of the finite certificate.

## 11. What is and is not machine-assisted

Analytic:

1. difference-body reduction and reconstruction;
2. at-most-one-interior-vertex perturbation;
3. MFCQ/KKT saturation;
4. gap localization and strong-concavity radius bound;
5. uniform twisted-Dirichlet operator estimate;
6. residual elimination inequalities;
7. best-code KKT uniqueness.

Finite certificate:

1. the one-dimensional feasible root and its perimeter bound;
2. interval enclosures of roots of unity;
3. exact coverage of all `2^31` normalized codes by meet-in-the-middle;
4. the partition of 96 survivors into three dihedral orbits;
5. code-specific operator-norm bounds and residual lower bounds;
6. the regular switch-root bound used for LICQ.

The verifier uses only the Python standard library (`Fraction`, integers,
`bisect`, and elementary data structures).  Its proof decisions do not use
binary floating point.  Decimal-looking values in the printed report are only
human-readable renderings after all assertions have passed.

## 12. Remaining external-validation requirements

This is a candidate proof of an open problem, not a literature-accepted
result.  Before public announcement it should be:

1. independently reimplemented, preferably in Arb/FLINT or Lean/Isabelle;
2. checked line by line by a convex-geometry expert, especially Lemmas 2.1,
   2.2 and the half-variable KKT normalization in Section 4;
3. compared against a separate exact code-orbit implementation;
4. supplemented by an interval-Newton enclosure of the unique KKT point if a
   certified decimal value is desired.

## References

- C. Bingane, *Maximal perimeter and maximal width of a convex small polygon*,
  2021.
- B. Mulansky and A. Potschka, *A zonogon approach for computing small polygons
  of maximum perimeter*, Mathematical Programming, 2025.
