# A computer-assisted proof candidate for Reinhardt's maximum-perimeter small hexadecagon

**Status (4 August 2026).** This document gives a complete internal proof chain for the case \(n=16\), conditional only on the successful execution of the accompanying exact-arithmetic verifier. The finite computation uses Python's standard library, exact rational arithmetic, and outward integer intervals; no floating-point comparison is used in a proof decision. The argument has not yet been independently refereed or reproduced in a second proof system, so it should be treated as a *computer-assisted proof candidate*, not yet as an established published theorem.

Verifier:

- `reinhardt_n16_verified_audit.py`
- SHA-256: `f179e7176877886e4fea606b6d77dc579c4e9463541d54e6cfdcced2e4fa35b3`

Recorded output:

- `reinhardt_n16_verified_audit_output.txt`
- SHA-256: `e716b3e891b104bddf5f8ea92e46a5c48156bdac2a9382a47214acaeebf90061`

The verifier certifies every numerical inequality explicitly invoked below and scans all \(2^{15}=32768\) normalized sign codes.

---

## 1. Statement

A **small polygon** is a planar convex polygon of diameter at most one. Write \(p(P)\) for its perimeter.

### Main theorem
Among all convex small hexadecagons, the maximum perimeter is attained by a unique congruence class (allowing translation, rotation, reflection, and cyclic relabeling). In the saturated difference-body representation its sign code is the dihedral class of

\[
 c_*=(+--+-++-)^2.
\]

Its perimeter is numerically

\[
 3.136547716486607386085967\ldots .
\]

The decimal value is included only to identify the optimizer; the global proof below does not depend on a decimal interval for the final stationary point.

---

## 2. Difference bodies and the first strict reduction

Let \(P\) be a convex small hexadecagon and put

\[
 Z=P-P=P+(-P).
\]

Then \(Z\) is centrally symmetric, \(Z\subseteq \overline B(0,1)\), and Minkowski additivity of planar perimeter gives

\[
 p(Z)=2p(P). \tag{2.1}
\]

The edge directions of a Minkowski sum are obtained by merging the edge-direction lists of its summands. Consequently, \(Z\) has at most \(32\) genuine edges.

We use the following certified lower bound, supplied by Bingane's explicit feasible construction \(C_{16}\):

\[
 p(C_{16})>L_0:=3.1365475. \tag{2.2}
\]

The exact verifier evaluates the published radical/trigonometric construction by nested radical intervals and proves (2.2).

### Lemma 2.1 — the extremal difference body has 32 strict vertices
Let \(P\) maximize perimeter. Then \(Z=P-P\) has exactly \(32\) genuine vertices and edges.

#### Proof
If a centrally symmetric polygon has fewer than \(32\) vertices, it has at most \(30\). An \(m\)-gon contained in the unit disk has perimeter at most

\[
 2m\sin\frac{\pi}{m},
\]

with equality only for the regular inscribed \(m\)-gon. Since \(m\sin(\pi/m)\) increases with \(m\), (2.1) yields

\[
 p(P)\le 30\sin\frac{\pi}{30}<3.1365475,
\]

where the last inequality is exactly certified. This contradicts (2.2). ∎

In particular, an extremal \(P\) has sixteen nonzero edges with pairwise distinct unoriented directions; otherwise the merged edge list of \(P+(-P)\) would have fewer than \(32\) members.

---

## 3. Exact reconstruction from a sign code

Let \(Z\) be a strictly convex centrally symmetric \(32\)-gon. Label one half of its vertices in counterclockwise order by

\[
 z_0,z_1,\ldots,z_{15},\qquad z_{16}=-z_0,
\]

and define its half-edge vectors

\[
 e_j=z_{j+1}-z_j,\qquad 0\le j\le15.
\]

The other sixteen edge vectors are \(-e_0,\ldots,-e_{15}\).

### Lemma 3.1 — reconstruction lemma
For \(c=(c_0,\ldots,c_{15})\in\{\pm1\}^{16}\), the condition

\[
 \sum_{j=0}^{15}c_je_j=0 \tag{3.1}
\]

is necessary and sufficient for the existence of a strictly convex hexadecagon \(P_c\) whose difference body is \(Z\) and whose edge set is

\[
 \{c_0e_0,\ldots,c_{15}e_{15}\}.
\]

The polygon \(P_c\) is unique up to translation.

#### Proof
Necessity follows from closure of the oriented edge list of \(P_c\).

Conversely, put \(f_j=c_je_j\). Because \(Z\) has \(32\) strict edges, the \(f_j\) are nonzero and have pairwise distinct directions. Sort them by polar angle. Equation (3.1) says that the sorted list closes. A cyclic list of nonzero vectors with strictly increasing directions and zero sum is the edge list of a strictly convex polygon, unique up to translation.

The edge list of \(-P_c\) is \(-f_0,\ldots,-f_{15}\). Hence the edge-direction merge for \(P_c+(-P_c)\) is exactly

\[
 \{f_j,-f_j:0\le j\le15\}
 =\{e_j,-e_j:0\le j\le15\},
\]

which is the full cyclic edge list of \(Z\). Two convex polygons with the same cyclic edge list differ by a translation. Both difference bodies are centrally symmetric about the origin, so the translation is zero and \(P_c-P_c=Z\). ∎

Expanding (3.1) by summation by parts gives the vertex form

\[
 \sum_{j=0}^{15}a_jz_j=0, \tag{3.2}
\]

where

\[
 a_0=-(c_0+c_{15}),\qquad
 a_j=c_{j-1}-c_j\quad(1\le j\le15). \tag{3.3}
\]

Thus every \(a_j\in\{-2,0,2\}\).

Lemma 3.1 is also the needed local-feasibility statement: any sufficiently small perturbation of the \(z_j\) that preserves central symmetry, strict cyclic order, disk containment, and (3.2) reconstructs a feasible small hexadecagon with the same code.

---

## 4. At most one interior vertex in a half difference body

For the half-edge representation, (2.1) gives

\[
 p(P)=\sum_{j=0}^{15}\lVert z_{j+1}-z_j\rVert. \tag{4.1}
\]

### Lemma 4.1 — two interior half-vertices are impossible
At a local maximizer satisfying Lemma 2.1, at most one of \(z_0,\ldots,z_{15}\) lies strictly inside the unit disk. If \(z_r\) is the unique interior half-vertex, then \(a_r\ne0\).

#### Proof
Suppose first that \(z_r,z_s\) are two distinct interior half-vertices.

If \(a_r=0\), set

\[
 z_r(t)=z_r+th
\]

and leave the other half-vertices fixed. If both coefficients are nonzero, set

\[
 z_r(t)=z_r+t a_sh,\qquad
 z_s(t)=z_s-t a_rh. \tag{4.2}
\]

If exactly one coefficient vanishes, vary the corresponding zero-coefficient vertex alone. In every case (3.2) is preserved.

Because the moved vertices are strictly inside the disk, both signs of sufficiently small \(t\) preserve disk containment. Strict convexity, nonzero edges, and cyclic vertex order are open conditions, so they too persist for \(|t|\) small. Lemma 3.1 therefore reconstructs a feasible polygon along the whole small two-sided interval.

Along this line, (4.1) is a sum of functions of the form

\[
 t\longmapsto\lVert u+tv\rVert,
\]

hence is convex. At least one edge changes. Choose \(h\) outside the finite set of directions parallel to the affected edges. For that edge, \(u\) and \(v\) are not parallel, so its norm is strictly convex. Thus the whole perimeter is strictly convex in \(t\). Strict convexity gives

\[
 p(P_0)<\max\{p(P_{-t}),p(P_t)\}
\]

for every sufficiently small \(t>0\), contradicting local maximality.

If a unique interior vertex \(z_r\) had \(a_r=0\), the same one-vertex perturbation would give the same contradiction. ∎

This proof includes the adjacent-vertex and endpoint cases: even if the variation of their common edge cancels, at least one outer incident edge has nonzero variation.

---

## 5. Quantitative near-regularity of every competitor

Let the full cyclic vertex list of \(Z\) be indexed modulo \(32\). For vertex \(z_j\), write

\[
 r_j=\lVert z_j\rVert,
\]

let \(\omega_j\in(0,\pi)\) be the width of its outward normal cone, let \(\eta_j\) be the midpoint direction of that cone, let \(\phi_j\) be the polar direction of \(z_j\), and choose

\[
 \delta_j=\eta_j-\phi_j\in[-\pi,\pi].
\]

Cauchy's perimeter formula, integrated separately over the normal cones, gives the exact identity

\[
 p(Z)=\sum_{j=0}^{31}2r_j\sin\frac{\omega_j}{2}\cos\delta_j,
 \qquad
 \sum_{j=0}^{31}\omega_j=2\pi. \tag{5.1}
\]

The universal inscribed-polygon bound from (5.1) and Jensen is

\[
 p(Z)\le U_Z:=64\sin\frac{\pi}{32}. \tag{5.2}
\]

For an optimizer, (2.2) and (2.1) imply

\[
 0\le U_Z-p(Z)<1.982\times10^{-6}. \tag{5.3}
\]

All numerical inequalities in the next lemma are exact verifier assertions.

### Lemma 5.1 — normal-cone, radial, and angular localization
Every vertex of an extremal \(Z\) satisfies

\[
 0.18<\omega_j<0.21, \tag{5.4}
\]

\[
 r_j>0.99998, \tag{5.5}
\]

\[
 |\delta_j|<0.005. \tag{5.6}
\]

Moreover, the projective angular distance between any two distinct half-vertices is greater than \(0.17\).

#### Proof
Fix \(\omega_j=t\). Dropping the factors \(r_k\cos\delta_k\le1\) and applying Jensen to the other 31 cone widths gives

\[
 p(Z)\le F_\omega(t)
 :=2\sin\frac t2+62\sin\frac{2\pi-t}{62}. \tag{5.7}
\]

The function increases up to \(2\pi/32\) and decreases afterwards. The verifier proves

\[
 F_\omega(0.18)<2L_0,
 \qquad
 F_\omega(0.21)<2L_0.
\]

Since \(p(Z)>2L_0\), (5.4) follows.

Rewrite the deficit as

\[
\begin{aligned}
 U_Z-p(Z)
 &=\left(U_Z-\sum_j2\sin\frac{\omega_j}{2}\right)\\
 &\quad+\sum_j2\sin\frac{\omega_j}{2}
       \left(1-r_j\cos\delta_j\right).
\end{aligned} \tag{5.8}
\]

Every term on the right is nonnegative. Hence

\[
 U_Z-p(Z)
 \ge2(1-r_j)\sin\frac{\omega_j}{2}.
\]

Equations (5.3), (5.4), and the certified inequality

\[
 2(1-0.99998)\sin(0.09)>1.982\times10^{-6}
\]

give (5.5).

Likewise, from

\[
 1-r_j\cos\delta_j
 =(1-r_j)+r_j(1-\cos\delta_j)
 \ge r_j(1-\cos\delta_j),
\]

we obtain

\[
 U_Z-p(Z)
 \ge2r_j\sin\frac{\omega_j}{2}(1-\cos\delta_j).
\]

The verifier checks that the right side at
\(r_j=0.99998\), \(\omega_j=0.18\), and \(|\delta_j|=0.005\)
exceeds the deficit bound (5.3). This proves (5.6).

If the consecutive normal-cone midpoint directions are lifted cyclically, then

\[
 \eta_{j+1}-\eta_j
 =\frac{\omega_j+\omega_{j+1}}2>0.18. \tag{5.9}
\]

Using (5.6), every consecutive radial angle gap is therefore greater than

\[
 0.18-2(0.005)=0.17. \tag{5.10}
\]

For two distinct half-vertices, their projective distance is a sum of one or more consecutive full-circle gaps, or the complementary sum to \(\pi\); either sum contains at least one gap from (5.10). ∎

---

## 6. Exclusion of the final unsaturated radius

### Lemma 6.1 — every difference-body vertex is on the unit circle
For an extremal hexadecagon,

\[
 \lVert z_j\rVert=1\qquad(0\le j\le31). \tag{6.1}
\]

#### Proof
By Lemma 4.1, it is enough to exclude a unique interior half-vertex \(z_r\). In this case \(a_r\ne0\), hence \(|a_r|=2\).

Consider the finite-dimensional program in variables \(z_0,\ldots,z_{15}\): maximize (4.1), subject to the two equality components of (3.2) and the disk inequalities

\[
 q_j(z)=\lVert z_j\rVert^2-1\le0.
\]

All strict convexity and ordering requirements hold on an open neighborhood and introduce no active constraints.

We first verify MFCQ. For every active vertex \(j\ne r\), choose \(d_j=-z_j\). Then

\[
 Dq_j(z)d_j=-2<0.
\]

Because \(a_r\ne0\), choose \(d_r\) uniquely so that

\[
 \sum_ja_jd_j=0.
\]

The equality Jacobian has rank two, again because its \(z_r\)-block is \(a_rI_2\). Thus MFCQ holds and the KKT equations are necessary.

Let \(g_j\) denote the gradient of (4.1) with respect to \(z_j\). The difference of the two unit incident edge tangents is the outward angle-bisector vector, so

\[
 g_j=2\sin\frac{\omega_j}{2}\,e^{i\eta_j}. \tag{6.2}
\]

At the interior vertex the disk multiplier vanishes. Therefore, for some equality multiplier \(\lambda\in\mathbb R^2\),

\[
 g_r=a_r\lambda. \tag{6.3}
\]

Consequently the projective direction of \(\lambda\) is \(\eta_r\), and

\[
 \lVert\lambda\rVert=\sin\frac{\omega_r}{2}. \tag{6.4}
\]

There is another index \(j\ne r\) with \(a_j\ne0\). Otherwise (3.2) would imply \(a_rz_r=0\), forcing \(z_r=0\); but the origin lies in the interior of the two-dimensional centrally symmetric polygon \(Z\), so it is not a vertex.

At this second index \(z_j\) lies on the unit circle. Project its KKT equation onto the tangent direction \(ie^{i\phi_j}\); the radial disk multiplier disappears. Using \(|a_j|=2\), (6.2), and (6.4), we get

\[
 \left|\sin(\eta_r-\phi_j)\right|
 =\frac{\sin(\omega_j/2)}{\sin(\omega_r/2)}
   |\sin\delta_j|. \tag{6.5}
\]

By Lemma 5.1 and exact interval evaluation,

\[
 \frac{\sin(0.105)}{\sin(0.09)}\sin(0.005)<0.00585. \tag{6.6}
\]

If \(d_{\mathbb{RP}^1}\) denotes projective angular distance, then for
\(0\le d\le\pi/2\), \(\sin d\ge2d/\pi\). Thus (6.5)--(6.6) imply

\[
 d_{\mathbb{RP}^1}(\eta_r,\phi_j)
 <\frac\pi2(0.00585).
\]

Together with \(|\eta_r-\phi_r|<0.005\), the verifier certifies

\[
 d_{\mathbb{RP}^1}(\phi_r,\phi_j)<0.0142<0.17,
\]

contradicting Lemma 5.1. Hence no interior vertex exists. ∎

This proves, for \(n=16\), the full-radius saturation property that the prior zonogon computation had assumed conjecturally.

---

## 7. The saturated angle model

Rotate so that

\[
 z_0=1,\qquad z_{16}=-1,
\]

and write

\[
 z_j=e^{i\phi_j},\qquad
 0=\phi_0<\phi_1<\cdots<\phi_{16}=\pi.
\]

Set

\[
 \alpha_j=\phi_{j+1}-\phi_j,
 \qquad
 \sum_{j=0}^{15}\alpha_j=\pi. \tag{7.1}
\]

Then

\[
 p(P)=\sum_{j=0}^{15}2\sin\frac{\alpha_j}{2}, \tag{7.2}
\]

and the code closure is

\[
 G_c(\phi):=
 \sum_{j=0}^{15}c_j
 \left(e^{i\phi_{j+1}}-e^{i\phi_j}\right)=0. \tag{7.3}
\]

Let

\[
 \alpha_0^*=\frac\pi{16},\qquad
 x_j=\alpha_j-\alpha_0^*,
 \qquad \sum_jx_j=0. \tag{7.4}
\]

### Lemma 7.1 — every global competitor is in a tiny regular neighborhood
Every saturated feasible configuration with perimeter greater than \(L_0\) satisfies

\[
 0.189<\alpha_j<0.204 \tag{7.5}
\]

and

\[
 \lVert x\rVert_2<0.0065. \tag{7.6}
\]

#### Proof
If one gap is fixed at \(t\), Jensen on the remaining fifteen gaps gives

\[
 p(P)\le F_\alpha(t)
 :=2\sin\frac t2+30\sin\frac{\pi-t}{30}. \tag{7.7}
\]

This function increases up to \(\pi/16\) and decreases afterwards. The verifier proves

\[
 F_\alpha(0.189)<L_0,
 \qquad
 F_\alpha(0.204)<L_0,
\]

which gives (7.5).

On this interval, for \(h(t)=2\sin(t/2)\),

\[
 h''(t)=-\frac12\sin\frac t2
 \le-\frac12\sin(0.189/2).
\]

Using \(\sum x_j=0\), strong concavity yields

\[
 32\sin\frac\pi{32}-p(P)
 \ge\frac14\sin(0.189/2)\lVert x\rVert_2^2. \tag{7.8}
\]

The verifier proves

\[
 32\sin\frac\pi{32}-L_0<9.91\times10^{-7}.
\]

and, from (7.8), the strict bound (7.6). ∎

---

## 8. Exact exclusion of all but one code class

Define cumulative angle perturbations

\[
 s_0=s_{16}=0,\qquad
 s_j=\sum_{k=0}^{j-1}x_k\quad(1\le j\le15), \tag{8.1}
\]

so that

\[
 \phi_j=\frac{j\pi}{16}+s_j.
\]

The discrete Dirichlet inequality is

\[
 \sum_{j=1}^{15}s_j^2
 \le\frac{1}{4\sin^2(\pi/32)}\sum_{j=0}^{15}x_j^2
 <26.1\lVert x\rVert_2^2. \tag{8.2}
\]

Let \(\xi=e^{i\pi/16}\). Taylor expansion of (7.3) at the regular point gives

\[
 G_c(x)=b_c+L_c(s)+R_c(s), \tag{8.3}
\]

where

\[
 b_c=\sum_{j=0}^{15}c_j(\xi^{j+1}-\xi^j), \tag{8.4}
\]

\[
 L_c(s)=i\sum_{j=1}^{15}(c_{j-1}-c_j)\xi^js_j, \tag{8.5}
\]

and, since \(|e^{it}-1-it|\le t^2/2\),

\[
 |R_c(s)|\le\sum_{j=1}^{15}s_j^2
 <26.1\lVert x\rVert_2^2. \tag{8.6}
\]

Let \(D\) be the \(15\times15\) Dirichlet path Laplacian, so that
\(\lVert x\rVert_2^2=s^TDs\). Write \(B_c\) for the real \(2\times15\) matrix representing (8.5), and define

\[
 \sigma_c^2=\lambda_{\max}(B_cD^{-1}B_c^T). \tag{8.7}
\]

Then \(|L_c(s)|\le\sigma_c\lVert x\rVert_2\). Consequently every code admitting a competitor with perimeter greater than \(L_0\) must satisfy

\[
 |b_c|
 \le0.0065\,\sigma_c+26.1(0.0065)^2. \tag{8.8}
\]

### Lemma 8.1 — finite exact code certificate
After normalizing \(c_0=+1\), exactly sixteen of the \(2^{15}=32768\) codes satisfy the necessary condition (8.8). They form precisely one dihedral orbit, represented by

\[
 c_*=(+--+-++-)^2. \tag{8.9}
\]

Every other code violates (8.8). The smallest certified violation margin is greater than

\[
 0.002111068044745996. \tag{8.10}
\]

#### Proof
This is the finite computer-assisted part. The accompanying verifier:

1. encloses \(\pi\) by exact rational Machin-series bounds;
2. constructs \(\sin(\pi/32)\), \(\cos(k\pi/16)\), and \(\sin(k\pi/16)\) using nested square-root integer intervals at scale \(10^{75}\);
3. uses the exact formula
   \[
   (D^{-1})_{jk}=\frac{\min(j,k)(16-\max(j,k))}{16};
   \]
4. computes a lower interval for \(|b_c|\) and an upper interval for \(\sigma_c\);
5. makes every keep/exclude decision using integers only;
6. scans all \(32768\) normalized codes and checks equality of the survivor set with the explicitly generated dihedral orbit of (8.9).

A fresh rerun reproduces the recorded output byte-for-byte and ends with `ALL ASSERTIONS PASSED`. ∎

---

## 9. Global uniqueness inside the surviving code

Fix the representative (8.9). Use variables

\[
 u=(\phi_1,\ldots,\phi_{15})
\]

with \(\phi_0=0\), \(\phi_{16}=\pi\), and minimize the negative perimeter

\[
 f(u)=-\sum_{j=0}^{15}2\sin\frac{\phi_{j+1}-\phi_j}{2}. \tag{9.1}
\]

The two real closure constraints are denoted by \(g(u)=0\). The nonzero coefficients
\(a_j=c_{j-1}-c_j\) occur at

\[
 S=\{1,3,4,5,7,8,9,11,12,13,15\},
 \qquad |S|=11. \tag{9.2}
\]

### Lemma 9.1 — uniform strong convexity
Throughout the high-perimeter box (7.5),

\[
 \nabla^2f(u)\succeq mI,
 \qquad m>0.0018. \tag{9.3}
\]

#### Proof
The Hessian is the weighted Dirichlet path Laplacian whose edge weights are

\[
 \frac12\sin\frac{\alpha_j}{2}.
\]

The smallest eigenvalue of the unweighted Dirichlet path Laplacian is
\(4\sin^2(\pi/32)\). Hence

\[
 m\ge2\sin(0.189/2)\sin^2(\pi/32)>0.0018,
\]

with the last inequality certified exactly. ∎

### Lemma 9.2 — uniform constraint qualification and multiplier bound
Throughout the same box,

\[
 \sigma_{\min}(Dg(u))>4.42, \tag{9.4}
\]

and every KKT multiplier \(y\in\mathbb R^2\) satisfies

\[
 \lVert y\rVert_2<0.000151. \tag{9.5}
\]

#### Proof
The nonzero columns of \(J=Dg(u)\) are

\[
 a_j(-\sin\phi_j,\cos\phi_j)^T,
 \qquad j\in S,
\]

with \(|a_j|=2\). Therefore

\[
 \lambda_{\min}(JJ^T)
 =2\left(11-\left|\sum_{j\in S}e^{2i\phi_j}\right|\right). \tag{9.6}
\]

At the regular point,

\[
 \sum_{j\in S}e^{2ij\pi/16}=-1. \tag{9.7}
\]

This follows directly by subtracting the complementary switch set from the full set of nontrivial sixteenth roots.

With \(\phi_j=j\pi/16+s_j\),

\[
 \left|\sum_{j\in S}e^{2i\phi_j}\right|
 \le1+2\sum_{j\in S}|s_j|
 \le1+2\sqrt{11}\lVert s\rVert_2. \tag{9.8}
\]

By (8.2) and (7.6),

\[
 \lVert s\rVert_2
 <\frac{0.0065}{2\sin(\pi/32)}<0.0332. \tag{9.9}
\]

The exact verifier combines (9.6)--(9.9) to prove (9.4).

For the gradient, let

\[
 q_j=\cos(\alpha_j/2)-\cos(\pi/32).
\]

The internal gradient is a first difference of the \(q_j\), whose operator norm is at most two. The mean-value theorem and (7.5) give

\[
 \lVert\nabla f(u)\rVert_2
 <\sin(0.102)(0.0065)<0.000663. \tag{9.10}
\]

At a KKT point, \(\nabla f+J^Ty=0\), so

\[
 \lVert y\rVert_2
 \le\frac{\lVert\nabla f\rVert_2}{\sigma_{\min}(J)}
 <\frac{0.000663}{4.42}<0.000151.
\]
∎

### Lemma 9.3 — at most one high-perimeter KKT point
The surviving code has at most one KKT point satisfying (7.5)--(7.6).

#### Proof
Suppose \((u,y_u)\) and \((v,y_v)\) are two feasible KKT pairs, and put \(d=u-v\). The line segment between \(u\) and \(v\) remains in the angle box, so Lemma 9.1 gives

\[
 d^T(\nabla f(u)-\nabla f(v))\ge m\lVert d\rVert_2^2. \tag{9.11}
\]

For the closure map, its second directional derivative has norm at most

\[
 2\sum_{j=1}^{15}d_j^2=2\lVert d\rVert_2^2,
\]

because \(|a_j|\le2\). Taylor's theorem, including its factor \(1/2\), therefore gives

\[
 \lVert Dg(u)d\rVert_2\le\lVert d\rVert_2^2,
 \qquad
 \lVert Dg(v)d\rVert_2\le\lVert d\rVert_2^2, \tag{9.12}
\]

where feasibility \(g(u)=g(v)=0\) was used.

Subtract the two stationarity equations and take the inner product with \(d\). Using (9.12),

\[
\begin{aligned}
 m\lVert d\rVert_2^2
 &\le d^T(\nabla f(u)-\nabla f(v))\\
 &\le(\lVert y_u\rVert_2+\lVert y_v\rVert_2)
      \lVert d\rVert_2^2\\
 &<0.000302\lVert d\rVert_2^2.
\end{aligned} \tag{9.13}
\]

But \(m>0.0018>0.000302\). Hence \(d=0\). ∎

---

## 10. Proof of the main theorem

The feasible set of convex hexadecagons of diameter at most one, modulo translation and allowing repeated limiting vertices, is compact, and perimeter is continuous. Hence a global maximizer exists. Bingane's construction makes its perimeter greater than \(L_0\).

Lemma 2.1 gives a strict \(32\)-vertex difference body. Lemma 4.1 and Lemma 6.1 show that every difference-body vertex is on the unit circle. Therefore the maximizer is represented by the saturated angle-code model of Section 7.

Lemma 7.1 puts it in the certified regular neighborhood. Lemma 8.1 forces its code into the single dihedral class of \(c_*\). The angle inequalities are strict, and Lemma 9.2 gives full row rank of the closure Jacobian, so the maximizer is a KKT point of the fixed-code problem. Lemma 9.3 says there is at most one such point in the entire region where a global optimizer can lie.

Existence of the original global maximizer supplies at least one such point. Therefore there is exactly one normalized fixed-code optimizer. The sixteen surviving sign strings are merely cyclic/reversal/sign representatives of the same difference body and reconstructed polygon; Lemma 3.1 removes translation ambiguity. This proves uniqueness up to congruence and relabeling. ∎

---

## 11. Numerical identification (not used for global validity)

A high-precision Newton solve for the unique point gives the half-circle gaps

\[
\begin{array}{rcl}
0.19831631349051794937,&&
0.19450334649319957491,\\
0.19450334649319957491,&&
0.19774551481013478084,\\
0.19499397164222383569,&&
0.19716378410929702498,\\
0.19716378410929702498,&&
0.19640626564702685354,\\
0.19640626564702685354,&&
0.19716378410929702498,\\
0.19716378410929702498,&&
0.19499397164222383569,\\
0.19774551481013478084,&&
0.19450334649319957491,\\
0.19450334649319957491,&&
0.19831631349051794937.
\end{array}
\]

They sum to \(\pi\), satisfy the closure equations for \(c_*\), and yield

\[
 p_{16}=3.13654771648660738608596703194\ldots .
\]

For a publication-quality numerical theorem, these coordinates should still be enclosed by an independent interval-Newton computation. That is not needed for the existence and global-uniqueness argument above, because existence comes from compactness and uniqueness from Lemma 9.3.

---

## 12. Adversarial audit and remaining external-validation items

The following previously dangerous points have been explicitly handled:

1. **Degenerate difference bodies:** excluded by the strict lower bound versus the best 30-gon in the disk.
2. **Perturbations leaving the original polygon problem:** prevented by Lemma 3.1.
3. **Adjacent moved vertices and endpoint \(z_0\):** included in Lemma 4.1; an outer incident edge still changes.
4. **A merely convex, locally constant perturbation:** a generic direction makes at least one norm term strictly convex.
5. **Use of KKT without constraint qualification:** MFCQ is constructed explicitly in Lemma 6.1.
6. **A hidden unsaturated radius:** contradicted quantitatively by projective angular separation.
7. **The earlier wrong Poincare constant:** the proof uses \(1/(4\sin^2(\pi/32))\), not \(1/(4\sin^2(\pi/16))\).
8. **Floating-point code decisions:** all are replaced by outward integer intervals.
9. **Objective normalization:** all gradients and Hessians here are derived from the full objective (9.1), avoiding the factor-two inconsistency in the displayed derivative formulas of the prior numerical paper.
10. **Local versus global fixed-code optimization:** Lemma 9.3 proves uniqueness throughout the entire high-perimeter region, not merely convergence of Newton's method near one point.

What remains before the result should be announced as established mathematics is external rather than a known logical gap:

- independent reproduction of the exact verifier, preferably in Arb/MPFI, Lean, Isabelle, or Coq;
- line-by-line expert review of Lemmas 3.1, 4.1, and 6.1;
- archival release of source, output, environment, and hashes;
- optional interval-Newton enclosure of the final decimal optimizer.

---

## References

1. C. Bingane, *Maximal perimeter and maximal width of a convex small polygon*, arXiv:2106.11831.
2. B. Mulansky and A. Potschka, *A zonogon approach for computing small convex polygons of maximum perimeter*, Mathematical Programming (2025), DOI 10.1007/s10107-025-02244-x.
3. K. Reinhardt, original work on extremal small polygons (1922).
