# Hexahedral elements Hexa8, Hexa20 and Hexa27: mathematical foundations of the implementation

This report details the mathematical construction of the shape functions, their derivatives, the Jacobian and the physical gradient for the trilinear (Hexa8), quadratic serendipity (Hexa20) and full triquadratic Lagrange (Hexa27) hexahedral elements. The shape functions are implemented in `ShapeFunctions.h`, while `computeGradientsAndJacobianHexa{8,20,27}` remains in `ArcaneFemFunctions.h`. The node numbering follows the VTK convention, consistent with the topology declared in `arcane/src/arcane/core/ItemTypeMng.cc` (`IT_Hexaedron8`, `IT_Hexaedron20`, `IT_Hexaedron27`).

---

## 1. The reference element and the coordinates (ξ, η, ζ)

### 1.1 Why a reference element?

A mesh contains hexahedra of arbitrary shapes and sizes (possibly distorted). Rather than building interpolation polynomials on each physical cell — which would be costly and unstable — the **isoparametric** finite element method defines the shape functions *once* on a fixed reference element, then transports all computations to the physical element through a geometric mapping.

The hexahedral reference element is the bi-unit cube:

$$
\hat{K} = [-1, 1]^3 = \{ (\xi, \eta, \zeta) \; : \; -1 \le \xi, \eta, \zeta \le 1 \}
$$

The three **reference coordinates** (also called *natural* or *parametric* coordinates) are:

| Symbol | Role | Values on the faces |
|---|---|---|
| ξ (xi)   | local abscissa, the "x" direction of the reference cube | ξ = −1 (left face), ξ = +1 (right face) |
| η (eta)  | local ordinate, the "y" direction | η = −1 (front face), η = +1 (back face) |
| ζ (zeta) | local elevation, the "z" direction | ζ = −1 (bottom face), ζ = +1 (top face) |

The center of the element is (ξ, η, ζ) = (0, 0, 0). The interval [−1, 1] (rather than [0, 1]) is chosen because it is the native interval of **Gauss–Legendre quadrature** (see §6) and because it makes the formulas symmetric.

### 1.2 The isoparametric mapping

Let $\mathbf{x}_a = (x_a, y_a, z_a)$, $a = 1, \dots, n$ be the physical coordinates of the $n$ nodes of the cell ($n \in \{8, 20, 27\}$). The geometry of the physical cell is described by the **same basis** as the solution (hence the name *iso*-parametric):

$$
\mathbf{x}(\xi, \eta, \zeta) \;=\; \sum_{a=1}^{n} N_a(\xi, \eta, \zeta)\, \mathbf{x}_a
$$

and the unknown field (e.g. temperature or displacement) is interpolated identically:

$$
u_h(\xi, \eta, \zeta) \;=\; \sum_{a=1}^{n} N_a(\xi, \eta, \zeta)\, u_a .
$$

The whole difficulty is therefore to construct the $N_a$, and then to know how to differentiate *in physical coordinates* a function defined *in reference coordinates* — that is the role of the Jacobian (§5).

---

## 2. Node positions — VTK convention

The VTK convention orders the nodes as follows: first the 8 corners, then the bottom edges, the top edges, then the vertical ones, then (Hexa27) the 6 face centers and the volume center. This is the order produced by `ItemTypeMng.cc` and assumed by every function in the code.

Throughout the figures below, node numbers are coloured by category — **black** for corners (0–7), **blue** for edge midpoints (8–19), **red** for face centres (20–25), **green** for the volume centre (26). Light dotted lines mark the three coordinate midplanes ξ = 0, η = 0, ζ = 0 on the six cube faces and the three interior axes through the centre, to help locate any node visually.

### 2.1 The 8 corners (common to all three elements)

![Hexa8 — 8 corner nodes at $(\pm 1, \pm 1, \pm 1)$ (VTK ordering).](https://github.com/user-attachments/assets/dbc4a1e4-3a18-4692-9589-1f2ba49214af)

| Corner | (ξ, η, ζ)     |     | Corner | (ξ, η, ζ)     |
|:------:|:-------------:|:---:|:------:|:-------------:|
|   0    | (−1, −1, −1)  |     |   4    | (−1, −1, +1)  |
|   1    | (+1, −1, −1)  |     |   5    | (+1, −1, +1)  |
|   2    | (+1, +1, −1)  |     |   6    | (+1, +1, +1)  |
|   3    | (−1, +1, −1)  |     |   7    | (−1, +1, +1)  |

Mnemonic rule: nodes 0–3 form the **bottom** face (ζ = −1) traversed counterclockwise when seen from above, and nodes 4–7 the **top** face (ζ = +1) in the same order (each top node sits "above" its bottom counterpart: 4 above 0, 5 above 1, etc.).

### 2.2 The 12 edge nodes (Hexa20 and Hexa27, indices 8–19)

Each edge node is the **midpoint** of its edge, so one of its coordinates is 0 and the other two are ±1.

![Hexa20 — 8 corners + 12 edge midpoints (indices 8–19), VTK ordering.](https://github.com/user-attachments/assets/55c36c3f-fe2a-4006-943e-4fab3bb358dd)

| Bottom edges (ζ = −1) | (ξ, η, ζ)     |     | Top edges (ζ = +1) | (ξ, η, ζ)     |     | Vertical edges (ζ = 0) | (ξ, η, ζ)     |
|:---------------------:|:-------------:|:---:|:------------------:|:-------------:|:---:|:----------------------:|:-------------:|
|          8            | ( 0, −1, −1)  |     |         12         | ( 0, −1, +1)  |     |          16            | (−1, −1,  0)  |
|          9            | (+1,  0, −1)  |     |         13         | (+1,  0, +1)  |     |          17            | (+1, −1,  0)  |
|         10            | ( 0, +1, −1)  |     |         14         | ( 0, +1, +1)  |     |          18            | (+1, +1,  0)  |
|         11            | (−1,  0, −1)  |     |         15         | (−1,  0, +1)  |     |          19            | (−1, +1,  0)  |

Order: first the 4 bottom edges (8–11), then the 4 top edges (12–15), then the 4 vertical edges (16–19).

The **Hexa20** stops here: 8 corners + 12 edge midpoints = 20 nodes.

### 2.3 The 6 face centers and the volume center (Hexa27 only, indices 20–26)

A face center has **two** zero coordinates; the volume center has all three.

![Hexa27 — full 3×3×3 grid: 8 corners + 12 edges + 6 face centres (20–25) + volume centre (26).](https://github.com/user-attachments/assets/b938fe9c-74c8-4cf7-bf7e-ddc64c1b6743)

| Face centre | (ξ, η, ζ)   | Face          |     | Face centre | (ξ, η, ζ)   | Face             |
|:-----------:|:-----------:|:--------------|:---:|:-----------:|:-----------:|:-----------------|
|     20      | (−1, 0, 0)  | ξ = −1 (left) |     |     24      | ( 0, 0, −1) | ζ = −1 (bottom)  |
|     21      | (+1, 0, 0)  | ξ = +1 (right)|     |     25      | ( 0, 0, +1) | ζ = +1 (top)     |
|     22      | ( 0,−1, 0)  | η = −1 (front)|     |     26      | ( 0, 0,  0) | *volume centre*  |
|     23      | ( 0,+1, 0)  | η = +1 (back) |     |             |             |                  |

Hexa27 total: 8 + 12 + 6 + 1 = 27 nodes — the complete 3 × 3 × 3 grid of points $\{-1, 0, +1\}^3$.

A word of caution: this face ordering is the VTK one. Other references use different orderings — the Bochum source in particular orders the face centres bottom, top, front, right, back, left (see §3.3 for the exact permutation). Numbering mismatches are the primary source of errors when implementing these elements, which is why the Kronecker property (§7) must always be tested *against the node table actually used by the code*.

---

## 3. Construction of the shape functions

### 3.1 Principle: Lagrange interpolation and the Kronecker property

We seek $n$ polynomials $N_a$ such that

$$
N_a(\boldsymbol{\xi}_b) = \delta_{ab} =
\begin{cases} 1 & \text{if } a = b\\ 0 & \text{otherwise} \end{cases}
$$

where $\boldsymbol{\xi}_b$ is the position of node $b$. This property (**Kronecker** delta) guarantees that the coefficients $u_a$ of the interpolation are exactly the nodal values of the field. The construction strategy differs from element to element.

### 3.2 Hexa8 — tensor product of linear Lagrange polynomials

**Official shape functions.** Following Felippa (*Advanced Finite Element Methods*, solid elements) and Šiljak (*Shape Functions generation, requirements, etc.*, Ruhr-Universität Bochum, 2009), whose corner numbering coincides with the VTK ordering (both references number the nodes $1{-}8$; VTK index $=$ reference index $-\,1$, used below), the eight shape functions are

$$
\begin{aligned}
N_0 &= \tfrac18(1-\xi)(1-\eta)(1-\zeta), &\quad N_4 &= \tfrac18(1-\xi)(1-\eta)(1+\zeta),\\
N_1 &= \tfrac18(1+\xi)(1-\eta)(1-\zeta), &\quad N_5 &= \tfrac18(1+\xi)(1-\eta)(1+\zeta),\\
N_2 &= \tfrac18(1+\xi)(1+\eta)(1-\zeta), &\quad N_6 &= \tfrac18(1+\xi)(1+\eta)(1+\zeta),\\
N_3 &= \tfrac18(1-\xi)(1+\eta)(1-\zeta), &\quad N_7 &= \tfrac18(1-\xi)(1+\eta)(1+\zeta),
\end{aligned}
$$

summarized in both references by the single expression

$$
\boxed{\;N^i(\xi,\eta,\zeta) = \tfrac{1}{8}\,(1 + \xi^i \xi)(1 + \eta^i \eta)(1 + \zeta^i \zeta)\;}
$$

where $(\xi^i, \eta^i, \zeta^i) \in \{-1,+1\}^3$ are the coordinates of node $i$.

**Construction.** On $[-1,1]$ with two nodes $\xi = \pm 1$, the two degree-1 Lagrange polynomials are

$$
\ell_-(\xi) = \frac{1 - \xi}{2}, \qquad \ell_+(\xi) = \frac{1 + \xi}{2},
$$

satisfying $\ell_-(-1) = 1$, $\ell_-(+1) = 0$ and vice versa (Kronecker property in 1D). The shape function of a corner node is the product of the three corresponding 1D polynomials; the factor $\tfrac18 = (\tfrac12)^3$ collects the three half-factors. The Kronecker property in 3D follows from the 1D one: if $b \ne a$, node $b$ differs from node $a$ in at least one direction and the corresponding factor vanishes.

**Correspondence with the code.** For node 0 at $(-1,-1,-1)$, the formula $N_0 = \tfrac18(1-\xi)(1-\eta)(1-\zeta)$ corresponds to the line `N(0) = 0.125 * one_minus_xi * one_minus_eta * one_minus_zeta`; the seven other lines follow the same pattern.

**Spanned space.** Expanding the products yields the basis
$$
\{1,\ \xi,\ \eta,\ \zeta,\ \xi\eta,\ \xi\zeta,\ \eta\zeta,\ \xi\eta\zeta\}:
$$
8 monomials for 8 nodes, the space $Q_1$. It is **linear in each direction separately** (hence *trilinear*) but contains **no squared term** — the element cannot represent a curvature such as $x^2$, which bounds its interpolation error at $\mathcal{O}(h^2)$ in the $L^2$ norm.

### 3.3 Hexa27 — tensor product of quadratic Lagrange polynomials

**Official shape functions.** The Hexa27 is the complete triquadratic Lagrange element $Q_2$ (DefElement, Hughes, Zienkiewicz–Taylor): each node $a$ lies on the grid $\{-1,0,+1\}^3$ and

$$
\boxed{\;N_a(\xi,\eta,\zeta) = \ell_{\xi_a}(\xi)\; \ell_{\eta_a}(\eta)\; \ell_{\zeta_a}(\zeta),\;}
$$

with the three degree-2 Lagrange polynomials on $\{-1, 0, +1\}$:

$$
\ell_{-1}(\xi) = \frac{\xi(\xi-1)}{2}, \qquad
\ell_{0}(\xi) = 1 - \xi^2, \qquad
\ell_{+1}(\xi) = \frac{\xi(\xi+1)}{2}.
$$

(Each follows from the general formula $\ell_i(\xi) = \prod_{j \ne i} \frac{\xi - \xi_j}{\xi_i - \xi_j}$; e.g. $\ell_{-1} = \frac{(\xi-0)(\xi-1)}{(-1-0)(-1-1)} = \frac{\xi(\xi-1)}{2}$.)

The per-family formulas below are constructed from the official shape functions of the Bochum reference (Šiljak, "Shape functions of Hexahedron with variable number of nodes (8–27 nodes)"), with all hierarchical corrections applied (see the equivalence paragraph below) and the node indices remapped to the VTK numbering used throughout this code — in particular the face-centre nodes, whose ordering differs between the two conventions. Since for a node coordinate $\xi^i = \pm 1$ one has $\ell_{\xi^i}(\xi) = \tfrac12\,\xi^i\xi\,(1 + \xi^i\xi)$, the 27 functions split into four families according to the node type (VTK indices; $(\xi^i, \eta^i, \zeta^i)$ are the coordinates of node $i$):

$$
\begin{aligned}
\text{corner nodes } i = 0,\dots,7: \quad
& N^i = \tfrac18\,\xi^i\xi\,\eta^i\eta\,\zeta^i\zeta\,(1+\xi^i\xi)(1+\eta^i\eta)(1+\zeta^i\zeta),\\[2pt]
\text{edge nodes } i = 8, 10, 12, 14\ (\xi^i = 0): \quad
& N^i = \tfrac14\,(1-\xi^2)\;\eta^i\eta\,\zeta^i\zeta\,(1+\eta^i\eta)(1+\zeta^i\zeta),\\[2pt]
\text{edge nodes } i = 9, 11, 13, 15\ (\eta^i = 0): \quad
& N^i = \tfrac14\,(1-\eta^2)\;\xi^i\xi\,\zeta^i\zeta\,(1+\xi^i\xi)(1+\zeta^i\zeta),\\[2pt]
\text{edge nodes } i = 16,\dots,19\ (\zeta^i = 0): \quad
& N^i = \tfrac14\,(1-\zeta^2)\;\xi^i\xi\,\eta^i\eta\,(1+\xi^i\xi)(1+\eta^i\eta),\\[2pt]
\text{face centres } i = 20, 21\ (\eta^i = \zeta^i = 0): \quad
& N^i = \tfrac12\,\xi^i\xi\,(1+\xi^i\xi)(1-\eta^2)(1-\zeta^2),\\[2pt]
\text{face centres } i = 22, 23\ (\xi^i = \zeta^i = 0): \quad
& N^i = \tfrac12\,\eta^i\eta\,(1+\eta^i\eta)(1-\xi^2)(1-\zeta^2),\\[2pt]
\text{face centres } i = 24, 25\ (\xi^i = \eta^i = 0): \quad
& N^i = \tfrac12\,\zeta^i\zeta\,(1+\zeta^i\zeta)(1-\xi^2)(1-\eta^2),\\[2pt]
\text{volume centre } i = 26: \quad
& N^{26} = (1-\xi^2)(1-\eta^2)(1-\zeta^2).
\end{aligned}
$$

**Examples** — one node of each family, as implemented:

$$
\begin{aligned}
N_0 &= \ell_{-1}(\xi)\,\ell_{-1}(\eta)\,\ell_{-1}(\zeta) = \tfrac18\,\xi(\xi-1)\,\eta(\eta-1)\,\zeta(\zeta-1) && \text{(corner, } (-1,-1,-1))\\
N_8 &= \ell_{0}(\xi)\,\ell_{-1}(\eta)\,\ell_{-1}(\zeta) = \tfrac14\,(1-\xi^2)\,\eta(\eta-1)\,\zeta(\zeta-1) && \text{(edge, } (0,-1,-1))\\
N_{20} &= \ell_{-1}(\xi)\,\ell_{0}(\eta)\,\ell_{0}(\zeta) = \tfrac12\,\xi(\xi-1)\,(1-\eta^2)\,(1-\zeta^2) && \text{(face centre, } (-1,0,0))\\
N_{26} &= \ell_{0}(\xi)\,\ell_{0}(\eta)\,\ell_{0}(\zeta) = (1-\xi^2)(1-\eta^2)(1-\zeta^2) && \text{(volume centre, } (0,0,0))
\end{aligned}
$$

The leading coefficients $\tfrac18, \tfrac14, \tfrac12, 1$ follow directly from the number of $\pm 1$ coordinates of the node — each contributes a factor $\tfrac12$ through the 1D polynomial $\tfrac12\,\xi(\xi \pm 1)$, while a zero coordinate uses $1-\xi^2$.

The centre node, a product of three factors $1-\xi^2$ each vanishing at $\xi = \pm 1$, is itself zero on the entire boundary of the element: it is the sole **interior-supported basis function** of the Hexa27, associated with the volume degree of freedom.

**Equivalence with the hierarchical form of the literature.** The Bochum reference presents the variable-node Hexa8–27 element in *hierarchical* form: trilinear corner functions $\tfrac18(1-\xi)(1-\eta)(1-\zeta)$, edge functions $\tfrac14(1-\xi^2)(1-\eta)(1-\zeta)$, face functions $\tfrac12(1-\xi^2)(1-\eta^2)(1-\zeta)$ and the interior bubble $(1-\xi^2)(1-\eta^2)(1-\zeta^2)$, followed by correction steps enforcing the Kronecker property at the higher-level nodes (corners: $-\tfrac12$ of the three adjacent edge functions, then $+\tfrac14$ of the three adjacent face functions, then $-\tfrac18$ of the bubble; edges: $-\tfrac12$ of the two adjacent face functions, then $+\tfrac14$ of the bubble; faces: $-\tfrac12$ of the bubble). With all 27 nodes present, the fully corrected functions coincide with the tensor-product form above, by virtue of the 1D identity

$$
\tfrac12(1-\xi) \;-\; \tfrac12(1-\xi^2) \;=\; \tfrac12\,\xi(\xi-1) \;=\; \ell_{-1}(\xi).
$$

When comparing node by node with these references, note that their corner and edge numbering 1–20 coincides with VTK 0–19 (offset by one), but their face-centre ordering differs: reference nodes 21…26 = bottom, top, front, right, back, left, whereas VTK 20…25 = left, right, front, back, bottom, top; the permutation reference $\to$ VTK is $\{21{\to}24,\ 22{\to}25,\ 23{\to}22,\ 24{\to}21,\ 25{\to}23,\ 26{\to}20,\ 27{\to}26\}$.

**Spanned space.** $Q_2$, all monomials $\xi^p \eta^q \zeta^r$ with $p,q,r \le 2$ — 27 monomials for 27 nodes. The element represents any quadratic field exactly (including $x^2 + y^2$), and its interpolation error is $\mathcal{O}(h^3)$ in the $L^2$ norm.

### 3.4 Hexa20 — the serendipity element

**Official shape functions.** Following Felippa ("The 20-node (serendipity) hexahedron") and Šiljak, whose numbering 1–20 coincides with VTK 0–19 (VTK indices used below), with $(\xi^i, \eta^i, \zeta^i)$ the coordinates of node $i$:

$$
\begin{aligned}
\text{corner nodes } i = 0,\dots,7: \quad
& N^i = \tfrac18(1+\xi^i\xi)(1+\eta^i\eta)(1+\zeta^i\zeta)\,(\xi^i\xi + \eta^i\eta + \zeta^i\zeta - 2),\\[2pt]
\text{midside nodes } i = 8, 10, 12, 14\ (\xi^i = 0): \quad
& N^i = \tfrac14(1-\xi^2)(1+\eta^i\eta)(1+\zeta^i\zeta),\\[2pt]
\text{midside nodes } i = 9, 11, 13, 15\ (\eta^i = 0): \quad
& N^i = \tfrac14(1-\eta^2)(1+\xi^i\xi)(1+\zeta^i\zeta),\\[2pt]
\text{midside nodes } i = 16,\dots,19\ (\zeta^i = 0): \quad
& N^i = \tfrac14(1-\zeta^2)(1+\xi^i\xi)(1+\eta^i\eta).
\end{aligned}
$$

**Examples** — one node of each family, as implemented:

$$
\begin{aligned}
N_0 &= \tfrac18(1-\xi)(1-\eta)(1-\zeta)(-\xi-\eta-\zeta-2) && \text{(corner, } (-1,-1,-1))\\
N_8 &= \tfrac14(1-\xi^2)(1-\eta)(1-\zeta) && \text{(midside, } (0,-1,-1))\\
N_9 &= \tfrac14(1+\xi)(1-\eta^2)(1-\zeta) && \text{(midside, } (+1,0,-1))\\
N_{16} &= \tfrac14(1-\xi)(1-\eta)(1-\zeta^2) && \text{(midside, } (-1,-1,0))
\end{aligned}
$$

The corner formula is the trilinear Hexa8 function corrected by the last factor $(\xi^i\xi + \eta^i\eta + \zeta^i\zeta - 2)$, which equals $+1$ at the corner itself ($1+1+1-2$) and $0$ at the midpoints of the three adjacent edges ($1+1+0-2$): it is what restores the Kronecker property once the midside nodes are introduced.

**Spanned space.** The 20 functions span all monomials of degree $\le 2$ in one direction and $\le 1$ in the others ($\xi^2$, $\xi^2\eta$, $\xi^2\eta\zeta$, …), but **not** the higher-order mixed monomials $\xi^2\eta^2$, $\xi^2\eta^2\zeta^2$, etc. The element contains the full $P_2$ (complete polynomials of degree 2), which is sufficient for $\mathcal{O}(h^3)$ convergence in $L^2$ — hence the name *serendipity*: the same asymptotic order as the Hexa27, achieved with 20 nodes instead of 27.



---

## 4. Derivatives of the shape functions

The derivatives $\partial N_a / \partial \xi$, $\partial N_a / \partial \eta$, $\partial N_a / \partial \zeta$ are obtained **analytically** (hand differentiation of the expressions of §3), and hard-coded in `computeGradientsAndJacobianHexa{8,20,27}` — never by finite differences, which would be both more expensive and less accurate.

The differentiation exploits the product structure: each $N_a$ is a product $f(\xi)\,g(\eta)\,h(\zeta)$ of three 1D polynomials (for the Hexa27 tensor-product functions and the Hexa20 edge functions; the Hexa20 corner functions carry an additional linear factor treated by the product rule), hence
$$
\frac{\partial N_a}{\partial \xi} = f'(\xi)\, g(\eta)\, h(\zeta),
$$

with the elementary 1D derivatives:

| $f(\xi)$ | $f'(\xi)$ |
|---|---|
| $\tfrac12(1 \pm \xi)$ | $\pm\tfrac12$ |
| $\tfrac12\,\xi(\xi - 1)$ | $\tfrac12(2\xi - 1)$ |
| $1 - \xi^2$ | $-2\xi$ |
| $\tfrac12\,\xi(\xi + 1)$ | $\tfrac12(2\xi + 1)$ |

**Examples from the code:**

- Hexa8, node 0: 
$$
N_0 = \tfrac18(1-\xi)(1-\eta)(1-\zeta) \Rightarrow \partial_\xi N_0 = -\tfrac18(1-\eta)(1-\zeta)
$$ 
— line `dN_dxi[0] = -0.125 * one_minus_eta * one_minus_zeta`.
- Hexa27, node 8: 
$$
N_8 = \tfrac14(1-\xi^2)\,\eta(\eta-1)\,\zeta(\zeta-1) \Rightarrow \partial_\xi N_8 = \tfrac14(-2\xi)\,\eta(\eta-1)\,\zeta(\zeta-1)
$$.
- Hexa20, corner 0: here $N_0$ is a product of 4 factors; the product rule gives, after grouping, 
$$
\partial_\xi N_0 = \tfrac18(1-\eta)(1-\zeta)(2\xi + \eta + \zeta + 1)
$$ 
— the $2\xi$ 
term betrays the presence of $\xi^2$ in the expanded $N_0$.


---

## 5. Jacobian and physical gradient

### 5.1 The Jacobian matrix

The isoparametric mapping $\mathbf{x}(\boldsymbol{\xi})$ of §1.2 has, at the evaluation point $(\xi, \eta, \zeta)$, the Jacobian matrix:

$$
J \;=\; \frac{\partial \mathbf{x}}{\partial \boldsymbol{\xi}} \;=\;
\begin{pmatrix}
\partial x/\partial \xi & \partial y/\partial \xi & \partial z/\partial \xi\\
\partial x/\partial \eta & \partial y/\partial \eta & \partial z/\partial \eta\\
\partial x/\partial \zeta & \partial y/\partial \zeta & \partial z/\partial \zeta
\end{pmatrix}.
$$

Differentiating $\mathbf{x} = \sum_a N_a \mathbf{x}_a$ term by term, each entry is a sum over the nodes:

$$
J_{00} = \sum_{a=1}^{n} \frac{\partial N_a}{\partial \xi}\, x_a, \qquad
J_{01} = \sum_{a=1}^{n} \frac{\partial N_a}{\partial \xi}\, y_a, \qquad
J_{10} = \sum_{a=1}^{n} \frac{\partial N_a}{\partial \eta}\, x_a, \quad \text{etc.}
$$

which is exactly the loop in the code:

```cpp
for (a = 0; a < n; ++a) {
  J[0][0] += dN_dxi[a]   * coord[a].x;   // ∂x/∂ξ
  J[0][1] += dN_dxi[a]   * coord[a].y;   // ∂y/∂ξ
  ...
  J[2][2] += dN_dzeta[a] * coord[a].z;   // ∂z/∂ζ
}
```

**Geometric interpretation.** $J$ describes how the reference cube is stretched, sheared and rotated to become the physical cell, *locally around the evaluation point*. For an undistorted Hexa8 (axis-aligned parallelepiped of dimensions $h_x \times h_y \times h_z$), $J = \mathrm{diag}(h_x/2,\, h_y/2,\, h_z/2)$ is constant; for a distorted cell or a quadratic element with curved edges, $J$ varies from one Gauss point to the next — which is why it is recomputed at every integration point.

### 5.2 The determinant

$$
\det J \;=\; \text{local volume dilation factor}, \qquad dV_{\text{physical}} = \det J \; d\xi\, d\eta\, d\zeta.
$$

It enters every integral (stiffness matrix, right-hand side, error norm). The code checks $\det J > 0$ and raises an error otherwise (`ARCANE_FATAL`): a negative or zero determinant signals an inverted or degenerate cell (misordered nodes, invalid mesh) — better to fail loudly than to assemble a wrong matrix.

### 5.3 From the reference gradient to the physical gradient

The $N_a$ are functions of $\boldsymbol{\xi}$, but the equation to solve (e.g. $-\Delta u = f$) involves derivatives **in physical space** $\partial/\partial x$. The chain rule links the two:

$$
\begin{pmatrix} \partial N_a/\partial \xi \\ \partial N_a/\partial \eta \\ \partial N_a/\partial \zeta \end{pmatrix}
= J \begin{pmatrix} \partial N_a/\partial x \\ \partial N_a/\partial y \\ \partial N_a/\partial z \end{pmatrix}
\qquad\Longrightarrow\qquad
\boxed{\;\nabla_{\!\mathbf{x}} N_a = J^{-1}\, \nabla_{\!\boldsymbol{\xi}} N_a\;}
$$

hence the second loop in the code:

```cpp
dN_dx[a] = invJ[0][0]*dN_dxi[a] + invJ[0][1]*dN_deta[a] + invJ[0][2]*dN_dzeta[a];
dN_dy[a] = invJ[1][0]*dN_dxi[a] + invJ[1][1]*dN_deta[a] + invJ[1][2]*dN_dzeta[a];
dN_dz[a] = invJ[2][0]*dN_dxi[a] + invJ[2][1]*dN_deta[a] + invJ[2][2]*dN_dzeta[a];
```

The inverse $J^{-1}$ is computed explicitly (3×3 cofactor formula, via `math::inverseMatrix(J, detJ)`) — for a 3×3 matrix this is both exact and cheaper than a factorization.

### 5.4 Gradient of a scalar field

Once the $\nabla_{\!\mathbf{x}} N_a$ are available, the gradient of the interpolated field is a simple linear combination of the nodal values:

$$
\nabla u_h \;=\; \sum_{a=1}^{n} u_a \, \nabla_{\!\mathbf{x}} N_a,
$$

which is what `computeGradientHexa8` implements (evaluated at the element center $(0,0,0)$). This is also the term appearing in the bilinear form of the stiffness matrix:

$$
K_{ab}^{e} = \int_{K^e} \nabla N_a \cdot \nabla N_b \; dV
= \int_{[-1,1]^3} \bigl(J^{-1}\nabla_{\!\boldsymbol{\xi}} N_a\bigr) \cdot \bigl(J^{-1}\nabla_{\!\boldsymbol{\xi}} N_b\bigr)\, \det J \;\, d\xi\, d\eta\, d\zeta.
$$

---

## 6. Numerical integration: Gauss–Legendre quadrature

The integral above is not computable analytically for a distorted cell (the integrand contains $J^{-1}$, rational in $\boldsymbol{\xi}$). It is replaced by Gauss–Legendre quadrature, tensorized in 3D:

$$
\int_{[-1,1]^3} f\, d\xi\, d\eta\, d\zeta \;\approx\; \sum_{i,j,k=1}^{p} w_i w_j w_k \; f(\xi_i, \eta_j, \zeta_k).
$$

A $p$-point rule in 1D integrates polynomials of degree $2p - 1$ exactly:

| Element | Degree of $N$ | Degree of $\nabla N \!\cdot\! \nabla N$ (per direction) | Usual minimal rule |
|---|---|---|---|
| Hexa8  | 1 per direction | 2 | $2^3 = 8$ points ($\xi_i = \pm 1/\sqrt{3}$, $w_i = 1$) |
| Hexa20 | 2 per direction | 4 | $3^3 = 27$ points ($\xi_i \in \{0, \pm\sqrt{3/5}\}$, $w_i \in \{8/9, 5/9\}$) |
| Hexa27 | 2 per direction | 4 | $3^3 = 27$ points |

(On a distorted cell the integrand is no longer polynomial and these rules become approximations — but of sufficient order to preserve the convergence rate of the element.)

---


## 8. References

- **Gaussian quadrature** — <https://en.wikipedia.org/wiki/Gaussian_quadrature> (Gauss–Legendre points and weights, $2p-1$ exactness).
- **Lagrange polynomials** — <https://en.wikipedia.org/wiki/Lagrange_polynomial> (general formula $\ell_i = \prod_{j\ne i}(\xi - \xi_j)/(\xi_i - \xi_j)$).
- **Kronecker delta** — <https://en.wikipedia.org/wiki/Kronecker_delta> ($\delta_{ab} = 1$ if $a = b$, $0$ otherwise; nodal-basis defining property).
- **Hexahedral elements (formal definitions, bases, numberings)** — DefElement: <https://defelement.org/elements/lagrange.html> (Lagrange $Q_k$) and <https://defelement.org/elements/serendipity.html> (serendipity family).
- **Hexa27 shape function** - Q2 lagrange on Hexahedron :
<https://defelement.org/elements/examples/hexahedron-lagrange-lobatto-2.html>
- **VTK node numbering convention** — VTK File Formats, cell types: <https://examples.vtk.org/site/Cxx/GeometricObjects/IsoparametricCellsDemo/> (`VTK_HEXAHEDRON`, `VTK_QUADRATIC_HEXAHEDRON`, `VTK_TRIQUADRATIC_HEXAHEDRON`).
- **Official shape-function formulas — Felippa** — C.A. Felippa, *Advanced Finite Element Methods (ASEN 6367)*, University of Colorado at Boulder, lecture notes (draft): chapters on solid (brick) elements — trilinear 8-node hexahedron, 20-node serendipity hexahedron, and variable-node elements.
- **Official shape-function formulas — Bochum** — E. Šiljak, *Shape Functions generation, requirements, etc.*, student presentation, Finite Elements in Structural Mechanics, Fakultät für Bauingenieurwesen, Ruhr-Universität Bochum, January 2009: 8-node and 20-node (serendipity) hexahedra, and the hexahedron with variable number of nodes (8–27). Note: this source gives the Hexa27 in *hierarchical* form (uncorrected building blocks + correction steps) with its own face-centre ordering; §3.3 details the equivalence with the tensor-product form used in the code and the renumbering to VTK.
- **Finite element method, isoparametric formulation** — T.J.R. Hughes, *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis*, Dover, 2000 (chap. 3: isoparametric elements, quadrature).
- **Method of manufactured solutions** — P.J. Roache, *Code Verification by the Method of Manufactured Solutions*, J. Fluids Eng. 124(1), 2002; and K. Salari, P. Knupp, *Code Verification by the Method of Manufactured Solutions*, Sandia Report SAND2000-1444, 2000.
