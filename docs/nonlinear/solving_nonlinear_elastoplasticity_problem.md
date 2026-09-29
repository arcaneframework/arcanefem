# Newton Solver and nonlinear assembly strategies in the ArcaneFEM
## 1. Problem setting and weak formulation

### 1.1 Strong and weak form

The quasi-static equilibrium problem on $\Omega$ with boundary
$\Gamma = \Gamma_D \cup \Gamma_N$ reads

$$
-\nabla \cdot  \sigma(\mathbf{u}) = \mathbf{f} \ \text{in } \Omega, \qquad
 \sigma\cdot\mathbf{n} = \mathbf{t} \ \text{on } \Gamma_N, \qquad
\mathbf{u} = \mathbf{g} \ \text{on } \Gamma_D.
$$

Its variational form is: finding $\mathbf{u}\in V$, $\mathbf{u}=\mathbf{g}$ on
$\Gamma_D$, such that

$$
\int_\Omega  \sigma(\mathbf{u}) \otimes  \varepsilon(\mathbf{v})\, d\Omega
= \int_\Omega \mathbf{f}\cdot\mathbf{v}\, d\Omega + \int_{\Gamma_N} \mathbf{t}\cdot\mathbf{v}\, d\Gamma
\qquad \forall\, \mathbf{v}\in V_0,
$$

with $\varepsilon_{ij}(\mathbf{v}) = \tfrac12(\partial v_i/\partial x_j + \partial v_j/\partial x_i)$.

Unlike linear elasticity, $\sigma$ is not $C\!:\!\varepsilon$ for a
fixed $C$, rather it is the output of a rate-independent elastoplastic constitutive
update with internal history variables, which is the source of the
nonlinearity resolved by Newton iteration.

Note that the implementation follows Voigt notation (plane strain, 2D), with
$$
\varepsilon = (\varepsilon_{xx}, \varepsilon_{yy}, \sqrt2\,\varepsilon_{xy})^\top.
$$

### 1.2 Linearisation and Newton method
In order to solve the nonlinear problem, it is linearised using the Newton method such that we end up with a discrete linear system to solve. 

The weak formulation that leads to such a system is given as,
find $\mathbf{du}\in V$, $\mathbf{u}=\mathbf{g}$ on
$\Gamma_D$, such that

$$
\int_\Omega \sigma(\mathbf{du}) \otimes \varepsilon(\mathbf{v})\, d\Omega
= \int_\Omega \mathbf{f} \cdot \mathbf{v}\, d\Omega + \int_{\Gamma_N} \mathbf{t}\cdot\mathbf{v}\, d\Gamma - \int_\Omega \sigma(\mathbf{U}) \odot \varepsilon(\mathbf{v})\, \qquad \forall\, \mathbf{v}\in V_0,
$$

with $U$ the nodal displacement vector, and more imporatntly

$$
 \sigma(\mathbf{du}) \otimes  \varepsilon(\mathbf{v}) = (C_{\mathrm{tang} :  \varepsilon(\mathbf{du})) \otimes  \varepsilon}(\mathbf{v})
$$
where $C_{\mathrm{tang}}$ is tangent material tensor for a given constitutive law. 


### 1.3 Discrete system: LHS and residual
With $U$ the nodal displacement vector, the discrete problem is: find $U$
such that

$$
R(du) = F_{\mathrm{ext}} + t_{\mathrm{ext}} - F_{\mathrm{int}}(U),
$$
that translates to linear system 
$$
A\delta u = b.
$$

The left-hand-side integrand reduces to the bilinear operator $A \delta u$ that requires the computation of a consistent $C_{\mathrm{tang}}$ based on the law as well as the element matrix based on the mesh type.  This bilinear form is assembled by `_assembleBilinearOperatorGlobal()`, `_assembleBilinearOperatorLocalVonMises()` or `_assembleBilinearOperatorLocalDruckerPrager()` depending on the strategy and law in use (§4).

The right-hand-side consists of $F_{\mathrm{ext}}$, $t_{\mathrm{ext}}$, and $F_{\mathrm{int}}(U)$, as well as the Dirichlet boundary conditions that require a different treatement. The assembly of the RHS is dispatched through `_assembleLinearOperator()`. It contains `_applyExternalBodyForce()`,  `_applyTraction`, and `_applyInternalBodyForce()`.

The quantity of interest for us is `_applyInternalBodyForce()`, which evaluates $\int_\Omega  \sigma(\mathbf{U}) \odot \varepsilon(\mathbf{v})$ and adds to $b$ element by element.

### 1.4 Dirichlet treatment
A new method `_applyDirichletNewton` is called to set Dirichlet BC on the discrete linear system obatined for the Newton method. Here, instead of setting the corresponding component of the $b$ vector to the Dirichlet value, a difference of Dirichlet value and the corresponding component displacement vector is used. This reflects the fact that in $R(du) = F_{\mathrm{ext}} + t_{\mathrm{ext}} - F_{\mathrm{int}}(U)$, we solve for $\delta u$ rather than displacement itself. 

---

## 2. The Newton solver

The solver runs one full Newton-Raphson iteration to equilibrium per load
step, driven by `_doStationarySolve()` -> `_solveNewton()`. Time itself
carries no physics: `compute()` simply advances `t += dt` until `t >= tmax`,
so each call to `_solveNewton()` corresponds to one quasi-static load
increment.

### 2.1 Kinematic bookkeeping

| Variable | Meaning                                                                                                                                                                                     |
| ----------| ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `m_U`    | Total displacement accumulated over all previously converged load steps. Updated only on Newton convergence, in `_updateTimeVariables()`.                                                   |
| `m_DUn`  | Displacement increment for the current load step, accumulated across Newton iterations. Reset to zero at the start of `_solveNewton()`, updated every iteration in `_incrementVariables()`. |
| `m_DUk`  | The Newton correction at the current iteration: the solution of $A\,\delta u = b$, fetched in `_updateNewtonIncrements()`.                                                                  |

The return-mapping kernels always take the **accumulated step increment**
`m_DUn` as their strain input, never the isolated Newton correction `m_DUk`.
Re-evaluating the constitutive update from the last **committed** state at
every iteration, using the full current increment each time, is what keeps
the algorithm path-independent within a load step and avoids compounding
error across iterations.

### 2.2 Newton loop

The algorithm is identical for the global and local strategies except for a
single step: how the tangent stiffness matrix is assembled. Everywhere else,
`_assembleBilinearOperator()` below stands for whichever concrete method is
selected once at `startInit()`.

```
_getMaterialParameters()
m_DUn = 0 ;  m_DUk = 0 ;  m_newton_iter = 0
_restoreConvergedState<Law>()

_assembleBilinearOperator()        // elastic tangent, iteration 0
_assembleLinearOperator()
m_residual_norm0 = ||R||

while (m_newton_iter < m_newton_max_iters and not converged) {
  m_newton_iter += 1

  _solve()                          // A . delta = b
  _updateNewtonIncrements()         // fetch delta into m_DUk
  _incrementVariables()             // m_DUn += m_DUk

  _assembleBilinearOperator()       // updated tangent stiffness
  _assembleLinearOperator()         // residual at the updated state
  _checkNewtonConvergence()
}

if (converged) {
  _updateTimeVariables()            // m_U += m_DUn
  _commitInternalVariables<Law>()   // *_old_gp = current state at Gauss-points
}
```

What `_assembleBilinearOperator()` does is where the two strategies diverge:

* **Global strategy.** `_assembleBilinearOperator()` expands to two
  sequential calls: `_integrateAndSaveConstitutiveLaw<Law>()`, a
  law-dependent sweep over every cell/Gauss point that evaluates the
  constitutive update and writes the resulting stress and tangent tensor
  into persistent mesh variables; followed by `_assembleBilinearOperatorGlobal()`,
  a **law-independent** sweep that only reads those stored values back to
  build $K_e$. At iteration 0, the update sweep is skipped and
  `_setGlobalElasticMaterialTensorAtGPs()` writes the elastic tensor
  directly.

* **Local strategy.** `_assembleBilinearOperator()` expands to a single,
  **law-dependent** call, `_assembleBilinearOperatorLocalVonMises()` or
  `_assembleBilinearOperatorLocalDruckerPrager()`, which evaluates the
  constitutive update and builds $K_e$ for each cell in one pass, with no
  separate update sweep and no persistent tangent storage. At iteration 0
  the same call is made with an `elastic_assembly` flag that bypasses the
  constitutive update entirely.

Section 4 details each branch.

### 2.3 Convergence check

`_checkNewtonConvergence()` uses a dual, OR-type criterion matching the
PETSc SNES convention: convergence is declared as soon as either the
relative increment or the relative residual falls below tolerance,

$$
\frac{\lVert \delta U \rVert_2}{\text{rtol}\cdot\lVert DU_n\rVert_2 + \text{atol}} \le 1
\qquad \text{or} \qquad
\frac{\lVert R(DU_n)\rVert_2}{\lVert R(DU_0)\rVert_2 + \varepsilon} \le \text{rtol}.
$$

Constrained DOFs are zeroed out of the residual before the norm is taken
(`_applyZeroRHSOnConstrainedDOFs`). For Drucker-Prager, the reference
residual norm $\lVert R(DU_0)\rVert$ is re-baselined after the first
iteration, since the first assembled RHS also carries the algebraic penalty
enforcing the imposed footing settlement and would otherwise dominate the
norm.

### 2.4 Commit on convergence

On convergence, `_commitInternalVariablesVonMises()` or
`_commitInternalVariablesDruckerPrager()` overwrite the `*_old_gp` history
arrays with the converged Gauss-point state (stress, accumulated plastic
strain, or plastic strain tensor), so the next load step's
`_restoreConvergedState<Law>()` starts from the correct committed state,
consistent with the standard time-discrete plasticity requirement that
irreversibility of the return mapping only applies across converged steps,
never within the Newton iterations of a single step.

---

## 3. Constitutive laws and internal state variables

Both laws expose the same interface to the module
(`_restoreConvergedState<Law>`, `_commitInternalVariables<Law>`, and either
`_integrateAndSaveConstitutiveLaw<Law>` or
`_assembleBilinearOperatorLocal<Law>`), and both route their physics through
a single reusable kernel that is called identically by the global update
sweep, the local element-matrix builder, and the GPU paths. This shared
kernel is the design choice that makes the global/local duality possible
without duplicating constitutive logic (§4).

### 3.1 Von Mises (J2) plasticity with linear isotropic hardening

**Gauss-point variables** (allocated in `_initConstitutiveLaw()`):

| Variable | Role |
|---|---|
| `m_sigma_gp` / `m_sigma_zz_gp` | Current stress $(\sigma_{xx},\sigma_{yy},\sigma_{xy})$ and out-of-plane $\sigma_{zz}$ (plane strain). |
| `m_sigma_old_gp` / `m_sigma_zz_old_gp` | Committed stress from the last converged load step. |
| `m_p_old_gp` | Committed accumulated equivalent plastic strain $\bar p$. |
| `m_dp_gp` | Plastic multiplier increment for the current state. |

The history variable is the committed **stress**, not the plastic strain:
this is the classical radial-return formulation, working in stress space
with a scalar isotropic-hardening variable.

**Elastic predictor.** With $C_{\mathrm{elas}}$ built from the Lame constants
$\lambda,\mu$,

$$
\sigma^{\mathrm{tr}} = \sigma_{\mathrm{old}} + C_{\mathrm{elas}} : \varepsilon(\mathbf{DU_n}), \qquad
\sigma_{zz}^{\mathrm{tr}} = \sigma_{zz,\mathrm{old}} + \lambda(\varepsilon_{xx}+\varepsilon_{yy}),
$$

and the trial equivalent stress $\sigma_{\mathrm{eq}}^{\mathrm{tr}} = \sqrt{\tfrac32\, s^{\mathrm{tr}}:s^{\mathrm{tr}}}$
including the plane-strain deviator $s_{zz}^{\mathrm{tr}}$.

**Plastic corrector.** The yield check $f = \sigma_{\mathrm{eq}}^{\mathrm{tr}} - \sigma_0 - H\bar p_{\mathrm{old}}$
is evaluated through a Macaulay-bracket switch (`yield_positive = <f>+`,
implemented as `(f + |f|)/2`) rather than a branch, which keeps the kernel
free of divergent control flow and portable to GPU. The plastic multiplier
is $\Delta p = \langle f\rangle_+ / (3\mu+H)$, the return direction is
$N = s^{\mathrm{tr}}/\sigma_{\mathrm{eq}}^{\mathrm{tr}}$, and with
$\beta = 3\mu\Delta p/\sigma_{\mathrm{eq}}^{\mathrm{tr}}$,

$$
\sigma = \sigma^{\mathrm{tr}} - \beta\, s^{\mathrm{tr}} .
$$

**Algorithmic tangent.** With $A = 3\mu\!\left(\dfrac{3\mu}{3\mu+H} - \beta\right)$,

$$
C_{\mathrm{tang}} = C_{\mathrm{elas}} - A\, N\otimes N - \tfrac{4}{3}\mu\beta\, I_{\mathrm{dev}},
$$

assembled term by term (diagonal, off-diagonal and shear entries) as the
standard consistent elastoplastic tangent for J2 plasticity with linear
isotropic hardening, restricted to plane strain.

### 3.2 Drucker-Prager

**Gauss-point variables:**

| Variable | Role |
|---|---|
| `m_sigma_gp` / `m_sigma_zz_gp` | Current stress. |
| `m_eps_p_gp` / `m_eps_p_zz_gp` | Current plastic strain tensor $\varepsilon^p$. |
| `m_eps_p_old_gp` / `m_eps_p_zz_old_gp` | Committed plastic strain from the last converged step. |

Here the history variable is the **plastic strain**, not the stress: the
formulation works from total strain, with the elastic trial strain built as
$\varepsilon^{\mathrm{tr}} = \varepsilon(\mathbf U + \mathbf{DU_n}) - \varepsilon^p_{\mathrm{old}}$,
using both the previously converged total displacement `m_U` and the current
step increment `m_DUn`. The friction and cohesion parameters are
pre-transformed once in `_getMaterialParameters()` to match the
Drucker-Prager cone to the **outer edge** of the Mohr-Coulomb hexagon in the
deviatoric plane, the standard "outer-tip" fit:

$$
\eta = \frac{3\tan\phi}{\sqrt{9+12\tan^2\phi}}, \qquad
\bar c = \frac{3\,c}{\sqrt{9+12\tan^2\phi}}.
$$

**Elastic predictor.** Deviator $e^{\mathrm{tr}}$ and mean strain of
$\varepsilon^{\mathrm{tr}}$, trial deviatoric-stress norm
$\rho^{\mathrm{tr}} = 2\mu\lVert e^{\mathrm{tr}}\rVert$, trial pressure
$p^{\mathrm{tr}} = K\operatorname{tr}\varepsilon^{\mathrm{tr}}$ ($K$: `bulk`).

**Two-surface return mapping.** The trial state is classified by

$$
\mathrm{criterion}_1 = \frac{\rho^{\mathrm{tr}}}{\sqrt2} + \eta\, p^{\mathrm{tr}} - \bar c, \qquad
\mathrm{criterion}_2 = \eta\, p^{\mathrm{tr}} - \frac{K\eta^2}{\mu\sqrt2}\rho^{\mathrm{tr}} - \bar c,
$$

giving, again through switch factors rather than branches,
`smooth_switch` for a return to the regular (differentiable) part of the
cone and `apex_switch` for a return to the cone's vertex, the classical
treatment of the Drucker-Prager singularity at the tip. On the smooth
branch, $\lambda_{\mathrm{smooth}} = \mathrm{criterion}_1/(\mu+K\eta^2)$
drives a radial correction along $N = e^{\mathrm{tr}}/\lVert e^{\mathrm{tr}}\rVert$
plus a volumetric term proportional to $\eta$; on the apex branch the stress
is projected directly to the apex, $\sigma = (\bar c/\eta)\mathbf 1$, matching
the "apex stress correction" of the standard Drucker-Prager return-mapping
algorithm.

**Algorithmic tangent.** On the smooth branch, with
$\kappa = 2\sqrt2\,\mu^2\lambda_{\mathrm{smooth}}/\rho^{\mathrm{tr}}$,

$$
C_{\mathrm{tang}} = C_{\mathrm{elas}} - \kappa\!\left(\tfrac23 I_{\mathrm{dev}} - N\otimes N\right) - \frac{1}{\mu+K\eta^2}\, \mathbf c \otimes \mathbf c,
$$

with $\mathbf c$ the friction-modified correction direction. On the apex
branch the whole tangent contribution is multiplied by `(1 - apex_switch)`,
i.e. set to zero rather than to a theoretically consistent apex tangent (the
apex tangent is singular in the classical sense). This is a pragmatic
simplification that trades convergence rate near the apex for robustness and
is worth flagging explicitly if this note is extended toward verification.

---

## 4. Assembling the tangent material tensor, history and internal state variables

The strategy is selected once at `startInit()` through
`m_gp_material_tensor_strategy` (`"global"` or `"local"`). It changes only
**where and how often** the constitutive kernel runs and **where its
$C_{\mathrm{tang}}$ output lives**; the physics kernel itself is identical in
both cases.

### 4.1 What is common to both strategies

Regardless of strategy, the Gauss-point **history/state** variables listed
in §3.1 and §3.2 are always allocated as mesh variables and persist across
the Newton loop and across load steps. What differs is whether the
**tangent material tensor** also gets a persistent mesh variable,
`m_C_tang_gp`:

```cpp
// startInit(), FemModule.cc
if (m_gp_material_tensor_strategy == "global")
  m_C_tang_gp.reshape({m_nGP, 3, 3});   // allocated only for the global strategy
```

For the local strategy, `m_C_tang_gp` is never allocated: the tangent tensor
for a given cell is a local (stack) variable computed and consumed inside a
single element-level function call.

### 4.2 Global strategy

Every Newton iteration, `_integrateAndSaveConstitutiveLawVonMises()` or
`_integrateAndSaveConstitutiveLawDruckerPrager()` sweeps all cells and
Gauss points, evaluates the constitutive update, and writes the resulting
stress, history variables and tangent tensor back into `m_sigma_gp` and
`m_C_tang_gp`. `_assembleBilinearOperatorGlobal()` then runs as a fully
separate, law-independent sweep that only reads `m_C_tang_gp` (and, through
`_applyInternalBodyForce()`, `m_sigma_gp`) to build $K_e$ and
$F_{\mathrm{int},e}$: no constitutive evaluation happens during this second
sweep. The constitutive law is therefore evaluated once per cell/Gauss point
per Newton iteration, decoupled from the assembly that consumes its output.

### 4.3 Local strategy

`_assembleBilinearOperatorLocalVonMises()` and
`_assembleBilinearOperatorLocalDruckerPrager()` replace both the update
sweep and the generic assembly of the global strategy with a single,
law-dependent pass. For each cell, the per-element callback (`Tria3`-specific
on the CPU side, since that is the only supported element type) reads the
history/state variables for that cell, evaluates the constitutive update to
get a local tangent tensor, writes the updated stress and history back to
the mesh variables, and immediately assembles $K_e$ from that local tensor,
all within the same function call. `m_C_tang_gp` is never written because it
is never allocated: the tangent tensor exists only for the lifetime of the
call. On the GPU/BSR path the same fusion holds, with an atomic-scattered
variant (`BSR`, one full element matrix per cell) and an atomic-free variant
(`AF-BSR`, one row block per cell/node pair, re-evaluating the constitutive
update per node-thread since the tangent is cheap relative to a shared-cache
scheme). `_assembleLinearOperator()` then reads `m_sigma_gp`, updated as a
side effect of the matrix assembly that precedes it in `_solveNewton()`.

### 4.4 Side-by-side comparison

| Aspect | Global strategy | Local strategy |
|---|---|---|
| `m_C_tang_gp` mesh variable | Allocated, persists across the Newton loop | Never allocated |
| Constitutive evaluations per Newton iteration | One dedicated sweep over all cells/GPs | Inline, once per cell, inside the element-matrix builder |
| Where $C_{\mathrm{tang}}$ lives | Mesh variable, read back during assembly | Local (stack) variable, consumed immediately |
| Stress update vs. matrix assembly | Decoupled: update sweep, then a pure-assembly sweep | Fused: one pass does both |
| Extra memory | $O(n_{\mathrm{cells}}\times n_{GP}\times 3\times3)$ for `m_C_tang_gp` | None beyond the always-present history variables |
| Iteration-0 elastic shortcut | `_setGlobalElasticMaterialTensorAtGPs()` writes `m_C_tang_gp` directly | `elastic_assembly = true` short-circuits the per-cell call |

Both strategies assemble the same discrete bilinear form
$K_e = \int_{\Omega_e} B^\top C_{\mathrm{tang}} B\, d\Omega_e$ and the same
internal-force vector $F_{\mathrm{int},e} = -\int_{\Omega_e} B^\top \sigma\, d\Omega_e$
per Newton iteration, and are numerically equivalent; they differ only in
whether $C_{\mathrm{tang}}$ is materialized mesh-wide or kept local to the
element-matrix computation.

---

## 5. Notes

* Both constitutive kernels avoid explicit branching in favor of
  Macaulay-bracket/switch constructions, which is what allows the same
  kernel to compile and run on both CPU and GPU.
* `_initConstitutiveLaw()` explicitly disables `m_use_gpu_functions` for
  Drucker-Prager, so its return mapping currently runs on the CPU only, even
  though GPU kernels exist in `DruckerPragerLaw.h`.
* The Drucker-Prager apex tangent is set to zero rather than linearized; see
  §3.2.
* 3D and non-`Tria3` elements are not yet supported for either law
  (`ARCANE_FATAL` guards throughout `_initConstitutiveLaw()` and the
  local-assembly dispatch functions). The weak form, Newton driver and the
  global/local dichotomy generalize beyond `Tria3`, but only the 2D linear
  triangle path is currently instantiated.
