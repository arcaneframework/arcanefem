# Adding a New Constitutive Law in ArcaneFEM

Setting up a new law takes two things: a `[MyLaw].cc` (and declarations in `[module].h`) file containing the law functions, and a few additions in the `module.cc` file that contains the Newton loop. 

Here, `[MyLaw]` is a placeholder to replace with the name of your law (e.g. `VonMises`, `DruckerPrager`).

---

## 1. Write `[MyLaw]Law.cc` with four functions

Example: `VonMisesLaw.cc`.

| Function | Role |
|---|---|
| `computeMaterialTensor[MyLaw]LawAtGpBase(...)` | Point-wise update at one Gauss point (`ARCCORE_HOST_DEVICE`): from the strain increment and the committed state, returns the new stress, the history increment and the consistent tangent. |
| `_integrateAndSaveConstitutiveLaw[MyLaw]()` | Called at every Newton iteration. Dispatches on dimension, element type and CPU/GPU, loops over cells and Gauss points, computes the strain increment, calls the point-wise function and stores the results in the global variables (stress, tangent material tensor, internal variables, history increment). |
| `_restoreConvergedState[MyLaw]()` | Called at the start of each time increment: copies the committed (`_old`) values into the current ones. |
| `_commitInternalVariables[MyLaw]()` | Called once the Newton solver has converged: updates the `_old` and history variables. |

---

## 2. Modify the Newton module file

**a) Initialize the variables, name and shape.** In `initConstitutiveLaw()`, add a branch for your law that declares and initializes every history and internal variable stored globally, with the shape `cells × number of Gauss points (× components)`. This includes the current and `_old` versions of the stress and history variables, and the tangent tensor.

**b) Register the law and its parameters.** Add the law name and its material parameters (options in the `.axl` file that are parsed in `initConstitutiveLaw()` and the derived parameters are then evaluated in `_getMaterialParameters()`.

**c) Add your functions in the Newton loop**, at three places:

```cpp
// 1. start of the increment
else if (m_constitutive_law == "[MyLaw]") _restoreConvergedState[MyLaw]();

// 2. inside the loop, before assembling the Jacobian and the residual
else if (m_constitutive_law == "[MyLaw]") _integrateAndSaveConstitutiveLaw[MyLaw]();

// 3. after Newton convergence
else if (m_constitutive_law == "[MyLaw]") _commitInternalVariables[MyLaw]();
```

The rest of the loop (solve, assembly, convergence check) does not change, since it only reads the stress tensor (as vector) and the tangent material tensor (as matrix) stored by your law.

> In the Newton loop as pasted, the `VonMises` commit block is missing its closing brace and there is no commit branch for `DruckerPrager`. Check the braces when adding your branch.
