### How are Dirichlet Boundary conditions imposed

To enforce Dirichlet boundary conditions in FEM many methods are available, some of the common ones include: 

|         **Method**         | **Reference**  | In ArcaneFEM |
| :------------------------: | :------------: | :----------: |
|    Weak penalty method     | Babuška, 1973a |     YES      |
|       Penalty method       | Babuška, 1973a |     YES      |
|      Row elimination       |                |     YES      |
|   Row/Column elimination   |                |     YES      |
| Lagrange multiplier method | Babuška, 1973b |      NO      |
|      Nitsche’s method      | Nitsche, 1971  |      NO      |

The first four are made available, so that one can choose the most appropriate method for  Dirichlet boundary condition implementation. The word 'appropriate'  in the preceding sentence is tricky, it depends on the physics, size of the problem, conditioning of the problem, used linear-solver, etc. So choose wisely. 

What remains common in these methods is these are applied after one assembles the FEM linear system. As such one could see these methods as a post-processing step to be applied to the assembled linear-systems  $\mathbf{Ax=b}$ before they are moved to the solving stage. More precisely, these methods alter the vector $\mathbf{b}$ and/or matrix $\mathbf{A}$.  

#### Weak penalty method and Penalty method ####

These methods are not exact rather weak sense of applying Dirichlet boundary condition, however these are the most straightforward way to apply the  Dirichlet boundary condition. These methods remain popular among practitioners and can be found in FEM packages such as MOOSE, FreeFEM, etc. Main benefits of these methods include, simplistic implementation,  conserved non-zero structure (sparsity) of  $\mathbf{A}$ as such matrix symmetry may be conserved, and avoiding insertion of zeros into $\mathbf{A}$ which is a tedious operation. 

In a nutshell,  to apply weak penalty method for a DOF $i$ which corresponds to a boundary node in the mesh the following needs to be changed: 

$a_{i,i} = a_{i,i} + \mathcal{P}$     and   $b_i = g_i \times \mathcal{P}$

here, $\mathcal{P}$  is the penalty term and $g_i$ is the given Dirichlet value. 

Similarly to apply penalty method for a DOF $i$ which corresponds to a boundary node in the mesh the following needs to be changed:

$a_{i,i} = \mathcal{P}$     and   $b_i = g_i \times \mathcal{P}$

here, $\mathcal{P}$  is the penalty term and $g_i$ is the given Dirichlet value.

The logic is simple as $\mathcal{P}$ is a large large value $\mathcal{P}>>1$, the diagonal contribution $a_{i,i}$ dominates the off-diagonal ones  (since $a_{i,j} << a_{i,i}$), this indeed means solution $b_i \approx g_i$. 

Generally, $\mathcal{P}$ is chosen to be large enough e.g, $1\times10^{31}$ so that  the resulting solution $\mathbf{x=A^{-1}b}$ is precise enough for a double precision calculation. This works most often, however sometimes this may result into $\mathbf{A}$ becoming  ill-conditioned, eventually  this might lead in failed solves and overall accuracy losses. So the user should select and tune the value of $\mathcal{P}$ in this case.

#### Row  elimination ####

To apply a Dirichlet condition $g_i$ on the Dirchlet DOF $i$ via row  elimination method two operations are performed. 

- First is on the matrix $\mathbf{A}$ we  zero out the $i$ th row and $i$ th column of the matrix $\mathbf{A}$  and impose $1$  at the diagonal $a_{i,i}=1$. For a FEM problem, if $N_{D}$ is the total number of Dirichlet DOFs and  if $\mathbf{A}$  has the size of  $N_{DOF}$, this method involves altering $\mathbf{A}$: 

$$\forall i=1, \ldots , N_D \quad ; \quad \forall j = 1 , \ldots , N_{DOF}$$

$$ a_{i,j}=0, \quad  a_{j,i}=0,  \quad a_{i,i} = 1$$

- Second operation is on the vector $\mathbf{b}$,  for  each Dirichlet condition to be imposed on the $i$ th DOF, we substitute the Dirichlet value $g_i$ to $i$ th component of the the vector $\mathbf{b}$.  Hence, this method involves:

$$\forall i=1, \ldots , N_D \quad \quad b_i = g_i$$

This method is more exact way of imposing Dirichlet boundary conditions , however the matrix $\mathbf{A}$ is not symmetric anymore (if it were before), the problem it is algorithmically more challenging  to implement as $\mathbf{A}$ in comparison to penalty methods.

#### Row column elimination ####

To apply a Dirichlet condition $g_i$ on the Dirchlet DOF $i$ via row column elimination method two operations are performed. 

- First is on the matrix $\mathbf{A}$ we  zero out the $i$ th row of the matrix $\mathbf{A}$  and impose $1$  at the diagonal $a_{i,i}=1$. For a FEM problem, if $N_{D}$ is the total number of Dirichlet DOFs and  if $\mathbf{A}$  has the size of  $N_{DOF}$, this method involves altering $\mathbf{A}$:

$$\forall i=1,\ldots ,N_D \quad ; \quad \forall j = 1 ,\ldots, N_{DOF}$$
  
$$a_{i,j}=0,  \quad a_{i,i} = 1$$

- Second operation is on the vector $\mathbf{b}$,  for  each Dirichlet condition to be imposed on the $i$ th DOF, we subtract the $i$ th column from the vector $\mathbf{b}$.  Hence, this method involves:

$$\forall i=1,\ldots ,N_D \quad ; \quad \forall j = 1 ,\ldots, N_{DOF}$$

$$b_j = b_j - a_{j,i} , \quad b_i = g_i$$

This method is more exact way of imposing Dirichlet boundary conditions, however it is algorithmically more challenging  to implement as $\mathbf{A}$ which is a sparse matrix is stored in compressed row format (CRS) then moving the  columns to the right-hand side becomes expensive for large systems with many Dirichlet DOFs.

## Referneces ##

- Babuška, I., 1973. The finite element method with penalty. *Mathematics of computation*, *27*(122), pp.221-228.

- Babuška, I., 1973. The finite element method with Lagrangian multipliers. *Numerische Mathematik*, *20*(3), pp.179-192.

- Nitsche, J., 1971, July. Über ein Variationsprinzip zur Lösung von  Dirichlet-Problemen bei Verwendung von Teilräumen, die keinen  Randbedingungen unterworfen sind. In *Abhandlungen aus dem mathematischen Seminar der Universität Hamburg* (Vol. 36, No. 1, pp. 9-15). Berlin/Heidelberg: Springer-Verlag.





