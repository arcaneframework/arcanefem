# linear elasticity

Here we deal with linear solid-mechanics governed by a system of PDE modeling the deformation of elastic bodies. The solver, here is a 2D unstructured mesh linear elasticity solver, which uses FEM to search for vector solution of displacement unknowns.



## Mathematics ##

#### Problem description ####

Under small elastic deformation, the steady  2D elastic deformation equation (linear elastic system) on  domain $\Omega$ reads



## The code ##



#### Post Process ####

For post processing the `ensight.case` file is outputted, which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).
