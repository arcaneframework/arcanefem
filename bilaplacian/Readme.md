# Notes on the solver #

The code here is a simple FEM code used to solve a bilaplacian problem on unstructured 2D mesh. 



## Theory ##

#### Problem description ####

The steady state 2D bilaplacian equation is solved for a closed square meshed domain $\Omega^h = (0,1)^2$ in order to know the vectorial unknowns $\{u_i(x,y)\}_{i=1}^2$, i.e, $u_1$ and $u_2$ within the domain. The system of equations read

$$\triangle u_1 + u_2  = 0  \quad \forall (x,y)\in\Omega^h $$

$$\triangle u_2  = f  \quad \forall (x,y)\in\Omega^h $$

There are no Neumann boundrary conditions attached to this problem, however to complete the problem description,  four first type (Dirichlet) boundary conditions are applied to this problem:

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Top}}\subset\partial \Omega^h,$

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Left}}\subset\partial \Omega^h,$

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Right}}\subset\partial \Omega^h,$

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Bot}}\subset\partial \Omega^h,$



#### Variational formulation



In this case  the variational formulation we use the subset of square integrable Sobolev functional space   $H^1_{0}(\Omega) \subset H^1{\Omega}$. The FEM formulation then reads:

search vectorial FEM trial function $(u^h_1,u^h_2)\in\left[H^1_0(\Omega^h)\right]^2$ satisfying

$$ \int_{\Omega^h}\nabla u^h_1 \nabla  v_2^h +  \int_{\Omega^h}\nabla u^h_2 \nabla  v_1^h + \int_{\Omega^h} u^h_2   v_2^h + \int_{\Omega^h}f v_1^h = 0 \quad \forall (v_1^h,v_2^h)\in \left[H^1_0(\Omega^h)\right]^2$$

given:

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Top}}\subset\partial \Omega^h,$

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Left}}\subset\partial \Omega^h,$

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Right}}\subset\partial \Omega^h,$

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Bot}}\subset\partial \Omega^h,$

$\int_{\Omega^h}f v_1^h=1\times10^5$



## The code ##



#### Post Process ####

For post processing the `ensight.case` file is outputted, which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).
