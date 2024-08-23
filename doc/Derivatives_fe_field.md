## Calculating Derivatives $dx(u)$ and $dy(u)$ Over a Triangular Finite Element


<img width="200" align="left" src="https://github.com/user-attachments/assets/124852e5-3169-4436-8fa9-b90bddfc9972"/>

To calculate the derivatives of a function $u$ with respect to $x$ and $y$ ( denoted as $\frac{\partial u}{\partial x}$, $\frac{\partial u}{\partial y}$, or more simply $dx(u)$, $dy(u)$ ) over a triangular finite element, we follow these steps:

### 1. Define the Triangle and the P1 Basis Functions

In a P1 (linear) finite element, the function $u$ is assumed to be linear over each triangle. 
Let the nodes of the triangle be $\mathbf{n}_0 = (x_0, y_0)$, $\mathbf{n}_1 = (x_1, y_1)$, and $\mathbf{n}_2 = (x_2, y_2)$,
with corresponding function values $u_0 = u(\mathbf{n}_0)$, $u_1 = u(\mathbf{n}_1)$, and $u_2 = u(\mathbf{n}_2)$.

The function $u(x, y)$ can be expressed as:

$$u(x, y) = a_0 + a_1 x + a_2 y$$

where $a_0$, $a_1$, and $a_2$ are coefficients to be determined.

### 2. Set Up the System of Equations

We solve for the coefficients $a_0$, $a_1$, and $a_2$ using the known values of $u$ at the triangle's vertices:

$$
\begin{aligned}
u_0 &= a_0 + a_1 x_0 + a_2 y_0  \\
u_1 &= a_0 + a_1 x_1 + a_2 y_1  \\
u_2 &= a_0 + a_1 x_2 + a_2 y_2 
\end{aligned}
$$

This system can be written in matrix form:

$$
\begin{pmatrix}
1 & x_0 & y_0  \\
1 & x_1 & y_1  \\
1 & x_2 & y_2 
\end{pmatrix}
\begin{pmatrix}
a_0 \\
a_1 \\
a_2 
\end{pmatrix}=
\begin{pmatrix}
u_0 \\
u_1 \\
u_2
\end{pmatrix}
$$


Solve this system to get $a_0$, $a_1$, and $a_2$.

### 3. Calculate the Derivatives $dx(u)$ and $dy(u)$

- The derivative $\frac{\partial u}{\partial x}$ is simply $a_1$. This is because the partial derivative of $u(x, y)$ with respect to $x$ is:

$$
\frac{\partial u}{\partial x} = a_1
$$

Therefore, $dx(u) = a_1$.

- The derivative $\frac{\partial u}{\partial y}$ is simply $a_2$. This is because the partial derivative of $u(x, y)$ with respect to $y$ is:

$$
\frac{\partial u}{\partial y} = a_2
$$

Therefore, $dy(u) = a_2$.

### 4. Shortcut Using the Area Formula (2D case)

For a 2D triangular element (i.e., $z = 0$), the coefficients $a_1$ and $a_2$ (which give the derivatives with respect to $x$ and $y$, respectively) can be directly computed using the area of the triangle and the coordinates of its vertices:

$$
a_1 = \frac{1}{2A} \left( u_0(y_1 - y_2) + u_1(y_2 - y_0) + u_2(y_0 - y_1) \right)
$$

$$
a_2 = \frac{1}{2A} \left( u_0(x_2 - x_1) + u_1(x_0 - x_2) + u_2(x_1 - x_0) \right)
$$

where $A$ is the area of the triangle, given by:

$$
A = \frac{1}{2} \left| x_0(y_1 - y_2) + x_1(y_2 - y_0) + x_2(y_0 - y_1) \right|
$$

This method directly gives the derivatives $dx(u)$ and $dy(u)$ without solving a system of equations.

## Implementation in Code: ##
```cpp
    static inline Real3 computeGradientTria3(Cell cell, const VariableNodeReal3& node_coord, const VariableNodeReal& u)
    {
      Real3 n0 = node_coord[cell.nodeId(0)];
      Real3 n1 = node_coord[cell.nodeId(1)];
      Real3 n2 = node_coord[cell.nodeId(2)];

      Real u0 = u[cell.nodeId(0)];
      Real u1 = u[cell.nodeId(1)];
      Real u2 = u[cell.nodeId(2)];

      Real A2 = ((n1.x - n0.x) * (n2.y - n0.y) - (n2.x - n0.x) * (n1.y - n0.y));

      return Real3 ( (u0*(n1.y - n2.y) + u1*(n2.y - n0.y) + u2*(n0.y - n1.y)) / A2 , (u0*(n2.x - n1.x) + u1*(n0.x - n2.x) + u2*(n1.x - n0.x)) / A2 , 0);
    }
```
