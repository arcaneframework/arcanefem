### How do we calculate gradients of shape functions for P1 Tetrahedron Finite Element ###

For a tetrahedron with nodes we assume four  nodes $m_0 , m_1. m_2, m_3$: 

<img width="100" align="left" src="https://github.com/user-attachments/assets/d4156c9d-1eb1-498d-8146-a678dc3eedd4"/>

 $m_0 = (x_0, y_0, z_0)$, $m_1 = (x_1, y_1, z_1)$, $m_2 = (x_2, y_2, z_2)$, and $m_3 = (x_3, y_3, z_3)$, 

The P1 finite element shape functions $\Phi_i$ for a tetrahedron can be expressed as:

$$
\Phi_i(x, y, z) = \frac{1}{6V} \left( a_i + b_i x + c_i y + d_i z \right)
$$

Where $V$ is the volume of the tetrahedron.

The gradients of the shape functions $\Phi_0$, $\Phi_1$, $\Phi_2$, $\Phi_3$ with respect to $x$, $y$, and $z$ can be computed as follows:

## Gradient with Respect to $x$ ($\frac{\partial \Phi_i}{\partial x}$)

The partial derivative with respect to $x$ is:

$$
\frac{\partial \Phi_i}{\partial x} = \frac{b_i}{6V}
$$

Here, the coefficients $b_i$ are:

$$
b_0 = \det \begin{pmatrix}
y_1 & z_1 \\
y_2 & z_2 \\
y_3 & z_3
\end{pmatrix}
$$

$$
b_1 = \det \begin{pmatrix}
y_0 & z_0 \\
y_2 & z_2 \\
y_3 & z_3
\end{pmatrix}
$$

$$
b_2 = \det \begin{pmatrix}
y_0 & z_0 \\
y_1 & z_1 \\
y_3 & z_3
\end{pmatrix}
$$

$$
b_3 = \det \begin{pmatrix}
y_0 & z_0 \\
y_1 & z_1 \\
y_2 & z_2
\end{pmatrix}
$$

## Gradient with Respect to $y$ ($\frac{\partial \Phi_i}{\partial y}$)

The partial derivative with respect to $y$ is:

$$
\frac{\partial \Phi_i}{\partial y} = \frac{c_i}{6V}
$$

Where the coefficients $c_i$ are:

$$
c_0 = \det \begin{pmatrix}
z_1 & x_1 \\
z_2 & x_2 \\
z_3 & x_3
\end{pmatrix}
$$

$$
c_1 = \det \begin{pmatrix}
z_0 & x_0 \\
z_2 & x_2 \\
z_3 & x_3
\end{pmatrix}
$$

$$
c_2 = \det \begin{pmatrix}
z_0 & x_0 \\
z_1 & x_1 \\
z_3 & x_3
\end{pmatrix}
$$

$$
c_3 = \det \begin{pmatrix}
z_0 & x_0 \\
z_1 & x_1 \\
z_2 & x_2
\end{pmatrix}
$$

## Gradient with Respect to $z$ ($\frac{\partial \Phi_i}{\partial z}$)

The partial derivative with respect to $z$ is:

$$
\frac{\partial \Phi_i}{\partial z} = \frac{d_i}{6V}
$$

Where the coefficients $d_i$ are:

$$
d_0 = \det \begin{pmatrix}
x_1 & y_1 \\
x_2 & y_2 \\
x_3 & y_3
\end{pmatrix}
$$

$$
d_1 = \det \begin{pmatrix}
x_0 & y_0 \\
x_2 & y_2 \\
x_3 & y_3
\end{pmatrix}
$$

$$
d_2 = \det \begin{pmatrix}
x_0 & y_0 \\
x_1 & y_1 \\
x_3 & y_3
\end{pmatrix}
$$

$$
d_3 = \det \begin{pmatrix}
x_0 & y_0 \\
x_1 & y_1 \\
x_2 & y_2
\end{pmatrix}
$$

## Implementation in Code for Gradient in $x$ ($\frac{\partial \Phi_i}{\partial x}$):

```cpp
    static inline  Real4 computeGradientXTetra4(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];
      Real3 vertex3 = node_coord[cell.nodeId(3)];

      Real3 v0 = vertex1 - vertex0;
      Real3 v1 = vertex2 - vertex0;
      Real3 v2 = vertex3 - vertex0;

      // 6 x Volume of tetrahedron
      Real V6 =  std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2)));

      Real4 dx;

      dx[0] = (vertex1.y * (vertex3.z - vertex2.z) + vertex2.y * (vertex1.z - vertex3.z) + vertex3.y * (vertex2.z - vertex1.z))/V6;
      dx[1] = (vertex0.y * (vertex2.z - vertex3.z) + vertex2.y * (vertex3.z - vertex0.z) + vertex3.y * (vertex0.z - vertex2.z))/V6;
      dx[2] = (vertex0.y * (vertex3.z - vertex1.z) + vertex1.y * (vertex0.z - vertex3.z) + vertex3.y * (vertex1.z - vertex0.z))/V6;
      dx[3] = (vertex0.y * (vertex1.z - vertex2.z) + vertex1.y * (vertex2.z - vertex0.z) + vertex2.y * (vertex0.z - vertex1.z))/V6; 

      return dx;
    }
```
