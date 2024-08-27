### Volume of a Tetrahedron ###
For a tetrahedron with nodes we assume four  nodes $\mathbf{n}_0 , \mathbf{n}_1, \mathbf{n}_2, \mathbf{n}_3$: 

<img width="100" align="left" src="https://github.com/user-attachments/assets/d4156c9d-1eb1-498d-8146-a678dc3eedd4"/>

 $\mathbf{n}_0 = (x_0, y_0, z_0)$, $\mathbf{n}_1 = (x_1, y_1, z_1)$, $\mathbf{n}_2 = (x_2, y_2, z_2)$, and $\mathbf{n}_3 = (x_3, y_3, z_3)$, 

 let us assume three vectors corresponding to the edges of the tetrahedron that start from vertex $\mathbf{n}_0$

$$\mathbf{v}_0​=\mathbf{n}_1​−\mathbf{n}_0​=(x_1​−x_0​,y_1​−y_0​,z_1​−z_0​)$$

$$\mathbf{v}_1​=\mathbf{n}_2​−\mathbf{n}_0​=(x_2​−x_0​,y_2​−y_0​,z_2​−z_0​)$$

$$\mathbf{v}_2​=\mathbf{n}_3​−\mathbf{n}_0​=(x_3​−x_0​,y_3−y_0​,z_3​−z_0​)$$

Now, we can think of the volume of a tetrahedron is $\frac{1}{6}$ of the volume of the parallelipiped formed by the vectors $\mathbf{v}_0, \mathbf{v}_1, \mathbf{v}_2$. 

Thus the volume $V$ of a tetrahedron is :

$$ V = \frac{1}{6} \| (\mathbf{v}_2 \times \mathbf{v}_1) \cdot \mathbf{v}_0 \|$$

### Implementation in Code ###

```cpp
    static inline Real computeVolumeTetra4(Cell cell, const VariableNodeReal3& node_coord)
    {
      Real3 vertex0 = node_coord[cell.nodeId(0)];
      Real3 vertex1 = node_coord[cell.nodeId(1)];
      Real3 vertex2 = node_coord[cell.nodeId(2)];
      Real3 vertex3 = node_coord[cell.nodeId(3)];

      Real3 v0 = vertex1 - vertex0;
      Real3 v1 = vertex2 - vertex0;
      Real3 v2 = vertex3 - vertex0;

      return std::abs(Arcane::math::dot(v0, Arcane::math::cross(v1, v2))) / 6.0;
    }
```
