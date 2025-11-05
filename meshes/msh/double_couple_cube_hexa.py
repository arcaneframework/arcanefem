'''
     -------------------------------------------------------------------
         This file is a part of ArcaneFEM (finite element tool in Arcane)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.fr
     Date     : 04-11-2025

     -------------------------------------------------------------------

     This file is distributed  in  the hope that it will be useful,
     but WITHOUT ANY WARRANTY; or without even the implied warranty
     of  FITNESS  FOR  A  PARTICULAR  PURPOSE.

     --------------------------------------------------------------------

     This file which produces a 3D cube mesh for testing soildynamics

     compile-run: python fyle.py
'''
import gmsh
import sys
import numpy as np

gmsh.initialize(sys.argv)
gmsh.model.add("cube_double_couple_3d")

# Cube dimensions
cube_size = 2.0
NN = 7  # Number of nodes per edge

# Calculate grid spacing
dcLen = cube_size / (NN - 1)

# Cube center
xC = cube_size / 2.0
yC = cube_size / 2.0
zC = cube_size / 2.0

# Create cube vertices
p1 = gmsh.model.geo.addPoint(0, 0, 0)
p2 = gmsh.model.geo.addPoint(cube_size, 0, 0)
p3 = gmsh.model.geo.addPoint(cube_size, cube_size, 0)
p4 = gmsh.model.geo.addPoint(0, cube_size, 0)
p5 = gmsh.model.geo.addPoint(0, 0, cube_size)
p6 = gmsh.model.geo.addPoint(cube_size, 0, cube_size)
p7 = gmsh.model.geo.addPoint(cube_size, cube_size, cube_size)
p8 = gmsh.model.geo.addPoint(0, cube_size, cube_size)

# Create lines for bottom face
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)

# Create lines for top face
l5 = gmsh.model.geo.addLine(p5, p6)
l6 = gmsh.model.geo.addLine(p6, p7)
l7 = gmsh.model.geo.addLine(p7, p8)
l8 = gmsh.model.geo.addLine(p8, p5)

# Create vertical lines
l9 = gmsh.model.geo.addLine(p1, p5)
l10 = gmsh.model.geo.addLine(p2, p6)
l11 = gmsh.model.geo.addLine(p3, p7)
l12 = gmsh.model.geo.addLine(p4, p8)

# Create curve loops for each face
cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])  # Bottom
cl2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])  # Top
cl3 = gmsh.model.geo.addCurveLoop([l1, l10, -l5, -l9])  # Front
cl4 = gmsh.model.geo.addCurveLoop([l2, l11, -l6, -l10])  # Right
cl5 = gmsh.model.geo.addCurveLoop([l3, l12, -l7, -l11])  # Back
cl6 = gmsh.model.geo.addCurveLoop([l4, l9, -l8, -l12])  # Left

# Create surfaces
s1 = gmsh.model.geo.addPlaneSurface([cl1])  # Bottom
s2 = gmsh.model.geo.addPlaneSurface([cl2])  # Top
s3 = gmsh.model.geo.addPlaneSurface([cl3])  # Front
s4 = gmsh.model.geo.addPlaneSurface([cl4])  # Right
s5 = gmsh.model.geo.addPlaneSurface([cl5])  # Back
s6 = gmsh.model.geo.addPlaneSurface([cl6])  # Left

# Create surface loop and volume
sl1 = gmsh.model.geo.addSurfaceLoop([s1, s2, s3, s4, s5, s6])
v1 = gmsh.model.geo.addVolume([sl1])

# Synchronize geometry
gmsh.model.geo.synchronize()

# Create physical groups for surfaces
gmsh.model.addPhysicalGroup(2, [s1, s2, s3, s4, s5, s6], 1)
gmsh.model.setPhysicalName(2, 1, "boundaries")

# Create physical group for 3D volume
gmsh.model.addPhysicalGroup(3, [v1], 1)
gmsh.model.setPhysicalName(3, 1, "volume")

# Set transfinite mesh for structured hexahedral elements
for curve in gmsh.model.getEntities(1):
    gmsh.model.mesh.setTransfiniteCurve(curve[1], NN)

for surface in gmsh.model.getEntities(2):
    gmsh.model.mesh.setTransfiniteSurface(surface[1])
    gmsh.model.mesh.setRecombine(surface[0], surface[1])

gmsh.model.mesh.setTransfiniteVolume(v1)

# Generate 3D mesh
print("Generating mesh...")
gmsh.model.mesh.generate(3)

# Procedure for physical points at target 
# Define target coordinates
targets = {
    'center': np.array([xC, yC, zC]),
    'north': np.array([xC, yC, zC + dcLen]),
    'south': np.array([xC, yC, zC - dcLen]),
    'east': np.array([xC + dcLen, yC, zC]),
    'west': np.array([xC - dcLen, yC, zC])
}

tol = 1e-8

# Get all mesh nodes
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
node_coords = node_coords.reshape(-1, 3)

# Find and tag each target point
entity_tag_counter = 1001
for name, target in targets.items():
    found_tag = None
    for tag, xyz in zip(node_tags, node_coords):
        if np.allclose(xyz, target, atol=tol):
            found_tag = int(tag)
            print(f"Found node {found_tag} at {xyz} for '{name}'")
            break

    if found_tag is None:
        print(f"Warning: No node found near {target} for '{name}'")
        continue

    # Create discrete entity with deterministic tag
    entity_dim = 0
    entity_tag = entity_tag_counter
    gmsh.model.addDiscreteEntity(entity_dim, entity_tag)
    gmsh.model.mesh.addNodes(entity_dim, entity_tag, [found_tag], target.tolist())

    # Create physical group with matching tag
    phys_tag = gmsh.model.addPhysicalGroup(entity_dim, [entity_tag], entity_tag)
    gmsh.model.setPhysicalName(entity_dim, phys_tag, name)

    entity_tag_counter += 1

# Set mesh format options
gmsh.option.setNumber("Mesh.SaveAll", 1)
gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
gmsh.option.setNumber("Mesh.Binary", 1)

# Write final mesh
gmsh.write("cube_double_couple_3d.hexa.msh")
print("Mesh saved to: cube_double_couple_3d.hexa.msh")
print("Physical points added successfully.")

gmsh.finalize()
