'''
     -------------------------------------------------------------------
         This file is a part of ArcaneFEM (finite element tool in Arcane)
     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.fr
     Date     : 07-11-2025

     -------------------------------------------------------------------

     This file is distributed  in  the hope that it will be useful,
     but WITHOUT ANY WARRANTY; or without even the implied warranty
     of  FITNESS  FOR  A  PARTICULAR  PURPOSE.

     --------------------------------------------------------------------

     This code produces a 3D french alps hexa mesh for testing soildynamics

     compile-run: python file.py
'''
import gmsh
import math
import os
import sys
import numpy as np

gmsh.initialize()

# --------------------------------------------------------------
# paramerters to control the mesh
# --------------------------------------------------------------

# depth of the top surface
depth_z = -1.0

# Numeber of nodes along each edge
NN = 50

# Top surface mesh file name
top_surface_mesh_file = 'top_fench_alps.msh'

# Define double-couple target points
targets = {
    'center': np.array([0.09153777977352227,    -0.2160065059623689,	-0.9001811733115577]),
    'north': np.array([0.1634282588162065,	    -0.2154647300608412,	-0.8986523346852678]),
    'south': np.array([0.02033057068205987,	    -0.2166168648590111,	-0.9042853291838218]),
    'east': np.array([0.0915768251709351,	    -0.2873443335212027,	-0.9004511926413347]),
    'west': np.array([0.09233912966822048,	    -0.1452554459451716,	-0.9025199403673613])
}

# tolerance for finding double-couple target points
tol = 1e-8

# --------------------------------------------------------------


path = os.path.dirname(os.path.abspath(__file__))

# load top surface
gmsh.merge(os.path.join(path, top_surface_mesh_file))

# classify the surface mesh according to given angle, and create discrete model
# entities (surfaces, curves and points) accordingly; curveAngle forces bounding
# curves to be split on sharp corners
gmsh.model.mesh.classifySurfaces(math.pi, curveAngle = math.pi / 10)

# create a geometry for the discrete curves and surfaces
gmsh.model.mesh.createGeometry()

# retrieve the surface, its boundary curves and corner points
s = gmsh.model.getEntities(2)
c = gmsh.model.getBoundary(s)

if (len(c) != 4):
    gmsh.logger.write('Should have 4 boundary curves!', level='error')

# get corner points and their coordinates
p = []
xyz = []
for e in c:
    pt = gmsh.model.getBoundary([e], combined=False)
    p.extend([pt[0][1]])
    xyz.extend(gmsh.model.getValue(0, pt[0][1], []))

# create other CAD entities to form one volume below the terrain surface; beware
# that only built-in CAD entities can be hybrid, i.e. have discrete entities on
# their boundary
p1 = gmsh.model.geo.addPoint(xyz[0], xyz[1], depth_z)
p2 = gmsh.model.geo.addPoint(xyz[3], xyz[4], depth_z)
p3 = gmsh.model.geo.addPoint(xyz[6], xyz[7], depth_z)
p4 = gmsh.model.geo.addPoint(xyz[9], xyz[10], depth_z)

c1 = gmsh.model.geo.addLine(p1, p2)
c2 = gmsh.model.geo.addLine(p2, p3)
c3 = gmsh.model.geo.addLine(p3, p4)
c4 = gmsh.model.geo.addLine(p4, p1)

c10 = gmsh.model.geo.addLine(p1, p[0])
c11 = gmsh.model.geo.addLine(p2, p[1])
c12 = gmsh.model.geo.addLine(p3, p[2])
c13 = gmsh.model.geo.addLine(p4, p[3])

ll1 = gmsh.model.geo.addCurveLoop([c1, c2, c3, c4])
s1 = gmsh.model.geo.addPlaneSurface([ll1])
ll3 = gmsh.model.geo.addCurveLoop([c1, c11, -c[0][1], -c10])
s3 = gmsh.model.geo.addPlaneSurface([ll3])
ll4 = gmsh.model.geo.addCurveLoop([c2, c12, -c[1][1], -c11])
s4 = gmsh.model.geo.addPlaneSurface([ll4])
ll5 = gmsh.model.geo.addCurveLoop([c3, c13, -c[2][1], -c12])
s5 = gmsh.model.geo.addPlaneSurface([ll5])
ll6 = gmsh.model.geo.addCurveLoop([c4, c10, -c[3][1], -c13])
s6 = gmsh.model.geo.addPlaneSurface([ll6])
sl1 = gmsh.model.geo.addSurfaceLoop([s1, s3, s4, s5, s6, s[0][1]])

# create the volume
v1 = gmsh.model.geo.addVolume([sl1])
tag = 2
gmsh.model.geo.addPhysicalGroup(3, [v1], tag)
gmsh.model.setPhysicalName(3, tag, "volume")

# Paraxial
tag = 1
gmsh.model.geo.addPhysicalGroup(2, [s1, s3, s4, s5, s6], tag)
gmsh.model.setPhysicalName(2, tag, "paraxial")

gmsh.model.geo.synchronize()

# Transfinite meshing
for c in gmsh.model.getEntities(1):
    gmsh.model.mesh.setTransfiniteCurve(c[1], NN)
for s in gmsh.model.getEntities(2):
    gmsh.model.mesh.setTransfiniteSurface(s[1])
    gmsh.model.mesh.setRecombine(s[0], s[1])
    gmsh.model.mesh.setSmoothing(s[0], s[1], 100)
gmsh.model.mesh.setTransfiniteVolume(v1)

# Generate a 3D mesh
gmsh.model.mesh.generate(3)

# --- get all nodes ---
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
node_coords = node_coords.reshape(-1, 3)

# --- find and tag each target point ---
entity_tag_counter = 1001
dc_point_not_found = False
for name, target in targets.items():
    found_tag = None
    for tag, xyz in zip(node_tags, node_coords):
        if np.allclose(xyz, target, atol=tol):
            found_tag = int(tag)
            print(f"Found node {found_tag} at {xyz} for '{name}'")
            break

    if found_tag is None:
        print(f"Warning: No node found near {target} for '{name}'")
        dc_point_not_found = True
        continue

    # create discrete entity with deterministic tag
    entity_dim = 0
    entity_tag = entity_tag_counter
    gmsh.model.addDiscreteEntity(entity_dim, entity_tag)  # create entity with tag 1001, 1002, ...
    gmsh.model.mesh.addNodes(entity_dim, entity_tag, [found_tag], target.tolist())

    # create a physical group with the same tag (so PhysicalNames will show 1001 etc.)
    phys_tag = gmsh.model.addPhysicalGroup(entity_dim, [entity_tag], entity_tag)  # pass desired phys tag
    gmsh.model.setPhysicalName(entity_dim, phys_tag, name)                        # give it a name

    entity_tag_counter += 1

if dc_point_not_found:
    print("\n\n********* ERROR ERROR ERROR ERROR ERROR ERROR *********\n")
    print("Some double-couple target points were NOT found in the mesh!")
    print("See Warning above. Adjust 'tol' if necessary.")
    print("\n********* EXITING WITHOUT WRITING MESH *********\n")
    exit(1)
else:
    print("Physical Points for couble couple added automatically.")

gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
gmsh.option.setNumber("Mesh.Binary", 1)

gmsh.write("french_alps_dc.hexa.msh")
gmsh.finalize()
