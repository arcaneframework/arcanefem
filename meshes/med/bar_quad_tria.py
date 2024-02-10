#!/usr/bin/env python

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()

###
### SHAPER component
###

from salome.shaper import model

model.begin()
partSet = model.moduleDocument()

### Create Part
Part_1 = model.addPart(partSet)
Part_1_doc = Part_1.document()
model.addParameter(Part_1_doc, "l", '0.5')

### Create Sketch
Sketch_1 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_1 = Sketch_1.addLine(0.5, 0, 0, 0)

### Create SketchProjection
SketchProjection_1 = Sketch_1.addProjection(model.selection("VERTEX", "PartSet/Origin"), False)
SketchPoint_1 = SketchProjection_1.createdFeature()
Sketch_1.setCoincident(SketchLine_1.endPoint(), SketchPoint_1.result())

### Create SketchLine
SketchLine_2 = Sketch_1.addLine(0, 0, 0, 0.1)

### Create SketchLine
SketchLine_3 = Sketch_1.addLine(0, 0.1, 0.5, 0.1)

### Create SketchLine
SketchLine_4 = Sketch_1.addLine(0.5, 0.1, 0.5, 0)
Sketch_1.setCoincident(SketchLine_4.endPoint(), SketchLine_1.startPoint())
Sketch_1.setCoincident(SketchLine_1.endPoint(), SketchLine_2.startPoint())
Sketch_1.setCoincident(SketchLine_2.endPoint(), SketchLine_3.startPoint())
Sketch_1.setCoincident(SketchLine_3.endPoint(), SketchLine_4.startPoint())
Sketch_1.setPerpendicular(SketchLine_1.result(), SketchLine_2.result())
Sketch_1.setPerpendicular(SketchLine_2.result(), SketchLine_3.result())
Sketch_1.setPerpendicular(SketchLine_3.result(), SketchLine_4.result())
Sketch_1.setLength(SketchLine_1.result(), "l")
Sketch_1.setLength(SketchLine_2.result(), "2*l/10")

### Create SketchLine
SketchLine_5 = Sketch_1.addLine(0.5, 0.1, 1, 0.1)
Sketch_1.setCoincident(SketchLine_3.endPoint(), SketchLine_5.startPoint())
Sketch_1.setHorizontal(SketchLine_5.result())

### Create SketchLine
SketchLine_6 = Sketch_1.addLine(1, 0.1, 1, 0)
Sketch_1.setCoincident(SketchLine_5.endPoint(), SketchLine_6.startPoint())

### Create SketchProjection
SketchProjection_2 = Sketch_1.addProjection(model.selection("EDGE", "PartSet/OX"), False)
SketchLine_7 = SketchProjection_2.createdFeature()
Sketch_1.setCoincident(SketchLine_6.endPoint(), SketchLine_7.result())
Sketch_1.setVertical(SketchLine_6.result())

### Create SketchLine
SketchLine_8 = Sketch_1.addLine(1, 0, 0.5, 0)
Sketch_1.setCoincident(SketchLine_6.endPoint(), SketchLine_8.startPoint())
Sketch_1.setCoincident(SketchLine_1.startPoint(), SketchLine_8.endPoint())
Sketch_1.setHorizontal(SketchLine_8.result())
Sketch_1.setEqual(SketchLine_1.result(), SketchLine_8.result())
model.do()

### Create Face
Face_1 = model.addFace(Part_1_doc, [model.selection("FACE", "Sketch_1/Face-SketchLine_4r-SketchLine_3r-SketchLine_2r-SketchLine_1r"), model.selection("FACE", "Sketch_1/Face-SketchLine_4f-SketchLine_8r-SketchLine_6r-SketchLine_5r")])

### Create Shell
Shell_1 = model.addShell(Part_1_doc, [model.selection("FACE", "Face_1_1"), model.selection("FACE", "Sketch_1/Face-SketchLine_4f-SketchLine_8r-SketchLine_6r-SketchLine_5r")])

### Create Group
Group_1 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Shell_1_1/Modified_Face&Face_1_1/Face_1_1"), model.selection("FACE", "Shell_1_1/Modified_Face&Sketch_1/Face-SketchLine_4f-SketchLine_8r-SketchLine_6r-SketchLine_5r")])
Group_1.setName("surface")
Group_1.result().setName("surface")

### Create Group
Group_2 = model.addGroup(Part_1_doc, "Edges", [model.selection("EDGE", "Face_1_1/Modified_Edge&Sketch_1/SketchLine_2")])
Group_2.setName("left")
Group_2.result().setName("left")

### Create Group
Group_3 = model.addGroup(Part_1_doc, "Edges", [model.selection("EDGE", "Sketch_1/SubEdge_2&Sketch_1/SketchLine_6")])
Group_3.setName("right")
Group_3.result().setName("right")

### Create Group
Group_4 = model.addGroup(Part_1_doc, "Edges", [model.selection("EDGE", "Shell_1_1/Modified_Edge&Sketch_1/SketchLine_3"), model.selection("EDGE", "Shell_1_1/Modified_Edge&Sketch_1/SketchLine_5")])
Group_4.setName("top")
Group_4.result().setName("top")

### Create Group
Group_5 = model.addGroup(Part_1_doc, "Edges", [model.selection("EDGE", "Shell_1_1/Modified_Edge&Sketch_1/SketchLine_1"), model.selection("EDGE", "Shell_1_1/Modified_Edge&Sketch_1/SketchLine_8")])
Group_5.setName("bottom")
Group_5.result().setName("bottom")

### Create Group
Group_6 = model.addGroup(Part_1_doc, "Edges", [model.selection("EDGE", "Shell_1_1/Modified_Edge&Sketch_1/SketchLine_4")])
Group_6.setName("interface")
Group_6.result().setName("interface")

### Create Group
Group_7 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Shell_1_1/Modified_Face&Sketch_1/Face-SketchLine_4f-SketchLine_8r-SketchLine_6r-SketchLine_5r")])
Group_7.setName("surface_right")
Group_7.result().setName("surface_right")

### Create Group
Group_8 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Shell_1_1/Modified_Face&Face_1_1/Face_1_1")])
Group_8.setName("surface_left")
Group_8.result().setName("surface_left")

model.end()

###
### SHAPERSTUDY component
###

model.publishToShaperStudy()
import SHAPERSTUDY
Face_1_2, = SHAPERSTUDY.shape(model.featureStringId(Face_1))
Shell_1_1, surface, left, right, top, bottom, interface, surface_right, surface_left, = SHAPERSTUDY.shape(model.featureStringId(Shell_1))
###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:

Mesh_1 = smesh.Mesh(Shell_1_1,'Mesh_1')
Regular_1D = Mesh_1.Segment()
Local_Length_1 = Regular_1D.LocalLength(0.03,None,1e-07)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
surface_1 = Mesh_1.GroupOnGeom(surface,'surface',SMESH.FACE)
left_1 = Mesh_1.GroupOnGeom(left,'left',SMESH.EDGE)
right_1 = Mesh_1.GroupOnGeom(right,'right',SMESH.EDGE)
top_1 = Mesh_1.GroupOnGeom(top,'top',SMESH.EDGE)
bottom_1 = Mesh_1.GroupOnGeom(bottom,'bottom',SMESH.EDGE)
interface_1 = Mesh_1.GroupOnGeom(interface,'interface',SMESH.EDGE)
isDone = Mesh_1.Compute()
[ surface_1, left_1, right_1, top_1, bottom_1, interface_1 ] = Mesh_1.GetGroups()
GMSH_2D = Mesh_1.Triangle(algo=smeshBuilder.GMSH_2D,geom=surface_left)
Gmsh_Parameters = GMSH_2D.Parameters()
Gmsh_Parameters.Set2DAlgo( 0 )
[ surface_1, left_1, right_1, top_1, bottom_1, interface_1 ] = Mesh_1.GetGroups()
Gmsh_Parameters.SetMinSize( 0.03 )
Gmsh_Parameters.SetMaxSize( 0.03 )
Gmsh_Parameters.SetIs2d( 1 )
[ surface_1, left_1, right_1, top_1, bottom_1, interface_1 ] = Mesh_1.GetGroups()
status = Mesh_1.RemoveHypothesis(GMSH_2D,surface_left)
NETGEN_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_2D,geom=surface_left)
status = Mesh_1.RemoveHypothesis(Gmsh_Parameters,surface_left)
isDone = Mesh_1.Compute()
[ surface_1, left_1, right_1, top_1, bottom_1, interface_1 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED( r'./bar_hybrid_quad_tria.med', 0, 41, 1, Mesh_1, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
Sub_mesh_1 = GMSH_2D.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(GMSH_2D.GetAlgorithm(), 'GMSH_2D')
smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
smesh.SetName(surface_1, 'surface')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(left_1, 'left')
smesh.SetName(right_1, 'right')
smesh.SetName(top_1, 'top')
smesh.SetName(bottom_1, 'bottom')
smesh.SetName(interface_1, 'interface')
smesh.SetName(Gmsh_Parameters, 'Gmsh Parameters')
smesh.SetName(Local_Length_1, 'Local Length_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
