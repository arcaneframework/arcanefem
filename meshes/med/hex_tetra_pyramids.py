#!/usr/bin/env python

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()

###
### SHAPER component
###

from SketchAPI import *

from salome.shaper import model

model.begin()
partSet = model.moduleDocument()

### Create Part
Part_1 = model.addPart(partSet)
Part_1_doc = Part_1.document()

### Create Sketch
Sketch_1 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_1 = Sketch_1.addLine(0.5, -0.1000000000000483, -2.334139327877319e-21, -0.1)

### Create SketchLine
SketchLine_2 = Sketch_1.addLine(-2.334139327877319e-21, -0.1, 0, 0)

### Create SketchLine
SketchLine_3 = Sketch_1.addLine(0, 0, 0.5, 0)

### Create SketchLine
SketchLine_4 = Sketch_1.addLine(0.5, 0, 0.5, -0.1000000000000483)
Sketch_1.setCoincident(SketchLine_4.endPoint(), SketchLine_1.startPoint())
Sketch_1.setCoincident(SketchLine_1.endPoint(), SketchLine_2.startPoint())
Sketch_1.setCoincident(SketchLine_2.endPoint(), SketchLine_3.startPoint())
Sketch_1.setCoincident(SketchLine_3.endPoint(), SketchLine_4.startPoint())
Sketch_1.setPerpendicular(SketchLine_1.result(), SketchLine_2.result())
Sketch_1.setPerpendicular(SketchLine_2.result(), SketchLine_3.result())
Sketch_1.setPerpendicular(SketchLine_3.result(), SketchLine_4.result())

### Create SketchProjection
SketchProjection_1 = Sketch_1.addProjection(model.selection("EDGE", "PartSet/OX"), False)
SketchLine_5 = SketchProjection_1.createdFeature()
Sketch_1.setCoincident(SketchLine_3.endPoint(), SketchLine_5.result())

### Create SketchLine
SketchLine_6 = Sketch_1.addLine(0.5, -0.1000000000000483, 1, -0.1000000000000483)
Sketch_1.setCoincident(SketchLine_1.startPoint(), SketchLine_6.startPoint())
Sketch_1.setHorizontal(SketchLine_6.result())

### Create SketchLine
SketchLine_7 = Sketch_1.addLine(1, -0.1000000000000483, 1, 0)
Sketch_1.setCoincident(SketchLine_6.endPoint(), SketchLine_7.startPoint())
Sketch_1.setCoincident(SketchLine_7.endPoint(), SketchLine_5.result())
Sketch_1.setVertical(SketchLine_7.result())

### Create SketchLine
SketchLine_8 = Sketch_1.addLine(0.5, 0, 1, 0)
Sketch_1.setCoincident(SketchLine_3.endPoint(), SketchLine_8.startPoint())
Sketch_1.setCoincident(SketchLine_7.endPoint(), SketchLine_8.endPoint())
Sketch_1.setCoincident(SketchLine_2.endPoint(), SketchAPI_Line(SketchLine_5).startPoint())
Sketch_1.setLength(SketchLine_2.result(), 0.1)
Sketch_1.setLength(SketchLine_8.result(), 0.5)
Sketch_1.setLength(SketchLine_3.result(), 0.5)
model.do()

### Create Face
Face_1 = model.addFace(Part_1_doc, [model.selection("FACE", "Sketch_1/Face-SketchLine_4f-SketchLine_6f-SketchLine_7f-SketchLine_8r"), model.selection("FACE", "Sketch_1/Face-SketchLine_4r-SketchLine_3r-SketchLine_2r-SketchLine_1r")])

### Create Shell
Shell_1 = model.addShell(Part_1_doc, [model.selection("FACE", "Face_1_2"), model.selection("FACE", "Face_1_1")])

### Create Extrusion
Extrusion_1 = model.addExtrusion(Part_1_doc, [model.selection("SHELL", "Shell_1_1")], model.selection(), 0.04, 0, "Faces|Wires")

### Create Group
Group_1 = model.addGroup(Part_1_doc, "Solids", [model.selection("SOLID", "Extrusion_1_1_1")])
Group_1.setName("quads_vol")
Group_1.result().setName("quads_vol")

### Create Group
Group_2 = model.addGroup(Part_1_doc, "Solids", [model.selection("SOLID", "Extrusion_1_1_2")])
Group_2.setName("Tria_vol")
Group_2.result().setName("Tria_vol")

### Create Group
Group_3 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Extrusion_1_1_1/Generated_Face&Sketch_1/SketchLine_2")])
Group_3.setName("Fixed")
Group_3.result().setName("Fixed")

### Create Group
Group_4 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Extrusion_1_1_2/Generated_Face&Sketch_1/SketchLine_7")])
Group_4.setName("right")
Group_4.result().setName("right")

### Create Group
Group_5_objects = [model.selection("FACE", "Extrusion_1_1_1/To_Face"),
                   model.selection("FACE", "Extrusion_1_1_2/To_Face"),
                   model.selection("FACE", "Extrusion_1_1_1/Generated_Face&Sketch_1/SketchLine_1"),
                   model.selection("FACE", "Extrusion_1_1_2/Generated_Face&Sketch_1/SketchLine_6"),
                   model.selection("FACE", "Extrusion_1_1_1/Generated_Face&Sketch_1/SketchLine_3"),
                   model.selection("FACE", "Extrusion_1_1_2/Generated_Face&Sketch_1/SketchLine_8"),
                   model.selection("FACE", "Extrusion_1_1_1/From_Face"),
                   model.selection("FACE", "Extrusion_1_1_2/From_Face")]
Group_5 = model.addGroup(Part_1_doc, "Faces", Group_5_objects)
Group_5.setName("top_bot")
Group_5.result().setName("top_bot")

### Create Group
Group_6 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Extrusion_1_1_1/Generated_Face&Sketch_1/SketchLine_4")])
Group_6.setName("interface")
Group_6.result().setName("interface")

model.end()

###
### SHAPERSTUDY component
###

model.publishToShaperStudy()
import SHAPERSTUDY
Extrusion_1_1, quads_vol, Tria_vol, Fixed, right, top_bot, interface, = SHAPERSTUDY.shape(model.featureStringId(Extrusion_1))
###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Extrusion_1_1,'Mesh_1')
Regular_1D = Mesh_1.Segment()
Local_Length_1 = Regular_1D.LocalLength(0.03,None,1e-07)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Hexa_3D = Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)
quads_vol_1 = Mesh_1.GroupOnGeom(quads_vol,'quads_vol',SMESH.VOLUME)
Tria_vol_1 = Mesh_1.GroupOnGeom(Tria_vol,'Tria_vol',SMESH.VOLUME)
Fixed_1 = Mesh_1.GroupOnGeom(Fixed,'Fixed',SMESH.FACE)
right_1 = Mesh_1.GroupOnGeom(right,'right',SMESH.FACE)
top_bot_1 = Mesh_1.GroupOnGeom(top_bot,'top_bot',SMESH.FACE)
interface_1 = Mesh_1.GroupOnGeom(interface,'interface',SMESH.FACE)
isDone = Mesh_1.Compute()
[ quads_vol_1, Tria_vol_1, Fixed_1, right_1, top_bot_1, interface_1 ] = Mesh_1.GetGroups()
NETGEN_3D = Mesh_1.Tetrahedron(geom=Tria_vol)
[ quads_vol_1, Tria_vol_1, Fixed_1, right_1, top_bot_1, interface_1 ] = Mesh_1.GetGroups()
NETGEN_3D_Parameters_1 = NETGEN_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( 1 )
NETGEN_3D_Parameters_1.SetMinSize( 1 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 2 )
NETGEN_3D_Parameters_1.SetElemSizeWeight( 0 )
NETGEN_3D_Parameters_1.SetCheckOverlapping( 2 )
NETGEN_3D_Parameters_1.SetCheckChartBoundary( 3 )
[ quads_vol_1, Tria_vol_1, Fixed_1, right_1, top_bot_1, interface_1 ] = Mesh_1.GetGroups()
isDone = Mesh_1.Compute()
[ quads_vol_1, Tria_vol_1, Fixed_1, right_1, top_bot_1, interface_1 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED(r'./hex_tetra_pyramics.med',auto_groups=1,version=1,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Sub_mesh_1 = NETGEN_3D.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Hexa_3D.GetAlgorithm(), 'Hexa_3D')
smesh.SetName(NETGEN_3D.GetAlgorithm(), 'NETGEN 3D')
smesh.SetName(Tria_vol_1, 'Tria_vol')
smesh.SetName(Fixed_1, 'Fixed')
smesh.SetName(top_bot_1, 'top_bot')
smesh.SetName(quads_vol_1, 'quads_vol')
smesh.SetName(right_1, 'right')
smesh.SetName(interface_1, 'interface')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(Local_Length_1, 'Local Length_1')

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
