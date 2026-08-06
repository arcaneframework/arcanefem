#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.15.0 with dump python functionality
###

import sys
from salome.kernel import salome

salome.salome_init()
from salome.kernel import salome_notebook
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

### Create Box
Box_1 = model.addBox(Part_1_doc, 1, 1, 1)

### Create Group
Group_1 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Box_1_1/Back")])
Group_1.setName("left")
Group_1.result().setName("left")

### Create Group
Group_2 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Box_1_1/Front")])
Group_2.setName("right")
Group_2.result().setName("right")

### Create Group
Group_3 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Box_1_1/Top")])
Group_3.setName("top")
Group_3.result().setName("top")

### Create Group
Group_4 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Box_1_1/Bottom")])
Group_4.setName("bottom")
Group_4.result().setName("bottom")

### Create Group
Group_5 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Box_1_1/Left")])
Group_5.setName("front")
Group_5.result().setName("front")

### Create Group
Group_6 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Box_1_1/Right")])
Group_6.setName("back")
Group_6.result().setName("back")

### Create Group
Group_7 = model.addGroup(Part_1_doc, "Solids", [model.selection("SOLID", "Box_1_1")])
Group_7.setName("volume")
Group_7.result().setName("volume")

model.end()

###
### SHAPERSTUDY component
###

model.publishToShaperStudy()
import SHAPERSTUDY
Box_1_1, left, right, top, bottom, front, back, volume, = SHAPERSTUDY.shape(model.featureStringId(Box_1))
###
### SMESH component
###

from salome.kernel import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Box_1_1,'Mesh_1')
Regular_1D = Mesh_1.Segment()
Number_of_Segments_1 = Regular_1D.NumberOfSegments(5,None,[])
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Hexa_3D = Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)
left_1 = Mesh_1.GroupOnGeom(left,'left',SMESH.FACE)
right_1 = Mesh_1.GroupOnGeom(right,'right',SMESH.FACE)
top_1 = Mesh_1.GroupOnGeom(top,'top',SMESH.FACE)
bottom_1 = Mesh_1.GroupOnGeom(bottom,'bottom',SMESH.FACE)
front_1 = Mesh_1.GroupOnGeom(front,'front',SMESH.FACE)
back_1 = Mesh_1.GroupOnGeom(back,'back',SMESH.FACE)
volume_1 = Mesh_1.GroupOnGeom(volume,'volume',SMESH.VOLUME)
[ left_1, right_1, top_1, bottom_1, front_1, back_1, volume_1 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportGMSHIO(r'unit_cube.hexa.msh', 'Gmsh 4.1 binary', Mesh_1)
  pass
except:
  print('ExportPartToGMSHIO() failed. Invalid file name?')
left_1.SetColor( SALOMEDS.Color( 1, 0.666667, 0 ))
right_1.SetColor( SALOMEDS.Color( 1, 0.666667, 0 ))
top_1.SetColor( SALOMEDS.Color( 1, 0.666667, 0 ))
bottom_1.SetColor( SALOMEDS.Color( 1, 0.666667, 0 ))
[ left_1, right_1, top_1, bottom_1, front_1, back_1, volume_1 ] = Mesh_1.GetGroups()
Number_of_Segments_1.SetNumberOfSegments( 5 )
Number_of_Segments_1.SetScaleFactor( 5 )
Number_of_Segments_1.SetReversedEdges( [] )
isDone = Mesh_1.Compute()
Mesh_1.CheckCompute()
[ left_1, right_1, top_1, bottom_1, front_1, back_1, volume_1 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportGMSHIO(r'unit_cube.hexa-non-cartesian.msh', 'Gmsh 4.1 binary', Mesh_1)
  pass
except:
  print('ExportPartToGMSHIO() failed. Invalid file name?')


## Set names of Mesh objects
smesh.SetName(left_1, 'left')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(bottom_1, 'bottom')
smesh.SetName(back_1, 'back')
smesh.SetName(front_1, 'front')
smesh.SetName(right_1, 'right')
smesh.SetName(Hexa_3D.GetAlgorithm(), 'Hexa_3D')
smesh.SetName(volume_1, 'volume')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(top_1, 'top')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
