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

### Create Sphere
Sphere_1 = model.addSphere(Part_1_doc, model.selection("VERTEX", "PartSet/Origin"), 1)

### Create Box
Box_1 = model.addBox(Part_1_doc, 1, 1, 1)

### Create Cut
Cut_1 = model.addCut(Part_1_doc, [model.selection("SOLID", "Sphere_1_1")], [model.selection("SOLID", "Box_1_1")], keepSubResults = True)

### Create Group
Group_1 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Cut_1_1/Modified_Face&Box_1_1/Bottom")])
Group_1.setName("horizontal")
Group_1.result().setName("horizontal")

### Create Group
Group_2 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Cut_1_1/Modified_Face&Sphere_1_1/Face_1")])
Group_2.setName("curved")
Group_2.result().setName("curved")

### Create Group
Group_3 = model.addGroup(Part_1_doc, "Solids", [model.selection("SOLID", "Cut_1_1")])
Group_3.setName("volume")
Group_3.result().setName("volume")

### Create Group
Group_4 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Cut_1_1/Modified_Face&Box_1_1/Back")])
Group_4.setName("verticalYZ")
Group_4.result().setName("verticalYZ")

### Create Group
Group_5 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Cut_1_1/Modified_Face&Box_1_1/Left")])
Group_5.setName("verticalXZ")
Group_5.result().setName("verticalXZ")

model.end()

###
### SHAPERSTUDY component
###

model.publishToShaperStudy()
import SHAPERSTUDY
Cut_1_1, horizontal, curved, volume, verticalYZ, verticalXZ, = SHAPERSTUDY.shape(model.featureStringId(Cut_1))
###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular 

Mesh_1 = smesh.Mesh(Cut_1_1,'Mesh_1')
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( 0.08 )
NETGEN_3D_Parameters_1.SetMinSize( 0.08 )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 2 )
NETGEN_3D_Parameters_1.SetChordalError( -1 )
NETGEN_3D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_1.SetCheckChartBoundary( 3 )
horizontal = Mesh_1.GroupOnGeom(horizontal,'horizontal',SMESH.FACE)
curved = Mesh_1.GroupOnGeom(curved,'curved',SMESH.FACE)
volume = Mesh_1.GroupOnGeom(volume,'volume',SMESH.VOLUME)
verticalYZ = Mesh_1.GroupOnGeom(verticalYZ,'verticalYZ',SMESH.FACE)
verticalXZ = Mesh_1.GroupOnGeom(verticalXZ,'verticalXZ',SMESH.FACE)
isDone = Mesh_1.Compute()
[ horizontal, curved, volume, verticalYZ, verticalXZ ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED( r'./sphere_cut.med', 0, 41, 1, Mesh_1, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
