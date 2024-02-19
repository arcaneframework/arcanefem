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
Sphere_1 = model.addSphere(Part_1_doc, model.selection("VERTEX", "PartSet/Origin"), 100)

### Create Box
Box_1 = model.addBox(Part_1_doc, 100, 100, 100)

### Create Cut
Cut_1 = model.addCut(Part_1_doc, [model.selection("SOLID", "Sphere_1_1")], [model.selection("SOLID", "Box_1_1")], keepSubResults = True)

### Create Group
Group_1_objects = [model.selection("FACE", "Cut_1_1/Modified_Face&Box_1_1/Bottom"),
                   model.selection("FACE", "Cut_1_1/Modified_Face&Box_1_1/Back"),
                   model.selection("FACE", "Cut_1_1/Modified_Face&Box_1_1/Left")]
Group_1 = model.addGroup(Part_1_doc, "Faces", Group_1_objects)
Group_1.setName("Cut")
Group_1.result().setName("Cut")

### Create Group
Group_2 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Cut_1_1/Modified_Face&Sphere_1_1/Face_1")])
Group_2.setName("sphere")
Group_2.result().setName("sphere")

### Create Group
Group_3 = model.addGroup(Part_1_doc, "Solids", [model.selection("SOLID", "Cut_1_1")])
Group_3.setName("volume")
Group_3.result().setName("volume")

model.end()

###
### SHAPERSTUDY component
###

model.publishToShaperStudy()
import SHAPERSTUDY
Cut_1_1, Cut, sphere, volume, = SHAPERSTUDY.shape(model.featureStringId(Cut_1))
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
NETGEN_3D_Parameters_1.SetMaxSize( 30 )
NETGEN_3D_Parameters_1.SetMinSize( 30 )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 2 )
NETGEN_3D_Parameters_1.SetChordalError( -1 )
NETGEN_3D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_1.SetCheckChartBoundary( 3 )
Cut_2 = Mesh_1.GroupOnGeom(Cut,'Cut',SMESH.FACE)
sphere_1 = Mesh_1.GroupOnGeom(sphere,'sphere',SMESH.FACE)
volume_1 = Mesh_1.GroupOnGeom(volume,'volume',SMESH.VOLUME)
isDone = Mesh_1.Compute()
[ Cut_2, sphere_1, volume_1 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED( r'./sphere_cut.med', 0, 41, 1, Mesh_1, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')


## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(Cut_2, 'Cut')
smesh.SetName(volume_1, 'volume')
smesh.SetName(sphere_1, 'sphere')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
