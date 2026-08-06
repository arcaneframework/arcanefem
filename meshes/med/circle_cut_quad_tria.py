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

from SketchAPI import *

from salome.shaper import model

model.begin()
partSet = model.moduleDocument()

### Create Part
Part_1 = model.addPart(partSet)
Part_1_doc = Part_1.document()

### Create Sketch
Sketch_1 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchProjection
SketchProjection_1 = Sketch_1.addProjection(model.selection("VERTEX", "PartSet/Origin"), False)
SketchPoint_1 = SketchProjection_1.createdFeature()

### Create SketchProjection
SketchProjection_2 = Sketch_1.addProjection(model.selection("EDGE", "PartSet/OY"), False)
SketchLine_1 = SketchProjection_2.createdFeature()

### Create SketchProjection
SketchProjection_3 = Sketch_1.addProjection(model.selection("EDGE", "PartSet/OX"), False)
SketchLine_2 = SketchProjection_3.createdFeature()

### Create SketchArc
SketchArc_1 = Sketch_1.addArc(0, 0, 0, 1, 1, 0, False)
Sketch_1.setCoincident(SketchPoint_1.result(), SketchArc_1.center(), True)
Sketch_1.setCoincident(SketchLine_1.result(), SketchArc_1.startPoint(), True)
Sketch_1.setCoincident(SketchLine_2.result(), SketchArc_1.endPoint(), True)

### Create SketchLine
SketchLine_3 = Sketch_1.addLine(0, 1, 0, 0)
Sketch_1.setCoincident(SketchArc_1.startPoint(), SketchLine_3.startPoint(), True)
Sketch_1.setCoincident(SketchAPI_Point(SketchPoint_1).coordinates(), SketchLine_3.endPoint(), True)

### Create SketchLine
SketchLine_4 = Sketch_1.addLine(0, 0, 1, 0)
Sketch_1.setCoincident(SketchLine_3.endPoint(), SketchLine_4.startPoint(), True)
Sketch_1.setCoincident(SketchArc_1.endPoint(), SketchLine_4.endPoint(), True)

### Create SketchLine
SketchLine_5 = Sketch_1.addLine(0, 0, -0.6427876096865393, -0.7660444431189781)
Sketch_1.setCoincident(SketchAPI_Point(SketchPoint_1).coordinates(), SketchLine_5.startPoint(), True)
Sketch_1.setCoincident(SketchLine_5.endPoint(), SketchArc_1.results()[1], True)
Sketch_1.setLength(SketchLine_4.result(), 1, True)

### Create SketchConstraintAngle
Sketch_1.setAngle(SketchLine_5.result(), SketchLine_3.result(), 140, type = "Direct", is_active = True)
model.do()

### Create Shell
Shell_1 = model.addShell(Part_1_doc, [model.selection("FACE", "Sketch_1/Face-SketchArc_1_2f-SketchLine_5r-SketchLine_3r"), model.selection("FACE", "Sketch_1/Face-SketchArc_1_2f-SketchLine_4r-SketchLine_5f")])

### Create Group
Group_1 = model.addGroup(Part_1_doc, "Edges", [model.selection("EDGE", "Sketch_1/SubEdge_1&Sketch_1/SketchLine_3")])
Group_1.setName("vertical")
Group_1.result().setName("vertical")

### Create Group
Group_2 = model.addGroup(Part_1_doc, "Edges", [model.selection("EDGE", "Sketch_1/SubEdge_2&Sketch_1/SketchLine_4")])
Group_2.setName("horizontal")
Group_2.result().setName("horizontal")

### Create Group
Group_3 = model.addGroup(Part_1_doc, "Edges", [model.selection("EDGE", "Sketch_1/SubEdge_1&Sketch_1/SketchArc_1_2"), model.selection("EDGE", "Sketch_1/SubEdge_2&Sketch_1/SketchArc_1_2")])
Group_3.setName("curved")
Group_3.result().setName("curved")

### Create Group
Group_4 = model.addGroup(Part_1_doc, "Edges", [model.selection("EDGE", "Sketch_1/SubEdge_2&Sketch_1/SketchLine_5")])
Group_4.setName("interface")
Group_4.result().setName("interface")

### Create Group
Group_5 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Sketch_1/Face-SketchArc_1_2f-SketchLine_5r-SketchLine_3r")])
Group_5.setName("surfaceQuad")
Group_5.result().setName("surfaceQuad")

### Create Group
Group_6 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Sketch_1/Face-SketchArc_1_2f-SketchLine_4r-SketchLine_5f")])
Group_6.setName("surfaceTria")
Group_6.result().setName("surfaceTria")

model.end()

###
### SHAPERSTUDY component
###

model.publishToShaperStudy()
import SHAPERSTUDY
Shell_1_1, vertical, horizontal, curved, interface, surfaceQuad, surfaceTria, = SHAPERSTUDY.shape(model.featureStringId(Shell_1))
###
### SMESH component
###

from salome.kernel import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Gmsh_Parameters = smesh.CreateHypothesis('GMSH_Parameters_2D', 'GMSHEngine')
Gmsh_Parameters.Set2DAlgo( 6 )
GMSH_2D = smesh.CreateHypothesis('GMSH_2D', 'GMSHEngine')
Gmsh_Parameters_1 = smesh.CreateHypothesis('GMSH_Parameters_2D', 'GMSHEngine')
Gmsh_Parameters_1.Set2DAlgo( 2 )
Gmsh_Parameters_1.SetMinSize( 0.2 )
Gmsh_Parameters_1.SetMaxSize( 0.2 )
Gmsh_Parameters_1.SetIs2d( 1 )
Gmsh_Parameters.SetMinSize( 0.1 )
Gmsh_Parameters.SetMaxSize( 0.1 )
Gmsh_Parameters.SetIs2d( 1 )
Mesh_1 = smesh.Mesh(Shell_1_1,'Mesh_1')
NETGEN_1D_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( 0.1 )
NETGEN_2D_Parameters_1.SetMinSize( 0.1 )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 2 )
NETGEN_2D_Parameters_1.SetChordalError( -1 )
NETGEN_2D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_1.SetWorstElemMeasure( 0 )
NETGEN_2D_Parameters_1.SetUseDelauney( 0 )
NETGEN_2D_Parameters_1.SetCheckChartBoundary( 0 )
vertical_1 = Mesh_1.GroupOnGeom(vertical,'vertical',SMESH.EDGE)
horizontal_1 = Mesh_1.GroupOnGeom(horizontal,'horizontal',SMESH.EDGE)
curved_1 = Mesh_1.GroupOnGeom(curved,'curved',SMESH.EDGE)
surfaceQuad_1 = Mesh_1.GroupOnGeom(surfaceQuad,'surfaceQuad',SMESH.FACE)
surfaceTria_1 = Mesh_1.GroupOnGeom(surfaceTria,'surfaceTria',SMESH.FACE)
Gmsh_Parameters_2 = smesh.CreateHypothesis('GMSH_Parameters_2D', 'GMSHEngine')
Gmsh_Parameters_2.SetMinSize( 0.1 )
Gmsh_Parameters_2.SetMaxSize( 0.1 )
status = Mesh_1.AddHypothesis(GMSH_2D,surfaceQuad)
status = Mesh_1.AddHypothesis(Gmsh_Parameters_2,surfaceQuad)
Gmsh_Parameters_2.Set2DAlgo( 6 )
Gmsh_Parameters_2.SetIs2d( 1 )
isDone = Mesh_1.Compute()
Mesh_1.CheckCompute()
[ vertical_1, horizontal_1, curved_1, smeshObj_1, surfaceQuad_1, surfaceTria_1 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportGMSHIO(r'circle_cut-quad-tria.msh', 'Gmsh 4.1 binary', Mesh_1)
  pass
except:
  print('ExportPartToGMSHIO() failed. Invalid file name?')
Sub_mesh_1 = Mesh_1.GetSubMesh( surfaceQuad, 'Sub-mesh_1' )

## some objects were removed
aStudyBuilder = salome.myStudy.NewBuilder()
SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(smeshObj_1))
if SO: aStudyBuilder.RemoveObjectWithChildren(SO)

## Set names of Mesh objects
smesh.SetName(Gmsh_Parameters_2, 'Gmsh Parameters')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
smesh.SetName(Gmsh_Parameters, 'Gmsh Parameters')
smesh.SetName(horizontal_1, 'horizontal')
smesh.SetName(surfaceTria_1, 'surfaceTria')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(curved_1, 'curved')
smesh.SetName(vertical_1, 'vertical')
smesh.SetName(surfaceQuad_1, 'surfaceQuad')
smesh.SetName(GMSH_2D, 'GMSH_2D')
smesh.SetName(Gmsh_Parameters_1, 'Gmsh Parameters')
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
