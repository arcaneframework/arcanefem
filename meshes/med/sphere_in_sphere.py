#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.14.0 with dump python functionality
###

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
Sphere_1 = model.addSphere(Part_1_doc, model.selection("VERTEX", "PartSet/Origin"), 0.1)

### Create Sphere
Sphere_2 = model.addSphere(Part_1_doc, model.selection("VERTEX", "PartSet/Origin"), 0.01)

### Create Cut
Cut_1 = model.addCut(Part_1_doc, [model.selection("COMPOUND", "all-in-Sphere_1")], [model.selection("SOLID", "Sphere_2_1")], keepSubResults = True)

### Create Group
Group_1 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Sphere_2_1/Face_1")])
Group_1.setName("inner")
Group_1.result().setName("inner")

### Create Group
Group_2 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Sphere_1_1/Face_1")])
Group_2.setName("outer")
Group_2.result().setName("outer")

### Create Group
Group_3 = model.addGroup(Part_1_doc, "Solids", [model.selection("SOLID", "Cut_1_1")])
Group_3.setName("vol")
Group_3.result().setName("vol")

model.end()

###
### SHAPERSTUDY component
###

model.publishToShaperStudy()
import SHAPERSTUDY
Cut_1_1, inner, outer, vol, = SHAPERSTUDY.shape(model.featureStringId(Cut_1))

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
