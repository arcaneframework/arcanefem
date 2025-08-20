#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.14.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/catA/mb258512/Install/arcanefem/arcanefem/meshes/med')

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

### Create Box
Box_2 = model.addBox(Part_1_doc, 0.75, 0.75, 0.75, 0.25, 0.25, 0.25)

### Create Cut
Cut_1 = model.addCut(Part_1_doc, [model.selection("SOLID", "Box_1_1")], [model.selection("SOLID", "Box_2_1")], keepSubResults = True)

### Create Group
Group_1 = model.addGroup(Part_1_doc, "Vertices", [model.selection("VERTEX", "[Box_2_1/Back][Box_2_1/Left][Box_2_1/Bottom]")])
Group_1.setName("center")
Group_1.result().setName("center")

### Create Group
Group_2 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Box_2_1/Bottom")])
Group_2.setName("horizontal")
Group_2.result().setName("horizontal")

### Create Group
Group_3 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Box_2_1/Back")])
Group_3.setName("verticalYZ")
Group_3.result().setName("verticalYZ")

### Create Group
Group_4 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Box_2_1/Left")])
Group_4.setName("verticalXZ")
Group_4.result().setName("verticalXZ")

### Create Group
Group_5 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Cut_1_1/Modified_Face&Box_1_1/Top")])
Group_5.setName("top")
Group_5.result().setName("top")

### Create Group
Group_6 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Box_1_1/Bottom")])
Group_6.setName("bottom")
Group_6.result().setName("bottom")

### Create Group
Group_7 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Cut_1_1/Modified_Face&Box_1_1/Front")])
Group_7.setName("left")
Group_7.result().setName("left")

### Create Group
Group_8 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Box_1_1/Back")])
Group_8.setName("right")
Group_8.result().setName("right")

### Create Group
Group_9 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Cut_1_1/Modified_Face&Box_1_1/Right")])
Group_9.setName("front")
Group_9.result().setName("front")

### Create Group
Group_10 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Box_1_1/Left")])
Group_10.setName("back")
Group_10.result().setName("back")

### Create Group
Group_11 = model.addGroup(Part_1_doc, "Solids", [model.selection("SOLID", "Cut_1_1")])
Group_11.setName("volume")
Group_11.result().setName("volume")

model.end()

###
### SHAPERSTUDY component
###

model.publishToShaperStudy()
import SHAPERSTUDY
Cut_1_1, center, horizontal, verticalYZ, verticalXZ, top, bottom, left, right, front, back, volume, = SHAPERSTUDY.shape(model.featureStringId(Cut_1))

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
