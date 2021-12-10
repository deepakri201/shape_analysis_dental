# find_distance_dental.py
#
#
# This script performs ICP (iterative closest point) registration between two sets of STL
# files. The mean absolute distance is calculated and displayed between the two meshes.
# A screenshot and the distance mesh in the vtk format are currently saved to disk.
#
# Notes for resources:
#     https://github.com/Kitware/VTK/blob/master/Examples/GUI/Python/ImplicitPlaneWidget.py
#     https://github.com/Kitware/VTK/blob/master/Examples/VisualizationAlgorithms/Python/ClipCow.py
#     https://gitlab.kitware.com/bxa/vtk/blob/a7996f418f2e1375aa91cbd98ab952a7162b4a26/Examples/Annotation/Python/annotatePick.py
#     http://vtk.1045678.n5.nabble.com/Picking-and-moving-a-vertex-td1233640.html
#     https://vtk.org/Wiki/VTK/Examples/Python/Visualization/AssignColorsCellFromLUT
#     https://github.com/Kitware/VTK/blob/master/Examples/GUI/Python/BoxWidget.py
#     https://python.hotexamples.com/examples/vtk/-/vtkPointPicker/python-vtkpointpicker-function-examples.html
#
#
# Deepa Krishnaswamy
# Brigham and Women's Hospital
# November 2021
#########################################################################
# -*- coding: utf-8 -*-

import numpy  as np
import vtk
import sys
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
print (str(vtk.VTK_MAJOR_VERSION) + '.' + str(vtk.VTK_MINOR_VERSION))
from PyQt5.QtWidgets import QApplication, QFileDialog

from vtk.util.colors import peacock, tomato
from vtk.util.colors import brown_ochre, tomato, banana, azure

import re

#########################################################################

def calculate_ICP(source, target):
    icp = vtk.vtkIterativeClosestPointTransform()
    icp.SetSource(source)
    icp.SetTarget(target)
    icp.GetLandmarkTransform().SetModeToRigidBody()
    #icp.DebugOn()
    # icp.SetMaximumNumberOfIterations(20)
    icp.SetMaximumNumberOfIterations(150)
    icp.StartByMatchingCentroidsOn()
    icp.Modified()
    icp.Update()

    icpTransformFilter = vtk.vtkTransformPolyDataFilter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        icpTransformFilter.SetInput(source)
    else:
        icpTransformFilter.SetInputData(source)

    icpTransformFilter.SetTransform(icp)
    icpTransformFilter.Update()

    transformedSource = icpTransformFilter.GetOutput()

    return transformedSource

def calculate_MAD(mesh_A, mesh_B):
    '''Calculates the mean absolute distance (MAD) between two meshes'''
    '''The distance is calculated between the cells'''

    ### calculate MAD metric ###
    distance_filter = vtk.vtkDistancePolyDataFilter()
    distance_filter.SetSignedDistance(0)
    # added
    distance_filter.ComputeSecondDistanceOff()
    distance_filter.SetComputeSecondDistance(0)
    # distance between ref and my mesh
    if vtk.VTK_MAJOR_VERSION <= 5:
        distance_filter.SetInput(0,mesh_A)
        distance_filter.SetInput(1,mesh_B)
    else:
        distance_filter.SetInputData(0,mesh_A)
        distance_filter.SetInputData(1,mesh_B)
    distance_filter.Update()
    # get distance
    # distances = vtk_to_numpy(distance_filter.GetOutput().GetPointData().GetScalars())
    distances = vtk_to_numpy(distance_filter.GetOutput().GetCellData().GetScalars())
    mean_distance = np.mean(distances)
    # get distance meshes
    distance_mesh = vtk.vtkPolyData()
    distance_mesh = distance_filter.GetOutput()

    return distance_mesh, distances, mean_distance

#########################################################################

app = QApplication(sys.argv)
options = QFileDialog.Options()
options |= QFileDialog.DontUseNativeDialog
fnames = QFileDialog.getOpenFileNames(None,"Select STL File","","STL Files (*.stl);;All Files (*)", "", options)[0]

if len(fnames) == 2:
    # Sort filenames so T1 is first, T2 is second
    fnames.sort(key=lambda f: int(re.sub('\D', '', f)))
    mesh1_filename = fnames[0]
    mesh2_filename = fnames[1]

else:
    sys.exit(app.exec_())

# Outputs
output_mesh_filename = 'distance_mesh.vtk'
output_jpg_filename = 'distance_mesh.jpg'
output_txt_filename = 'mean_absolute_distance_in_mm.txt'

print ('mesh1_filename: ' + str(mesh1_filename))
print ('mesh2_filename: ' + str(mesh2_filename))
#########################################################################

# read mesh
reader1 = vtk.vtkSTLReader()
reader1.SetFileName(mesh1_filename)
reader1.Update()
mesh1 = reader1.GetOutput()

# read gt mesh
reader2 = vtk.vtkSTLReader()
reader2.SetFileName(mesh2_filename)
reader2.Update()
mesh2 = reader2.GetOutput()

mapper1, mapper2, mapper1t = [vtk.vtkPolyDataMapper() for i in range(3)]
actor1, actor2, actor1t = [vtk.vtkActor() for i in range(3)]

mapper1.SetInputData(mesh1)
mapper2.SetInputData(mesh2)

actor1.SetMapper(mapper1)
actor2.SetMapper(mapper2)

actor1.GetProperty().SetColor(0.8,0.4,0.2) # orange = T1
actor2.GetProperty().SetColor(0.4,0.8,0.2) # green = T2

# calculate MAD

mesh1t = calculate_ICP(mesh1, mesh2)
distance_mesh, distances, mean_distance = calculate_MAD(mesh2, mesh1t)

# print maximum distance
print ('max distance: ' + str(np.max(distances)))

distance_meshb, distancesb, mean_distanceb = calculate_MAD(mesh1t, mesh2)

dist_min = distance_mesh.GetPointData().GetScalars().GetRange()[0]
dist_minb = distance_meshb.GetPointData().GetScalars().GetRange()[0]

dist_max = distance_mesh.GetPointData().GetScalars().GetRange()[1]
dist_maxb = distance_meshb.GetPointData().GetScalars().GetRange()[1]

# save distance mesh
polydata_writer = vtk.vtkPolyDataWriter()
polydata_writer.SetFileName(output_mesh_filename)
polydata_writer.SetInputData(distance_mesh)
polydata_writer.Write()

# save the MAD to a text file
with open(output_txt_filename, 'w') as f:
    f.write(str(np.round(mean_distance,4)))

# Create renderers
renWin = vtk.vtkRenderWindow()
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

mapper1t.SetInputData(mesh1t)
actor1t.SetMapper(mapper1t)
actor1t.GetProperty().SetColor(0.8,0.4,0.2)

mapper, mapperb = [vtk.vtkPolyDataMapper() for i in range(2)]
actor, actorb = [vtk.vtkActor() for i in range(2)]
ren1, ren1b, ren2, ren3 = [vtk.vtkRenderer() for i in range(4)]

if vtk.VTK_MAJOR_VERSION <= 5:
    mapper.SetInput(distance_mesh)
    mapperb.SetInput(distance_meshb)
else:
    mapper.SetInputData(distance_mesh)
    mapperb.SetInputData(distance_meshb)

# Set range between 0 and 3
# mapper.SetScalarRange(0, 3)
# mapper.SetScalarRange(0, 1)
# mapper.SetScalarRange(0, 1.6)
# mapper.SetScalarRange(0, 2)
mapper.SetScalarRange(0, 2.4)

# Create own look up table with set colors
hueLut = vtk.vtkLookupTable()
# hueLut.SetNumberOfTableValues(6)
hueLut.SetNumberOfTableValues(5)
hueLut.Build()
nc = vtk.vtkNamedColors()
# hueLut.SetTableValue(0,nc.GetColor4d("Blue"))
# hueLut.SetTableValue(1,nc.GetColor4d("Purple"))
# hueLut.SetTableValue(2,nc.GetColor4d("Green"))
# hueLut.SetTableValue(3,nc.GetColor4d("Yellow"))
# hueLut.SetTableValue(4,nc.GetColor4d("Orange"))
# hueLut.SetTableValue(5,nc.GetColor4d("Red"))
hueLut.SetTableValue(0,nc.GetColor4d("Blue"))
hueLut.SetTableValue(1,nc.GetColor4d("Green"))
hueLut.SetTableValue(2,nc.GetColor4d("Yellow"))
hueLut.SetTableValue(3,nc.GetColor4d("Orange"))
hueLut.SetTableValue(4,nc.GetColor4d("Red"))

mapper.SetLookupTable(hueLut)

# Set the scalar bar
scalarBar = vtk.vtkScalarBarActor()
scalarBar.SetLookupTable(hueLut)
scalarBar.SetTitle("mean distance")
# scalarBar.SetMaximumNumberOfColors(6); # 0.5 between each
# scalarBar.SetNumberOfLabels(7);
scalarBar.SetMaximumNumberOfColors(5) # 0.2 between each
scalarBar.SetNumberOfLabels(6)

# Set range between 0 and 3 for b
# mapperb.SetScalarRange(0, 3)
# mapperb.SetScalarRange(0, 1)
# mapperb.SetScalarRange(0, 1.6)
# mapperb.SetScalarRange(0, 2)
mapperb.SetScalarRange(0, 2.4)

# Create own look up table with set colors for b
hueLutb = vtk.vtkLookupTable()
# hueLutb.SetNumberOfTableValues(6)
hueLutb.SetNumberOfTableValues(5)
hueLutb.Build()
nc = vtk.vtkNamedColors()
# hueLutb.SetTableValue(0,nc.GetColor4d("Blue"))
# hueLutb.SetTableValue(1,nc.GetColor4d("Purple"))
# hueLutb.SetTableValue(2,nc.GetColor4d("Green"))
# hueLutb.SetTableValue(3,nc.GetColor4d("Yellow"))
# hueLutb.SetTableValue(4,nc.GetColor4d("Orange"))
# hueLutb.SetTableValue(5,nc.GetColor4d("Red"))
hueLutb.SetTableValue(0,nc.GetColor4d("Blue"))
hueLutb.SetTableValue(1,nc.GetColor4d("Green"))
hueLutb.SetTableValue(2,nc.GetColor4d("Yellow"))
hueLutb.SetTableValue(3,nc.GetColor4d("Orange"))
hueLutb.SetTableValue(4,nc.GetColor4d("Red"))

mapperb.SetLookupTable(hueLutb)

# Set the scalar bar b
scalarBarb = vtk.vtkScalarBarActor()
scalarBarb.SetLookupTable(hueLutb)
scalarBarb.SetTitle("mean distance")
# scalarBarb.SetMaximumNumberOfColors(6); # 0.5 between each
# scalarBarb.SetNumberOfLabels(7);
scalarBarb.SetMaximumNumberOfColors(5); # 0.2 between each
scalarBarb.SetNumberOfLabels(6);

# Add both actors to renderer
ren1.AddActor2D(scalarBar)
ren1b.AddActor2D(scalarBarb)

#####################################################
### Set up cutting box - 3 planes and the clipper ###
#####################################################

normals = vtk.vtkPolyDataNormals()
normals.SetInputData(distance_mesh)
normals.Update()

points = vtk_to_numpy(distance_mesh.GetPoints().GetData())
bounds_min = np.min(points,axis=0)
bounds_max = np.max(points,axis=0)
planes = vtk.vtkPlanes()
planes.SetBounds(bounds_min[0], bounds_max[0], bounds_min[1], bounds_max[1], bounds_min[2], bounds_max[2])

clipper = vtk.vtkClipPolyData()
clipper.SetInputData(distance_mesh)
clipper.SetClipFunction(planes)
clipper.InsideOutOn()
clipper.GenerateClippedOutputOn()

clipMapper = vtk.vtkPolyDataMapper()
clipMapper.SetInputConnection(clipper.GetOutputPort())
clipMapper.ScalarVisibilityOn()
# clipMapper.SetScalarRange(0, 3)
# clipMapper.SetScalarRange(0, 1)
# clipMapper.SetScalarRange(0, 1.6)
# clipMapper.SetScalarRange(0, 2)
clipMapper.SetScalarRange(0, 2.4)
clipMapper.SetLookupTable(hueLut)

clipActor = vtk.vtkActor()
clipActor.SetMapper(clipMapper)

# create a MAD annotation
corner_MAD = vtk.vtkCornerAnnotation()
corner_MAD.SetMinimumFontSize(16)
corner_MAD.SetMaximumFontSize(16)
corner_MAD.SetText(7, "MAD: " + str(np.round(mean_distance,4)))

boxWidget = vtk.vtkBoxWidget()
boxWidget.SetCurrentRenderer(ren1)
boxWidget.SetInteractor(iren)
boxWidget.SetPlaceFactor(1.25)

def myCallback(obj, event):
    '''Callback function for the plane'''

    global clipActor, planes
    obj.GetPlanes(planes)
    clipActor.VisibilityOn()

    ### Update the distance mesh here by the new clip plane ###
    # Don't need to update scalar bars since they are set always
    # from 0 to 3.
    # Clip mesh2 and mesh1t using the same clip plane and
    # calculate a new distance mesh

    global mesh2
    global mesh1t
    global distance_mesh

    clipper_mesh2 = vtk.vtkClipPolyData()
    clipper_mesh2.SetInputData(mesh2)
    clipper_mesh2.SetClipFunction(planes)
    # clipper_mesh2.InsideOutOn()
    clipper_mesh2.GenerateClippedOutputOn()
    clipper_mesh2.Update()
    mesh2_clipped = vtk.vtkPolyData()
    mesh2_clipped = clipper_mesh2.GetClippedOutput()

    clipper_mesh1t = vtk.vtkClipPolyData()
    clipper_mesh1t.SetInputData(mesh1t)
    clipper_mesh1t.SetClipFunction(planes)
    # clipper_mesh1t.InsideOutOn()
    clipper_mesh1t.GenerateClippedOutputOn()
    clipper_mesh1t.Update()
    mesh1t_clipped = vtk.vtkPolyData()
    mesh1t_clipped = clipper_mesh1t.GetClippedOutput()

    distance_mesh, distances, mean_distance = calculate_MAD(mesh2_clipped, mesh1t_clipped)

    global corner_MAD
    corner_MAD.SetMinimumFontSize(14)
    corner_MAD.SetMaximumFontSize(14)
    corner_MAD.SetText(7, "MAD: " + str(np.round(mean_distance,4)))

    # save distance mesh that has been clipped
    polydata_writer = vtk.vtkPolyDataWriter()
    polydata_writer.SetFileName(output_mesh_filename)
    polydata_writer.SetInputData(distance_mesh)
    polydata_writer.Write()

    # save the MAD to a text file
    with open(output_txt_filename, 'w') as f:
        f.write(str(np.round(mean_distance,4)))



# Interaction event
boxWidget.SetInputData(distance_mesh)
boxWidget.PlaceWidget()
boxWidget.On()
boxWidget.AddObserver("EndInteractionEvent", myCallback) # should be faster compared to "InteractionEvent"

##############################################################
### Set up the clicking of a point to display the distance ###
##############################################################

# # Create the point picker
# picker = vtk.vtkPointPicker()
# # create a text actor
# txt1 = vtk.vtkTextActor()
#
# def annotate_pick(obj, event):
#     '''Callback function for displaying the scalar value of a clicked point'''
#
#     global ren1
#     global iren
#     global txt_actor
#     global distance_mesh
#
#     pointPicker = vtk.vtkPointPicker()
#     eventPosition = obj.GetEventPosition()
#     result = picker.Pick(float(eventPosition[0]),float(eventPosition[1]),0.0,ren1)
#     pickId = picker.GetPointId()
#
#     # Get scalar that is associated with this pickId
#     scalars = vtk_to_numpy(distance_mesh.GetPointData().GetScalars())
#     dist_value = scalars[pickId]
#
#     # Set as input to text actor
#     txt1.SetInput("Picked point distance value: " + str(np.round(dist_value,4)))


# Create the cell picker
picker = vtk.vtkCellPicker()
# create a text actor
txt1 = vtk.vtkTextActor()

def annotate_pick(obj, event):
    '''Callback function for displaying the scalar value of a clicked point'''

    global ren1
    global iren
    global txt_actor
    global distance_mesh

    # pointPicker = vtk.vtkPointPicker()
    cellPicker = vtk.vtkCellPicker()
    eventPosition = obj.GetEventPosition()
    result = picker.Pick(float(eventPosition[0]),float(eventPosition[1]),0.0,ren1)
    pickId = picker.GetCellId()

    # Get scalar that is associated with this pickId
    scalars = vtk_to_numpy(distance_mesh.GetCellData().GetScalars())
    dist_value = scalars[pickId]

    # Set as input to text actor
    txt1.SetInput("Picked point distance value: " + str(np.round(dist_value,4)))


iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
iren.AddObserver("RightButtonPressEvent", annotate_pick)

# Text actor properties
txtprop=txt1.GetTextProperty()
txtprop.SetFontFamilyToArial()
txtprop.SetFontSize(14)
txtprop.SetColor(1,1,1)
txt1.SetTextProperty(txtprop)
txt1.SetPosition2(20,30)

corner1 = vtk.vtkCornerAnnotation()
corner1.SetMinimumFontSize(14)
corner1.SetMaximumFontSize(14)
corner1.SetText(2, "Mean absolute distance \n Rotate with left mouse \n Move object or box planes by middle mouse button \n Turn on/off box planes by pressing i \n Move planes with left mouse \n Select cell on condyle with right mouse")

corner2 = vtk.vtkCornerAnnotation()
corner2.SetMinimumFontSize(14)
corner2.SetMaximumFontSize(14)
# corner2.SetText(2,"Before alignment \n T1 (orange) and T2 (green)")
corner2.SetText(2,"Before alignment \n scan 1 (orange) and scan 2 (green)")

corner3 = vtk.vtkCornerAnnotation()
corner3.SetMinimumFontSize(14)
corner3.SetMaximumFontSize(14)
# corner3.SetText(2,"After alignment \n T1 (orange) and T2 (green)")
corner3.SetText(2,"After alignment \n scan 1 (orange) and scan 2 (green)")

#################
### Rendering ###
#################

# actor.SetMapper(mapper)
# ren1.AddActor(actor) # This is the original distance mesh

ren1.AddActor(clipActor)
ren1.AddActor2D(txt1)
ren1.AddActor2D(corner1)
ren1.AddActor2D(corner_MAD)
ren1.SetViewport(0.5, 0.5, 1.0, 1.0)

actorb.SetMapper(mapperb)
ren1b.AddActor(actorb)
ren1b.SetViewport(0.5, 0.0, 1.0, 0.5)

ren2.AddActor(actor1)
ren2.AddActor(actor2)
ren2.AddActor2D(corner2)
ren2.SetViewport(0.0, 0.5, 0.5, 1.0)

ren3.AddActor(actor1t)
ren3.AddActor(actor2)
ren3.AddActor2D(corner3)
ren3.SetViewport(0.0, 0.0, 0.5, 0.5)

ren1b.SetBackground(0.05, 0.05, 0.05)
ren3.SetBackground(0.1, 0.1, 0.1)
ren2.SetBackground(0.15, 0.15, 0.15)

renWin.AddRenderer(ren1)
renWin.AddRenderer(ren1b)
renWin.AddRenderer(ren2)
renWin.AddRenderer(ren3)
renWin.Render()

renWin.SetSize(1920,1080)
renWin.Render()

# write out screenshot
w2if = vtk.vtkWindowToImageFilter()
w2if.SetInput(renWin)
w2if.Update()

writer = vtk.vtkJPEGWriter()
writer.SetFileName(output_jpg_filename)
writer.SetInputConnection(w2if.GetOutputPort())
writer.Write()

# iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

# added
# iren.SetPicker(picker)

iren.Initialize()

iren.Start()

# sys.exit(app.exec_())
