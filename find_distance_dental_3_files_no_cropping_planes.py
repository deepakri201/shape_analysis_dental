# find_distance_dental_3_files_no_cropping_planes.py
# 
# This script performs ICP (iterative closest point) registration between two sets of STL
# files, the T1, T2 and the impression. T1 is registered to T2, and the impression is also 
# registered to T2. 
# 
# The mean absolute distance is calculated and displayed the three sets of meshes: 
#   1. Between T2 and T1 registered to T2 
#   2. Between T2 and the Impression registered to T2 
#   3. Between T1 registered to T2, and the Impression registered to T2 
# 
# The distance mesh in the vtk format is saved for each comparison above. A screenshot 
# is also saved to the same folder. 
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
# February 2023
#########################################################################
# -*- coding: utf-8 -*-

import os 
import sys
import csv 
import numpy  as np
import vtk
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

### Select the stl files ### 

app = QApplication(sys.argv)
options = QFileDialog.Options()
options |= QFileDialog.DontUseNativeDialog
fnames = QFileDialog.getOpenFileNames(None,"Select STL File","","STL Files (*.stl);;All Files (*)", "", options)[0]

if len(fnames) == 3:
    # Sort filenames so T1 is first, T2 is second, then impression
    fnames.sort(key=lambda f: int(re.sub('\D', '', f)))
    mesh1_filename = fnames[0]
    mesh2_filename = fnames[1]
    mesh3_filename = fnames[2] 

else:
    sys.exit(app.exec_())

#########################################################################

output_directory = os.path.join(os.getcwd(), 'results')
if not os.path.isdir(output_directory):
  os.mkdir(output_directory)

output_mesh_filename_T2_to_T1t = os.path.join(output_directory, 'distance_mesh_T2_to_T1t.vtk')
output_mesh_filename_T2_to_impt = os.path.join(output_directory, 'distance_mesh_T2_to_impt.vtk')
output_mesh_filename_T1t_to_impt = os.path.join(output_directory, 'distance_mesh_T1t_to_impt.vtk')
# Write out to one text file 
output_csv_filename = os.path.join(output_directory, 'mean_absolute_distance_in_mm.csv')
# output_txt_filename_T2_to_T1t = 'mean_absolute_distance_in_mm_T2_to_T1t.txt'
# output_txt_filename_T2_to_impt = 'mean_absolute_distance_in_mm_T2_to_impt.txt'
# output_txt_filename_T1t_to_impt = 'mean_absolute_distance_in_mm_T1t_to_impt.txt'

output_jpg_filename = os.path.join(output_directory,'distance_mesh.jpg')

# # Outputs
# output_mesh_filename_1 = 'distance_mesh_T1.vtk'
# output_txt_filename_1 = 'mean_absolute_distance_in_mm_T1.txt'
# output_jpg_filename = 'distance_mesh.jpg'
#
# output_mesh_filename_3 = 'distance_mesh_impression.vtk'
# output_txt_filename_3 = 'mean_absolute_distance_in_mm_impression.txt'
#
# print ('mesh1_filename: ' + str(mesh1_filename))
# print ('mesh2_filename: ' + str(mesh2_filename))
# print ('mesh3_filename: ' + str(mesh3_filename))
#########################################################################

# read mesh T1
reader1 = vtk.vtkSTLReader()
reader1.SetFileName(mesh1_filename)
reader1.Update()
mesh1 = reader1.GetOutput()

# read gt mesh T2
reader2 = vtk.vtkSTLReader()
reader2.SetFileName(mesh2_filename)
reader2.Update()
mesh2 = reader2.GetOutput()

# read 3rd mesh impression
reader3 = vtk.vtkSTLReader()
reader3.SetFileName(mesh3_filename)
reader3.Update()
mesh3 = reader3.GetOutput()

#########################################################################
### Set the mappers and actors ### 

mapper1, mapper2, mapper1t, mapper3, mapper3t  = [vtk.vtkPolyDataMapper() for i in range(5)]
actor1, actor2, actor1t, actor3, actor3t = [vtk.vtkActor() for i in range(5)]

mapper1.SetInputData(mesh1)
mapper2.SetInputData(mesh2)
mapper3.SetInputData(mesh3)

actor1.SetMapper(mapper1)
actor2.SetMapper(mapper2)
actor3.SetMapper(mapper3)

actor1.GetProperty().SetColor(0.8,0.4,0.2) # orange = T1
actor2.GetProperty().SetColor(0.4,0.8,0.2) # green = T2
actor3.GetProperty().SetColor(1.0,0.2,0.4) # red = impression 

#########################################################################
### calculate MAD ###

mesh1t = calculate_ICP(mesh1, mesh2)
distance_mesh_1t2, distances_1t2, mean_distance_1t2 = calculate_MAD(mesh2, mesh1t)

mesh3t = calculate_ICP(mesh3, mesh2)
distance_mesh_3t2, distances_3t2, mean_distance_3t2 = calculate_MAD(mesh2, mesh3t)

distance_mesh_1t3t, distances_1t3t, mean_distance_1t3t = calculate_MAD(mesh1t, mesh3t)

### get min and max distances and save distance meshes ### 

# print maximum distance
# print ('max distance: ' + str(np.max(distances)))

# distance_meshb, distancesb, mean_distanceb = calculate_MAD(mesh1t, mesh2)
# dist_min = distance_mesh.GetPointData().GetScalars().GetRange()[0]
# dist_minb = distance_meshb.GetPointData().GetScalars().GetRange()[0]
# dist_max = distance_mesh.GetPointData().GetScalars().GetRange()[1]
# dist_maxb = distance_meshb.GetPointData().GetScalars().GetRange()[1]

dist_min_1t2 = distance_mesh_1t2.GetPointData().GetScalars().GetRange()[0] 
dist_max_1t2 = distance_mesh_1t2.GetPointData().GetScalars().GetRange()[1] 

# save distance mesh 1t2 
polydata_writer = vtk.vtkPolyDataWriter()
polydata_writer.SetFileName(output_mesh_filename_T2_to_T1t)
polydata_writer.SetInputData(distance_mesh_1t2) 
polydata_writer.Write()

dist_min_3t2 = distance_mesh_3t2.GetPointData().GetScalars().GetRange()[0] 
dist_max_3t2 = distance_mesh_3t2.GetPointData().GetScalars().GetRange()[1] 

# save distance mesh 3t2 
polydata_writer = vtk.vtkPolyDataWriter()
polydata_writer.SetFileName(output_mesh_filename_T2_to_impt)
polydata_writer.SetInputData(distance_mesh_3t2)
polydata_writer.Write()

dist_min_1t3t = distance_mesh_1t3t.GetPointData().GetScalars().GetRange()[0] 
dist_max_1t3t = distance_mesh_1t3t.GetPointData().GetScalars().GetRange()[1] 

# save distance mesh 1t3t
polydata_writer = vtk.vtkPolyDataWriter()
polydata_writer.SetFileName(output_mesh_filename_T1t_to_impt)
polydata_writer.SetInputData(distance_mesh_1t3t)
polydata_writer.Write()

# save the MAD to a text file
# with open(output_txt_filename_T2_to_T1t, 'w') as f:
#     f.write(str(np.round(mean_distance_1t2,4)))
#
# # save the MAD to a text file
# with open(output_txt_filename_T2_to_impt, 'w') as f:
#     f.write(str(np.round(mean_distance_3t2,4)))
#
# # save the MAD to a text file
# with open(output_txt_filename_T1t_to_impt, 'w') as f:
#     f.write(str(np.round(mean_distance_1t3t,4)))
#

with open(output_csv_filename, 'w', newline='') as file:
  writer = csv.writer(file)
  writer.writerow(["Distance (mm) between T2 and T1 transformed to T2", 
                   "Distance (mm) between T2 and Impression transformed to T2",
                   "Distance (mm) between T1 transformed to T2, and Impression transformed to T2"])
  writer.writerow([str(np.round(mean_distance_1t2,4)), 
                   str(np.round(mean_distance_3t2,4)),
                   str(np.round(mean_distance_1t3t,4))])
  file.close() 


#########################################################################
   
# Create renderers
renWin = vtk.vtkRenderWindow()
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

mapper1t.SetInputData(mesh1t)
actor1t.SetMapper(mapper1t)
actor1t.GetProperty().SetColor(0.8,0.4,0.2)

mapperA, mapperB, mapperC = [vtk.vtkPolyDataMapper() for i in range(3)]
actorA, actorB, actorC = [vtk.vtkActor() for i in range(3)]
renA, renB, renC, renD, renE = [vtk.vtkRenderer() for i in range(5)]

mapper3t.SetInputData(mesh3t)
actor3t.SetMapper(mapper3t)
actor3t.GetProperty().SetColor(1.0,0.2,0.4)

mapper3 = vtk.vtkPolyDataMapper()
actor3 = vtk.vtkActor()

mapper3.SetInputData(mesh3)
actor3.SetMapper(mapper3)
actor3.GetProperty().SetColor(1.0,0.2,0.4)


if vtk.VTK_MAJOR_VERSION <= 5:
  mapperA.SetInput(distance_mesh_1t2)
  mapperB.SetInput(distance_mesh_3t2)
  mapperC.SetInput(distance_mesh_1t3t)

else:
    mapperA.SetInputData(distance_mesh_1t2)
    mapperB.SetInputData(distance_mesh_3t2)
    mapperC.SetInputData(distance_mesh_1t3t)

    
    
# Set range between 0 and 3
mapperA.SetScalarRange(0, 2)
mapperB.SetScalarRange(0, 2)
mapperC.SetScalarRange(0, 2)
# Create own look up table with set colors
hueLut = vtk.vtkLookupTable()
hueLut.SetNumberOfTableValues(5)
hueLut.Build()
nc = vtk.vtkNamedColors()
hueLut.SetTableValue(0,nc.GetColor4d("Blue"))
hueLut.SetTableValue(1,nc.GetColor4d("Green"))
hueLut.SetTableValue(2,nc.GetColor4d("Yellow"))
hueLut.SetTableValue(3,nc.GetColor4d("Orange"))
hueLut.SetTableValue(4,nc.GetColor4d("Red"))
mapperA.SetLookupTable(hueLut)
mapperB.SetLookupTable(hueLut)
mapperC.SetLookupTable(hueLut)
# Set the scalar bar
scalarBarA = vtk.vtkScalarBarActor()
scalarBarA.SetLookupTable(hueLut)
scalarBarA.SetTitle("mean distance")
scalarBarA.SetMaximumNumberOfColors(5) # 0.2 between each
scalarBarA.SetNumberOfLabels(6)   
scalarBarB = vtk.vtkScalarBarActor()
scalarBarB.SetLookupTable(hueLut)
scalarBarB.SetTitle("mean distance")
scalarBarB.SetMaximumNumberOfColors(5) # 0.2 between each
scalarBarB.SetNumberOfLabels(6)  
scalarBarC = vtk.vtkScalarBarActor()
scalarBarC.SetLookupTable(hueLut)
scalarBarC.SetTitle("mean distance")
scalarBarC.SetMaximumNumberOfColors(5) # 0.2 between each
scalarBarC.SetNumberOfLabels(6)   

# scalarBarA.SizeTitle(12)
# scalarBarB.SizeTitle(12) 
# scalarBarC.SizeTitle(12)
scalarBarA.SetUnconstrainedFontSize(12)
scalarBarB.SetUnconstrainedFontSize(12) 
scalarBarC.SetUnconstrainedFontSize(12)

actorA.SetMapper(mapperA)
actorB.SetMapper(mapperB)
actorC.SetMapper(mapperC)

# cornerA = vtk.vtkCornerAnnotation()
# cornerA.SetMinimumFontSize(14)
# cornerA.SetMaximumFontSize(14)
# # textProperty = cornerA.GetTextProperty()
# # textProperty.SetColor((0,0,0))
# # corner2.SetText(2,"Before alignment \n T1 (orange) and T2 (green)")
# cornerA.SetText(2,"Before alignment \n scan 1 (orange) and scan 2 (green) and scan 3 impression (red)")
#
# cornerB = vtk.vtkCornerAnnotation()
# cornerB.SetMinimumFontSize(14)
# cornerB.SetMaximumFontSize(14)
# # textProperty = cornerB.GetTextProperty()
# # textProperty.SetColor((0,0,0))
# # corner3.SetText(2,"After alignment \n T1 (orange) and T2 (green)")
# cornerB.SetText(2,"After alignment \n scan 1 (orange) and scan 2 (green) and scan 3 impression (red)")
#
# cornerC = vtk.vtkCornerAnnotation()
# cornerC.SetMinimumFontSize(14)
# cornerC.SetMaximumFontSize(14)
# # textProperty = cornerC.GetTextProperty()
# # textProperty.SetColor((0,0,0))
# cornerC.SetText(2, "Rotate with left mouse \n Move object or box planes by middle mouse button \n Turn on/off box planes by pressing i \n Move planes with left mouse")
#
# cornerD = vtk.vtkCornerAnnotation()
# cornerD.SetMinimumFontSize(14)
# cornerD.SetMaximumFontSize(14)
# # textProperty = cornerD.GetTextProperty()
# # textProperty.SetColor((0,0,0))
# cornerD.SetText(2, "Displaying T1 transformed to T2, and T2")
#
# cornerE = vtk.vtkCornerAnnotation()
# cornerE.SetMinimumFontSize(14)
# cornerE.SetMaximumFontSize(14)
# # textProperty = cornerE.GetTextProperty()
# # textProperty.SetColor((0,0,0))
# cornerE.SetText(2, "Displaying Impression transformed to T2, and T2")
#
# cornerF = vtk.vtkCornerAnnotation()
# cornerF.SetMinimumFontSize(14)
# cornerF.SetMaximumFontSize(14)
# textProperty = cornerF.GetTextProperty()
# textProperty.SetColor((0,0,0))
# cornerF.SetText(2, "Displaying T1 transformed to T2, and impression transformed to T2")

#################
### Rendering ###
#################

### A - before alignment ###
# renA.AddActor2D(cornerA)
# renA.AddViewProp(cornerA)
txtA = vtk.vtkTextActor()
txtA.SetInput("Before alignment \n scan 1 (orange) and scan 2 (green) and scan 3 impression (red)")
renA.AddActor(txtA)
renA.AddActor(actor1)
renA.AddActor(actor2)
renA.AddActor(actor3)
# renA.SetViewport(0.0, 0.5, 0.33, 1)
renA.SetViewport(0.0, 0.5, 0.5, 1)

### B - after alignment ### 
# renB.AddActor2D(cornerB)
# renB.AddViewProp(cornerB)
txtB = vtk.vtkTextActor()
txtB.SetInput("After alignment \n scan 1 (orange) and scan 2 (green) and scan 3 impression (red)")
renB.AddActor(txtB)
renB.AddActor(actor1t)
renB.AddActor(actor2)
renB.AddActor(actor3t)
# renB.SetViewport(0.0, 0.0, 0.33, 0.5)
renB.SetViewport(0.0, 0.0, 0.5, 0.5)

# ### C - middle to clip ### --> change to clip actor later 
# # renC.AddActor(cornerC)
# # renC.AddViewProp(cornerC)
# txtC = vtk.vtkTextActor()
# txtC.SetInput("Rotate with left mouse \n Move object or box planes by middle mouse button \n Turn on/off box planes by pressing i \n Move planes with left mouse")
# renC.AddActor(txtC)
# renC.AddActor(actor1t)
# renC.AddActor(actor2)
# renC.AddActor(actor3t)
# renC.AddActor(clipActor)
# renC.SetViewport(0.33, 0.0, 0.67, 1)

### C - T1 transformed to T2, and T2 ### 
# renC.AddActor(cornerD)
# renC.AddViewProp(cornerD)
txtC = vtk.vtkTextActor()
txtC.SetInput("Displaying distance between T1 transformed to T2, and T2")
renC.AddActor(txtC)
renC.AddActor(actorA)
renC.AddActor2D(scalarBarA)

# create a MAD annotation
corner_MAD = vtk.vtkCornerAnnotation()
corner_MAD.SetMinimumFontSize(12)
corner_MAD.SetMaximumFontSize(12)
corner_MAD.SetText(7, "MAD: " + str(np.round(mean_distance_1t2,4)))
renC.AddActor2D(corner_MAD)

# renC.SetViewport(0.67, 0.67, 1, 1)
renC.SetViewport(0.5, 0.67, 1, 1)

### D - Impression transformed to T2, and T2 ###
# renD.AddActor(cornerD)
# renD.AddViewProp(cornerD)
txtD = vtk.vtkTextActor()
txtD.SetInput("Displaying distance between Impression transformed to T2, and T2")
renD.AddActor(txtD)
renD.AddActor(actorB)
renD.AddActor2D(scalarBarB)

# create a MAD annotation
corner_MAD = vtk.vtkCornerAnnotation()
corner_MAD.SetMinimumFontSize(12)
corner_MAD.SetMaximumFontSize(12)
corner_MAD.SetText(7, "MAD: " + str(np.round(mean_distance_3t2,4)))
renD.AddActor2D(corner_MAD)

# renD.SetViewport(0.67, 0.33, 1, 0.67)
renD.SetViewport(0.5, 0.33, 1, 0.67)

### E - T1 transformed to T2, and Impression transformed to T2 ###
# renE.AddActor(cornerE) 
# renE.AddViewProp(cornerE)
txtE = vtk.vtkTextActor()
txtE.SetInput("Displaying distance between T1 transformed to T2, and Impression transformed to T2")
renE.AddActor(txtE) 
renE.AddActor(actorC)
renE.AddActor2D(scalarBarC)

# create a MAD annotation
corner_MAD = vtk.vtkCornerAnnotation()
corner_MAD.SetMinimumFontSize(12)
corner_MAD.SetMaximumFontSize(12)
corner_MAD.SetText(7, "MAD: " + str(np.round(mean_distance_1t3t,4)))
renE.AddActor2D(corner_MAD)

# renF.SetViewport(0.67, 0.0, 1.0, 0.33)
renE.SetViewport(0.5, 0.0, 1.0, 0.33)

# renA.SetBackground(0.05, 0.05, 0.05)
# renB.SetBackground(0.1, 0.1, 0.1)
# renC.SetBackground(0.15, 0.15, 0.15)
# renD.SetBackground(0.1, 0.1, 0.1)
# renE.SetBackground(0.1, 0.1, 0.1)
# renF.SetBackground(0.1, 0.1, 0.1)
renA.SetBackground(0.05, 0.05, 0.05)
renB.SetBackground(0.1, 0.1, 0.1)
renC.SetBackground(0.15, 0.15, 0.15)
renD.SetBackground(0.05, 0.05, 0.05)
renE.SetBackground(0.1, 0.1, 0.1)

renWin.AddRenderer(renA)
renWin.AddRenderer(renB)
renWin.AddRenderer(renC)
renWin.AddRenderer(renD)
renWin.AddRenderer(renE)
# renWin.AddRenderer(renF)
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

iren.Initialize()

iren.Start()

# sys.exit(app.exec_())
