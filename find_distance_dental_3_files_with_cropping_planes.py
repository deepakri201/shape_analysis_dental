# find_distance_dental_3_files.py
# 
# This script performs ICP (iterative closest point) registration between two sets of STL
# files, the T1, T2 and the impression. T1 is reigstered to T2, and the impression is also 
# registered to T2. 
# The mean absolute distance is calculated and displayed between the two sets of meshes.
# A screenshot and the distance mesh in the vtk format are currently saved to disk for each
# registration. 
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
# September 2022
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

output_mesh_filename_T2_to_T1t = 'distance_mesh_T2_to_T1t.vtk'
output_txt_filename_T2_to_T1t = 'mean_absolute_distance_in_mm_T2_to_T1t.txt'
output_mesh_filename_T2_to_impt = 'distance_mesh_T2_to_impt.vtk'
output_txt_filename_T2_to_impt = 'mean_absolute_distance_in_mm_T2_to_impt.txt'
output_mesh_filename_T1t_to_impt = 'distance_mesh_T1t_to_impt.vtk'
output_txt_filename_T1t_to_impt = 'mean_absolute_distance_in_mm_T1t_to_impt.txt'

output_jpg_filename = 'distance_mesh.jpg'

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

# distance_meshb_3, distancesb_3, mean_distanceb_3 = calculate_MAD(mesh3t, mesh2)
# dist_min_3 = distance_mesh_3.GetPointData().GetScalars().GetRange()[0]
# dist_minb_3 = distance_meshb_3.GetPointData().GetScalars().GetRange()[0]
# dist_max_3 = distance_mesh_3.GetPointData().GetScalars().GetRange()[1]
# dist_maxb_3 = distance_meshb_3.GetPointData().GetScalars().GetRange()[1]

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
with open(output_txt_filename_T2_to_T1t, 'w') as f:
    f.write(str(np.round(mean_distance_1t2,4)))
    
# save the MAD to a text file
with open(output_txt_filename_T2_to_impt, 'w') as f:
    f.write(str(np.round(mean_distance_3t2,4)))
    
# save the MAD to a text file
with open(output_txt_filename_T1t_to_impt, 'w') as f:
    f.write(str(np.round(mean_distance_1t3t,4)))
    
#########################################################################
   

# Create renderers
renWin = vtk.vtkRenderWindow()
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

mapper1t.SetInputData(mesh1t)
actor1t.SetMapper(mapper1t)
actor1t.GetProperty().SetColor(0.8,0.4,0.2)

# mapper, mapperb = [vtk.vtkPolyDataMapper() for i in range(2)]
# actor, actorb = [vtk.vtkActor() for i in range(2)]
# ren1, ren1b, ren2, ren3 = [vtk.vtkRenderer() for i in range(4)]
# mapper = vtk.vtkPolyDataMapper()
# actor = vtk.vtkActor()
mapperA, mapperB, mapperC = [vtk.vtkPolyDataMapper() for i in range(3)]
actorA, actorB, actorC = [vtk.vtkActor() for i in range(3)]
renA, renB, renC, renD, renE, renF = [vtk.vtkRenderer() for i in range(6)]

mapper3t.SetInputData(mesh3t)
actor3t.SetMapper(mapper3t)
actor3t.GetProperty().SetColor(1.0,0.2,0.4)

# mapper3, mapper3b = [vtk.vtkPolyDataMapper() for i in range(2)]
# actor3, actor3b = [vtk.vtkActor() for i in range(2)]
mapper3 = vtk.vtkPolyDataMapper()
actor3 = vtk.vtkActor()

mapper3.SetInputData(mesh3)
actor3.SetMapper(mapper3)
actor3.GetProperty().SetColor(1.0,0.2,0.4)

# mapper3b.SetInputData(mesh3)
# actor3b.SetMapper(mapper3b)
# actor3b.GetProperty().SetColor(1.0,0.2,0.4)


if vtk.VTK_MAJOR_VERSION <= 5:
  mapperA.SetInput(distance_mesh_1t2)
  mapperB.SetInput(distance_mesh_3t2)
  mapperC.SetInput(distance_mesh_1t3t)
    # mapper.SetInput(distance_mesh_1t2)
    # mapperb.SetInput(distance_mesh_3)
    
    # mapperb.SetInput(distance_meshb)
    # mapper3.SetInput(distance_mesh_3)
    # mapper3b.SetInput(distance_meshb_3)
else:
    mapperA.SetInputData(distance_mesh_1t2)
    mapperB.SetInputData(distance_mesh_3t2)
    mapperC.SetInputData(distance_mesh_1t3t)
    # mapper.SetInputData(distance_mesh_1t2)
    # mapperb.SetInputData(distance_mesh_3)
    
    # mapperb.SetInputData(distance_meshb)
    # mapper3.SetInputData(distance_mesh_3)
    # mapper3b.SetInputData(distance_meshb_3)
    
    
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

actorA.SetMapper(mapperA)
actorB.SetMapper(mapperB)
actorC.SetMapper(mapperC)

# # Add both actors to renderer
# renD.AddActor2D(scalarBarA)
# renE.AddActor2D(scalarBarB)
# renF.AddActor2D(scalarBarC)

# # Set range between 0 and 3
# # mapper.SetScalarRange(0, 3)
# # mapper.SetScalarRange(0, 1)
# # mapper.SetScalarRange(0, 1.6)
# mapper.SetScalarRange(0, 2)
# # mapper.SetScalarRange(0, 2.4)
#
# # Create own look up table with set colors
# hueLut = vtk.vtkLookupTable()
# # hueLut.SetNumberOfTableValues(6)
# hueLut.SetNumberOfTableValues(5)
# hueLut.Build()
# nc = vtk.vtkNamedColors()
# # hueLut.SetTableValue(0,nc.GetColor4d("Blue"))
# # hueLut.SetTableValue(1,nc.GetColor4d("Purple"))
# # hueLut.SetTableValue(2,nc.GetColor4d("Green"))
# # hueLut.SetTableValue(3,nc.GetColor4d("Yellow"))
# # hueLut.SetTableValue(4,nc.GetColor4d("Orange"))
# # hueLut.SetTableValue(5,nc.GetColor4d("Red"))
# hueLut.SetTableValue(0,nc.GetColor4d("Blue"))
# hueLut.SetTableValue(1,nc.GetColor4d("Green"))
# hueLut.SetTableValue(2,nc.GetColor4d("Yellow"))
# hueLut.SetTableValue(3,nc.GetColor4d("Orange"))
# hueLut.SetTableValue(4,nc.GetColor4d("Red"))
#
# mapper.SetLookupTable(hueLut)
#
# # Set the scalar bar
# scalarBar = vtk.vtkScalarBarActor()
# scalarBar.SetLookupTable(hueLut)
# scalarBar.SetTitle("mean distance")
# # scalarBar.SetMaximumNumberOfColors(6); # 0.5 between each
# # scalarBar.SetNumberOfLabels(7);
# scalarBar.SetMaximumNumberOfColors(5) # 0.2 between each
# scalarBar.SetNumberOfLabels(6)




# # Set range between 0 and 3 for b
# # mapperb.SetScalarRange(0, 3)
# # mapperb.SetScalarRange(0, 1)
# # mapperb.SetScalarRange(0, 1.6)
# mapperb.SetScalarRange(0, 2)
# # mapperb.SetScalarRange(0, 2.4)
#
# # Create own look up table with set colors for b
# hueLutb = vtk.vtkLookupTable()
# # hueLutb.SetNumberOfTableValues(6)
# hueLutb.SetNumberOfTableValues(5)
# hueLutb.Build()
# nc = vtk.vtkNamedColors()
# # hueLutb.SetTableValue(0,nc.GetColor4d("Blue"))
# # hueLutb.SetTableValue(1,nc.GetColor4d("Purple"))
# # hueLutb.SetTableValue(2,nc.GetColor4d("Green"))
# # hueLutb.SetTableValue(3,nc.GetColor4d("Yellow"))
# # hueLutb.SetTableValue(4,nc.GetColor4d("Orange"))
# # hueLutb.SetTableValue(5,nc.GetColor4d("Red"))
# hueLutb.SetTableValue(0,nc.GetColor4d("Blue"))
# hueLutb.SetTableValue(1,nc.GetColor4d("Green"))
# hueLutb.SetTableValue(2,nc.GetColor4d("Yellow"))
# hueLutb.SetTableValue(3,nc.GetColor4d("Orange"))
# hueLutb.SetTableValue(4,nc.GetColor4d("Red"))
#
# mapperb.SetLookupTable(hueLutb)
#
# # Set the scalar bar b
# scalarBarb = vtk.vtkScalarBarActor()
# scalarBarb.SetLookupTable(hueLutb)
# scalarBarb.SetTitle("mean distance")
# # scalarBarb.SetMaximumNumberOfColors(6); # 0.5 between each
# # scalarBarb.SetNumberOfLabels(7);
# scalarBarb.SetMaximumNumberOfColors(5); # 0.2 between each
# scalarBarb.SetNumberOfLabels(6);

# Add both actors to renderer
# ren1.AddActor2D(scalarBar)
# ren1b.AddActor2D(scalarBarb)


################################# Clipping planes ##################
### For mesh 2 and mesh1t ### 

normals = vtk.vtkPolyDataNormals()
normals.SetInputData(mesh2)
normals.Update()

points = vtk_to_numpy(mesh2.GetPoints().GetData())
bounds_min = np.min(points,axis=0)
bounds_max = np.max(points,axis=0)
planes = vtk.vtkPlanes()
planes.SetBounds(bounds_min[0], bounds_max[0], bounds_min[1], bounds_max[1], bounds_min[2], bounds_max[2])

clipper = vtk.vtkClipPolyData()
clipper.SetInputData(mesh2)
clipper.SetClipFunction(planes)
clipper.InsideOutOn()
# clipper.GenerateClippedOutputOn()
clipper.SetValue(0)
clipper.Update()

clipMapper = vtk.vtkPolyDataMapper()
clipMapper.SetInputConnection(clipper.GetOutputPort())
clipMapper.ScalarVisibilityOff()
# clipMapper.ScalarVisibilityOn()
# clipMapper.SetScalarRange(0, 3)
# clipMapper.SetScalarRange(0, 1)
# clipMapper.SetScalarRange(0, 1.6)
# clipMapper.SetScalarRange(0, 2)
# clipMapper.SetScalarRange(0, 2.4)
# clipMapper.SetLookupTable(hueLut)

clipActor = vtk.vtkActor()
clipActor.SetMapper(clipMapper)

# # create a MAD annotation
# corner_MAD = vtk.vtkCornerAnnotation()
# corner_MAD.SetMinimumFontSize(16)
# corner_MAD.SetMaximumFontSize(16)
# corner_MAD.SetText(7, "MAD: " + str(np.round(mean_distance,4)))

boxWidget = vtk.vtkBoxWidget()
boxWidget.SetCurrentRenderer(renC)
boxWidget.SetInteractor(iren)
boxWidget.SetPlaceFactor(1.25)

def myCallback_central_mesh2(obj, event):
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
    # global mesh1t
    # global distance_mesh

    clipper_mesh2 = vtk.vtkClipPolyData()
    clipper_mesh2.SetInputData(mesh2)
    clipper_mesh2.SetClipFunction(planes)
    # clipper_mesh2.InsideOutOn()
    
    # clipper_mesh2.GenerateClippedOutputOn()
    clipper_mesh2.SetValue(0)

    
    clipper_mesh2.Update()
    mesh2_clipped = vtk.vtkPolyData()
    mesh2_clipped = clipper_mesh2.GetClippedOutput()

    # clipper_mesh1t = vtk.vtkClipPolyData()
    # clipper_mesh1t.SetInputData(mesh1t)
    # clipper_mesh1t.SetClipFunction(planes)
    # # clipper_mesh1t.InsideOutOn()
    # clipper_mesh1t.GenerateClippedOutputOn()
    # clipper_mesh1t.Update()
    # mesh1t_clipped = vtk.vtkPolyData()
    # mesh1t_clipped = clipper_mesh1t.GetClippedOutput()
    #
    # distance_mesh, distances, mean_distance = calculate_MAD(mesh2_clipped, mesh1t_clipped)
    #
    # global corner_MAD
    # corner_MAD.SetMinimumFontSize(14)
    # corner_MAD.SetMaximumFontSize(14)
    # corner_MAD.SetText(7, "MAD: " + str(np.round(mean_distance,4)))
    #
    # # save distance mesh that has been clipped
    # polydata_writer = vtk.vtkPolyDataWriter()
    # polydata_writer.SetFileName(output_mesh_filename_1)
    # polydata_writer.SetInputData(distance_mesh)
    # polydata_writer.Write()
    #
    # # save the MAD to a text file
    # with open(output_txt_filename_1, 'w') as f:
    #     f.write(str(np.round(mean_distance,4)))


# Interaction event
boxWidget.SetInputData(mesh2)
boxWidget.PlaceWidget()
boxWidget.On()
boxWidget.AddObserver("EndInteractionEvent", myCallback_central_mesh2) # should be faster compared to "InteractionEvent"


#####################################################
### Set up cutting box - 3 planes and the clipper ###
#####################################################

# ### For mesh 2 and mesh1t ### 
#
# normals = vtk.vtkPolyDataNormals()
# normals.SetInputData(distance_mesh)
# normals.Update()
#
# points = vtk_to_numpy(distance_mesh.GetPoints().GetData())
# bounds_min = np.min(points,axis=0)
# bounds_max = np.max(points,axis=0)
# planes = vtk.vtkPlanes()
# planes.SetBounds(bounds_min[0], bounds_max[0], bounds_min[1], bounds_max[1], bounds_min[2], bounds_max[2])
#
# clipper = vtk.vtkClipPolyData()
# clipper.SetInputData(distance_mesh)
# clipper.SetClipFunction(planes)
# clipper.InsideOutOn()
# clipper.GenerateClippedOutputOn()
#
# clipMapper = vtk.vtkPolyDataMapper()
# clipMapper.SetInputConnection(clipper.GetOutputPort())
# clipMapper.ScalarVisibilityOn()
# # clipMapper.SetScalarRange(0, 3)
# # clipMapper.SetScalarRange(0, 1)
# # clipMapper.SetScalarRange(0, 1.6)
# clipMapper.SetScalarRange(0, 2)
# # clipMapper.SetScalarRange(0, 2.4)
# clipMapper.SetLookupTable(hueLut)
#
# clipActor = vtk.vtkActor()
# clipActor.SetMapper(clipMapper)
#
# # create a MAD annotation
# corner_MAD = vtk.vtkCornerAnnotation()
# corner_MAD.SetMinimumFontSize(16)
# corner_MAD.SetMaximumFontSize(16)
# corner_MAD.SetText(7, "MAD: " + str(np.round(mean_distance,4)))
#
# boxWidget = vtk.vtkBoxWidget()
# boxWidget.SetCurrentRenderer(ren1)
# boxWidget.SetInteractor(iren)
# boxWidget.SetPlaceFactor(1.25)
#
# ### For mesh2 and mesh3t ###
#
# normals_3 = vtk.vtkPolyDataNormals()
# normals_3.SetInputData(distance_mesh_3)
# normals_3.Update()
#
# points_3 = vtk_to_numpy(distance_mesh_3.GetPoints().GetData())
# bounds_min_3 = np.min(points_3,axis=0)
# bounds_max_3 = np.max(points_3,axis=0)
# planes_3 = vtk.vtkPlanes()
# planes_3.SetBounds(bounds_min_3[0], bounds_max_3[0], bounds_min_3[1], bounds_max_3[1], bounds_min_3[2], bounds_max_3[2])
#
# clipper_3 = vtk.vtkClipPolyData()
# clipper_3.SetInputData(distance_mesh_3)
# clipper_3.SetClipFunction(planes_3)
# clipper_3.InsideOutOn()
# clipper_3.GenerateClippedOutputOn()
#
# clipMapper_3 = vtk.vtkPolyDataMapper()
# clipMapper_3.SetInputConnection(clipper_3.GetOutputPort())
# clipMapper_3.ScalarVisibilityOn()
# # clipMapper.SetScalarRange(0, 3)
# # clipMapper.SetScalarRange(0, 1)
# # clipMapper.SetScalarRange(0, 1.6)
# clipMapper_3.SetScalarRange(0, 2)
# # clipMapper.SetScalarRange(0, 2.4)
# clipMapper_3.SetLookupTable(hueLut)
#
# clipActor_3 = vtk.vtkActor()
# clipActor_3.SetMapper(clipMapper_3)
#
# # create a MAD annotation
# corner_MAD_3 = vtk.vtkCornerAnnotation()
# corner_MAD_3.SetMinimumFontSize(16)
# corner_MAD_3.SetMaximumFontSize(16)
# corner_MAD_3.SetText(7, "MAD: " + str(np.round(mean_distance_3,4)))
#
# boxWidget_3 = vtk.vtkBoxWidget()
# boxWidget_3.SetCurrentRenderer(ren1b)
# boxWidget_3.SetInteractor(iren)
# boxWidget_3.SetPlaceFactor(1.25)

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
    polydata_writer.SetFileName(output_mesh_filename_1)
    polydata_writer.SetInputData(distance_mesh)
    polydata_writer.Write()

    # save the MAD to a text file
    with open(output_txt_filename_1, 'w') as f:
        f.write(str(np.round(mean_distance,4)))
        
def myCallback_3(obj, event):
    '''Callback function for the plane'''

    global clipActor_3, planes_3
    obj.GetPlanes(planes_3)
    clipActor_3.VisibilityOn()

    ### Update the distance mesh here by the new clip plane ###
    # Don't need to update scalar bars since they are set always
    # from 0 to 3.
    # Clip mesh2 and mesh1t using the same clip plane and
    # calculate a new distance mesh

    global mesh2
    global mesh3t
    global distance_mesh_3

    clipper_mesh2_3 = vtk.vtkClipPolyData()
    clipper_mesh2_3.SetInputData(mesh2)
    clipper_mesh2_3.SetClipFunction(planes_3)
    # clipper_mesh2.InsideOutOn()
    clipper_mesh2_3.GenerateClippedOutputOn()
    clipper_mesh2_3.Update()
    mesh2_clipped_3 = vtk.vtkPolyData()
    mesh2_clipped_3 = clipper_mesh2_3.GetClippedOutput()

    clipper_mesh3t = vtk.vtkClipPolyData()
    clipper_mesh3t.SetInputData(mesh3t)
    clipper_mesh3t.SetClipFunction(planes_3)
    # clipper_mesh1t.InsideOutOn()
    clipper_mesh3t.GenerateClippedOutputOn()
    clipper_mesh3t.Update()
    mesh3t_clipped = vtk.vtkPolyData()
    mesh3t_clipped = clipper_mesh3t.GetClippedOutput()

    distance_mesh_3, distances_3, mean_distance_3 = calculate_MAD(mesh2_clipped, mesh3t_clipped)

    global corner_MAD_3
    corner_MAD_3.SetMinimumFontSize(14)
    corner_MAD_3.SetMaximumFontSize(14)
    corner_MAD_3.SetText(7, "MAD: " + str(np.round(mean_distance_3,4)))

    # save distance mesh that has been clipped
    polydata_writer = vtk.vtkPolyDataWriter()
    polydata_writer.SetFileName(output_mesh_filename_3)
    polydata_writer.SetInputData(distance_mesh_3)
    polydata_writer.Write()

    # save the MAD to a text file
    with open(output_txt_filename_3, 'w') as f:
        f.write(str(np.round(mean_distance_3,4)))



# # Interaction event
# boxWidget.SetInputData(distance_mesh)
# boxWidget.PlaceWidget()
# boxWidget.On()
# boxWidget.AddObserver("EndInteractionEvent", myCallback) # should be faster compared to "InteractionEvent"

# # Interaction event
# boxWidget_3.SetInputData(distance_mesh_3)
# boxWidget_3.PlaceWidget()
# boxWidget_3.On()
# boxWidget_3.AddObserver("EndInteractionEvent", myCallback_3) # should be faster compared to "InteractionEvent"

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
# picker = vtk.vtkCellPicker()
# create a text actor
# txt1 = vtk.vtkTextActor()
# txt3 = vtk.vtkTextActor()

# def annotate_pick(obj, event):
#     '''Callback function for displaying the scalar value of a clicked point'''
#
#     global ren1
#     global iren
#     global txt_actor
#     global distance_mesh
#
#     # pointPicker = vtk.vtkPointPicker()
#     cellPicker = vtk.vtkCellPicker()
#     eventPosition = obj.GetEventPosition()
#     result = picker.Pick(float(eventPosition[0]),float(eventPosition[1]),0.0,ren1)
#     pickId = picker.GetCellId()
#
#     # Get scalar that is associated with this pickId
#     scalars = vtk_to_numpy(distance_mesh.GetCellData().GetScalars())
#     dist_value = scalars[pickId]
#
#     # Set as input to text actor
#     txt1.SetInput("Picked point distance value: " + str(np.round(dist_value,4)))
#
#
# iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
# iren.AddObserver("RightButtonPressEvent", annotate_pick)


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


# # Text actor properties
# txtprop=txt1.GetTextProperty()
# txtprop.SetFontFamilyToArial()
# txtprop.SetFontSize(14)
# txtprop.SetColor(1,1,1)
# txt1.SetTextProperty(txtprop)
# txt1.SetPosition2(20,30)

# corner1 = vtk.vtkCornerAnnotation()
# corner1.SetMinimumFontSize(14)
# corner1.SetMaximumFontSize(14)
# corner1.SetText(2, "Mean absolute distance \n Rotate with left mouse \n Move object or box planes by middle mouse button \n Turn on/off box planes by pressing i \n Move planes with left mouse \n Select cell on condyle with right mouse")

# corner2 = vtk.vtkCornerAnnotation()
# corner2.SetMinimumFontSize(14)
# corner2.SetMaximumFontSize(14)
# # corner2.SetText(2,"Before alignment \n T1 (orange) and T2 (green)")
# corner2.SetText(2,"Before alignment \n scan 1 (orange) and scan 2 (green) and scan 3 impression (red)")

# corner3 = vtk.vtkCornerAnnotation()
# corner3.SetMinimumFontSize(14)
# corner3.SetMaximumFontSize(14)
# # corner3.SetText(2,"After alignment \n T1 (orange) and T2 (green)")
# corner3.SetText(2,"After alignment \n scan 1 (orange) and scan 2 (green) and scan 3 impression (red)")

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
renA.SetViewport(0.0, 0.5, 0.33, 1)

### B - after alignment ### 
# renB.AddActor2D(cornerB)
# renB.AddViewProp(cornerB)
txtB = vtk.vtkTextActor()
txtB.SetInput("After alignment \n scan 1 (orange) and scan 2 (green) and scan 3 impression (red)")
renB.AddActor(txtB)
renB.AddActor(actor1t)
renB.AddActor(actor2)
renB.AddActor(actor3t)
renB.SetViewport(0.0, 0.0, 0.33, 0.5)

### C - middle to clip ### --> change to clip actor later 
# renC.AddActor(cornerC)
# renC.AddViewProp(cornerC)
txtC = vtk.vtkTextActor()
txtC.SetInput("Rotate with left mouse \n Move object or box planes by middle mouse button \n Turn on/off box planes by pressing i \n Move planes with left mouse")
renC.AddActor(txtC)
renC.AddActor(actor1t)
renC.AddActor(actor2)
renC.AddActor(actor3t)
renC.AddActor(clipActor)
renC.SetViewport(0.33, 0.0, 0.67, 1)

### D - T1 transformed to T2, and T2 ### 
# renD.AddActor(cornerD)
# renD.AddViewProp(cornerD)
txtD = vtk.vtkTextActor()
txtD.SetInput("Displaying T1 transformed to T2, and T2")
renD.AddActor(txtD)
renD.AddActor(actorA)
renD.AddActor2D(scalarBarA)

# create a MAD annotation
corner_MAD = vtk.vtkCornerAnnotation()
corner_MAD.SetMinimumFontSize(12)
corner_MAD.SetMaximumFontSize(12)
corner_MAD.SetText(7, "MAD: " + str(np.round(mean_distance_1t2,4)))
renD.AddActor2D(corner_MAD)

renD.SetViewport(0.67, 0.67, 1, 1)

### E - Impression transformed to T2, and T2 ###
# renE.AddActor(cornerE)
# renE.AddViewProp(cornerE)
txtE = vtk.vtkTextActor()
txtE.SetInput("Displaying Impression transformed to T2, and T2")
renE.AddActor(txtE)
renE.AddActor(actorB)
renE.AddActor2D(scalarBarB)

# create a MAD annotation
corner_MAD = vtk.vtkCornerAnnotation()
corner_MAD.SetMinimumFontSize(12)
corner_MAD.SetMaximumFontSize(12)
corner_MAD.SetText(7, "MAD: " + str(np.round(mean_distance_3t2,4)))
renE.AddActor2D(corner_MAD)

renE.SetViewport(0.67, 0.33, 1, 0.67)

### F - T1 transformed to T2, and Impression transformed to T2 ###
# renF.AddActor(cornerF) 
# renF.AddViewProp(cornerF)
txtF = vtk.vtkTextActor()
txtF.SetInput("Displaying T1 transformed to T2, and impression transformed to T2")
renF.AddActor(txtF) 
renF.AddActor(actorC)
renF.AddActor2D(scalarBarC)

# create a MAD annotation
corner_MAD = vtk.vtkCornerAnnotation()
corner_MAD.SetMinimumFontSize(12)
corner_MAD.SetMaximumFontSize(12)
corner_MAD.SetText(7, "MAD: " + str(np.round(mean_distance_1t3t,4)))
renF.AddActor2D(corner_MAD)

renF.SetViewport(0.67, 0.0, 1.0, 0.33)

renA.SetBackground(0.05, 0.05, 0.05)
renB.SetBackground(0.1, 0.1, 0.1)
renC.SetBackground(0.15, 0.15, 0.15)
renD.SetBackground(0.1, 0.1, 0.1)
renE.SetBackground(0.1, 0.1, 0.1)
renF.SetBackground(0.1, 0.1, 0.1)

renWin.AddRenderer(renA)
renWin.AddRenderer(renB)
renWin.AddRenderer(renC)
renWin.AddRenderer(renD)
renWin.AddRenderer(renE)
renWin.AddRenderer(renF)
renWin.Render()

renWin.SetSize(1920,1080)
renWin.Render()


# ### upper right - ren1 - Clip 1t and 2 ### 
# ren1.AddActor(clipActor)
# ren1.AddActor2D(txt1)
# ren1.AddActor2D(corner1)
# ren1.AddActor2D(corner_MAD)
# ren1.SetViewport(0.5, 0.5, 1.0, 1.0)

### (previously) Same as above 1 and 2 from other side ### 
# actorb.SetMapper(mapperb)
# ren1b.AddActor(actorb)

# actor3b.SetMapper(mapper3b)
# ren1b.AddActor(actor3b)

# ### lower right - ren1b - Clip 2 and 3t ###
# ren1b.AddActor(clipActor_3) 
# ren1b.AddActor(txt3)
# ren1b.AddActor2D(corner1)
# ren1b.AddActor2D(corner_MAD_3)

# ren1b.AddActor(actor2)
# ren1b.AddActor(actor3t)
# ren1b.SetViewport(0.5, 0.0, 1.0, 0.5)

# #### upper left - ren2 - The original ### 
# ren2.AddActor(actor1)
# ren2.AddActor(actor2)
# ren2.AddActor(actor3)
# ren2.AddActor2D(corner2)
# ren2.SetViewport(0.0, 0.5, 0.5, 1.0)

# ### lower left - ren3 - The transformed ### 
# ren3.AddActor(actor1t)
# ren3.AddActor(actor2)
# ren3.AddActor(actor3t)
# ren3.AddActor2D(corner3)
# ren3.SetViewport(0.0, 0.0, 0.5, 0.5)

# ren1b.SetBackground(0.05, 0.05, 0.05)
# ren3.SetBackground(0.1, 0.1, 0.1)
# ren2.SetBackground(0.15, 0.15, 0.15)

# renWin.AddRenderer(ren1)
# renWin.AddRenderer(ren1b)
# renWin.AddRenderer(ren2)
# renWin.AddRenderer(ren3)
# renWin.Render()
#
# renWin.SetSize(1920,1080)
# renWin.Render()

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
