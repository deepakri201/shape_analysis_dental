This repository contains code for comparing T1, T2 and impression dental scans. The T1 and impression are both registered to T2 space. 
The mean absolute distance (MAD) is calculated between the three sets of meshes: 
1. T2 and T1 registered to T2 
2. T2 and Impression registered to T2 
3. T1 regsitered to T2 and Impression registered to T2 
The output directory ("results") contains the following: 
1. distance_mesh_T2_to_T1t.vtk - MAD mesh between T2 and T1 registered to T2
2. distance_mesh_T2_to_impt.vtk - MAD mesh between T2 and Impression registered to T2 
3. distance_mesh_T1t_to_impt.vtk - MAD mesh between T1 registered to T2 and Impression registered to T2 
4. mean_absolute_distance_in_mm.csv - csv file containing the distances in mm for the above three comparisons 
3. distance_mesh - Screenshot of the meshes 

The function find_distance_dental_3_files_no_cropping_planes.py is the one that should be run. The one with cropping planes is a WIP. 
