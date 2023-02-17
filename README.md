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

To run on a mac without git: 
1. Download and unzip the zip file 
2. Open a terminal in the directory above 
3. Type in the terminal and run python3 -m venv shape_analysis_dental_env. This creates a python environment that we will use without affecting other installations on your computer.
4. To activate your environment you created, run source shape_analysis_dental_env/bin/activate. You should now see something like (shape_analysis_dental_env) at the beginning of your terminal prompt. 
5. To install the necessary packages, run pip install vtk, and pip install PyQt5
6. To the run python code, type python find_distance_dental_3_files_no_cropping_planes.py
