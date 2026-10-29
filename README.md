# Files used in the study: 
Tectonic memory and reactant depletion regulate serpentinization-driven hydrogen generation

Yuan Wang, Zhonglan Liu, Gianreto Manatschal, Shuming Feng, Wenliang Xu


Model runs in this study used dealii 9.2.0.

The software version used for this study can be found at: 

ASPECT 2.3

        https://github.com/Djneu/aspect/tree/fault_analysis
	
FastScape commit 18f2588

(elimated)

        https://github.com/fastscape-lem/fastscapelib-fortran 
	
(backup)

        https://github.com/Jlu-geodynamic/fastscape_backup_commit18f2588

All source code used in this study can be found in the backup.


# Additional ASPECT plugins statement:

1."comp" and "rift" were modified from the plugins used in the study: 

Neuharth, D. et al. Evolution of Rift Systems and Their Fault Networks in Response to Surface Processes. Tectonics 41, e2021TC007166 (2022). https://doi.org/https://doi.org/10.1029/2021TC007166


We change some initial setting to fit our study.
Note: The initial topography function of this part is under development, it is not used in this study. 

2."grain" were modified from the sorce code. 

# Additional Perple_X plugins statement:

the input file 20260310_basalt.dat was modified from the study: 

Klein, B. Z. & Behn, M. D. The Impact of Mélange Formation on the Rheology and H2O Budget of the Subduction Zone Interface, With Implications for Diapir Nucleation. Geophysical Research Letters 52, e2024GL113049 (2025). https://doi.org/https://doi.org/10.1029/2024GL113049


# To run the model and visualize the results:

1.Create a build directory for fastscape and compile it with an added flag for creating a shared library.

        cmake -DBUILD_FASTSCAPELIB_SHARED=ON /path/to/fastscapemake
        make

2.Compile ASPECT with a flag FASTSCAPE_DIR pointing to the fastscape build folder with the shared library

        cmake -DFASTSCAPE_DIR=/path/to/fastscape/build -DDEAL_II_DIR=/path/to/dealii -DCMAKE_BUILD_TYPE=Release /path/to/aspect
        make

3. Using ASPECT to compile plugin. 
   Note: the path of the '#include <files>' in .h and .cc should be modified.

        cmake -DAspect_DIR=/path/to/ASPECT/build ..
        make

4. Change the path of the plugin in the prm files and run.

        mpirun -n [ ] /path/to/aspect/build /path/to/prm/file

    Note 1: the results of Model M4 were calculated in paraview based on '00_Model_M4_serp_calculate.pvsm' and the output files of Model M1.
    Note 2: we recommend to change the input 'set Basalt material file name' as the absolute path of '20260310_basalt.tab'.
    Note 3: each model has 5 prm files, they should run in order of stage. 

5. Use Paraview to visualize results by the .pvsm file.

