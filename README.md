# BrainFEM
Complete finite element code in Cython, extended for assigning 'mechanisms' to individual elements or groups of elements. Initially created to model the effect of molecular mechanisms in the axon.\
Note, this code is based on Python 2.7

#-------------------------------------#\
Directory structure\
#-------------------------------------#\
pyxd/                     Contains all .pyx and .pxd files for a standard FEM code, such as objects for Nodes, Elements, SPC, MPC, Solver, etc.\
pyxdX/                    Contains all .pyx and .pxd files that extend the standard FEM code to adding molecular mechanisms.\
All other directories     Each directory represents a different model. These directories contain files with model description, model parameters, and post processing files.\


#-------------------------------------#\
Getting started (for MacOS)\
#-------------------------------------#
1) Download and install anaconda: https://www.anaconda.com/download/#macos\
2) Download and install paraview: https://www.paraview.org/download/\
3) Install vtk by running from the terminal: conda install -c anaconda vtk \
4) Test that Cython is working correctly by running Main.py in the TestCython directory. From terminal, run: python Main.py\
5) Test that standard FEM is working by running Main_BarTest_B3.py in the BarTest_B3 directory. Note, this is a simple finite element model of a single bar element that is rotated and stretched.\
6) Test that the extension of FEM is working by running Main_DyneinModel_B3.py in DyneinModel_B3 directory. This is a simple model of dynein molecules pushing a microtubules forward.\

#-------------------------------------#\
Basic info\
#-------------------------------------#\
Every time a model is ran, a subdirectory called RES_*** . The subdirectory contains a Input/, Output/, and PostProcess/ directory.\
The Input/ directory contains copies of all input files that can later be used to reproduce the results.\
The Output/ directory contains a log.txt file that contains all terminal display output. It further contains of optional, user defined files that are model dependent and can be coded in the virtual WriteStepOutput(self, Solver s) function.\
The PostProcess/ directory is intended to save all Paraview output files and other output files, such as figures or videos.
