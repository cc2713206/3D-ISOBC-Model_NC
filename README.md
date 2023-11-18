1. System requirements
Our 3D ISOBC model V1.0 was developed on Windows 10 Professional Edition x64. The source codes, including three programs of main.f90, param.f90, and sub.f90, were written in Fortran. So Visual Studio and Fortran Compiler need to be installed first for the development of the 3D ISOBC model.
2. Installation guide
First, install Visual Studio, and its latest edition (Professional 2022) can be freely downloaded from the website: https://visualstudio.microsoft.com/.
Second, install Fortran Compiler, including two software packages: Intel® oneAPI Base Toolkit and Intel® HPC Toolkit, both of which can be freely downloaded from the website: https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#base-kit. Install Intel® oneAPI Base Toolkit first, then install the other one.
3. Demo
The present three programs are for the reference run (i.e., H140W8), and they can be directly run together in Visual Studio (see the next statement). The output information is given in sub.f90, and all the corresponding output data can be directly obtained after about 2 days run time in the path you set.
4. Instructions for use
First, run Visual Studio; second, create a project based on Fortran; third, name the project, set its path, and create it; forth, load the three programs (i.e., main.f90, param.f90, and sub.f90) in the Source Files of your project; Fifth, set the Release solution, compile the codes, and run it without debug.
# 3D-ISOBC-Model_NC
