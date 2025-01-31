# Influence of Jupiter-like Planets on Centaur Transits to the Inner Solar System
This repository is for the course "Project Computational Science" at the University of Amsterdam.
For a quick overview of the project we refer to (poster pdf location). For a more detailled description of
the project, we refer to (report-file location).

#### Structure and brief description
- __Archive\\__ > contains code among other things which are not used anymore in the final version
- __Data\\__ > contains all data that is (generated and) used in our project
- __Figures\\__ > contains all figures and videos that are generated in our project
- __Include\\__ > contains various .hpp file(s) with helper functions used elsewhere
- __vtkGenerating\\__ > contains the code for generating the .vtp files for the animation as well as a CMakeLists.txt
- _centaurClassesCSVgenerating.cpp_ > generates the class data for the Centaurs in .csv files for use in vtkGenerating\\solar_system_animateVTK.cpp 
- _energyTesting.cpp_ > calculates the total Centaur energy over time and saves it as a .csv file
- _hypothesisTesting.cpp_ > calculates the inner, outer and impact Centaur asteroids over time in our simulation
for different Jupiter masses and saves them separately as csv files
- _hypothesisTesting.ipynb_ > creates visualisations for our hypothesis testing
- _planetVisualisation.ipynb_ > creates visualisations for checking accuracy of planets
- _processCentaurs.ipynb_ > clones Centaurs and creates visualisations of the cloning process
- _solar\_system.hpp_ > the c++ header file that contains the basis classes for running the simulation

### Dependencies
- C++
- g++ (or another c++ compiler)
- CMake (only for vtk-visualisation)
- Python with the packages:
    - Numpy
    - Pandas
    - Matplotlib
    - Pyorb (for kepler orbits)
- Jupiter notebook

### g++ compilation command (or use tasks.json)
> g++.exe -fdiagnostics-color=always -O3 -march=native -std=c++17 {filepath_to_compile} -o {filepath_to_save_.exe} -g

## google drive:
For large files like the made animations using ParaView, and the data used for these animations we refer to our google drive
> https://drive.google.com/drive/folders/1oEzvrlVWRl6ETvyb8Z7ytmVTRoND0eae?usp=sharing

### planet- & celestialData source
NASA Jet Propulsion Laboratory, "JPL Horizons System," California Institute of Technology, Pasadena, CA, 2023. [Online]. Available: https://ssd.jpl.nasa.gov/horizons/. [Accessed: January. 2025].
