# Influence of Jupiter-like Planets on Centaur Transits to the Inner Solar System
This repository is for the course "Project Computational Science" at the University of Amsterdam.
For a quick overview of the project we refer to (poster pdf location). For a more detailled description of
the project, we refer to (report-file location).

#### Structure and brief description
__Archive\\__ > contains code among other things which are not used anymore in the final version <br>
__Data\\__ > contains all data that is (generated and) used in our project <br>
__Figures\\__ > contains all figures and videos that are generated in our project <br>
__Include\\__ > contains various .hpp file(s) with helper functions used elsewhere <br>

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

### planet- & celestialData source
NASA Jet Propulsion Laboratory, "JPL Horizons System," California Institute of Technology, Pasadena, CA, 2023. [Online]. Available: https://ssd.jpl.nasa.gov/horizons/. [Accessed: January. 2025].
