/* This file is used to create the .vtp files for the animations of the solar system and centaurs.
 * The .vtp files are created for the biased simulations with MGj = 0.25, 1.0, 1.75 and the unbiased simulation with MGj = 1.0.
 * The centaur classes are determined by the centaurClassesCSVgenerating.cpp file and are used for the colors of the celestial objects.
 * 
 * NOTE:
 * The generated .vtp files are stored in seperate folders in the Data folder and can be visualized using the ParaView software.
 * 
 * IMPORTANT:
 * For this file to run you need to have the vtk library installed and added to your path. This way
 * cmake can find the necessary library files and link it to the project. Also for execution of the compiled
 * code it is necessary to have the vtk dll files also added to your path.
 * Of course having cmake installed is also required.
 */

#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <algorithm>

#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkStringArray.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataWriter.h>

#include "..\Include\util.hpp"
#include "..\solar_system.hpp"

// Relative Path to Data folder
const std::string PATH_TO_DATA = "..\\..\\..\\Data\\";

// Function to save the positions of celestial objects and centaurs to .vtp files
void saveObjectPositionsVTK(SolarSystem sim, 
                            std::vector<std::vector<std::array<double, 3>>> celestialPosData, 
                            std::vector<std::vector<std::array<double, 3>>> centaurPosData,
                            std::vector<int> centaurClasses,
                            std::string dirName) {
    // Create the VTK directory
    std::string vtkDir = PATH_TO_DATA + dirName;
    clearAndCreateDirectory(vtkDir);

    // Determine number of digits for file naming
    int digits = std::to_string(centaurPosData[0].size()).length();
    
    for (size_t t = 0; t < centaurPosData[0].size(); ++t) {
        // Create VTK objects
        auto points = vtkSmartPointer<vtkPoints>::New();
        auto radii = vtkSmartPointer<vtkFloatArray>::New();
        auto colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
        auto labels = vtkSmartPointer<vtkStringArray>::New();
        
        // Initialize arrays
        radii->SetName("radius");
        colors->SetName("colors");
        colors->SetNumberOfComponents(3); // RGB
        labels->SetName("labels");
        
        int celColors[7][3] = {
            {255, 204, 51}, //Sun color
            {79,76,176}, //Earth color
            {193,68,14}, //Mars color
            {201,144,57}, //Jupiter color
            {226,191,125}, //Saturn color
            {147,205,241}, //Uranus color
            {61,94,249} //Neptune color
        };
        // Logarithmically scaled radii of celestial objects between 4 and 10, for visualization purposes
        // in ParaViewS
        double celestialRadii[7] = {
            10.0,       //Sun
            4.7110485,   //Earth
            4.0,        //Mars
            7.4100869,   //Jupiter
            7.2041347,   //Saturn
            6.2676178,   //Uranus
            6.2342535    //Neptune
        };

        // Add celestial object data
        for (int i = 0; i < celestialPosData.size(); i++) {
            points->InsertNextPoint(celestialPosData[i][t].data());
            radii->InsertNextValue(celestialRadii[i]);
            colors->InsertNextTuple3(celColors[i][0], celColors[i][1], celColors[i][2]);
        }

        // Add centaur data
        for (int i = 0; i < centaurPosData.size(); i++) {
            points->InsertNextPoint(centaurPosData[i][t].data());
            if (centaurClasses[i] == 0) {
                colors->InsertNextTuple3(0, 255, 0); // Normal Centaur, color: green
                radii->InsertNextValue(3);
            } else if (centaurClasses[i] == 1) {
                colors->InsertNextTuple3(255, 0, 0); // Inner Centaur, color: red
                radii->InsertNextValue(3);
            } else if (centaurClasses[i] == 2) {
                colors->InsertNextTuple3(255, 255, 0); // Outer Centaur, color: yellow
                radii->InsertNextValue(3);
            } else {
                
            }
        }
        
        // Create polydata and attach arrays
        auto polyData = vtkSmartPointer<vtkPolyData>::New();
        polyData->SetPoints(points);
        polyData->GetPointData()->AddArray(radii);
        polyData->GetPointData()->AddArray(colors);
        
        // Write to VTK files with leading zeros to make ParaView understand it is
        // a time series
        std::ostringstream filename;
        filename << vtkDir << "\\out_" 
                << std::setw(digits) 
                << std::setfill('0') << t << ".vtp";
        auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName(filename.str().c_str());
        writer->SetInputData(polyData);
        writer->Write();
    }
}

// Runs the simulation and saves the position of celestial objects and centaurs of every second time step
// calls saveObjectPositionsVTK to save these positions to .vtp files
void runSimulation(int nCentaurs, int lowerbound, int upperbound, std::string dirName, std::vector<int> centaurClasses,
                   double jupiterMassMultiplier, double bias = 0.6) {
    
    SolarSystem sim(PATH_TO_DATA, nCentaurs, bias);
    sim.dt = 5;
    sim.celestials[3].mass *= jupiterMassMultiplier;

    std::vector<std::vector<std::array<double, 3>>> celestialPosData(sim.nCelestials);
    std::vector<std::vector<std::array<double, 3>>> centaurPosData(sim.nCentaurs);

    for (int i = 0; i < upperbound; i++)
    {   
        if (i > lowerbound) {
            if (i % 2 == 0) { //save every second timestep to have 10 days per file
                for (int j = 0; j < sim.nCelestials; j++) {
                    celestialPosData[j].push_back(sim.celestials[j].position);
                }
                for (int j = 0; j < sim.nCentaurs; j++) {
                    centaurPosData[j].push_back(sim.centaurs[j].position);
                }
            }
        }
        sim.performTimestep();
    }
    saveObjectPositionsVTK(sim, celestialPosData, centaurPosData, centaurClasses, dirName);
}


int main()
{   
    // -------------------------------------------------
    // Create the .vtp files for the biased animations with MGj = 0.25, 1.0, 1.75
    // -------------------------------------------------
    for (double factor=0.25; factor<1.76; factor+=0.75) {
        std::string readName = "vtkCentaurClasses/classes" + std::to_string(factor).substr(0, 4);
        std::string writeName = "vtkdata" + std::to_string(factor).substr(0, 4);

        std::vector<int> classes =  read1DVector(PATH_TO_DATA + readName + ".csv");
        auto startTime = std::chrono::high_resolution_clock::now();
        runSimulation(24375, 0, 10000, writeName, classes, factor);
        auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);
        std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;
    }

    // -------------------------------------------------
    // Create the .vtp files for the unbiased animation with MGj = 1.0
    // -------------------------------------------------
    std::string readName = "vtkCentaurClasses\\classes_nobias";
    std::string writeName = "vtkdataNo_bias";

    std::vector<int> classes =  read1DVector(PATH_TO_DATA + readName + ".csv");
    auto startTime = std::chrono::high_resolution_clock::now();
    runSimulation(24375, 0, 10000, writeName, classes, 1.0, 1.0);
    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);
    std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;

    return 0;
}