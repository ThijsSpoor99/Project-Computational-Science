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

const std::string PATH_TO_DATA = "..\\..\\..\\Data\\";

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
            {255, 204, 51}, //sun color
            {79,76,176}, //earth color
            {193,68,14}, //mars color
            {201,144,57}, //jupiter color
            {226,191,125}, //saturn color
            {147,205,241}, //uranus color
            {61,94,249} //neptune color
        };
        //Logarithmically scaled radii of celestial objects between 0.4 and 1, for visualization purposes
        double celestialRadii[7] = {
            10.0,       // Sun
            4.7110485,   // Earth
            4.0,        // Mars
            7.4100869,   // Jupiter
            7.2041347,   // Saturn
            6.2676178,   // Uranus
            6.2342535    // Neptune
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
                colors->InsertNextTuple3(0, 255, 0); // color: green
                radii->InsertNextValue(3); //set Radii, smaller as mars
            } else if (centaurClasses[i] == 1) {
                colors->InsertNextTuple3(255, 0, 0); // color: red
                radii->InsertNextValue(3); //set Radii, same size as Jupiter
            } else if (centaurClasses[i] == 2) {
                colors->InsertNextTuple3(255, 255, 0); // color: yellow
                radii->InsertNextValue(3); //set Radii, same size as Jupiter
            } else {
                
            }
        }
        
        // Create polydata and attach arrays
        auto polyData = vtkSmartPointer<vtkPolyData>::New();
        polyData->SetPoints(points);
        polyData->GetPointData()->AddArray(radii);
        polyData->GetPointData()->AddArray(colors);
        
        // Write to VTK file
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

void runSimulation(int nCentaurs, int lowerbound, int upperbound, std::string dirName, std::vector<int> centaurClasses,
                   double jupiterMassMultiplier) {
    
    SolarSystem sim(PATH_TO_DATA, nCentaurs, 0.8);
    sim.dt = 10;
    sim.celestials[3].mass *= jupiterMassMultiplier;

    std::vector<std::vector<std::array<double, 3>>> celestialPosData(sim.nCelestials);
    std::vector<std::vector<std::array<double, 3>>> centaurPosData(sim.nCentaurs);

    for (int i = 0; i < upperbound; i++)
    {   
        if (i > lowerbound) {
            for (int j = 0; j < sim.nCelestials; j++) {
                celestialPosData[j].push_back(sim.celestials[j].position);
            }
            for (int j = 0; j < sim.nCentaurs; j++) {
                centaurPosData[j].push_back(sim.centaurs[j].position);
            }
        }
        sim.performTimestep();
    }
    saveObjectPositionsVTK(sim, celestialPosData, centaurPosData, centaurClasses, dirName);
}

//24375
int main()
{   
    
    for (double factor=0.25; factor<1.76; factor+=0.75) {
        std::string readName = "vtkSim/classes" + std::to_string(factor).substr(0, 4);
        std::string writeName = "vtkdata" + std::to_string(factor).substr(0, 4);
        std::vector<int> classes =  read1DVector(PATH_TO_DATA + readName + ".csv");
        auto startTime = std::chrono::high_resolution_clock::now();
        runSimulation(24375, 0, 5000, writeName, classes, factor);
        auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);
        std::cout << "Total simulation time 1: " << duration.count() << " ms" << std::endl;
    }

    return 0;
}