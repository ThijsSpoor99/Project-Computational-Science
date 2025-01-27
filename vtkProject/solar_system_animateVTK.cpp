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
                            std::vector<std::vector<std::array<double, 3>>> centaurPosData) {
    // Create the VTK directory
    std::string vtkDir = PATH_TO_DATA + "vtpData";
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
            1.,  // Sun 
            0.47110485,    // Earth
            0.4,    // Mars
            0.74100869,   // Jupiter
            0.72041347,   // Saturn
            0.62676178,   // Uranus
            0.62342535    // Neptune
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
            if (sim.centaurs[i].inner) {
                colors->InsertNextTuple3(255, 0, 0); // color: red
                radii->InsertNextValue(celestialRadii[3]); //set Radii, same size as Jupiter
            } else if (sim.centaurs[i].outer) {
                colors->InsertNextTuple3(255, 255, 0); // color: yellow
                radii->InsertNextValue(celestialRadii[3]); //set Radii, same size as Jupiter
            } else {
                colors->InsertNextTuple3(0, 255, 0); // color: green
                radii->InsertNextValue(celestialRadii[2]); //set Radii, same size as Mars
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

int main()
{
    int SimulationTime = 1e6;
    //int nCentaurs = 66884;

    int nCentaurs = 1000;

    SolarSystem sim(PATH_TO_DATA, nCentaurs);
    
    std::vector<std::vector<std::array<double, 3>>> celestialPosData;
    std::vector<std::vector<std::array<double, 3>>> centaurPosData;

    celestialPosData.resize(sim.nCelestials);
    centaurPosData.resize(sim.nCentaurs);

    auto startTime = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < SimulationTime; i++)
    {
        if (i % 10000 == 0) {
            std::cout << "Timestep: " << i << std::endl;
        }
        if (i < 500) {
            for (int j = 0; j < sim.nCelestials; j++) {
                celestialPosData[j].push_back(sim.celestials[j].position);
            }
            for (int j = 0; j < sim.nCentaurs; j++) {
                centaurPosData[j].push_back(sim.centaurs[j].position);
            }
        }
        sim.performTimestep();
    }

    saveObjectPositionsVTK(sim, celestialPosData, centaurPosData);

    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);

    std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;

    return 0;
}