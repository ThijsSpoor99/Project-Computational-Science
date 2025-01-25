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

void saveObjectPositionsVTK(SolarSystem sim) {
    // Create the VTK directory
    std::string vtkDir = PATH_TO_DATA + "vtpData";
    clearAndCreateDirectory(vtkDir);

    // Determine number of digits for file naming
    int digits = std::to_string(sim.maxTimesteps).length();
    
    for (int t = 0; t < sim.maxTimesteps; ++t) {
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
        //radii of celestial objects in km
        double celestialRadii[7] = {
            696340.0,  // Sun 
            6371.0,    // Earth
            3389.5,    // Mars
            69911.0,   // Jupiter
            58232.0,   // Saturn
            25362.0,   // Uranus
            24622.0    // Neptune
        };

        double minRadius = *std::min_element(celestialRadii, celestialRadii + 7);
        double maxRadius = *std::max_element(celestialRadii, celestialRadii + 7);

        // minmax scaling of celestial radii
        for (int i = 0; i < 7; ++i) {
            celestialRadii[i] = (celestialRadii[i] - minRadius) / (maxRadius - minRadius);
        }

        // Add celestial object data
        for (int i = 0; i < sim.nCelestials; i++) {
            points->InsertNextPoint(sim.celestials[i].positionHistory[t].data());
            radii->InsertNextValue(celestialRadii[i]);
            colors->InsertNextTuple3(celColors[i][0], celColors[i][1], celColors[i][2]);
            labels->InsertNextValue(sim.celestials[i].name);
        }

        for (const auto& centaur : sim.centaurs) {
            if (centaur.positionHistory.size() < sim.maxTimesteps) {
                std::cerr << "Centaur '" << centaur.name << "' has insufficient positionHistory size: " 
                        << centaur.positionHistory.size() << ", expected: " << sim.maxTimesteps << ", alive" << centaur.exist << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        // Add centaur data
        for (const auto& centaur : sim.centaurs) {
            points->InsertNextPoint(centaur.positionHistory[t].data());
            radii->InsertNextValue(celestialRadii[2]); //set Radii, same size as smallest planet
            if (centaur.inner) {
                colors->InsertNextTuple3(255, 0, 0); // color: red
            } else if (centaur.outer) {
                colors->InsertNextTuple3(255, 255, 0); // color: yellow
            } else {
                colors->InsertNextTuple3(0, 255, 0); // color: green
            }
            labels->InsertNextValue(centaur.name);
        }
        
        // Create polydata and attach arrays
        auto polyData = vtkSmartPointer<vtkPolyData>::New();
        polyData->SetPoints(points);
        polyData->GetPointData()->AddArray(radii);
        polyData->GetPointData()->AddArray(colors);
        polyData->GetPointData()->AddArray(labels);
        
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
    int SimulationTime = 1500;
    int timeStepsSaved = 1000;
    int nCentaurs = 66884;

    //int nCentaurs = 10000;

    SolarSystem sim(PATH_TO_DATA, timeStepsSaved, nCentaurs);
    
    auto startTime = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < SimulationTime; i++)
    {
        if (i % 1000 == 0) {
            std::cout << "Timestep: " << i << std::endl;
        }
        sim.performTimestep();
    }

    saveObjectPositionsVTK(sim);

    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);

    std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;

    return 0;
}