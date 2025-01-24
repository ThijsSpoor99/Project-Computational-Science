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

void saveObjectPositionsVTK(SolarSystem system, int t, int fname_digits) {
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

    double celestialRadii[7] = {
        696340.0,  // Sun radius in km
        6371.0,    // Earth radius in km
        3389.5,    // Mars radius in km
        69911.0,   // Jupiter radius in km
        58232.0,   // Saturn radius in km
        25362.0,   // Uranus radius in km
        24622.0    // Neptune radius in km
    };

    // Add celestial object data
    for (int i = 0; i < system.nCelestials; i++) {
        points->InsertNextPoint(system.celestials[i].position.data());
        radii->InsertNextValue(celestialRadii[i]);
        colors->InsertNextTuple3(celColors[i][0], celColors[i][1], celColors[i][2]); // Example color: cyan
        labels->InsertNextValue(system.celestials[i].name);
    }

    // Add centaur data
    for (const auto& centaur : system.centaurs) {
        points->InsertNextPoint(centaur.position.data());
        radii->InsertNextValue(1000); // Example small radius for centaurs
        colors->InsertNextTuple3(255, 0, 0); // Example color: red
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
    filename << PATH_TO_DATA + "vtpData\\" << "out_" 
            << std::setw(fname_digits) 
            << std::setfill('0') << t << ".vtp";

    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.str().c_str());
    writer->SetInputData(polyData);
    writer->Write();
}

int main()
{
    SolarSystem sim(PATH_TO_DATA);
    clearAndCreateDirectory(PATH_TO_DATA + "vtpData");
    int timeSteps = 10;
    int fname_digits = std::to_string(timeSteps).length();

    auto startTime = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < timeSteps; i++)
    {
        sim.performTimestep();
        saveObjectPositionsVTK(sim, i, fname_digits);
    }
    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);

    std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;

    return 0;
}