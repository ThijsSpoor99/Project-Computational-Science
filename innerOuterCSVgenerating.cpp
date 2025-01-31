/* This file will generate the class data for the centaurs in the simulation. 
 * The centaurs are classified as either inner, outer or neither. 
 * The class data is saved in a csv file. The simulation is run for three different values of the mass ratio MGj = 0.25, 1.0, 1.75. 
 * The simulation is also run for an unbiased case with MGj = 1.0. The class data is saved in the Data folder.
 * 
 * The resulting class data will be used for the colors of celestial objects in the .vtp files.
 * These .vtp files are generated in the vtkGenerating/solar_system_animateVTK.cpp file.
 */


#include "solar_system.hpp"

int main()
{
    SolarSystem sim("Data\\", 24375, 1.0);
    sim.dt = 5;

    int nSteps = int(4e6);
    std::cout << "Beginning simulation" << std::endl;


    // -------------------------------------------------
    // Get the class data for biased simulations with MGj = 0.25, 1.0, 1.75
    // -------------------------------------------------
    for (double factor=0.25; factor<1.76; factor+=0.75) {
        sim.resetSimulation(); // to be sure
        sim.celestials[3].GM *= factor; 
        
        std::cout << "Running simulation with factor: " << factor << std::endl;
        auto startTime = std::chrono::high_resolution_clock::now();
        for (int t = 1; t<nSteps; t++) {
            if (t % 100000 == 0) {
                std::cout << "Step: " << t << std::endl;
            }
            sim.performTimestep();
        }

        auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);
        std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;

        std::vector<int> classes(sim.nCentaurs);
        for (int i=0; i<sim.nCentaurs; i++) {
            if (sim.centaurs[i].inner) {
                classes[i] = 1;
            } else if (sim.centaurs[i].outer) {
                classes[i] = 2;
            } else {
                classes[i] = 0;
            }
        }
        std::string classesFilepath = "Data/vtkCentaurClasses/classes" + std::to_string(factor).substr(0, 4) + ".csv";
        save1DIntVector(classes, classesFilepath);
    }

    // -------------------------------------------------
    // Get the class data for an unbiased simulation with MGj = 1
    // -------------------------------------------------
    sim.resetSimulation(); // to be sure
    sim.celestials[3].GM *= 1;
    sim.bias = 1.0;

    std::cout << "Running unbiassed simulation" << std::endl;
    auto startTime = std::chrono::high_resolution_clock::now();
    for (int t = 1; t<nSteps; t++) {
        if (t % 100000 == 0) {
            std::cout << "Step: " << t << std::endl;
        }
        sim.performTimestep();
    }
    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);
    std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;
    std::vector<int> classes(sim.nCentaurs);
    for (int i=0; i<sim.nCentaurs; i++) {
        if (sim.centaurs[i].inner) {
            classes[i] = 1;
        } else if (sim.centaurs[i].outer) {
            classes[i] = 2;
        } else {
            classes[i] = 0;
        }
    }
    save1DIntVector(classes, "Data/vtkCentaurClasses/classes_nobias.csv");

    return 0;
}