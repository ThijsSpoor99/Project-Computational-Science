#include "solar_system.hpp"

int main()
{
    SolarSystem sim("Data\\", 24375, 0.8);
    sim.dt = 10;

    int nSteps = int(1e6);
    std::cout << "Beginning simulation" << std::endl;

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
            if (sim.centaurs[i].exist) {
                classes[i] = 0;
            } else if (sim.centaurs[i].inner) {
                classes[i] = 1;
            } else {
                classes[i] = 2;
            }
        }
        std::string classesFilepath = "Data/vtkSim/classes" + std::to_string(factor).substr(0, 4) + ".csv";

        save1DIntVector(classes, classesFilepath);
    }

    return 0;
}