#include "solar_system.hpp"

int main()
{
    SolarSystem sim("Data\\", 24375, 1.0);
    sim.dt = 5;

    int nSteps = int(1e5);
    std::cout << "Beginning simulation" << std::endl;

    for (double factor=1.0; factor<1.1; factor+=1.5) {
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
        save1DIntVector(classes, "Data/vtkCentaurClasses/classes_nobias.csv");
    }

    return 0;
}