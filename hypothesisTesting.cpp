#include "solar_system.hpp"

int main()
{
    SolarSystem sim("Data\\", 24375, 0.8);
    sim.dt = 10;

    int nSteps = int(1e6);
    std::vector<int> nInner(nSteps);
    std::vector<int> nOuter(nSteps);

    for (double factor=0.0; factor<2.1; factor+=0.25) {
        sim.resetSimulation(); // to be sure
        sim.celestials[3].GM *= factor; 
        
        nInner[0] = 0;
        nOuter[0] = 0;
        
        auto startTime = std::chrono::high_resolution_clock::now();
        for (int t = 1; t<nSteps; t++) {
            sim.performTimestep();
            nInner[t] = sim.nInner;
            nOuter[t] = sim.nOuter;
        }

        auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);
        std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;

        std::string innerFilepath = "Data/hypothesis/inner" + std::to_string(factor).substr(0, 4) + ".csv";
        std::string outerFilepath = "Data/hypothesis/outer" + std::to_string(factor).substr(0, 4)  + ".csv";

        save1DIntVector(nInner, innerFilepath);
        save1DIntVector(nOuter, outerFilepath);
    }

    return 0;
}