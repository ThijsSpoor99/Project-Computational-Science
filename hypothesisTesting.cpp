#include "solar_system.hpp"

int main()
{
    SolarSystem sim("Data\\", 24375, 0.8);
    
    int nSteps = int(1e6);
    std::vector<int> nInner(nSteps);
    std::vector<int> nOuter(nSteps);

    sim.celestials[3].GM *= 1; 
    sim.dt = 10;

    nInner[0] = 0;
    nOuter[0] = 0;
    
    auto startTime = std::chrono::high_resolution_clock::now();
    for (int t = 1; t<nSteps; t++) {
        sim.performTimestep();
        nInner[t] = sim.nInner;
        nOuter[t] = sim.nOuter;

        if (t % int(1e4) == 0) {
            std::cout << 100 * float(t) / nSteps << std::endl;
        }
    }

    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);
    std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;

    save1DIntVector(nInner, "Data/hypothesis/innerTest.csv");
    save1DIntVector(nOuter, "Data/hypothesis/outerTest.csv");

    return 0;
}