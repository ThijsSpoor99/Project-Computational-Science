/* 
This program calculates the inner, outer and impact Centaur asteroids over time in our simulation
for different Jupiter masses and saves them separately as csv files.
*/

#include "solar_system.hpp"

int main()
{
    // create simulation class
    SolarSystem sim("Data\\", 24375, 0.6);
    sim.dt = 5;

    int nSteps = int(4e6);
    std::vector<int> nInner(nSteps);
    std::vector<int> nOuter(nSteps);
    std::vector<int> nImpacts(nSteps);

    // perform simulation for every Jupiter mass
    for (double factor = 0.0; factor < 2.1; factor += 0.20)
    {
        sim.resetSimulation();
        sim.celestials[3].GM *= factor;

        // set t=0
        nInner[0] = 0;
        nOuter[0] = 0;

        auto startTime = std::chrono::high_resolution_clock::now();
        for (int step = 1; step < nSteps; step++)
        {
            sim.performTimestep();
            nInner[step] = sim.nInner;
            nOuter[step] = sim.nOuter;
            nImpacts[step] = sim.nImpacts;
        }

        auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);
        std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;

        std::string innerFilepath = "Data/hypothesis/inner" + std::to_string(factor).substr(0, 4) + ".csv";
        std::string outerFilepath = "Data/hypothesis/outer" + std::to_string(factor).substr(0, 4) + ".csv";
        std::string impactFilepath = "Data/hypothesis/impact" + std::to_string(factor).substr(0, 4) + ".csv";

        save1DIntVector(nInner, innerFilepath);
        save1DIntVector(nOuter, outerFilepath);
        save1DIntVector(nImpacts, impactFilepath);
    }

    return 0;
}