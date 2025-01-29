#include "solar_system.hpp"

int main()
{
    SolarSystem sim("Data/", 24375, 0.6);

    std::array<int, 3> dtArray = {5, 100};

    for (int i = 0; i < 2; i++)
    {
        sim.resetSimulation();
        std::vector<double> centaurEnergy = {};
        centaurEnergy.push_back(sim.centaurEnergy);

        sim.dt = dtArray[i];

        auto startTime = std::chrono::high_resolution_clock::now();
        for (int t = dtArray[i]; t < 1e7; t += dtArray[i])
        {
            sim.performTimestep();
            sim.calcSystemEnergy();
            centaurEnergy.push_back(sim.centaurEnergy);
        }

        auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stopTime - startTime);

        std::cout << "Simulation time: " << duration.count() << " seconds" << std::endl;

        std::string filepath = "Data/energyTesting/centaurEnergy" + std::to_string(dtArray[i]) + "dt.csv";
        save1DVector(centaurEnergy, filepath);
    }

    return 0;
}