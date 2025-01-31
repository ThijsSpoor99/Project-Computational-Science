/*
This program calculates the total Centaur energy over time and saves it as a csv file.
*/

#include "solar_system.hpp"

int main()
{
    SolarSystem sim("Data/", 24375, 0.6);

    // dt's to calc energy for
    std::array<int, 3> dtArray = {5};

    for (size_t i = 0; i < dtArray.size(); i++)
    {
        sim.resetSimulation();
        std::vector<double> centaurEnergy = {};
        centaurEnergy.push_back(sim.centaurEnergy);

        sim.dt = dtArray[i];

        auto startTime = std::chrono::high_resolution_clock::now();
        for (size_t t = dtArray[i]; t < size_t(2e7); t += dtArray[i])
        {
            sim.performTimestep();
            sim.calcSystemEnergy();
            centaurEnergy.push_back(sim.centaurEnergy);
        }

        auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stopTime - startTime);

        std::cout << "Simulation time: " << duration.count() << " seconds" << std::endl;

        std::string filepath = "Data/energyTesting/centaurEnergy_" + std::to_string(dtArray[i]) + "dt.csv";
        save1DVector(centaurEnergy, filepath);
    }

    return 0;
}