#include "solar_system.hpp"

int main()
{
    SolarSystem sim;
    std::vector<double> centaurEnergy = {};
    sim.dt = 50;

    for (double t = 1; t < 1e6; t += 500)
    {
        sim.performTimestep();
        sim.calcSystemEnergy();
        centaurEnergy.push_back(sim.centaurEnergy);
    }

    save1DVector(centaurEnergy, "Data/centaurEnergy50.csv");

    return 0;
}