#include "solar_system.hpp"

int main() {
    SolarSystem sim;
    std::vector<double> centaurEnergy = {};
    sim.dt = 1;

    for (double t=1; t<1e6; t += 1) {
        sim.performTimestep();
        sim.calcSystemEnergy();
        centaurEnergy.push_back(sim.totalEnergy);
    }    

    save1DVector(centaurEnergy , "Data/totalEnergy1.csv");

    return 0;
}