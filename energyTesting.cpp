#include "solar_system.hpp"

int main()
{
    SolarSystem sim("Data/", 24375, 0.8);
    sim.dt = 10;
    double tempError;

    sim.resetSimulation();
    std::vector<double> centaurEnergyError = {};

    std::vector<double> centaurEnergy1 = {};
    std::vector<double> centaurEnergy2 = {};

    std::vector<double> centaurKineticEnergy(sim.nCentaurs);
    std::vector<double> centaurPotentialEnergy(sim.nCentaurs);
    for (int i; i<sim.nCentaurs; i++) {
            sim.calcCentaurEnergy(i);
            centaurKineticEnergy[i] = sim.centaurs[i].kineticEnergy;
            centaurPotentialEnergy[i] = sim.centaurs[i].potentialEnergy;
    }
    
    for (int t = 1; t < 1e6; t += 10)
    {   
        sim.performTimestep();
        
        tempError = 0.0;
        for (int i; i<sim.nCentaurs; i++) {
            sim.calcCentaurEnergy(i);
            tempError += centaurKineticEnergy[i] - sim.centaurs[i].kineticEnergy + centaurPotentialEnergy[i] - sim.centaurs[i].potentialEnergy;
            centaurKineticEnergy[i] = sim.centaurs[i].kineticEnergy;
            centaurPotentialEnergy[i] = sim.centaurs[i].potentialEnergy;
        }

        centaurEnergy1.push_back(sim.centaurs[22662].totalEnergy);
        centaurEnergy2.push_back(sim.centaurs[10201].totalEnergy);

        centaurEnergyError.push_back(tempError);

    }

    save1DVector(centaurEnergyError, "Data/centaurEnergyError10dt.csv");
    save1DVector(centaurEnergy1, "Data/centaurEnergy1.csv");
    save1DVector(centaurEnergy2, "Data/centaurEnergy2.csv");


    return 0;
}