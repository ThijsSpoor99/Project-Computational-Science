// more types
#include <array>
#include <vector>
#include <string>

// more streams
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

// math stuff
#include <cmath>
#include <chrono>

#include "Include/util.hpp"

class Celestial {
public:

    // create class variables
    std::string name;
    std::vector<std::array<double, 3>> position;
    std::array<double, 3> velocity;
    double mass;
    double GM;

    // constructor
    Celestial() = default;
    Celestial(const std::string inputName,
              const std::array<double, 3> &startPos,
              const std::array<double, 3> &startVel,  
              double inputMass,
              const size_t nSteps
              ) {
        name = inputName;
        position.resize(nSteps);
        position[0] = distToAU3D(startPos);
        velocity = velToAstro3D(startVel);
        mass = inputMass;
        GM = GMtoAstro(Constants::G * mass);
    }
};

class SolarSystem {
public:

    // set simulation variables
    const int nPlanets = 8;
    int dt = 1;
    int nSteps = 365*200;

    // create objects
    Celestial Sun;
    std::array<Celestial, 8> planets;

    // input data
    std::vector<std::vector<std::string>> planetData;

    // calculation variables
    std::array<double, 3> rVec = {0.0, 0.0, 0.0}; 
    std::array<double, 3> aVec = {0.0, 0.0, 0.0};
    double GMrCubed = 0.0;
    
    // constructor
    SolarSystem() {
        Sun = Celestial("Sun", {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.989e30, nSteps);

        planetData = readCSV("Data//planetData.csv");

        for (int i = 0; i < nPlanets; i++) {
            planets[i] = 
            Celestial(
                planetData[i][0],
                {std::stod(planetData[i][1]), std::stod(planetData[i][2]), std::stod(planetData[i][3])},
                {std::stod(planetData[i][4]), std::stod(planetData[i][5]), std::stod(planetData[i][6])},
                std::stod(planetData[i][7]),
                nSteps
            );
        }
    }

    void updatePlanets(int &t) {
        // update planet i
        for (int i=0; i<nPlanets; i++){
            aVec = {0.0, 0.0, 0.0};

            // calc gravitational pull from planet j
            for (int j=0; j<nPlanets; j++) {
                if (i != j) {
                    rVec = subtract3D(planets[i].position[t-1], planets[j].position[t-1]);
                    GMrCubed = planets[j].GM / calcCube(rVec);

                    aVec[0] -= GMrCubed * rVec[0];
                    aVec[1] -= GMrCubed * rVec[1];
                    aVec[2] -= GMrCubed * rVec[2];
                }
            }

            // calc gravitational pull from SUN
            // Stationary sun, thus pos=rVec
            rVec = planets[i].position[t-1];
            GMrCubed = Sun.GM / calcCube(rVec);
            aVec[0] -= GMrCubed * rVec[0];
            aVec[1] -= GMrCubed * rVec[1];
            aVec[2] -= GMrCubed * rVec[2];

            // update velocity of planet i
            planets[i].velocity[0] = planets[i].velocity[0] + dt * aVec[0];
            planets[i].velocity[1] = planets[i].velocity[1] + dt * aVec[1];
            planets[i].velocity[2] = planets[i].velocity[2] + dt * aVec[2];

            // update position of planet i
            planets[i].position[t][0] = planets[i].position[t-1][0] + dt * planets[i].velocity[0];
            planets[i].position[t][1] = planets[i].position[t-1][1] + dt * planets[i].velocity[1];
            planets[i].position[t][2] = planets[i].position[t-1][2] + dt * planets[i].velocity[2];
        }
    }

    void simulate() {
        std::cout << "Simulating..." << std::endl;
        for (int t=1; t < nSteps; t++) {
            updatePlanets(t);
        }
        std::cout << "Simulation complete" << std::endl;

    }
};

int main() {
    SolarSystem sim;
    sim.simulate();

    for (int i=0; i<sim.nPlanets; i++) {
        saveToCSV(sim.planets[i].position, "Data\\planetPositions\\" + sim.planets[i].name + ".csv");
    }
    return 0;
}