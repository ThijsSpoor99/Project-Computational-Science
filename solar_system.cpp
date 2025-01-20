#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cmath>
#include <chrono>

#include "Include\util.hpp"

class Celestial
{
public:
    std::string name;
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    std::array<double, 3> acceleration; // for Velocity Verlet
    double mass;
    double GM;
    double radiusCubed;

    Celestial() = default;
    Celestial(const std::string &inputName,
              const std::array<double, 3> &startPos,
              const std::array<double, 3> &startVel,
              double inputMass,
              double inputRadius)
        : name(inputName),
          mass(inputMass),
          GM(Constants::G * inputMass),
          velocity(startVel),
          radiusCubed(10)     
    {
        position = startPos;
        acceleration = {0.0, 0.0, 0.0};
        GM = GMtoAstro(GM);
    }
};

class Centaur
{
public:
    std::string name;
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    std::array<double, 3> acceleration; // for Velocity Verlet
    bool exist = true;

    Centaur(const std::string name,
            const std::array<double, 3> &startPos,
            const std::array<double, 3> &startVel)
        : name(name), velocity(startVel)
    {
        position = startPos;
        acceleration = {0.0, 0.0, 0.0};
    }
};

class SolarSystem
{
public:
    int nCelestials;
    const int nCentaurs = 1000;
    const int dt = 1;        // 1 day
    const long nSteps = 1e6; // 1e6

    Celestial Sun;
    std::vector<Celestial> celestials;
    std::vector<Centaur> centaurs;

    std::vector<std::vector<std::string>> celestialData;
    std::vector<std::vector<std::string>> centaurData;

    // temporary vectors used in computations
    std::array<double, 3> rVec;
    std::array<double, 3> aVec;

    double vSquared = 0.0;
    double rNorm = 0.0;
    double rCubed = 0.0;
    double GMrCubed = 0.0;
    double tempEnergy = 0.0;
    std::vector<double> celestialEnergy;
    std::vector<double> centaurEnergy;

    SolarSystem()
        :celestialData(readCSV("Data\\celestialDataReduced.csv")), centaurData(readCSV("Data\\CentaursCartesian.csv"))
    {
        // create celestials
        nCelestials = celestialData.size();
        for (int i = 0; i < nCelestials; i++)
        {
            celestials.push_back(Celestial(
                celestialData[i][0],
                {std::stod(celestialData[i][1]),
                 std::stod(celestialData[i][2]),
                 std::stod(celestialData[i][3])},
                {std::stod(celestialData[i][4]),
                 std::stod(celestialData[i][5]),
                 std::stod(celestialData[i][6])},
                std::stod(celestialData[i][7]),
                std::stod(celestialData[i][8])));
        }

        // create centaurs
        for (int i = 0; i < nCentaurs; i++)
        {
            centaurs.push_back(Centaur(
                centaurData[i][0],
                {std::stod(centaurData[i][1]),
                 std::stod(centaurData[i][2]),
                 std::stod(centaurData[i][3])},
                {std::stod(centaurData[i][4]),
                 std::stod(centaurData[i][5]),
                 std::stod(centaurData[i][6])}));
        }
    }

    // -------------------------------------------------
    // Compute the net acceleration for each planet at index t
    // -------------------------------------------------
    void computeCelestialsAcceleration(int &t)
    {
        // recalculate acceleration for each celestial
        for (int i = 0; i < nCelestials; i++)
        {
            celestials[i].acceleration = {0.0, 0.0, 0.0};
            // gravitational pull from other celestials
            for (int j = 0; j < nCelestials; j++)
            {
                if (j == i) {
                    continue;
                } else {
                    rVec = subtract3D(celestials[i].position, celestials[j].position);
                    rCubed = calcCube(rVec);
                    GMrCubed = celestials[j].GM / rCubed;
                    celestials[i].acceleration[0] -= GMrCubed * rVec[0];
                    celestials[i].acceleration[1] -= GMrCubed * rVec[1];
                    celestials[i].acceleration[2] -= GMrCubed * rVec[2];
                }
            }
        }
    }

    // -------------------------------------------------
    // Compute the net acceleration for each centaur at index t
    // -------------------------------------------------
    void computeCentaursAcceleration(int &t)
    {
        // recalculate acceleration for each celestial
        for (int i = 0; i < nCentaurs; i++)
        {
            centaurs[i].acceleration = {0.0, 0.0, 0.0};
            // acceleration from celestials
            for (int j = 0; j < nCelestials; j++)
            {
                rVec = subtract3D(centaurs[i].position, celestials[j].position);
                rCubed = calcCube(rVec);
                GMrCubed = celestials[j].GM / rCubed;
                centaurs[i].acceleration[0] -= GMrCubed * rVec[0];
                centaurs[i].acceleration[1] -= GMrCubed * rVec[1];
                centaurs[i].acceleration[2] -= GMrCubed * rVec[2];
            }
        }
    }

    // -------------------------------------------------
    // Velocity Verlet update for celestials
    // -------------------------------------------------
    void updateCelestials(int &t)
    {
        // half-step velocity + position
        for (int i = 0; i < nCelestials; i++)
        {
            // half-step velocity
            celestials[i].velocity[0] += 0.5 * dt * celestials[i].acceleration[0];
            celestials[i].velocity[1] += 0.5 * dt * celestials[i].acceleration[1];
            celestials[i].velocity[2] += 0.5 * dt * celestials[i].acceleration[2];

            // position update
            celestials[i].position[0] = celestials[i].position[0] + dt * celestials[i].velocity[0];
            celestials[i].position[1] = celestials[i].position[1] + dt * celestials[i].velocity[1];
            celestials[i].position[2] = celestials[i].position[2] + dt * celestials[i].velocity[2];
        }

        // recompute acceleration at new positions
        for (int i = 0; i < nCelestials; i++)
        {
            celestials[i].acceleration = {0.0, 0.0, 0.0};
            // gravitational pull from other celestials
            for (int j = 0; j < nCelestials; j++)
            {
                if (j == i) {
                    continue;
                } else{
                    rVec = subtract3D(celestials[i].position, celestials[j].position);
                    rCubed = calcCube(rVec);
                    GMrCubed = celestials[j].GM / rCubed;
                    celestials[i].acceleration[0] -= GMrCubed * rVec[0];
                    celestials[i].acceleration[1] -= GMrCubed * rVec[1];
                    celestials[i].acceleration[2] -= GMrCubed * rVec[2];
                }
            }

            celestials[i].velocity[0] += 0.5 * dt * celestials[i].acceleration[0];
            celestials[i].velocity[1] += 0.5 * dt * celestials[i].acceleration[1];
            celestials[i].velocity[2] += 0.5 * dt * celestials[i].acceleration[2];
        }
    }

    // -------------------------------------------------
    // Velocity Verlet update for centaurs
    // -------------------------------------------------
    void updateCentaurs(int &t)
    {
        // unlike the celestial class, the acceleration of a centaur
        // does not depend on the position of the other centaurs
        // thus the calculation can be done in one loop

        // 1) Half-step velocity + position
        for (int i = 0; i < nCentaurs; i++)
        {   
            if (centaurs[i].exist) {

                centaurs[i].velocity[0] += 0.5 * dt * centaurs[i].acceleration[0];
                centaurs[i].velocity[1] += 0.5 * dt * centaurs[i].acceleration[1];
                centaurs[i].velocity[2] += 0.5 * dt * centaurs[i].acceleration[2];

                centaurs[i].position[0] = centaurs[i].position[0] + dt * centaurs[i].velocity[0];
                centaurs[i].position[1] = centaurs[i].position[1] + dt * centaurs[i].velocity[1];
                centaurs[i].position[2] = centaurs[i].position[2] + dt * centaurs[i].velocity[2];

                // acceleration from celestials
                centaurs[i].acceleration = {0.0, 0.0, 0.0};
                for (int j = 0; j < nCelestials; j++)
                {
                    rVec = subtract3D(centaurs[i].position, celestials[j].position);
                    rCubed = calcCube(rVec);

                    if (rCubed < celestials[j].radiusCubed) {
                        centaurs[i].exist = false;
                        std::cout << "impact" << std::endl;
                        break;
                    } else {
                        GMrCubed = celestials[j].GM / rCubed;
                        centaurs[i].acceleration[0] -= GMrCubed * rVec[0];
                        centaurs[i].acceleration[1] -= GMrCubed * rVec[1];
                        centaurs[i].acceleration[2] -= GMrCubed * rVec[2];
                    }
                }

                centaurs[i].velocity[0] += 0.5 * dt * centaurs[i].acceleration[0];
                centaurs[i].velocity[1] += 0.5 * dt * centaurs[i].acceleration[1];
                centaurs[i].velocity[2] += 0.5 * dt * centaurs[i].acceleration[2];

            } else {
                std::cout << "bad" << std::endl;
            }
        }
    }

    // -------------------------------------------------
    // Compute the total (kinetic + potential) energy in the system
    // -------------------------------------------------
    void calcSystemEnergy(int &t)
    {
        // celestials
        tempEnergy = 0.0;
        for (int i = 0; i < nCelestials; i++)
        {
            // kinetic: 1/2 m v^2
            vSquared = calcSquare(celestials[i].velocity);
            tempEnergy += 0.5 * celestials[i].mass * vSquared;

            // potential from planet-planet
            for (int j = i + 1; j < nCelestials; j++)
            {
                rVec = subtract3D(celestials[i].position, celestials[j].position);
                rNorm = calcNorm(rVec);
                tempEnergy -= celestials[i].GM * celestials[j].mass / rNorm;
            }
        }
        celestialEnergy.push_back(tempEnergy);
        std::cout << "Celestial energy at step " << t << ": " << tempEnergy << std::endl;

        // centaurs: assume negligible mass => treat as m=1 for KE
        tempEnergy = 0.0;
        for (int i = 0; i < nCentaurs; i++)
        {
            // kinetic
            vSquared = calcSquare(centaurs[i].velocity);
            tempEnergy += 0.5 * vSquared;

            // potential from centaur-planet
            for (int j = 0; j < nCelestials; j++)
            {
                rVec = subtract3D(centaurs[i].position, celestials[j].position);
                rNorm = calcNorm(rVec);
                tempEnergy -= celestials[j].GM / rNorm;
            }
        }
        centaurEnergy.push_back(tempEnergy);
        std::cout << "Centaur energy at step " << t << ": " << tempEnergy << std::endl;
    }

    // -------------------------------------------------
    // Run the sim
    // -------------------------------------------------
    void simulate()
    {
        // compute initial accelerations at t=0
        int t = 0;
        computeCelestialsAcceleration(t);
        computeCentaursAcceleration(t);

        // check energy at t=0
        calcSystemEnergy(t);

        std::cout << "Simulating..." << std::endl;
        for (int t = 1; t < nSteps; t++)
        {
            updateCelestials(t);
            updateCentaurs(t);

            // check energy every 100000 steps
            if (t % 100000 == 0)
            {
                calcSystemEnergy(t);
            }
        }
    }

};

int main()
{
    SolarSystem sim;

    // measure runtime
    auto startTime = std::chrono::high_resolution_clock::now();
    sim.simulate();
    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);

    std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;

    return 0;
}