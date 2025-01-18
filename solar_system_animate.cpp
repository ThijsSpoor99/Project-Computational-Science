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
    std::vector<std::array<double, 3>> position;
    std::array<double, 3> velocity;
    std::array<double, 3> acceleration; // for Velocity Verlet
    double mass;
    double GM;

    Celestial(const std::string inputName,
              const std::array<double, 3> &startPos,
              const std::array<double, 3> &startVel,
              double inputMass,
              const size_t nSteps)
        : name(inputName),
          mass(inputMass),
          GM(Constants::G * inputMass),
          velocity(startVel)
    {
        position.resize(nSteps);
        position[0] = startPos;
        acceleration = {0.0, 0.0, 0.0};
        GM = GMtoAstro(GM);
    }
};

class Centaur
{
public:
    std::string name;
    std::vector<std::array<double, 3>> position;
    std::array<double, 3> velocity;
    std::array<double, 3> acceleration; // for Velocity Verlet

    Centaur(const std::string name,
            const std::array<double, 3> &startPos,
            const std::array<double, 3> &startVel,
            const size_t nSteps)
        : name(name), velocity(startVel)
    {
        position.resize(nSteps);
        position[0] = startPos;
        acceleration = {0.0, 0.0, 0.0};
    }
};

class SolarSystem
{
public:
    const int nPlanets = 8;
    const int nCentaurs = 10;
    const int dt = 1;        // 1 day
    const long nSteps = 1e7; // 1e6

    Celestial Sun;
    std::vector<Celestial> planets;
    std::vector<Centaur> centaurs;

    std::vector<std::vector<std::string>> planetData;
    std::vector<std::vector<std::string>> centaurData;

    // temporary vectors used in computations
    std::array<double, 3> rVec;
    std::array<double, 3> aVec;

    double vSquared = 0.0;
    double rNorm = 0.0;
    double rCubed = 0.0;
    double GMrCubed = 0.0;
    double tempEnergy = 0.0;
    std::vector<double> systemEnergy;

    SolarSystem()
        : Sun(Celestial("Sun",
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        1.989e30,
                        nSteps))
          // read data
          ,
          planetData(readCSV("Data\\planetData.csv")), centaurData(readCSV("Data\\CentaursCartesian.csv"))
    {
        // create planets
        for (int i = 0; i < nPlanets; i++)
        {
            planets.push_back(Celestial(
                planetData[i][0],
                {std::stod(planetData[i][1]),
                 std::stod(planetData[i][2]),
                 std::stod(planetData[i][3])},
                {std::stod(planetData[i][4]),
                 std::stod(planetData[i][5]),
                 std::stod(planetData[i][6])},
                std::stod(planetData[i][7]),
                nSteps));
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
                 std::stod(centaurData[i][6])},
                nSteps));
        }
    }

    // -------------------------------------------------
    // Compute the net acceleration for each planet at index t
    // -------------------------------------------------
    void computePlanetsAcceleration(int t)
    {
        for (int i = 0; i < nPlanets; i++)
        {
            std::array<double, 3> aLocal = {0.0, 0.0, 0.0};

            // planet-planet gravitational pull
            for (int j = 0; j < nPlanets; j++)
            {
                if (j == i)
                    continue;
                rVec = subtract3D(planets[i].position[t], planets[j].position[t]);
                rCubed = calcCube(rVec);
                GMrCubed = planets[j].GM / rCubed;
                aLocal[0] -= GMrCubed * rVec[0];
                aLocal[1] -= GMrCubed * rVec[1];
                aLocal[2] -= GMrCubed * rVec[2];
            }

            // planet-Sun gravitational pull (Sun is at origin)
            rVec = planets[i].position[t];
            rCubed = calcCube(rVec);
            GMrCubed = Sun.GM / rCubed;
            aLocal[0] -= GMrCubed * rVec[0];
            aLocal[1] -= GMrCubed * rVec[1];
            aLocal[2] -= GMrCubed * rVec[2];

            planets[i].acceleration = aLocal;
        }
    }

    // -------------------------------------------------
    // Compute the net acceleration for each centaur at index t
    // -------------------------------------------------
    void computeCentaursAcceleration(int t)
    {
        for (int i = 0; i < nCentaurs; i++)
        {
            std::array<double, 3> aLocal = {0.0, 0.0, 0.0};

            // from planets
            for (int j = 0; j < nPlanets; j++)
            {
                rVec = subtract3D(centaurs[i].position[t], planets[j].position[t]);
                rCubed = calcCube(rVec);
                GMrCubed = planets[j].GM / rCubed;
                aLocal[0] -= GMrCubed * rVec[0];
                aLocal[1] -= GMrCubed * rVec[1];
                aLocal[2] -= GMrCubed * rVec[2];
            }
            // from Sun
            rVec = centaurs[i].position[t];
            rCubed = calcCube(rVec);
            GMrCubed = Sun.GM / rCubed;
            aLocal[0] -= GMrCubed * rVec[0];
            aLocal[1] -= GMrCubed * rVec[1];
            aLocal[2] -= GMrCubed * rVec[2];

            centaurs[i].acceleration = aLocal;
        }
    }

    // -------------------------------------------------
    // Velocity Verlet update for planets
    // -------------------------------------------------
    void updatePlanets(int t)
    {
        // half-step velocity + position
        for (int i = 0; i < nPlanets; i++)
        {
            // half-step velocity
            planets[i].velocity[0] += 0.5 * dt * planets[i].acceleration[0];
            planets[i].velocity[1] += 0.5 * dt * planets[i].acceleration[1];
            planets[i].velocity[2] += 0.5 * dt * planets[i].acceleration[2];

            // position update
            planets[i].position[t][0] = planets[i].position[t - 1][0] + dt * planets[i].velocity[0];
            planets[i].position[t][1] = planets[i].position[t - 1][1] + dt * planets[i].velocity[1];
            planets[i].position[t][2] = planets[i].position[t - 1][2] + dt * planets[i].velocity[2];
        }

        // recompute acceleration at new positions
        std::vector<std::array<double, 3>> newA(nPlanets, {0.0, 0.0, 0.0});
        for (int i = 0; i < nPlanets; i++)
        {
            std::array<double, 3> aLocal = {0.0, 0.0, 0.0};

            // planet-planet
            for (int j = 0; j < nPlanets; j++)
            {
                if (j == i)
                    continue;
                rVec = subtract3D(planets[i].position[t], planets[j].position[t]);
                rCubed = calcCube(rVec);
                GMrCubed = planets[j].GM / rCubed;
                aLocal[0] -= GMrCubed * rVec[0];
                aLocal[1] -= GMrCubed * rVec[1];
                aLocal[2] -= GMrCubed * rVec[2];
            }
            // planet-Sun
            rVec = planets[i].position[t];
            rCubed = calcCube(rVec);
            GMrCubed = Sun.GM / rCubed;
            aLocal[0] -= GMrCubed * rVec[0];
            aLocal[1] -= GMrCubed * rVec[1];
            aLocal[2] -= GMrCubed * rVec[2];

            newA[i] = aLocal;
        }

        // second half-step velocity update
        for (int i = 0; i < nPlanets; i++)
        {
            planets[i].velocity[0] += 0.5 * dt * newA[i][0];
            planets[i].velocity[1] += 0.5 * dt * newA[i][1];
            planets[i].velocity[2] += 0.5 * dt * newA[i][2];

            planets[i].acceleration = newA[i]; // store for next iteration
        }
    }

    // -------------------------------------------------
    // Velocity Verlet update for centaurs
    // -------------------------------------------------
    void updateCentaurs(int t)
    {
        // 1) Half-step velocity + position
        for (int i = 0; i < nCentaurs; i++)
        {
            centaurs[i].velocity[0] += 0.5 * dt * centaurs[i].acceleration[0];
            centaurs[i].velocity[1] += 0.5 * dt * centaurs[i].acceleration[1];
            centaurs[i].velocity[2] += 0.5 * dt * centaurs[i].acceleration[2];

            centaurs[i].position[t][0] = centaurs[i].position[t - 1][0] + dt * centaurs[i].velocity[0];
            centaurs[i].position[t][1] = centaurs[i].position[t - 1][1] + dt * centaurs[i].velocity[1];
            centaurs[i].position[t][2] = centaurs[i].position[t - 1][2] + dt * centaurs[i].velocity[2];
        }

        // compute new acceleration
        std::vector<std::array<double, 3>> newA(nCentaurs, {0.0, 0.0, 0.0});
        for (int i = 0; i < nCentaurs; i++)
        {
            std::array<double, 3> aLocal = {0.0, 0.0, 0.0};

            // from planets
            for (int j = 0; j < nPlanets; j++)
            {
                rVec = subtract3D(centaurs[i].position[t], planets[j].position[t]);
                rCubed = calcCube(rVec);
                GMrCubed = planets[j].GM / rCubed;
                aLocal[0] -= GMrCubed * rVec[0];
                aLocal[1] -= GMrCubed * rVec[1];
                aLocal[2] -= GMrCubed * rVec[2];
            }
            // from Sun
            rVec = centaurs[i].position[t];
            rCubed = calcCube(rVec);
            GMrCubed = Sun.GM / rCubed;
            aLocal[0] -= GMrCubed * rVec[0];
            aLocal[1] -= GMrCubed * rVec[1];
            aLocal[2] -= GMrCubed * rVec[2];

            newA[i] = aLocal;
        }

        // second half-step velocity update
        for (int i = 0; i < nCentaurs; i++)
        {
            centaurs[i].velocity[0] += 0.5 * dt * newA[i][0];
            centaurs[i].velocity[1] += 0.5 * dt * newA[i][1];
            centaurs[i].velocity[2] += 0.5 * dt * newA[i][2];

            centaurs[i].acceleration = newA[i]; // store for next iteration
        }
    }

    // -------------------------------------------------
    // Compute the total (kinetic + potential) energy in the system
    // -------------------------------------------------
    void calcSystemEnergy(int t)
    {
        tempEnergy = 0.0;

        // planets
        for (int i = 0; i < nPlanets; i++)
        {
            // kinetic: 1/2 m v^2
            vSquared = calcSquare(planets[i].velocity);
            tempEnergy += 0.5 * planets[i].mass * vSquared;

            // potential from planet-planet
            for (int j = i + 1; j < nPlanets; j++)
            {
                rVec = subtract3D(planets[i].position[t], planets[j].position[t]);
                rNorm = calcNorm(rVec);
                tempEnergy -= planets[i].GM * planets[j].mass / rNorm;
            }
            // potential from planet-sun
            rNorm = calcNorm(planets[i].position[t]);
            tempEnergy -= Sun.GM * planets[i].mass / rNorm;
        }

        // centaurs: assume negligible mass => treat as m=1 for KE
        for (int i = 0; i < nCentaurs; i++)
        {
            // kinetic
            vSquared = calcSquare(centaurs[i].velocity);
            tempEnergy += 0.5 * vSquared;

            // potential from centaur-planet
            for (int j = 0; j < nPlanets; j++)
            {
                rVec = subtract3D(centaurs[i].position[t], planets[j].position[t]);
                rNorm = calcNorm(rVec);
                tempEnergy -= planets[j].GM / rNorm;
            }
            // potential from centaur-sun
            rNorm = calcNorm(centaurs[i].position[t]);
            tempEnergy -= Sun.GM / rNorm;
        }

        std::cout << "Energy at step " << t << ": " << tempEnergy << std::endl;
        systemEnergy.push_back(tempEnergy);
    }

    // -------------------------------------------------
    // Run the sim
    // -------------------------------------------------
    void simulate()
    {
        // compute initial accelerations at t=0
        computePlanetsAcceleration(0);
        computeCentaursAcceleration(0);

        // check energy at t=0
        calcSystemEnergy(0);

        std::cout << "Simulating..." << std::endl;
        for (int t = 1; t < nSteps; t++)
        {
            updatePlanets(t);
            updateCentaurs(t);

            // check energy every 100000 steps
            if (t % 100000 == 0)
            {
                calcSystemEnergy(t);
            }
        }
    }

    void savePlanetPositions() {
        for (int i=0; i<nPlanets; i++) {
            saveToCSV(planets[i].position, "Data\\planetPositions\\" + planets[i].name + ".csv");
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

    sim.savePlanetPositions();

    std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;

    return 0;
}