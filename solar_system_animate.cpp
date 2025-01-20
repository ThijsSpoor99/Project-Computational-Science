#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cmath>
#include <chrono>
#include <iomanip>

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

    Celestial() = default;
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
    int nCelestials;
    const int nCentaurs = 30;
    const int dt = 1;        // 1 day
    const long nSteps = 365e3; // 1e6

    Celestial Sun;
    std::vector<Celestial> planets;
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
    std::vector<double> systemEnergy;

    SolarSystem()
        : celestialData(readCSV("Data\\celestialDataReduced.csv")), centaurData(readCSV("Data\\CentaursCartesian.csv"))
    {
        // create planets
        nCelestials = celestialData.size();
        for (int i = 0; i < nCelestials; i++)
        {
            planets.push_back(Celestial(
                celestialData[i][0],
                {std::stod(celestialData[i][1]),
                 std::stod(celestialData[i][2]),
                 std::stod(celestialData[i][3])},
                {std::stod(celestialData[i][4]),
                 std::stod(celestialData[i][5]),
                 std::stod(celestialData[i][6])},
                std::stod(celestialData[i][7]),
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
        for (int i = 0; i < nCelestials; i++)
        {
            std::array<double, 3> aLocal = {0.0, 0.0, 0.0};

            // planet-planet gravitational pull
            for (int j = 0; j < nCelestials; j++)
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
            for (int j = 0; j < nCelestials; j++)
            {
                rVec = subtract3D(centaurs[i].position[t], planets[j].position[t]);
                rCubed = calcCube(rVec);
                GMrCubed = planets[j].GM / rCubed;
                aLocal[0] -= GMrCubed * rVec[0];
                aLocal[1] -= GMrCubed * rVec[1];
                aLocal[2] -= GMrCubed * rVec[2];
            }

            centaurs[i].acceleration = aLocal;
        }
    }

    // -------------------------------------------------
    // Velocity Verlet update for planets
    // -------------------------------------------------
    void updatePlanets(int t)
    {
        // half-step velocity + position
        for (int i = 0; i < nCelestials; i++)
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
        std::vector<std::array<double, 3>> newA(nCelestials, {0.0, 0.0, 0.0});
        for (int i = 0; i < nCelestials; i++)
        {
            std::array<double, 3> aLocal = {0.0, 0.0, 0.0};

            // planet-planet
            for (int j = 0; j < nCelestials; j++)
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

            newA[i] = aLocal;
        }

        // second half-step velocity update
        for (int i = 0; i < nCelestials; i++)
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
            for (int j = 0; j < nCelestials; j++)
            {
                rVec = subtract3D(centaurs[i].position[t], planets[j].position[t]);
                rCubed = calcCube(rVec);
                GMrCubed = planets[j].GM / rCubed;
                aLocal[0] -= GMrCubed * rVec[0];
                aLocal[1] -= GMrCubed * rVec[1];
                aLocal[2] -= GMrCubed * rVec[2];
            }

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
        for (int i = 0; i < nCelestials; i++)
        {
            // kinetic: 1/2 m v^2
            vSquared = calcSquare(planets[i].velocity);
            tempEnergy += 0.5 * planets[i].mass * vSquared;

            // potential from planet-planet
            for (int j = i + 1; j < nCelestials; j++)
            {
                rVec = subtract3D(planets[i].position[t], planets[j].position[t]);
                rNorm = calcNorm(rVec);
                tempEnergy -= planets[i].GM * planets[j].mass / rNorm;
            }
        }

        // centaurs: assume negligible mass => treat as m=1 for KE
        for (int i = 0; i < nCentaurs; i++)
        {
            // kinetic
            vSquared = calcSquare(centaurs[i].velocity);
            tempEnergy += 0.5 * vSquared;

            // potential from centaur-planet
            for (int j = 0; j < nCelestials; j++)
            {
                rVec = subtract3D(centaurs[i].position[t], planets[j].position[t]);
                rNorm = calcNorm(rVec);
                tempEnergy -= planets[j].GM / rNorm;
            }
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
        for (int i=0; i<nCelestials; i++) {
            saveToCSV(planets[i].position, "Data\\planetPositions\\" + planets[i].name + ".csv");
        }
    }

    void saveObjectPositionsVTK(int nSteps) {
        //determine number of digits in nSteps for file name padding
        int originalSteps = nSteps;
        int digits = 0; while (nSteps != 0) { nSteps /= 10; digits++; }
        clearAndCreateDirectory("Data\\vtkData");
        
        for (int i = 0; i < originalSteps; ++i) {
            std::vector<std::array<double, 3>> stepVec;
            
            //add xyz of planets at timestep i to stepVec
            for (int ip = 0; ip < nCelestials; ++ip) {
                stepVec.push_back(planets[ip].position[i]);
            }

            //add xyz of asteroids at timestep i to stepVec
            for (int ic = 0; ic < nCentaurs; ++ic) {
                stepVec.push_back(centaurs[ic].position[i]);
            }

            std::ostringstream filename;
            filename << "Data\\vtkData\\out_"
                    << std::setfill('0')  // fill with zeros
                    << std::setw(digits)  // width based on time steps
                    << i
                    << ".vtk";

            // Save the data to a VTK file for the current timestep
            saveToVTK(stepVec, filename.str());
             
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
    
    //If you want to save the positions of the first 1000 time steps
    //beware though 1000 files will be created
    //sim.saveObjectPositionsVTK(1000);

    std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;

    return 0;
}