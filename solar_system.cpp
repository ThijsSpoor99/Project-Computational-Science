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
    std::array<double, 3> acceleration;
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
          radiusCubed(inputRadius)
    {
        position = startPos;
        acceleration = {0.0, 0.0, 0.0};
        GM = convertGM(GM);
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

    Centaur(const std::string &name,
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
    const int nCentaurs = 903;
    const int dt = 1;        // 1 day
    const long nSteps = 1e6; // 1e6

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

    // vectors for storing energy at each checkpoint
    std::vector<double> celestialEnergy;
    std::vector<double> centaurEnergy;

    // store individual centaur energies
    std::vector<double> centaurIndividualEnergies;

    int nImpacts = 0;
    int nInner = 0;

    SolarSystem()
        : celestialData(readCSV("Data\\celestialData.csv")),
          centaurData(readCSV("Data\\centaurData.csv"))
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

    void reset() 
    {
        nImpacts = 0;
        nInner = 0;
        // create celestials
        celestials = {};
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
        centaurs = {};
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
                if (j == i)
                    continue;
                rVec = subtract3D(celestials[i].position, celestials[j].position);
                rCubed = calcCube(rVec);
                GMrCubed = celestials[j].GM / rCubed;
                celestials[i].acceleration[0] -= GMrCubed * rVec[0];
                celestials[i].acceleration[1] -= GMrCubed * rVec[1];
                celestials[i].acceleration[2] -= GMrCubed * rVec[2];
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
            if (!centaurs[i].exist)
                continue; // skip "destroyed" centaurs
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
            celestials[i].position[0] += dt * celestials[i].velocity[0];
            celestials[i].position[1] += dt * celestials[i].velocity[1];
            celestials[i].position[2] += dt * celestials[i].velocity[2];
        }

        // recompute acceleration at new positions
        for (int i = 0; i < nCelestials; i++)
        {
            celestials[i].acceleration = {0.0, 0.0, 0.0};
            // gravitational pull from other celestials
            for (int j = 0; j < nCelestials; j++)
            {
                if (j == i)
                    continue;
                rVec = subtract3D(celestials[i].position, celestials[j].position);
                rCubed = calcCube(rVec);
                GMrCubed = celestials[j].GM / rCubed;
                celestials[i].acceleration[0] -= GMrCubed * rVec[0];
                celestials[i].acceleration[1] -= GMrCubed * rVec[1];
                celestials[i].acceleration[2] -= GMrCubed * rVec[2];
            }

            // complete the velocity step
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
        // half-step velocity + position
        for (int i = 0; i < nCentaurs; i++)
        {
            rVec = subtract3D(celestials[0].position, centaurs[i].position);
            rNorm = calcNorm(rVec);

            if (!centaurs[i].exist)
                continue;

            if (rNorm < 1.5e8) {
                centaurs[i].exist = false;
                nInner += 1;
            }

            centaurs[i].velocity[0] += 0.5 * dt * centaurs[i].acceleration[0];
            centaurs[i].velocity[1] += 0.5 * dt * centaurs[i].acceleration[1];
            centaurs[i].velocity[2] += 0.5 * dt * centaurs[i].acceleration[2];

            centaurs[i].position[0] += dt * centaurs[i].velocity[0];
            centaurs[i].position[1] += dt * centaurs[i].velocity[1];
            centaurs[i].position[2] += dt * centaurs[i].velocity[2];

            // acceleration from celestials at new position
            centaurs[i].acceleration = {0.0, 0.0, 0.0};
            for (int j = 0; j < nCelestials; j++)
            {
                rVec = subtract3D(centaurs[i].position, celestials[j].position);
                rCubed = calcCube(rVec);

                // check impact
                if (rCubed < celestials[j].radiusCubed)
                {
                    centaurs[i].exist = false;
                    nImpacts += 1;
                    std::cout << nImpacts << " impact" << std::endl;
                    // remain in the loop to preserve sums, but no more velocity update
                }
                else
                {
                    GMrCubed = celestials[j].GM / rCubed;
                    centaurs[i].acceleration[0] -= GMrCubed * rVec[0];
                    centaurs[i].acceleration[1] -= GMrCubed * rVec[1];
                    centaurs[i].acceleration[2] -= GMrCubed * rVec[2];
                }
            }

            // complete velocity step if still exists
            if (centaurs[i].exist)
            {
                centaurs[i].velocity[0] += 0.5 * dt * centaurs[i].acceleration[0];
                centaurs[i].velocity[1] += 0.5 * dt * centaurs[i].acceleration[1];
                centaurs[i].velocity[2] += 0.5 * dt * centaurs[i].acceleration[2];
            }
        }
    }

    // -------------------------------------------------
    // Compute the total (kinetic + potential) energy
    // for celestials AND each centaur individually.
    // -------------------------------------------------
    void calcSystemEnergy(int &t)
    {
        //
        // 1) Celestials
        //
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

        // Centaurs
        //  - also track each centaur's energy individually in a separate vector.
        double centaurKE = 0.0;
        double centaurPE = 0.0;
        tempEnergy = 0.0;

        // resize to hold each centaur's energy
        centaurIndividualEnergies.resize(nCentaurs, 0.0);

        for (int i = 0; i < nCentaurs; i++)
        {
            if (!centaurs[i].exist)
            {
                // store 0 if "destroyed"
                centaurIndividualEnergies[i] = 0.0;
                continue;
            }

            // Kinetic = 1/2 * v^2 (assuming m=1 for each centaur)
            vSquared = calcSquare(centaurs[i].velocity);
            double thisKE = 0.5 * vSquared;
            centaurKE += thisKE;

            // Potential from centaur-planet
            double thisPE = 0.0;
            for (int j = 0; j < nCelestials; j++)
            {
                rVec = subtract3D(centaurs[i].position, celestials[j].position);
                rNorm = calcNorm(rVec);
                thisPE -= celestials[j].GM / rNorm;
            }
            centaurPE += thisPE;

            // store this centaur's total energy
            centaurIndividualEnergies[i] = thisKE + thisPE;
        }

        // aggregated total
        tempEnergy = centaurKE + centaurPE;
        centaurEnergy.push_back(tempEnergy);

        std::cout << "   (Kinetic = " << centaurKE << ", Potential = " << centaurPE << ")" << std::endl;
        std::cout << "Centaur energy at step " << t << ": " << tempEnergy << std::endl;
    }

    // -------------------------------------------------
    // Run the sim
    // -------------------------------------------------
    void simulate()
    {
        // open CSV file for writing
        std::ofstream outFile("energy_output_verlet.csv");

        outFile << "t,celestialEnergy";
        for (int i = 0; i < nCentaurs; i++)
        {
            outFile << ",centaur_" << i;
        }
        outFile << "\n";

        // compute initial accelerations at t=0
        int t = 0;
        computeCelestialsAcceleration(t);
        computeCentaursAcceleration(t);

        // check energy at t=0
        calcSystemEnergy(t);

        // write energies at t=0 (including each centaur)
        outFile << t << "," << celestialEnergy.back();
        for (int i = 0; i < nCentaurs; i++)
        {
            outFile << "," << centaurIndividualEnergies[i];
        }
        outFile << "\n";

        std::cout << "Simulating..." << std::endl;
        for (t = 1; t < nSteps; t++)
        {
            updateCelestials(t);
            updateCentaurs(t);

            // check energy every 100000 steps
            if (t % 100000 == 0)
            {
                calcSystemEnergy(t);

                // Write new row to CSV
                outFile << t << "," << celestialEnergy.back();
                for (int i = 0; i < nCentaurs; i++)
                {
                    outFile << "," << centaurIndividualEnergies[i];
                }
                outFile << "\n";
            }
        }

        // Close the output file
        outFile.close();
    }
};

int main()
{
    SolarSystem sim;
    for (double n = 0; n < 2.1; n+=0.2) {     
        sim.celestials[3].GM *= n;

        // measure runtime
        auto startTime = std::chrono::high_resolution_clock::now();
        sim.simulate();
        auto stopTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);

        std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;
        std::cout << "Inner asteroids: " << sim.nInner << " (n=" << n << ")" <<std::endl;

        sim.reset();
    }

    return 0;
}