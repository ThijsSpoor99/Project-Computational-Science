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

#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkStringArray.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataWriter.h>

#include "..\Include\util.hpp"

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
        GM = convertGM(GM);
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
    const int nCentaurs = 100;
    const int dt = 1;        // 1 day
    const long nSteps = 1e5; // 1e6

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
        : celestialData(readCSV("..\\..\\..\\Data\\celestialData.csv")), centaurData(readCSV("..\\..\\..\\Data\\centaurData.csv"))
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

    void saveObjectPositionsVTK(int nSteps) {
        // Create the VTK directory
        std::string vtkDir = "..\\..\\..\\Data\\vtpData";
        clearAndCreateDirectory(vtkDir);

        // Determine number of digits for file naming
        int digits = std::to_string(nSteps).length();

        for (int t = 0; t < nSteps; ++t) {
            // Create VTK objects
            auto points = vtkSmartPointer<vtkPoints>::New();
            auto radii = vtkSmartPointer<vtkFloatArray>::New();
            auto colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
            auto labels = vtkSmartPointer<vtkStringArray>::New();

            // Initialize arrays
            radii->SetName("radius");
            colors->SetName("colors");
            colors->SetNumberOfComponents(3); // RGB
            labels->SetName("labels");

            int celColors[7][3] = {
                {255, 204, 51}, //sun color
                {79,76,176}, //earth color
                {193,68,14}, //mars color
                {201,144,57}, //jupiter color
                {226,191,125}, //saturn color
                {147,205,241}, //uranus color
                {61,94,249} //neptune color
            };

            double celestialRadii[7] = {
                696340.0,  // Sun radius in km
                6371.0,    // Earth radius in km
                3389.5,    // Mars radius in km
                69911.0,   // Jupiter radius in km
                58232.0,   // Saturn radius in km
                25362.0,   // Uranus radius in km
                24622.0    // Neptune radius in km
            };

            // Add celestial object data
            for (int i = 0; i < nCelestials; i++) {
                points->InsertNextPoint(planets[i].position[t].data());
                radii->InsertNextValue(celestialRadii[i]);
                colors->InsertNextTuple3(celColors[i][0], celColors[i][1], celColors[i][2]); // Example color: cyan
                labels->InsertNextValue(planets[i].name);
            }

            // Add centaur data
            for (const auto& centaur : centaurs) {
                points->InsertNextPoint(centaur.position[t].data());
                radii->InsertNextValue(1000); // Example small radius for centaurs
                colors->InsertNextTuple3(255, 0, 0); // Example color: red
                labels->InsertNextValue(centaur.name);
            }

            // Create polydata and attach arrays
            auto polyData = vtkSmartPointer<vtkPolyData>::New();
            polyData->SetPoints(points);
            polyData->GetPointData()->AddArray(radii);
            polyData->GetPointData()->AddArray(colors);
            polyData->GetPointData()->AddArray(labels);

            // Write to VTK file
            std::ostringstream filename;
            filename << vtkDir << "\\out_" 
                    << std::setw(digits) 
                    << std::setfill('0') << t << ".vtp";

            auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
            writer->SetFileName(filename.str().c_str());
            writer->SetInputData(polyData);
            writer->Write();
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

    //sim.savePlanetPositions();
    
    //If you want to save the positions of the first 1000 time steps
    //beware though 1000 files will be created
    sim.saveObjectPositionsVTK(3000);

    std::cout << "Total simulation time: " << duration.count() << " ms" << std::endl;

    return 0;
}