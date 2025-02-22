#ifndef solarSystem_HPP
#define solarSystem_HPP

#include "Include\util.hpp"

// -------------------------------------------------
// The class for sun and planets
// -------------------------------------------------
class Celestial
{
public:
    std::string name;
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    std::array<double, 3> acceleration = {0.0, 0.0, 0.0};
    double mass;
    double GM;
    double radius;

    double kineticEnergy = 0.0;
    double potentialEnergy = 0.0;
    double totalEnergy = 0.0;

    Celestial() = default;
    Celestial(const std::string &inputName,
              const std::array<double, 3> &startPos,
              const std::array<double, 3> &startVel,
              double inputMass,
              double inputRadius) {
        name = inputName;
        position = startPos;
        velocity = startVel;
        mass = inputMass;
        GM = Constants::G * inputMass;
        GM = convertGM(GM);
        radius = 1e6;
    }
};

// -------------------------------------------------
// The class for the Centaur asteroids
// -------------------------------------------------
class Centaur
{
public:
    std::string name;
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    std::array<double, 3> acceleration = {0.0, 0.0, 0.0}; // for Velocity Verlet
    bool exist = true;
    bool inner = false;
    bool outer = false;

    double kineticEnergy = 0.0;
    double potentialEnergy = 0.0;
    double totalEnergy = 0.0;

    double initialEnergy = 0.0;

    Centaur() = default;
    Centaur(const std::string &inputName,
            const std::array<double, 3> &startPos,
            const std::array<double, 3> &startVel) {
        name = inputName;
        position = startPos;
        velocity = startVel;
    }
};

// -------------------------------------------------
// The simulation class. Only stores the information for the current timestep. 
// Simulation itself is preferably called in a different c++ file which also saves the information of all steps (if necessary).
// -------------------------------------------------
class SolarSystem
{
public:
    std::string pathToData;
    int nCelestials;
    int nCentaurs;
    int dt = 1;
    double bias;

    // vector to hold simulated objects
    std::vector<Celestial> celestials;
    std::vector<Centaur> centaurs;
    
    // vector to store imported initial state data
    std::vector<std::vector<std::string>> celestialData;
    std::vector<std::vector<std::string>> centaurData;

    // temporary vectors used in computations
    std::array<double, 3> rVec;
    std::array<double, 3> aVec;

    // temporary doubles used in computations
    double vSquared = 0.0;
    double rNorm = 0.0;
    double rCubed = 0.0;
    double GMrCubed = 0.0;
    double tempEnergy = 0.0;

    // total energy of the (sub)system
    // in current model all 3 should be constant
    double celestialEnergy = 0.0;
    double centaurEnergy = 0.0;
    double totalEnergy = 0.0;

    // track hypothesis information
    int nImpacts = 0;
    int nInner = 0;
    int nOuter = 0;

    // boundary in km
    double innerBoundary = 2 * Constants::AU * 1.0e-3;
    double outerBoundary = 1000.0 * Constants::AU * 1.0e-3;

    // constructor
    SolarSystem() = default;
    SolarSystem(const std::string& inputPath = "Data\\", int inputNCentaurs = 24375, double inputBias = 0.8) 
        : pathToData(inputPath), nCentaurs(inputNCentaurs), bias(inputBias) {

        // create celestials
        celestialData = readCSV(pathToData + "celestialData.csv");
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

        // calculate acceleration at t=0
        for (int i = 0; i < nCelestials; i++) {
            computeCelestialAcceleration(i);
        }

        // create centaurs
        centaurData = readCSV(pathToData + "centaurData.csv");

        if (centaurData.size() < nCentaurs) {
            throw std::runtime_error("Not enough data in centaurData.csv to initialize all centaurs.");
        }

        for (int i = 0; i < nCentaurs; i++)
        {
            centaurs.push_back(Centaur(
                centaurData[i][0],
                {std::stod(centaurData[i][1]),
                 std::stod(centaurData[i][2]),
                 std::stod(centaurData[i][3])},
                {std::stod(centaurData[i][4]) * bias,
                 std::stod(centaurData[i][5]) * bias,
                 std::stod(centaurData[i][6]) * bias}));

            // calculate acceleration at t=0
            computeCentaurAcceleration(i);
        }

        // calculate all energies at t=0
        calcSystemEnergy();
        
        // store initial energy for impacted (removed) asteroids
        // close encounters with planets do not conserve energy with the used integrator
        for (int i = 0; i < nCentaurs; i++) {
            centaurs[i].initialEnergy = centaurs[i].totalEnergy;
        }
    }

    // -------------------------------------------------
    // Reset the simulation to the initial state (similar to constructor)
    // -------------------------------------------------
    void resetSimulation() {

        // reset celestials
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

        // calculate acceleration at t=0
        for (int i = 0; i < nCelestials; i++) {
            computeCelestialAcceleration(i);
        }

        // reset centaurs
        centaurs = {};
        for (int i = 0; i < nCentaurs; i++)
        {
            centaurs.push_back(Centaur(
                centaurData[i][0],
                {std::stod(centaurData[i][1]),
                 std::stod(centaurData[i][2]),
                 std::stod(centaurData[i][3])},
                {std::stod(centaurData[i][4]) * bias,
                 std::stod(centaurData[i][5]) * bias,
                 std::stod(centaurData[i][6]) * bias}));
            
            // calculate acceleration at t=0
            computeCentaurAcceleration(i);
        }

        // calculate all energies at t=0
        calcSystemEnergy();
        for (int i = 0; i < nCentaurs; i++) {
            centaurs[i].initialEnergy = centaurs[i].totalEnergy;
        }

        // reset hypothesis counters
        nInner = 0;
        nOuter = 0;
        nImpacts = 0;

    }

    // -------------------------------------------------
    // Compute the net acceleration for the given celestial index
    // -------------------------------------------------
    void computeCelestialAcceleration(int &i) {
        // acceleration from other celestials
        celestials[i].acceleration = {0.0, 0.0, 0.0};
        for (int j=0; j<nCelestials; j++) {
            
            // celestial does not pull on itself
            if (i==j) {
                continue;
            }

            // acceleration from other celestials
            rVec = subtract3D(celestials[i].position, celestials[j].position);
            rCubed = calcCube(rVec);
            GMrCubed = celestials[j].GM / rCubed;
            celestials[i].acceleration[0] -= GMrCubed * rVec[0];
            celestials[i].acceleration[1] -= GMrCubed * rVec[1];
            celestials[i].acceleration[2] -= GMrCubed * rVec[2];
        }
    }

    // -------------------------------------------------
    // Compute the net acceleration for the given centaur index
    // -------------------------------------------------
    void computeCentaurAcceleration(int &i) {

        // acceleration from celestials
        centaurs[i].acceleration = {0.0, 0.0, 0.0};
        for (int j=0; j<nCelestials; j++) {
            rVec = subtract3D(centaurs[i].position, celestials[j].position);
            rNorm = calcNorm(rVec);

            // check for impact
            if (rNorm < celestials[j].radius) {
                centaurs[i].exist = false;
                nImpacts += 1;
                centaurs[i].totalEnergy = centaurs[i].initialEnergy;
                std::cout << "Centaur (i=" << i << ") impact with " << celestials[j].name << std::endl;
            }

            // update acceleration
            GMrCubed = celestials[j].GM / (rNorm * rNorm * rNorm);
            centaurs[i].acceleration[0] -= GMrCubed * rVec[0];
            centaurs[i].acceleration[1] -= GMrCubed * rVec[1];
            centaurs[i].acceleration[2] -= GMrCubed * rVec[2];
        }
    }

    // -------------------------------------------------
    // Velocity Verlet update for celestials
    // -------------------------------------------------
    void updateCelestials() {
        for (int i=0; i<nCelestials; i++)
        {
            // first half velocity step
            celestials[i].velocity[0] += 0.5 * dt * celestials[i].acceleration[0];
            celestials[i].velocity[1] += 0.5 * dt * celestials[i].acceleration[1];
            celestials[i].velocity[2] += 0.5 * dt * celestials[i].acceleration[2];

            // full position step
            celestials[i].position[0] += dt * celestials[i].velocity[0];
            celestials[i].position[1] += dt * celestials[i].velocity[1];
            celestials[i].position[2] += dt * celestials[i].velocity[2];
        }

        for (int i=0; i<nCelestials; i++) {

            // calculate acceleration after full position step
            computeCelestialAcceleration(i);

            // second half velocity step
            celestials[i].velocity[0] += 0.5 * dt * celestials[i].acceleration[0];
            celestials[i].velocity[1] += 0.5 * dt * celestials[i].acceleration[1];
            celestials[i].velocity[2] += 0.5 * dt * celestials[i].acceleration[2];
        }
    }

    // -------------------------------------------------
    // Velocity Verlet update for centaurs
    // -------------------------------------------------
    void updateCentaurs() {
        for (int i=0; i<nCentaurs; i++) {

            // check if centaur still needs to be updated
            if (!centaurs[i].exist) {
                continue;
            }

            // check if centaur has entered inner region of solar system
            rVec = subtract3D(centaurs[i].position, celestials[0].position);
            rNorm = calcNorm(rVec);
            if (rNorm < innerBoundary) {
                centaurs[i].inner = true;
                centaurs[i].exist = false;
                nInner += 1;
            }

            // check if centaur has 'left' the solar system
            if (rNorm > outerBoundary) {
                centaurs[i].outer = true;
                centaurs[i].exist = false;
                nOuter += 1;
            }

            // first half velocity step
            centaurs[i].velocity[0] += 0.5 * dt * centaurs[i].acceleration[0];
            centaurs[i].velocity[1] += 0.5 * dt * centaurs[i].acceleration[1];
            centaurs[i].velocity[2] += 0.5 * dt * centaurs[i].acceleration[2];

            // full position step
            centaurs[i].position[0] += dt * centaurs[i].velocity[0];
            centaurs[i].position[1] += dt * centaurs[i].velocity[1];
            centaurs[i].position[2] += dt * centaurs[i].velocity[2];

            // calculate acceleration after full position step
            // does not depend on other centaurs
            // thus can be done within same loop
            computeCentaurAcceleration(i);

            // second half velocity step
            centaurs[i].velocity[0] += 0.5 * dt * centaurs[i].acceleration[0];
            centaurs[i].velocity[1] += 0.5 * dt * centaurs[i].acceleration[1];
            centaurs[i].velocity[2] += 0.5 * dt * centaurs[i].acceleration[2];
        }
    }

    // -------------------------------------------------
    // Compute the energy for the given celestial index
    // -------------------------------------------------
    void calcCelestialEnergy(int &i) {

        // kinetic
        celestials[i].kineticEnergy = 0.5 * celestials[i].mass * calcSquare(celestials[i].velocity);
        
        // potential
        celestials[i].potentialEnergy = 0.0;
        for (int j=0; j<nCelestials; j++) {

            // no potential with itself
            if (i==j) {
                continue;
            }

            // potential from other celestials
            rVec = subtract3D(celestials[i].position, celestials[j].position);
            rNorm = calcNorm(rVec);
            celestials[i].potentialEnergy -= celestials[i].mass * celestials[j].GM / rNorm;
        }

        // total = kinetic + potential
        celestials[i].totalEnergy = celestials[i].kineticEnergy + celestials[i].potentialEnergy;
    }

    // -------------------------------------------------
    // Compute the energy for the given centaur index
    // -------------------------------------------------
    void calcCentaurEnergy(int &i) {

        // kinetic
        // centaurs have neglible mass (m=1 for easy calculation)    
        centaurs[i].kineticEnergy = 0.5 * calcSquare(centaurs[i].velocity);
        
        // potential from the celestials
        centaurs[i].potentialEnergy = 0.0;
        for (int j=0; j<nCelestials; j++) {
            rVec = subtract3D(centaurs[i].position, celestials[j].position);
            rNorm = calcNorm(rVec);
            centaurs[i].potentialEnergy -= celestials[j].GM / rNorm;
        }

        // total = kinetic + potential
        centaurs[i].totalEnergy = centaurs[i].kineticEnergy + centaurs[i].potentialEnergy;
    }

    // -------------------------------------------------
    // Compute the energy of the system (celestials, centaurs & total)
    // -------------------------------------------------
    void calcSystemEnergy() {

        // energy from celestials
        celestialEnergy = 0.0;
        for (int i=0; i<nCelestials; i++) {
            calcCelestialEnergy(i);
            celestialEnergy += celestials[i].totalEnergy;
        }

        // energy from centaurs
        centaurEnergy = 0.0;
        for (int i=0; i<nCentaurs; i++) {
            if (!centaurs[i].exist) {
                centaurEnergy += centaurs[i].totalEnergy;
                continue;
            }
            calcCentaurEnergy(i);
            centaurEnergy += centaurs[i].totalEnergy;
        }

        // as the centaurs do not have an effect on the celestials (neglible mass)
        // we expect both celestial and centaur energy to be conserved
        // as well as the total
        totalEnergy = celestialEnergy + centaurEnergy;
    } 

    // -------------------------------------------------
    // Updates the position, velocity and acceleration of celestials & centaurs.
    // -------------------------------------------------
    void performTimestep() {
        updateCelestials();
        updateCentaurs();
    }

};

#endif // solarSystem_HPP