// more types
#include <array>
#include <vector>
#include <string>

// more streams
#include <iostream>
#include <fstream>
#include <sstream>
#include "Include\json.hpp"
#include "Include/json.hpp"
using json = nlohmann::json;
#include <filesystem>

// math stuff
#include <cmath>

// universal constants
namespace Constants {
    constexpr double G = 6.67430e-11;
}

// how is this readable, look at the length man
std::vector<std::vector<std::string>> getPlanetData() {
    // open file
    std::ifstream file("Data/planetData.csv");
    std::string line;
    
    // initialize 2D vector to store data
    std::vector<std::vector<std::string>> planetData;

    // skip first row
    std::getline(file, line);

    // continue untill all rows read
    while (std::getline(file, line)) {

        // convert row to stringstream to split up row
        std::stringstream ss(line);

        // temp vector for split up row
        std::string value;
        std::vector<std::string> row;
        
        // continue untill full row is read
        while (std::getline(ss, value, ',')) {

            // add value to temp row vector
            row.push_back(value); 
        }
        
        // add row to 2D vector
        planetData.push_back(row);
    }

    // always close the file :)
    file.close();

    return planetData;
}

std::array<double, 3> acceleration(std::array<double, 3> rVec, double& GM) {
    std::array<double, 3> aVec;

    double rSquared = rVec[0] * rVec[0] + rVec[1] * rVec[1] + rVec[2] * rVec[2];
    double rCubed = rSquared * sqrt(rSquared);

    aVec[0] = -GM * rVec[0] / rCubed;
    aVec[1] = -GM * rVec[1] / rCubed;
    aVec[2] = -GM * rVec[2] / rCubed;

    return aVec;
}

class Celestial {
public:

    // create publicly accessible objects
    std::vector<std::array<double, 3>> position;
    std::vector<std::array<double, 3>> velocity;
    double mass;
    double GM;

    std::string name;

    // constructor
    // & means reference
    Celestial(std::string inputName,
            const std::array<double, 3>& startPos, 
            const std::array<double, 3>& startVel, 
            double inputMass,
            size_t nSteps
            )

            // constructor list
            : 
            name(inputName), 
            mass(inputMass), 
            GM(Constants::G * inputMass)
                
    {
        // resize arrays
        position.resize(nSteps);
        velocity.resize(nSteps);

        // set starting conditions
        position[0] = startPos;
        velocity[0] = startVel;
    }

    // convert to json style
    json toJson() const {
        return json{{"name", name}, 
                    {"position", position},
                    {"velocity", velocity},
                    {"mass", mass}}; 
    }
};

class SolarSystem {
public:

    // set variables that should not change
    const int nPlanets = 8;
    const int dt = 86400;
    const int nSteps = 365*1000;

    // vector allows dynamic resizing and memory allocation (but slower)
    Celestial Sun;
    std::vector<Celestial> planets;
    std::vector<Celestial> centaurs;

    std::vector<std::vector<std::string>> planetData = getPlanetData();

    // constructor
    SolarSystem() 
        :Sun(Celestial("Sun", {0, 0, 0}, {0, 0, 0}, 1.989e30, nSteps)) 
    {
        for (int i = 0; i < nPlanets; i++) {
            planets.push_back(Celestial(
                planetData[i][0],
                {std::stod(planetData[i][1]), std::stod(planetData[i][2]), std::stod(planetData[i][3])},
                {std::stod(planetData[i][4]), std::stod(planetData[i][5]), std::stod(planetData[i][6])},
                std::stod(planetData[i][7]),
                nSteps
                )   
            );
        }
    }

    void updatePlanets(int t) {
        for (int i = 0; i < nPlanets; i++) {
            for (int j = 0; j < nPlanets; j++) {
                if (i != j) {
                    std::array<double, 3> rVec;
                    rVec[0] = planets[i].position[0][0] - planets[j].position[0][0]; 
                    rVec[1] = planets[i].position[0][1] - planets[j].position[0][1]; 
                    rVec[2] = planets[i].position[0][2] - planets[j].position[0][2]; 

                    std::array<double, 3> aVec = acceleration(rVec, planets[j].GM);

                    planets[i].velocity[t][0] = planets[i].velocity[t-1][0] + dt * aVec[0];
                    planets[i].velocity[t][1] = planets[i].velocity[t-1][1] + dt * aVec[1];
                    planets[i].velocity[t][2] = planets[i].velocity[t-1][2] + dt * aVec[2];

                    planets[i].position[t][0] = planets[i].position[t-1][0] + dt * planets[i].velocity[t][0];
                    planets[i].position[t][1] = planets[i].position[t-1][1] + dt * planets[i].velocity[t][1];
                    planets[i].position[t][2] = planets[i].position[t-1][2] + dt * planets[i].velocity[t][2];
                }
        }
    }
    }

    void simulate() {
        for (int t=1; t < nSteps + 1; t++) {
            updatePlanets(t);
        }
    }

    void savePlanetSimulation() {
        std::vector<json> planetsJson; 
        std::cout << nPlanets << std::endl;

        for (int i = 0; i < nPlanets; i++) {
            std::cout << "Processing planet " << i + 1 << "..." << std::endl;
            planetsJson.push_back(planets[i].toJson());
        }

        std::ofstream file("Data\\planets.json");
        if (file.is_open()) {
            json outputJson = planetsJson;
            file << outputJson.dump(4);
            file.close();
            std::cout << "Planets saved to planets.json" << std::endl;
        } else {
            std::cerr << "Error opening file!" << std::endl;
    }
    }
};

int main() {
    SolarSystem sim;

    sim.simulate();

    sim.savePlanetSimulation();

    return 0;
}