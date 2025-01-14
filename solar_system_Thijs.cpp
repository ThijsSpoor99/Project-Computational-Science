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

// universal constants
namespace Constants {
    constexpr double G = 6.67430e-11;
}

// returns v1 - v2
std::array<double, 3> subtract3D(const std::array<double, 3> &v1, const std::array<double, 3> &v2) {
    std::array<double, 3> result;
    result[0] = v1[0] - v2[0];
    result[1] = v1[1] - v2[1];
    result[2] = v1[2] - v2[2];
    return result;
}

// returns v1 + v2
std::array<double, 3> add3D(const std::array<double, 3> &v1, const std::array<double, 3> &v2) {
    std::array<double, 3> result;
    result[0] = v1[0] + v2[0];
    result[1] = v1[1] + v2[1];
    result[2] = v1[2] + v2[2];
    return result;
}

// returns |v|
double calcNorm(const std::array<double, 3> &v) {
    double norm = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return norm;
}

// returns |v|^2
double calcSquare(const std::array<double, 3> &v) {
    double square = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    return square;
}

// returns |v|^3 
double calcCube(const std::array<double, 3> &v) {
    double square = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    double cube = square * sqrt(square);
    return cube;
}

// save a vector as .csv
void saveVector(const std::vector<double> &data, const std::string& filename) {
    const std::string directory = "Data\\";
    const std::string filepath = directory + filename;
    std::ofstream file(filepath);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filepath << " for writing." << std::endl;
        return;
    }

    for (size_t i = 0; i < data.size(); ++i) {
        file << data[i];
        if (i != data.size() - 1) {
            file << ","; 
        }
    }

    file.close();
    std::cout << "Vector written to " << filepath << " successfully." << std::endl;
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

class Celestial {
public:

    // create class attributes
    std::string name;
    std::vector<std::array<double, 3>> position;
    std::vector<std::array<double, 3>> velocity;
    double mass;
    double GM;

    // constructor
    // & means reference
    Celestial(std::string inputName,
              const std::array<double, 3> &startPos, 
              const std::array<double, 3> &startVel, 
              double inputMass,
              size_t nSteps
              )

              // constructor list
              : 
              name(inputName), 
              mass(inputMass), 
              GM(Constants::G * inputMass)
                
        {
        // resize vectors
        position.resize(nSteps);
        velocity.resize(nSteps);

        // set starting conditions
        position[0] = startPos;
        velocity[0] = startVel;
        }
};

class SolarSystem {
public:

    // set unchanging constants
    const int nPlanets = 8;
    const int dt = 86400;
    const long nSteps = 10000000;

    // vector allows dynamic resizing and memory allocation (but slower)
    Celestial Sun;
    std::vector<Celestial> planets;
    std::vector<Celestial> centaurs;

    // input data
    std::vector<std::vector<std::string>> planetData;

    // allocate memory for calculations
    std::array<double, 3> rVec;
    std::array<double, 3> aVec;
    double vSquared = 0.0;
    double rNorm = 0.0;
    double rCubed = 0.0;
    double GMrCubed = 0.0;
    double GMrSquared = 0.0;
    double tempEnergy = 0.0;
    std::vector<double> systemEnergy;

    // constructor
    SolarSystem() 
        :
        Sun(Celestial("Sun", {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.989e30, nSteps)),
        planetData(getPlanetData()) 

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

    void updatePlanets(int &t) {
        // calc aVec for each planet
        for (int i = 0; i < nPlanets; i++) {

            // calc gravitational pull from other planets
            aVec = {0.0, 0.0, 0.0};
            for (int j = 0; j < nPlanets; j++) {
                if (i != j) {
                    rVec = subtract3D(planets[i].position[t-1], planets[j].position[t-1]);
                    rCubed = calcCube(rVec);
                    GMrCubed = planets[j].GM / rCubed;

                    aVec[0] -= GMrCubed * rVec[0];
                    aVec[1] -= GMrCubed * rVec[1];
                    aVec[2] -= GMrCubed * rVec[2];
                    }
                }

            // calc gravitational pull from SUN
            // Stationary sun, thus pos=rVec
            rVec = planets[i].position[t-1];
            rCubed = calcCube(rVec);
            GMrCubed = Sun.GM / rCubed;
            aVec[0] -= GMrCubed * rVec[0];
            aVec[1] -= GMrCubed * rVec[1];
            aVec[2] -= GMrCubed * rVec[2];

            planets[i].velocity[t][0] = planets[i].velocity[t-1][0] + dt * aVec[0];
            planets[i].velocity[t][1] = planets[i].velocity[t-1][1] + dt * aVec[1];
            planets[i].velocity[t][2] = planets[i].velocity[t-1][2] + dt * aVec[2];

            planets[i].position[t][0] = planets[i].position[t-1][0] + dt * planets[i].velocity[t][0];
            planets[i].position[t][1] = planets[i].position[t-1][1] + dt * planets[i].velocity[t][1];
            planets[i].position[t][2] = planets[i].position[t-1][2] + dt * planets[i].velocity[t][2];
            }   
    }

    void calcSystemEnergy(int &t) {
        tempEnergy = 0.0;
        for (int i = 0; i < nPlanets; i++) {

            // kinetic energy
            vSquared = calcSquare(planets[i].velocity[t]);
            tempEnergy += 0.5 * planets[i].mass * vSquared;

            // potential from planet-planet
            // gravity is pair-wise (avoid double counting)
            for (int j = i + 1; j < nPlanets; j++) {
                rVec = subtract3D(planets[i].position[t], planets[j].position[t]);
                rNorm = calcNorm(rVec);
                tempEnergy -= planets[i].GM * planets[j].mass / rNorm;
            }
            
            // potential from planet-sun
            rNorm = calcNorm(planets[i].position[t]);
            tempEnergy -= Sun.GM * planets[i].mass / rNorm;

        }
        std::cout << tempEnergy << std::endl;
        systemEnergy.push_back(tempEnergy);
    }

    void simulate() {
        int t = 0;
        calcSystemEnergy(t);
        std::cout << "Simulating..." << std::endl;
        for (int t=1; t < nSteps; t++) {
            updatePlanets(t);

            if (t % 100000 == 0) {
                calcSystemEnergy(t);
            }
        }
    }
};

int main() {
    SolarSystem sim;

    // time simulation
    auto startTime = std::chrono::high_resolution_clock::now();
    sim.simulate();
    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);
    std::cout << duration.count() << " ms" << std::endl;

    saveVector(sim.systemEnergy, "systemEnergy.csv");

    // sim.savePlanetSimulation();

    return 0;
}