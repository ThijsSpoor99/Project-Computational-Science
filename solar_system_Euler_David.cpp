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

// check if point (3D vec) is in sphere (3D vec wiht radius r)
bool isPointInSphere(const std::array<double, 3> pVec,
                     const std::array<double, 3> cVec, 
                     double r) {
    
    // Calculate difference vector
    std::array<double, 3> dVec = subtract3D(pVec, cVec);

    // Calculate the squared distance from the point to the sphere's center
    double dSquared = calcSquare(dVec);
                             
    // Check if the squared distance is less than the squared radius
    return dSquared < std::pow(r, 2);
}

// imports csv data as 2D vector
std::vector<std::vector<std::string>> readCSV(const std::string filepath) {
    // open file
    std::ifstream file(filepath);
    std::string line;
    
    // initialize 2D vector to store data
    std::vector<std::vector<std::string>> data;

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
        data.push_back(row);
    }

    // always close the file :)
    file.close();

    return data;
}

class Celestial {
public:

    // create class variables
    std::string name;
    std::vector<std::array<double, 3>> position;
    std::array<double, 3> velocity;
    double mass;
    double GM;

    // constructor
    // & means reference
    Celestial(const std::string inputName,
              const std::array<double, 3> &startPos,
              const std::array<double, 3> &startVel,  
              double inputMass,
              const size_t nSteps
              )

              // constructor list
              : 
              name(inputName), 
              mass(inputMass), 
              GM(Constants::G * inputMass),
              velocity(startVel)

        {
        // resize vectors
        position.resize(nSteps);

        // set starting conditions
        position[0] = startPos;
        }
};

class Centaur {
public:

    // create class variables
    std::string name;
    std::vector<std::array<double, 3>> position;
    std::array<double, 3> velocity;

    // constructor
    Centaur(const std::string name,
             const std::array<double, 3> &startPos,
             const std::array<double, 3> &startVel,
             const size_t nSteps)
             :
             name(name),
             velocity(startVel)
        {
        position.resize(nSteps);
        position[0] = startPos;
        }
             
};

class SolarSystem {
public:

    // set unchanging constants
    const int nPlanets = 8;
    const int nCentaurs = 50;
    const int dt = 86400; // 1 day
    const long nSteps = 1e6;

    // vector allows dynamic resizing and memory allocation (but slower)
    Celestial Sun;
    std::vector<Celestial> planets;
    std::vector<Centaur> centaurs;

    // input data
    std::vector<std::vector<std::string>> planetData;
    std::vector<std::vector<std::string>> centaurData;

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
    double earthHits = 0;

    // constructor
    SolarSystem() 
        :
        Sun(Celestial("Sun", {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.989e30, nSteps)),
        planetData(readCSV("Data\\planetData.csv")),
        centaurData(readCSV("Data\\CentaursCartesian.csv")) 

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

        for (int i = 0; i < nCentaurs; i++) {
            centaurs.push_back(Centaur(
                centaurData[i][0],
                {std::stod(centaurData[i][1]), std::stod(centaurData[i][2]), std::stod(centaurData[i][3])},
                {std::stod(centaurData[i][4]), std::stod(centaurData[i][5]), std::stod(centaurData[i][6])},
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

            planets[i].velocity[0] = planets[i].velocity[0] + dt * aVec[0];
            planets[i].velocity[1] = planets[i].velocity[1] + dt * aVec[1];
            planets[i].velocity[2] = planets[i].velocity[2] + dt * aVec[2];

            planets[i].position[t][0] = planets[i].position[t-1][0] + dt * planets[i].velocity[0];
            planets[i].position[t][1] = planets[i].position[t-1][1] + dt * planets[i].velocity[1];
            planets[i].position[t][2] = planets[i].position[t-1][2] + dt * planets[i].velocity[2];
            }   
    }

    void updateCentaurs(int &t) {
        // calc aVec for each centaur
        for (int i = 0; i < nCentaurs; i++) {
            aVec = {0.0, 0.0, 0.0};

            // calc gravitional pull from planets
            for (int j = 0; j <nPlanets; j++) {
                rVec = subtract3D(centaurs[i].position[t-1], planets[j].position[t-1]);
                rCubed = calcCube(rVec);
                GMrCubed = planets[j].GM / rCubed;

                aVec[0] -= GMrCubed * rVec[0];
                aVec[1] -= GMrCubed * rVec[1];
                aVec[2] -= GMrCubed * rVec[2];
            }

            // calc gravitational pull from SUN
            // Stationary sun, thus pos=rVec
            rVec = centaurs[i].position[t-1];
            rCubed = calcCube(rVec);
            GMrCubed = Sun.GM / rCubed;
            aVec[0] -= GMrCubed * rVec[0];
            aVec[1] -= GMrCubed * rVec[1];
            aVec[2] -= GMrCubed * rVec[2];

            centaurs[i].velocity[0] = centaurs[i].velocity[0] + dt * aVec[0];
            centaurs[i].velocity[1] = centaurs[i].velocity[1] + dt * aVec[1];
            centaurs[i].velocity[2] = centaurs[i].velocity[2] + dt * aVec[2];

            centaurs[i].position[t][0] = centaurs[i].position[t-1][0] + dt * centaurs[i].velocity[0];
            centaurs[i].position[t][1] = centaurs[i].position[t-1][1] + dt * centaurs[i].velocity[1];
            centaurs[i].position[t][2] = centaurs[i].position[t-1][2] + dt * centaurs[i].velocity[2];
        }
    }
    void checkPlanethHits(int &t, const std::string &pName) {
        //get planet and its position at time t
        std::array<double, 3> pPos;
        bool planetFound = false;

        for (int i = 0; i < nPlanets; ++i) {
            if (planets[i].name == pName) {
                pPos = planets[i].position[t];
                planetFound = true;
                break;
            }
        }

        // if no planet with name pName was found
        if (!planetFound) {
            std::cerr << "Planet " << pName << " not found.\n";
            return;
        }
        
        
        //for all Centaurs, check if planet has been hit
        for (int i = 0; i < nCentaurs; ++i) {
            if (isPointInSphere(centaurs[i].position[t-1], pPos, 1e9)) {
                earthHits += 1;
            }
        }
    }

    void calcSystemEnergy(int &t) {
        tempEnergy = 0.0;
        for (int i = 0; i < nPlanets; i++) {

            // kinetic energy
            vSquared = calcSquare(planets[i].velocity);
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

        // calc centaur energy
        // we assume neglilible mass (m=1 for easy computing)
        for (int i = 0; i < nCentaurs; i++) {

            // kinetic energy
            vSquared = calcSquare(centaurs[i].velocity);
            tempEnergy += 0.5 * vSquared;

            // potential from centaur-planet
            for (int j = 1; j < nPlanets; j++) {
                rVec = subtract3D(centaurs[i].position[t], planets[j].position[t]);
                rNorm = calcNorm(rVec);
                tempEnergy -= planets[j].GM / rNorm;
            }
            
            // potential from centaur-sun
            rNorm = calcNorm(centaurs[i].position[t]);
            tempEnergy -= Sun.GM / rNorm;
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
            updateCentaurs(t);
            checkPlanethHits(t, std::string("Earth"));

            if (t % 100000 == 0) {
                calcSystemEnergy(t);
            }
        }
        std::cout << "Earth hits " << earthHits << std::endl;
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

    return 0;
}