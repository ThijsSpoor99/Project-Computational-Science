// solar_system_project_Joel.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream> // for file streams
#include <sstream> // for string streams
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <cassert>
#include <chrono>

// ------------------------------------------- Constants -------------------------------------------

constexpr double G = 6.674e-11; // gravitational constant

// ------------------------------------------- Helper functions -------------------------------------------
inline double norm3D(const std::array<double, 3> &v)
{
    // Euclidean norm of a 3D vector
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// subtract two 3D vectors
inline std::array<double, 3> operator-(const std::array<double, 3> &a,
                                       const std::array<double, 3> &b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

// add two 3D vectors
inline std::array<double, 3> operator+(const std::array<double, 3> &a,
                                       const std::array<double, 3> &b)
{
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}
// add two 3D vectors in place
inline std::array<double, 3> &operator+=(std::array<double, 3> &a,
                                         const std::array<double, 3> &b)
{
    a[0] += b[0];
    a[1] += b[1];
    a[2] += b[2];
    return a;
}

// multiply 3D vector by a scalar
inline std::array<double, 3> operator*(double dt, const std::array<double, 3> &v)
{
    return {dt * v[0], dt * v[1], dt * v[2]};
}

// Acceleration due to gravity from ANOTHER body of mass M at position "otherPos"
// exerted on an object at CURRENT position "posObj"
std::array<double, 3> gravitationalAcceleration(
    const std::array<double, 3> &posObj,
    const std::array<double, 3> &otherPos,
    double otherMass)
{
    std::array<double, 3> rVec = posObj - otherPos;
    double r = norm3D(rVec);
    if (r <= 0.0)
    {
        // Avoid division by zero
        std::cerr << "Warning: zero distance in gravitationalAcceleration.\n";
        return {0.0, 0.0, 0.0};
    }
    double factor = -G * otherMass / (r * r * r);
    return {factor * rVec[0], factor * rVec[1], factor * rVec[2]};
}

// ------------------------------------------- Celestial class -------------------------------------------

class Celestial
{
public:
    Celestial(size_t Nsteps)
        : N(Nsteps)
    {
        positions.resize(N, {0.0, 0.0, 0.0});
        velocities.resize(N, {0.0, 0.0, 0.0});
    }

    // Set initial conditions for the body -> like __init__ in Python
    void setInitialConditions(const std::array<double, 3> &startPos,
                              const std::array<double, 3> &startVel,
                              double massVal)
    {
        positions[0] = startPos;
        velocities[0] = startVel;
        mass = massVal;
    }

    // Update the i-th step, given a list of other bodies that influence this one
    void update(const std::vector<Celestial *> &otherBodies, size_t i, double dt)
    {
        // Start from velocity at step i-1
        velocities[i] = velocities[i - 1];
        // Accumulate gravitational accelerations from each "other" body
        for (auto *other : otherBodies)
        {
            std::array<double, 3> acc = gravitationalAcceleration(
                positions[i - 1], other->positions[i - 1], other->mass);
            // Update velocity using Euler's method
            velocities[i][0] += dt * acc[0];
            velocities[i][1] += dt * acc[1];
            velocities[i][2] += dt * acc[2];
        }
        // Update position
        positions[i] = positions[i - 1] + (dt * velocities[i]);
    }

    std::vector<std::array<double, 3>> positions;  // positions[i] -> x,y,z
    std::vector<std::array<double, 3>> velocities; // velocities[i] -> vx,vy,vz
    double mass{0.0};
    size_t N{0};
};

// ------------------------------------------- CSV loading helper -------------------------------------------

struct PlanetData
{
    std::string name;
    double x, y, z;
    double vx, vy, vz;
    double M;
};

// simple CSV loader -> expects columns: x,y,z,Vx,Vy,Vz,M
// and the row label is stored as "planet name" in the first column.
std::vector<PlanetData> loadPlanetData(const std::string &filename)
{
    std::vector<PlanetData> data;
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Could not open file: " << filename << std::endl;
        return data;
    }

    // CSV has a header, so we discard
    std::string header;
    if (std::getline(file, header))
    {
    }

    // read rows
    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        PlanetData pd;
        std::string planetName;
        std::getline(ss, planetName, ','); // read up to comma
        pd.name = planetName;

        std::string val;
        // read x
        std::getline(ss, val, ',');
        pd.x = std::stod(val);
        // read y
        std::getline(ss, val, ',');
        pd.y = std::stod(val);
        // read z
        std::getline(ss, val, ',');
        pd.z = std::stod(val);
        // read vx
        std::getline(ss, val, ',');
        pd.vx = std::stod(val);
        // read vy
        std::getline(ss, val, ',');
        pd.vy = std::stod(val);
        // read vz
        std::getline(ss, val, ',');
        pd.vz = std::stod(val);
        // read M
        std::getline(ss, val, ',');
        pd.M = std::stod(val);

        data.push_back(pd);
    }
    return data;
}

// ------------------------------------------- SolarSystem class -------------------------------------------

class SolarSystem
{
public:
    // pass in how many steps (N) and timestep (dt).
    SolarSystem(size_t Nsteps, double dtSec, const std::string &csvFile)
        : N(Nsteps), dt(dtSec)
    {
        // load planet data
        std::vector<PlanetData> planetInfo = loadPlanetData(csvFile);

        // for each planet, create Celestial object
        for (const auto &p : planetInfo)
        {
            Celestial c(N);
            // set initial conditions
            std::array<double, 3> startPos = {p.x, p.y, p.z};
            std::array<double, 3> startVel = {p.vx, p.vy, p.vz};
            c.setInitialConditions(startPos, startVel, p.M);
            // add planet to list of bodies
            bodies.push_back(c);
            // keep track of planet names
            names.push_back(p.name);
        }

        // create Sun manually
        Celestial sun(N);
        sun.setInitialConditions({0.0, 0.0, 0.0},
                                 {0.0, 0.0, 0.0},
                                 1.989e30);
        bodies.push_back(sun);
        names.push_back("Sun");
    }

    // ------------------------------------------- update the positions of all bodies in the system -------------------------------------------

    void simulate()
    {
        // list of pointers to pass to update() calls so each body can see the others
        for (size_t i = 1; i < N; ++i)
        {
            // For each body EXCEPT the last one (the Sun)
            for (size_t b = 0; b < bodies.size() - 1; ++b)
            {
                // construct a vector of other bodies
                std::vector<Celestial *> others;
                others.reserve(bodies.size() - 1);
                // skip the b-th body
                for (size_t ob = 0; ob < bodies.size(); ++ob)
                {
                    if (ob == b)
                        continue; // skip self
                    others.push_back(&bodies[ob]);
                }
                // update the b-th body using the "others"
                bodies[b].update(others, i, dt);
            }
        }
    }

    // write results to CSV to plot in Python
    void writeResults(const std::string &outputFile) const
    {
        std::ofstream out(outputFile);
        if (!out.is_open())
        {
            std::cerr << "Failed to open " << outputFile << " for writing.\n";
            return;
        }
        // write header
        out << "Step,Planet,x,y,z\n";
        for (size_t i = 0; i < N; ++i)
        {
            for (size_t b = 0; b < bodies.size(); ++b)
            {
                out << i << ","
                    << names[b] << ","
                    << bodies[b].positions[i][0] << ","
                    << bodies[b].positions[i][1] << ","
                    << bodies[b].positions[i][2] << "\n";
            }
        }
        out.close();
    }

private:
    size_t N;
    double dt;
    std::vector<Celestial> bodies;
    std::vector<std::string> names;
};

// ------------------------------------------- main() -------------------------------------------

int main()
{
    // 100 years with 1-day time steps
    size_t N = 10000000;
    double dayInSeconds = 60.0 * 60.0 * 24.0;
    double dt = dayInSeconds;

    // path to planetData csv file
    std::string csvFile = "../Data/planetData.csv";

    // create SolarSystem
    SolarSystem sim(N, dt, csvFile);

    // Run simulation
    auto startT = std::chrono::high_resolution_clock::now();
    sim.simulate();
    auto endT = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endT - startT);
    std::cout << duration.count() << " ms" << std::endl;

    // write data to plot in Python
    // sim.writeResults("solar_system_output.csv");

    // std::cout << "Simulation completed! Results saved to solar_system_output.csv\n";
    return 0;
}