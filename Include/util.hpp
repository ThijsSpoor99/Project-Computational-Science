#ifndef UTIL_HPP
#define UTIL_HPP

#include <array>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <iomanip> 

// Universal constants
namespace Constants {
    constexpr double G = 6.67430e-11;
    constexpr double AU = 1.495978707e11;
    constexpr double day = 8.64e4;
}

// Function definitions

// Converts GM from m^3/s^2 to AU^3/day^2
inline double GMtoAstro(const double &GM) {
    double factor = (8.6400e4 * 8.6400e4) / (Constants::AU * Constants::AU * Constants::AU);
    return factor * GM;
}

// Converts GM from m^3/s^2 to km^3/day^2
inline double convertGM(const double &GM) {
    double factor = (Constants::day * Constants::day) / (1e9);
    return factor * GM; 
}

// Converts a 3D vector from meters to astronomical units
inline std::array<double, 3> distToAU3D(const std::array<double, 3> &v) {
    return {v[0] / Constants::AU, v[1] / Constants::AU, v[2] / Constants::AU};
}

// Converts a 3D velocity vector from m/s to AU/day
inline std::array<double, 3> velToAstro3D(const std::array<double, 3> &v) {
    return {Constants::day * v[0] / Constants::AU, Constants::day * v[1] / Constants::AU, Constants::day * v[2] / Constants::AU};
}

// Returns the difference of two 3D vectors
inline std::array<double, 3> subtract3D(const std::array<double, 3> &v1, const std::array<double, 3> &v2) {
    return {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
}

// Returns the sum of two 3D vectors
inline std::array<double, 3> add3D(const std::array<double, 3> &v1, const std::array<double, 3> &v2) {
    return {v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]};
}

// Returns the norm (magnitude) of a 3D vector
inline double calcNorm(const std::array<double, 3> &v) {
    return std::hypot(v[0], v[1], v[2]);
}

// Returns the square of the norm of a 3D vector
inline double calcSquare(const std::array<double, 3> &v) {
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

// Returns the cube of the norm of a 3D vector
inline double calcCube(const std::array<double, 3> &v) {
    double square = calcSquare(v);
    return square * std::sqrt(square);
}

inline // imports csv data as 2D vector
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

inline void saveToCSV(const std::vector<std::array<double, 3>>& data, const std::string& filepath) {
    std::ofstream file(filepath);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file with path: " << filepath << std::endl;
        return;
    }

    // Loop through the rows and columns to write the data
    for (size_t i = 0; i < data.size(); i+=1) {
        for (size_t j = 0; j < data[i].size(); j+=1) {
            file << data[i][j];
            if (j != data[i].size() - 1) {
                file << ",";  // Add comma between elements in a row
            }
        }
        file << "\n";  // New line after each row
    }

    file.close();
    std::cout << "Data successfully written to " << filepath << std::endl;
}

inline void save1DVector(const std::vector<double> &data, const std::string &filepath) {
    std::ofstream file(filepath);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file with path: " << filepath << std::endl;
        return;
    }

    // Use a string stream for efficient string construction
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(15);

    for (size_t i = 0; i < data.size(); ++i) {
        oss << data[i];
        if (i != data.size() - 1) {
            oss << ",";
        }
    }
    
    file << oss.str() << std::endl; // Write the constructed string to the file

    file.close();
    std::cout << "Vector successfully written to " << filepath << std::endl;
}

inline void save1DIntVector(const std::vector<int> &data, const std::string &filepath) {
    std::ofstream file(filepath);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file with path: " << filepath << std::endl;
        return;
    }

    // Use a string stream for efficient string construction
    std::ostringstream oss;

    for (size_t i = 0; i < data.size(); ++i) {
        oss << data[i];
        if (i != data.size() - 1) {
            oss << ",";
        }
    }
    
    file << oss.str() << std::endl; // Write the constructed string to the file

    file.close();
    std::cout << "Vector successfully written to " << filepath << std::endl;
}

inline void save2DVector(const std::vector<std::vector<double>>& data, const std::string& filepath) {
    std::ofstream file(filepath);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file with path: " << filepath << std::endl;
        return;
    }

    // Loop through the rows and columns to write the data
    for (size_t i = 0; i < data.size(); i+=1) {
        for (size_t j = 0; j < data[i].size(); j+=1) {
            file << data[i][j];
            if (j != data[i].size() - 1) {
                file << ",";  // Add comma between elements in a row
            }
        }
        file << "\n";  // New line after each row
    }

    file.close();
    std::cout << "Data successfully written to " << filepath << std::endl;
}

/* Function to write a single `.vtk` file, input is a vector which contains tripples (xyz) of all the
 * solar system objects
 *
 *  data: Vector of triples which are the xyz coordinates of all the celestial at a given time step
 *  filepath: string which states the name of the file to which it should be written
 * 
 * Side effects: creates a file which contains the xyz coordinates of the celestial objects at a givne
 * time step.
 */
inline void saveToVTK(const std::vector<std::array<double, 3>>& data, const std::string& filepath) {
    std::ofstream vtkFile(filepath, std::ios::binary);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: Could not open file " << filepath << std::endl;
        return;
    }

    //legacy vtk header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Solar System Simulation\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET POLYDATA\n";

    // Write points
    vtkFile << "POINTS " << data.size() << " float\n";
    for (int i = 0; i < data.size(); ++i) {
        vtkFile << data[i][0] << " " << data[i][1]  << " " << data[i][2]  << "\n";
    }

    vtkFile.close();
}

//util function to (re)create a folder and clearing the previous contents
inline void clearAndCreateDirectory(const std::string& path) {
    try {
        // If directory exists, remove it and its contents
        if (std::filesystem::exists(path)) {
            std::filesystem::remove_all(path);
            std::cout << "Existing directory deleted: " << path << std::endl;
        }

        // Create the directory again
        std::filesystem::create_directories(path);
        std::cout << "Directory recreated: " << path << std::endl;

    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error with filesystem operation: " << e.what() << std::endl;
    }
}

#endif // UTILS_HPP