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
#include "solar_system.hpp"

int main(void) {
    SolarSystem solarSystem("Data\\", 0, 20000);
    //solarSystem.celestials[3].GM *= 2;

    // run the simulation
    for (int t=0; t<1e6; t++) {
        if (t % 100000 == 0) {
            std::cout << "t = " << t << std::endl;
        }
        solarSystem.performTimestep();
    }

    std::cout << "Number of inner centaurs: " << solarSystem.nInner << std::endl;
    std::cout << "Number of outer centaurs: " << solarSystem.nOuter << std::endl;

    return 0;
}