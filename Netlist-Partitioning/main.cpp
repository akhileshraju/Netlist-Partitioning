//
//  main.cpp
//  Netlist-Partitioning
//
//  Created by Akhilesh Raju on 12/2/16.
//  Copyright © 2016 Akhilesh Raju. All rights reserved.
//

/***********************************************************************
 * Simulated Annealing Algorithm
 * Execute using this ---> "./main <benchmark_file_name> <output_file_name>"
 *********************************************************************/

#include <iostream>
#include "annealing.hpp"

/***************************
 *      Main function
 **************************/

int main(int argc, char *argv[]) {
    std::string inputFilePath;
    std::string outputFilePath;
    if(argc == 2) {
        inputFilePath = argv[1];
        outputFilePath = argv[2];
    }
    else {
        inputFilePath = "/Users/akhileshraju/Coding/Xcode/partitioningTool/Netlist-Partitioning/Netlist-Partitioning/Benchmarks/B1";
        outputFilePath = "/Users/akhileshraju/Coding/Xcode/partitioningTool/Netlist-Partitioning/Netlist-Partitioning/Benchmarks/R1";
    }
    simulatedAnnealing sa(inputFilePath);
    sa.performAnnealing(outputFilePath);
    std::cout << "\nOutput file \"" << outputFilePath << "\" generated in working folder\n\n";
    return 0;
}
