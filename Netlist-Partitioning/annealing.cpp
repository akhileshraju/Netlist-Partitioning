//
//  annealing.cpp
//  Netlist-Partitioning
//
//  Created by Akhilesh Raju on 12/2/16.
//  Copyright © 2016 Akhilesh Raju. All rights reserved.
//

#include "annealing.hpp"

/*******************************************************************
 *                      Private member functions
 ******************************************************************/

void simulatedAnnealing::setParameters() {
    if (noOfCells <= 50) {
        alpha = 0.98;
        T = 5;
        Tf = 0.01;
        steps = 1000;
    }
    else if (noOfCells == 500) {
        alpha = 0.979;
        T = 80;
        Tf = 0.001;
        steps = 5000;
    }
    else if (noOfCells == 4500) {
        alpha = 0.975;
        T = 10;
        Tf = 0.000001;
        steps = 25000;
    }
    else if (noOfCells == 10000) {
        alpha = 0.978;
        T = 10;
        Tf = 0.0000000001;
        steps = 40000;
    }
    else if (noOfCells == 25000) {
        alpha = 0.98;
        T = 10;
        Tf = 0.00000000000001;
        steps = 65000;
    }
    else {
        alpha = 0.98;
        T = 15;
        Tf = 0.0000000000000001;
        steps = 100000;
    }
}

void simulatedAnnealing::addToAdjList(unsigned int data[2]) {
    if (adjList[data[0]] == NULL) {
        adjList[data[0]] = new node;
        adjList[data[0]]->value = data[1];
        adjList[data[0]]->weight = 1;
        adjList[data[0]]->next = NULL;
    }
    else {
        node *newNode = adjList[data[0]];
        while (newNode->next != NULL)
            newNode = newNode->next;
        newNode->next = new node;
        newNode->next->value = data[1];
        newNode->next->weight = 1;
        newNode->next->next = NULL;
    }
    if (adjList[data[1]] == NULL) {
        adjList[data[1]] = new node;
        adjList[data[1]]->value = data[0];
        adjList[data[1]]->weight = 1;
        adjList[data[1]]->next = NULL;
    }
    else {
        node *newNode = adjList[data[1]];
        while (newNode->next != NULL)
            newNode = newNode->next;
        newNode->next = new node;
        newNode->next->value = data[0];
        newNode->next->weight = 1;
        newNode->next->next = NULL;
    }
}

void simulatedAnnealing::updateEValue(unsigned int nodeNumber) {
    if (nodeNumber == 0) {
        for (unsigned int iter = 0; iter < partition1.size(); iter++) {
            node *current = adjList[partition1[iter]];
            while (current->next != NULL) {
                if (find(partition2.begin(), partition2.end(), current->value)!= partition2.end())
                    eValues[partition1[iter]]++;
                current = current->next;
            }
            if (find(partition2.begin(), partition2.end(), current->value)!= partition2.end())
                eValues[partition1[iter]]++;
        }
        for (unsigned int iter = 0; iter < partition2.size(); iter++) {
            node *current = adjList[partition2[iter]];
            while (current->next != NULL) {
                if (find(partition1.begin(), partition1.end(), current->value)!= partition1.end())
                    eValues[partition2[iter]]++;
                current = current->next;
            }
            if (find(partition1.begin(), partition1.end(), current->value)!= partition1.end())
                eValues[partition2[iter]]++;
        }
    }
    else {
        eValues[nodeNumber] = 0;
        if (find(partition1.begin(), partition1.end(), nodeNumber)!= partition1.end()) {
            node *current = adjList[nodeNumber];
            while (current->next != NULL) {
                if (find(partition2.begin(), partition2.end(), current->value)!= partition2.end())
                    eValues[nodeNumber]++;
                current = current->next;
            }
            if (find(partition2.begin(), partition2.end(), current->value)!= partition2.end())
                eValues[nodeNumber]++;
        }
        else {
            node *current = adjList[nodeNumber];
            while (current->next != NULL) {
                if (find(partition1.begin(), partition1.end(), current->value)!= partition1.end())
                    eValues[nodeNumber]++;
                current = current->next;
            }
            if (find(partition1.begin(), partition1.end(), current->value)!= partition1.end())
                eValues[nodeNumber]++;
        }
    }
}

void simulatedAnnealing::updateIValue(unsigned int nodeNumber) {
    if (nodeNumber == 0) {
        for (unsigned int iter = 0; iter < partition1.size(); iter++) {
            node *current = adjList[partition1[iter]];
            while (current->next != NULL) {
                if (find(partition1.begin(), partition1.end(), current->value)!= partition1.end())
                    iValues[partition1[iter]]++;
                current = current->next;
            }
            if (find(partition1.begin(), partition1.end(), current->value)!= partition1.end())
                iValues[partition1[iter]]++;
        }
        for (unsigned int iter = 0; iter < partition2.size(); iter++) {
            node *current = adjList[partition2[iter]];
            while (current->next != NULL) {
                if (find(partition2.begin(), partition2.end(), current->value)!= partition2.end())
                    iValues[partition2[iter]]++;
                current = current->next;
            }
            if (find(partition2.begin(), partition2.end(), current->value)!= partition2.end())
                iValues[partition2[iter]]++;
        }
    }
    else {
        iValues[nodeNumber] = 0;
        if(find(partition1.begin(), partition1.end(), nodeNumber)!= partition1.end()) {
            node *current = adjList[nodeNumber];
            while (current->next != NULL) {
                if (find(partition1.begin(), partition1.end(), current->value)!= partition1.end())
                    iValues[nodeNumber]++;
                current = current->next;
            }
            if (find(partition1.begin(), partition1.end(), current->value)!= partition1.end())
                iValues[nodeNumber]++;
        }
        else {
            node *current = adjList[nodeNumber];
            while (current->next != NULL) {
                if (find(partition2.begin(), partition2.end(), current->value)!= partition2.end())
                    iValues[nodeNumber]++;
                current = current->next;
            }
            if (find(partition2.begin(), partition2.end(), current->value)!= partition2.end())
                iValues[nodeNumber]++;
        }
    }
}

void simulatedAnnealing::updateDValues(unsigned int nodeNumber) {
    if (nodeNumber == 0) {
        updateIValue(0);
        updateEValue(0);
        for (unsigned int iter = 0; iter < noOfCells; iter++)
            dValues[iter] = eValues[iter] - iValues[iter];
    }
    else {
        updateIValue(nodeNumber);
        updateEValue(nodeNumber);
        dValues[nodeNumber] = eValues[nodeNumber] - iValues[nodeNumber];
    }
}

void simulatedAnnealing::perturb() {
    index1 = rand() % (noOfCells / 2) ;
    index2 = rand() % (noOfCells / 2) ;
}

void simulatedAnnealing::calculateGain() {
    unsigned int connectivity = 0;
    perturb();
    node *current = adjList[partition1[index1]];
    while(current->next != NULL) {
        if (current->value == partition2[index2])
            connectivity++;
        current = current->next;
    }
    if(current->value == partition2[index2])
        connectivity++;
    gain = dValues[partition1[index1]] + dValues[partition2[index2]] - 2 * connectivity;
}

// Returns 'true' if move is to be accepted, 'false' if to be rejected
bool simulatedAnnealing::acceptMove() {
    double random, boltz;
    if (gain < 0)
        return 1;
    else {
        random = ((double) rand() / ((double) RAND_MAX + 1) * 1);
        boltz = pow(e, (-(gain / T)));
        if (random < boltz)
            return 1;
        else
            return 0;
    }
}

void simulatedAnnealing::computeFinalCost() {
    for(unsigned int iter = 0; iter < partition1.size(); iter++)
        finalCost += eValues[partition1[iter]];
}

void simulatedAnnealing::writeResult(std::string fileName) {
    std::ofstream fileHandle (fileName);
    if (!fileHandle.is_open()) {
        std::cout << "\n!! Could not open benchmark file !!\n";
        exit(0);
    }
    fileHandle << finalCost;
    fileHandle << std::endl;
    for(unsigned int iter = 0; iter < partition1.size(); iter++)
        fileHandle << partition1[iter] << ' ';
    fileHandle << std::endl;
    for(unsigned int iter = 0; iter < partition2.size(); iter++)
        fileHandle << partition2[iter] << ' ';
}

/*******************************************************************
 *                      Public member functions
 ******************************************************************/

simulatedAnnealing::simulatedAnnealing(std::string fileName) {
    unsigned int data[2];
    std::ifstream benchmark(fileName);
    if (!benchmark.is_open()) {
        std::cout << "\n!! Could not open benchmark file !!\n";
        exit(0);
    }
    // Reads the first two lines of the files to get the number of cells and nets
    benchmark >> noOfCells;
    benchmark >> noOfNets;
    // Resizing the vectors to the desired size based on data read from the file
    adjList.resize(noOfCells + 1, NULL);
    eValues.resize(noOfCells + 1, 0);
    iValues.resize(noOfCells + 1, 0);
    dValues.resize(noOfCells + 1, 0);
    while (!benchmark.eof()) {
        benchmark >> data[0];
        benchmark >> data[1];
        // Sending the read data to the function which inserts into the adjacency list
        addToAdjList(data);
    }
    // Creating the two partitions. At the end of this loop
    // partition1 contains cells from 1....noOfCells/2
    // partition2 contains cells from noOfCells/2 + 1....noOfCells
    for (unsigned int iter = 1; iter <= noOfCells / 2; iter++) {
        partition1.push_back(iter);
        partition2.push_back(iter + noOfCells / 2);
    }
    // Passing zero as argument tells the function to compute the D values for all the nodes
    // The updateIvalues() and updateEvalues() are called from updateDvalues()
    updateDValues(0);
    // Sets the initial, freezing temperatures, the reduction factor alpha and the number of steps
    setParameters();
}

void simulatedAnnealing::performAnnealing(std::string fileName) {
    while(T > Tf) {
        for(unsigned int iterator = 0; iterator < steps; iterator++) {
            calculateGain();
            gain *= -1;
            if(acceptMove()) {
                unsigned int temp = partition1[index1];
                partition1[index1] = partition2[index2];
                partition2[index2] = temp;
                unsigned int node1 = partition1[index1];
                unsigned int node2 = partition2[index2];
                sort(partition1.begin(), partition1.end());
                sort(partition2.begin(), partition2.end());
                updateDValues(node1);
                updateDValues(node2);
            }
        }
        T *= alpha;
    }
    computeFinalCost();
    writeResult(fileName);
}

void simulatedAnnealing::printAdjList() {
    std::cout << std::endl << std::endl;
    for (unsigned int iter = 0; iter < adjList.size(); iter++) {
        if (adjList[iter] == NULL)
            continue;
        else {
            std::cout << "\nNode " << iter << " -> ";
            node *traverse = adjList[iter];
            while (traverse->next != NULL) {
                std::cout << '(' << traverse->value << ',' << traverse->weight << ')';
                traverse = traverse->next;
            }
            std::cout << '(' << traverse->value << ',' << traverse->weight << ')';
        }
    }
}

void simulatedAnnealing::printPartitions() {
    std::cout << std::endl << std::endl;
    for (unsigned int iter = 0; iter < noOfCells / 2; iter++) {
        std::cout << partition1[iter] << "\t\t" << partition2[iter] << std::endl;
    }
}

void simulatedAnnealing::printEvalues() {
    std::cout << "\n\nE Values " << std::endl << "-----------------------\n";
    for (unsigned int iter = 1; iter < eValues.size(); iter++)
        std::cout << "Node " << iter << ": " << eValues[iter] << std::endl;
}

void simulatedAnnealing::printIvalues() {
    std::cout << "\n\nI Values " << std::endl << "-----------------------\n";
    for (unsigned int iter = 1; iter < iValues.size(); iter++)
        std::cout << "Node " << iter << ": " << iValues[iter] << std::endl;
}

void simulatedAnnealing::printDvalues() {
    std::cout << "\n\nD Values " << std::endl << "-----------------------\n";
    for (unsigned int iter = 1; iter < dValues.size(); iter++)
        std::cout << "Node " << iter << ": " << dValues[iter] << std::endl;
}
