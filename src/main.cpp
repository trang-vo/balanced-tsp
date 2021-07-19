#include "tsplib_reader.h" // TsplibReader class
#include "utils.h"
#include "fairoptimization.h"

#include <fstream>
#include <iostream>
#include <map>
#include <time.h>
#include <cmath>
#include <string>
#include <vector>
#include <experimental/filesystem>
#include <sstream>
#include <stdlib.h>
#include <time.h>

#include <ilcplex/ilocplex.h>
#include <lemon/list_graph.h>
#include <lemon/gomory_hu.h>
#include <lemon/lgf_reader.h>

namespace fs = std::experimental::filesystem;
using namespace fairtsp;

void runOneFile(std::string fileName) {
    std::ifstream file("data/tsplib/brazil58.tsp");
    std::cout << "Start solve problem" << std::endl;
    if (file.is_open())
    {
        TsplibReader reader(file, 100);
        if (reader.checkInput)
        {
            std::cout << "start solve model" << std::endl;
            FairOpt model(reader.name, reader.graph, reader.edgeMap);
            model.solve(0);
        };
    };
};

int main()
{
    std::string filePaths = "listProblems.txt";
    writePathsInFolder2File("data/tsplib", "tsp", filePaths);
    std::vector<std::string> listProb;
    std::ifstream inputFilePaths("listProblems.txt");
    getListProblem(inputFilePaths, listProb);
    // std::ofstream outputResults("resultsMaxminDirectedproblem.txt");

    for (int i = 0; i < listProb.size(); ++i) {
        std::ifstream file(listProb[i]);
        if (file.is_open()) {
            TsplibReader reader(file, 100);
            if (reader.checkInput) {
                std::cout << "solving the problem " << listProb[i] << std::endl;
                for (uint version = 0u; version < 3; ++version) {
                    FairOpt model(reader.name, reader.graph, reader.edgeMap);
                    model.solve(version);
                };
            };
        };
    };

    return 0;
};