#include "utils.h"

#include <map>
#include <fstream>
#include <iostream>
#include <experimental/filesystem>

#include <lemon/list_graph.h>
#include <lemon/lgf_writer.h>
#include <lemon/connectivity.h>

namespace fs = std::experimental::filesystem;

namespace fairtsp 
{
    bool getFilesInFolder(std::string folder, std::string extenstion, std::vector<std::string>& listFiles)
    {
        try
        {
            for (const auto& entry : fs::directory_iterator(folder))
            {
                std::string fileName = entry.path();
                std::vector<std::string> partFile = splitSentence(fileName.substr(folder.size() + 1), '.');
                if (partFile.back() == extenstion)
                    listFiles.push_back(fileName);
            };  
        }
        catch(const std::exception& e)
        {
            std::cout << "Error in getFilesInFolder" << std::endl;
            std::cerr << e.what() << '\n';
            return false;
        };

        return true;
    };

    bool writeVectorToFile(std::vector<std::string> vectors, std::string fileName)
    {
        std::ofstream file(fileName);
        try
        {
            for (std::string element : vectors)
            {
                file << element << std::endl;
            };
        }
        catch (const std::exception &e)
        {
            std::cout << "Error in writeVectorToFile" << std::endl;
            std::cerr << e.what() << '\n';
            return false;
        };

        return true;        
    };

    void writePathsInFolder2File(std::string folderName, std::string extension, std::string fileName) {
        std::vector<std::string> listFile;
        getFilesInFolder(folderName, extension, listFile);
        writeVectorToFile(listFile, fileName);
    };

    std::map<std::string, int> getNumberCuts(IloCplex model)
    {
        std::map<std::string, int> numberCuts;
        numberCuts["CutCover"] = model.getNcuts(IloCplex::CutType::CutCover);
        numberCuts["CutGubCover"] = model.getNcuts(IloCplex::CutType::CutGubCover);
        numberCuts["CutFrac"] = model.getNcuts(IloCplex::CutType::CutFrac);
        numberCuts["CutMir"] = model.getNcuts(IloCplex::CutType::CutMir);
        numberCuts["CutFlowPath"] = model.getNcuts(IloCplex::CutType::CutFlowPath);
        numberCuts["CutDisj"] = model.getNcuts(IloCplex::CutType::CutDisj);
        numberCuts["CutImplBd"] = model.getNcuts(IloCplex::CutType::CutImplBd);
        numberCuts["CutZeroHalf"] = model.getNcuts(IloCplex::CutType::CutZeroHalf);
        numberCuts["CutLocalCover"] = model.getNcuts(IloCplex::CutType::CutLocalCover);
        numberCuts["CutTighten"] = model.getNcuts(IloCplex::CutType::CutTighten);
        numberCuts["CutObjDisj"] = model.getNcuts(IloCplex::CutType::CutObjDisj);
        numberCuts["CutUser"] = model.getNcuts(IloCplex::CutType::CutUser);
        numberCuts["CutTable"] = model.getNcuts(IloCplex::CutType::CutTable);
        numberCuts["CutSolnPool"] = model.getNcuts(IloCplex::CutType::CutSolnPool);

        return numberCuts;
    };

    bool getListProblem(std::ifstream& file, std::vector<std::string> &listProblem)
    {
        std::string line;
        try
        {
            while (file)
            {
                getline(file, line);
                trim(line);
                if (line.size() > 0)
                    listProblem.push_back(line);
            };
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
            return false;
        };

        return true;
    };

    std::vector<std::string> splitSentence(const std::string& str, char delimiter)
    {
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(str);
        while (std::getline(tokenStream, token, delimiter))
        {
            trim(token);
            if (token.size() > 0)
                tokens.push_back(token);
        };
        return tokens;
    };

    void generateRandomGraph(std::string fileName, int numNodes)
    {
        std::ofstream output(fileName);
        // write properties
        output << "DIMENSION : " << numNodes << std::endl;
        output << "WEIGHT_TYPE : 1" << std::endl;

        // generate weights for instance
        for (int i = 0; i < numNodes; ++i)
        {
            for (int j = i + 1; j < numNodes; ++j)
            {
                float weight;
                if (rand() % 6 == 0)
                    weight = (float) rand() / RAND_MAX;
                else
                {
                    weight = 0;
                };
                output << weight << " ";
            };
            output << std::endl;
        };

        output << "EOF";
        output.close();
    };

    std::vector<std::string> readWeightLine(std::ifstream& fin)
    {
        std::string line;
        std::vector<std::string> weights;
        while (fin)
        {
            getline(fin, line);
            trim(line);
            if (line != "EOF" & line != "DISPLAY_DATA_SECTION")
            {
                std::vector<std::string> weightLine = splitSentence(line, ' ');
                weights.insert(weights.end(), weightLine.begin(), weightLine.end());
            }
            else
            {
                break;
            };        
        };

        return weights;
    };

    std::vector<int> getMinCutSet(const lemon::ListGraph &graph, std::vector<int> edgesUsed)
    {
        std::vector<int> setWithout1;
        for (int i = 0; i < edgesUsed.size(); ++i)
        {
            lemon::ListGraph::Edge edge = graph.edgeFromId(edgesUsed[i]);
            int endpoint1 = graph.id(graph.u(edge));
            int endpoint2 = graph.id(graph.v(edge));

            if (endpoint1 == 0)
                setWithout1.push_back(endpoint2 + 1);
            else if (endpoint2 == 0)
            {
                setWithout1.push_back(endpoint1 + 1);
            };
        };

        return setWithout1;
    };

    void sortEdges(std::vector<int> indexEdgesUsed, std::vector<int> valueEdgesUsed, std::vector<int>& sortedEdges)
    {
        std::map<int, std::vector<int>> linker;
        for (int i = 0; i < indexEdgesUsed.size(); ++i)
        {
            linker[valueEdgesUsed[i]].push_back(indexEdgesUsed[i]);
        };

        for (std::map<int, std::vector<int>>::iterator it = linker.begin(); it != linker.end(); ++it)
        {
            for (int j = 0; j < it->second.size(); ++j)
            {
                sortedEdges.push_back(it->second[j]);
            };
        };
    };

    bool checkBiConnected(const lemon::ListGraph &graph, const lemon::ListGraph::EdgeMap<double>& edgeVarValue) {
        typedef lemon::ListGraph Graph;
        Graph supportGraph;
        lemon::GraphCopy<Graph, Graph> cg(graph, supportGraph);
        Graph::NodeMap<Graph::Node> nr(graph);
        cg.nodeRef(nr);
        cg.run();

        for (Graph::EdgeIt e(graph); e != lemon::INVALID; ++e) {
            if (std::abs(edgeVarValue[(Graph::Edge)e]) > TOLERANCE) {
                Graph::Node endpoint1 = graph.u((Graph::Edge)e);
                Graph::Node endpoint2 = graph.v((Graph::Edge)e);
                supportGraph.addEdge(endpoint1, endpoint2);
            };
        };

        return lemon::biEdgeConnected(supportGraph);
    };

    void generateGGI(std::string folder, int numNodes) {
        std::string fileName = folder + "weight_ggi_" + std::to_string(numNodes);
        std::ofstream output(fileName);
        for (int i = 1; i <= numNodes; ++i) {
            output << (float) (2 * (numNodes - i) + 1) / (float) (numNodes * numNodes);
            if (i < numNodes) {
                output << ",";
            };
        };
        output.close();
    };
};