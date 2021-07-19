#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include <functional>
#include <map>
#include <fstream>

#include <lemon/list_graph.h>
#include <lemon/lgf_writer.h>

// Magic tricks to have CPLEX behave well:
#ifndef IL_STD
#define IL_STD
#endif
#include <cstring>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
// End magic tricks

namespace fairtsp
{
    const double TOLERANCE = .00001;
    bool getFilesInFolder(std::string folder, std::string extenstion, std::vector<std::string> &listFile);
    bool writeVectorToFile(std::vector<std::string> vectors, std::string fileName);
    void writePathsInFolder2File(std::string folderName, std::string extension, std::string fileName);

    std::map<std::string, int> getNumberCuts(IloCplex model);
    bool getListProblem(std::ifstream& file, std::vector<std::string> &listProblem);

    std::vector<std::string> splitSentence(const std::string& str, char delimiter);
    static inline void ltrim(std::string &s) {s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {return !std::isspace(ch);}));};
    static inline void rtrim(std::string &s) {s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {return !std::isspace(ch);}).base(), s.end());}
    static inline void trim(std::string &s) {
        ltrim(s);
        rtrim(s);
    };

    template <class K, class V>
    std::vector<K> getKeysMap(std::map<K, V> dict)
    {
        std::vector<K> keys;
        for (const std::tuple<K,V>& element : dict)
            keys.push_back(std::get<0>(element));

        return keys;
    };

    std::vector<std::string> readWeightLine(std::ifstream &fin);

    void generateRandomGraph(std::string fileName, int numNodes);
    std::vector<int> getMinCutSet(const lemon::ListGraph &graph, std::vector<int> edgesUsed);
    void getCombinations(std::vector<std::vector<int>> &result ,int array[], int combi[], int start, int end, int Nnum, int index);

    void sortEdges(std::vector<int> indexEdgesUsed, std::vector<int> valueEdgesUsed, std::vector<int>& sortedEdges);
    bool checkBiConnected(const lemon::ListGraph &graph, const lemon::ListGraph::EdgeMap<double>& edgeVarValue);

    void generateGGI(std::string folder, int numNodes);
    bool writeVectorToFile(std::vector<std::string> vectors, std::string fileName);
    void writePathsInFolder2File(std::string folderName, std::string extension, std::string fileName);

    template <class WeightType>
    class RandomGraph
    {
        void generateGraph()
        {
            for (int i = 0; i < this->dimension; ++i)
            {
                this->graph.addNode();
            };

            // generate completed graph
            for (lemon::ListGraph::NodeIt it(this->graph); it != lemon::INVALID; ++it)
            {
                lemon::ListGraph::NodeIt jt = it;
                ++jt;
                while (jt != lemon::INVALID)
                {
                    this->graph.addEdge(it, jt);
                    ++jt;
                };
            };
        };

    public:
        int dimension;
        int weightType;
        lemon::ListGraph graph;
        lemon::ListGraph::EdgeMap<WeightType> edgeMap;
        RandomGraph() : edgeMap(graph){};
        void loadData(std::string fileName)
        {
            std::ifstream input(fileName);
            std::string line;

            getline(input, line);
            sscanf(line.c_str(), "DIMENSION : %d", &this->dimension);
            getline(input, line);
            sscanf(line.c_str(), "WEIGHT_TYPE : %d", &this->weightType);

            generateGraph();

            std::vector<std::string> weights = readWeightLine(input);

            int id = 0;
            for (int i = 0; i < this->dimension; ++i)
            {
                for (int j = i + 1; j < this->dimension; ++j)
                {
                    lemon::ListGraph::Edge edge = lemon::findEdge(this->graph, this->graph.nodeFromId(i), this->graph.nodeFromId(j));
                    if (this->weightType == 0)
                        this->edgeMap[edge] = std::stoi(weights[id]);
                    else if (this->weightType == 1)
                    {
                        this->edgeMap[edge] = std::stod(weights[id]);
                    };

                    if (j == i + 1 & i % 2 == 0)
                    {
                        this->edgeMap[edge] = 999999;
                    };
                    ++id;

                    if (i % 2 == 1 & j % 2 == 1)
                    {
                        this->edgeMap[edge] = 999999;
                    };

                    if (i == 0 & j % 2 == 0)
                    {
                        this->edgeMap[edge] = 999999;
                    }
                };
            };
        };
    };

    template <class WT>
    void generateGraphFile(const char* fileName, int numNodes, int maxVal = 100, int percentEdges = 20, int bigValIndex = 1)
    {
        lemon::ListGraph graph;
        lemon::ListGraph::EdgeMap<WT> edgeMap(graph);

        for (int i = 0; i < numNodes; ++i)
        {
            graph.addNode();
        };

        // generate completed graph
        std::vector<lemon::ListGraph::Edge> edges;
        for (lemon::ListGraph::NodeIt it(graph); it != lemon::INVALID; ++it)
        {
            lemon::ListGraph::NodeIt jt = it;        
            ++jt;
            while (jt != lemon::INVALID)
            {
                lemon::ListGraph::Edge edge = graph.addEdge(it, jt);
                edges.push_back(edge);
                ++jt;
            };
        };

        // generate big edge costs
        for (int i = 0; i < numNodes - 2; ++i)
        {
            lemon::ListGraph::Edge edge = lemon::findEdge(graph, graph.nodeFromId(i), graph.nodeFromId(i + 2));
            edgeMap[edge] = 999999;
        };
        
        // generate random edge costs
        std::random_shuffle(edges.begin(), edges.end());
        int numRealEdges = (int) (lemon::countEdges(graph) * percentEdges / 100);
        while (numRealEdges >= 0)
        {
            if (std::is_same<int, WT>::value)
            {
                edgeMap[edges[numRealEdges]] = (int) rand() % maxVal;
            }
            else if (std::is_same<WT, float>::value)
            {
                edgeMap[edges[numRealEdges]] = (float) rand() / RAND_MAX;
            }
            else
            {
                edgeMap[edges[numRealEdges]] = 1;
            }
            numRealEdges--;
        };

        lemon::GraphWriter<lemon::ListGraph> writerGraph(graph, fileName);
        lemon::IdMap<lemon::ListGraph, lemon::ListGraph::Node> nodeLabelMap(graph);
        writerGraph.nodeMap("label", nodeLabelMap);
        writerGraph.edgeMap("edge_costs", edgeMap);
        writerGraph.run();
    };

}; // namespace tsp

#endif