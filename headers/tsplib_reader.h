#ifndef TSPLIBREADER_H
#define TSPLIBREADER_H

#include "utils.h"

#include <cmath>
#include <string>

#include <lemon/list_graph.h>
#include <lemon/adaptors.h>


namespace fairtsp {
    class TsplibReader
    {
    private:
        const char delimiter = ':';
        int limitDimension;

        bool readProblem(std::ifstream &fin);
        bool checkKeywords(std::string part, std::string value);

        lemon::ListGraph::EdgeMap<std::vector<float>> getNodeMap(lemon::ListGraph g);

        void processNodeCoord(std::ifstream& fin);
        void processEdgeWeight(std::ifstream& fin);

        // functions process node coord
        int metric(std::vector<float> x, std::vector<float> y);
        std::vector<float> convert2Geographical(std::vector<float> x);
        int distanceGeographical(std::vector<float> x, std::vector<float> y);
        int distancePseudoEuclidean(std::vector<float> x, std::vector<float> y);

        // functions process edge weight
        void readFullMatrix(std::vector<std::string> weights);
        void readUpperRow(std::vector<std::string> weights);
        void readLowerRow(std::vector<std::string> weights);
        void readUpperDiagRow(std::vector<std::string> weights);
        void readLowerDiagRow(std::vector<std::string> weights);
        void readUpperCol(std::vector<std::string> weights);
        void readLowerCol(std::vector<std::string> weights);
        void readUpperDiagCol(std::vector<std::string> weights);
        void readLowerDiagCol(std::vector<std::string> weights);

    public:    
        std::string name;
        std::string type;
        std::string comment;
        int dimension;
        int capacity;
        std::string edge_weight_type;
        std::string edge_weight_format;
        std::string edge_data_format;
        std::string node_coord_type;
        std::string display_data_type;
        std::vector<int> solution;  
        bool checkInput;
        bool checkOutput;
          
        lemon::ListGraph graph;
        std::vector<lemon::ListGraph::Edge> edges;
        lemon::ListGraph::NodeMap<std::vector<float>> nodeMap;
        lemon::ListGraph::EdgeMap<int> edgeMap;

        lemon::ListGraph sparseGraph;
        lemon::ListGraph::EdgeMap<int> edgeMapSparse;

        explicit TsplibReader(std::ifstream& inputFile, int limitDimension);
        explicit TsplibReader(std::ifstream& inputFile, std::ifstream& outputFile, int limitDimension);

        void getSparseGraph(int percent);
    };
};

#endif