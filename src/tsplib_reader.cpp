#include "tsplib_reader.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm> 

namespace fairtsp {
    using lemon::ListGraph;
    using std::vector;

    inline int nint(float x) { return int(x + 0.5); };

    TsplibReader::TsplibReader(std::ifstream &inputFile, int limitDimension) : name(), type(), comment(), dimension(), capacity(),
    edge_weight_type(), edge_weight_format(), edge_data_format(), node_coord_type(), display_data_type(), limitDimension(limitDimension),
    graph(), edges(), nodeMap(this->graph), edgeMap(this->graph), sparseGraph(), edgeMapSparse(this->sparseGraph)
    {
        if (!readProblem(inputFile)) 
        {
            std::cout << "Input file is invalid" << std::endl;
            this->checkInput = false;
        }
        else
        {
            this->checkInput = true;
            std::cout << "Imported problem" << std::endl;
            std::cout << "Name: " << this->name << std::endl;
            std::cout << "Dimension: " << this->dimension << std::endl;
        };  
        inputFile.close();      
    };

    TsplibReader::TsplibReader(std::ifstream &inputFile, std::ifstream &outputFile, int limitDimension) : name(), type(), comment(), dimension(), capacity(),
    edge_weight_type(), edge_weight_format(), edge_data_format(), node_coord_type(), display_data_type(), limitDimension(limitDimension),
    graph(), edges(), nodeMap(this->graph), edgeMap(this->graph), sparseGraph(), edgeMapSparse(this->sparseGraph)
    {
        if (!readProblem(inputFile)) 
        {
            std::cout << "Input file is invalid" << std::endl;
            this->checkInput = false;
        }
        else
        {
            this->checkInput = true;
            std::cout << "Imported problem" << std::endl;
            std::cout << "Name: " << this->name << std::endl;
            std::cout << "Dimension: " << this->dimension << std::endl;

            std::string line;
            try
            {
                while (outputFile)
                {
                    int value;
                    getline(outputFile, line);
                    if (sscanf(line.c_str(), "%d", &value) == 1)
                    {
                        if (value != -1)
                            this->solution.push_back(value);
                        else
                        {
                            this->solution.push_back(1);
                        };                
                    };
                };
                this->checkOutput = true;
            }
            catch (const std::exception& e)
            {
                this->checkOutput = false;
                std::cerr << "error " << std::endl;
            };
        };        

        inputFile.close();
        outputFile.close();
    };

    bool TsplibReader::readProblem(std::ifstream &fin) 
    {
        std::string line;
        while (fin) 
        {
            getline(fin, line);
            if (line == "EOF") {
                break;
            };

            line.erase(remove(line.begin(), line.end(), '\r'), line.end());

            if (line.find(delimiter) != line.npos) {
                std::vector<std::string> keywords = splitSentence(line, this->delimiter);
                std::string part = keywords[0];
                std::string value = keywords[1];

                if (!checkKeywords(part, value)) {
                    return false;
                };
            };

            if (this->dimension != 0) {
                if (this->dimension > this->limitDimension)
                {
                    std::cout << "Dimension exceed " << this->dimension << std::endl;
                    return false;
                };
                
                trim(line);
                if (line == "NODE_COORD_SECTION")
                {
                    processNodeCoord(fin);
                    break;
                };

                if (line == "EDGE_WEIGHT_SECTION")
                {
                    std::cout << "Process edges weight" << std::endl;
                    processEdgeWeight(fin);
                    break;
                };
            };
        };

        return true;
    };

    bool TsplibReader::checkKeywords(std::string part, std::string value) {
        if (part == "NAME") {
            this->name = value;
        } 
        else if (part == "TYPE") {
            if (value == "TSP" || value == "ATSP") {
                this->type = value;
            } else {
                std::cout << part << " not supported" << std::endl;
                return false;
            };
        }
        else if (part == "COMMENT") {
            this->comment = value;
        }
        else if (part == "DIMENSION") {
            this->dimension = stoi(value);
        }
        else if (part == "CAPACITY") {
            this->capacity = stoi(value);
        }
        else if (part == "EDGE_WEIGHT_TYPE") {
            if (value == "EXPLICIT" ||
                value == "EUC_2D" || value == "EUC_3D" ||
                value == "MAX_2D" || value == "MAX_3D" ||
                value == "MAN_2D" || value == "MAN_3D" ||
                value == "CEIL_2D" || value == "GEO" ||
                value == "XRAY1" || value == "XRAY2" ||
                value == "ATT" || value == "SPECIAL") {
                this->edge_weight_type = value;
                } else {
                std::cout << part << " not supported" << std::endl;
                return false;
                };
        }
        else if (part == "EDGE_WEIGHT_FORMAT") {
            if (value == "FUNCTION" || value == "FULL_MATRIX" ||
                value == "UPPER_ROW" || value == "LOWER_ROW" ||
                value == "UPPER_DIAG_ROW" || value == "LOWER_DIAG_ROW" ||
                value == "UPPER_COL" || value == "LOWER_COL" ||
                value == "UPPER_DIAG_COL" || value == "LOWER_DIAG_COL") {
                this->edge_weight_format = value;
                } else {
                std::cout << part << " not supported" << std::endl;
                return false;
                };
        }
        else if (part == "EDGE_DATA_FORMAT") {
            if (value == "EDGE_LIST" || value == "ADJ_LIST") {
                this->edge_data_format = value;
            } else {
                std::cout << part << " not supported" << std::endl;
                return false;
            };
        }
        else if (part == "NODE_COORD_TYPE") {
            if (value == "TWOD_COORDS" || value == "THREED_COORDS" || value == "NO_COORDS") {
                this->node_coord_type = value;
            } else {
                std::cout << part << " not supported" << std::endl;
                return false;
            };
        }
        else if (part == "DISPLAY_DATA_TYPE") {
            this->display_data_type = value;
        }
        else {
            std::cout << part << " not supported" << std::endl;
            return false;
        };

        return true;
    };

    void TsplibReader::processNodeCoord(std::ifstream& fin){
        int vertex;
        vector<float> vertex_coord;
        ListGraph::Node node;
        ListGraph::Edge edge;

        float x, y, z;
        std::string line;
        while (fin) {
            getline(fin, line);

            if (fin) {
                if (sscanf(line.c_str(), "%d %f %f %f", &vertex, &x, &y, &z) == 4) {
                    vertex_coord = {x, y, z};
                }
                else if (sscanf(line.c_str(), "%d %f %f %f", &vertex, &x, &y, &z) == 3) {
                    vertex_coord = {x, y};
                }
                else
                    break;

                node = this->graph.addNode();
                this->nodeMap[node] = vertex_coord;
            };
        };

        lemon::graphCopy(this->graph, this->sparseGraph).run();

        // generate completed graph
        for (ListGraph::NodeIt it(this->graph); it != lemon::INVALID; ++it)
        {
            ListGraph::NodeIt jt = it;
            ++jt;
            while (jt != lemon::INVALID)
            {
                edge = this->graph.addEdge(it, jt);
                this->edges.push_back(edge);
                edgeMap[edge] = this->metric(nodeMap[it], nodeMap[jt]);
                ++jt;
            };
        };
    };


    void TsplibReader::processEdgeWeight(std::ifstream& fin)
    {
        for (int i = 0; i < this->dimension; ++i)
        {
            this->graph.addNode();
        };

        lemon::graphCopy(this->graph, this->sparseGraph).run();

        // generate completed graph
        for (ListGraph::NodeIt it(this->graph); it != lemon::INVALID; ++it)
        {
            ListGraph::NodeIt jt = it;
            ++jt;
            while (jt != lemon::INVALID)
            {
                lemon::ListGraph::Edge edge = this->graph.addEdge(it, jt);
                edges.push_back(edge);
                ++jt;
            };
        };

        std::vector<std::string> weights = readWeightLine(fin);

        if (this->edge_weight_format == "FULL_MATRIX")
        {
            readFullMatrix(weights);
        }
        else if (this->edge_weight_format == "UPPER_ROW")
        {
            readUpperRow(weights);
        }
        else if (this->edge_weight_format == "LOWER_ROW")
        {
            readLowerRow(weights);
        }
        else if (this->edge_weight_format == "UPPER_DIAG_ROW")
        {
            readUpperDiagRow(weights);
        }
        else if (this->edge_weight_format == "LOWER_DIAG_ROW")
        {
            readLowerDiagRow(weights);
        };
    };

    int TsplibReader::metric(vector<float> x, vector<float> y) {
        if (this->edge_weight_type == "EUC_2D") {
            return nint(std::sqrt(std::pow(x[0] - y[0], 2) + std::pow(x[1] - y[1], 2)));
        }
        else if (this->edge_weight_type == "EUC_3D") {
            return nint(std::sqrt(std::pow(x[0] - y[0], 2) + std::pow(x[1] - y[1], 2) + std::pow(x[2] - y[2], 2)));
        }
        else if (this->edge_weight_type == "MAX_2D") {
            return nint(std::max(std::abs(x[0] - y[0]), std::abs(x[1] - y[1])));
        }
        else if (this->edge_weight_type == "MAX_3D") {
            return nint(std::max(std::max(std::abs(x[0] - y[0]), std::abs(x[1] - y[1])), std::abs(x[2] - y[2])));
        }
        else if (this->edge_weight_type == "MAN_2D") {
            return nint(std::abs(x[0] - y[0]) + std::abs(x[1] - y[1]));
        }
        else if (this->edge_weight_type == "MAN_3D") {
            return nint(std::abs(x[0] - y[0]) + std::abs(x[1] - y[1]) + std::abs(x[2] - y[2]));
        }
        else if (this->edge_weight_type == "GEO") {
            return distanceGeographical(x, y); 
        }
        else if (this->edge_weight_type == "ATT") {
            return distancePseudoEuclidean(x, y); 
        }
        else if (this->edge_weight_type == "CEIL_2D") {
            return std::ceil(std::sqrt(std::pow(x[0] - y[0], 2) + std::pow(x[1] - y[1], 2)));
        };

        return 0;
    };

    vector<float> TsplibReader::convert2Geographical(vector<float> x) {
        vector<float> geo_coord;
        const float PI = 3.141592;
        
        int deg = int(x[0]);
        float min = x[0] - deg;
        float latitude = PI * (deg + 5.0 * min / 3.0) / 180.0;
        geo_coord.push_back(latitude);

        deg = int(x[1]);
        min = x[1] - deg;
        float longitude = PI * (deg + 5.0 * min / 3.0) / 180.0;
        geo_coord.push_back(longitude);

        return geo_coord;
    };

    int TsplibReader::distanceGeographical(vector<float> x, vector<float> y) {
        const float RRR = 6378.388;
        vector<float> xg = convert2Geographical(x);
        vector<float> yg = convert2Geographical(y);

        float q1 = std::cos(xg[1] - yg[1]);
        float q2 = std::cos(xg[0] - yg[0]);
        float q3 = std::cos(xg[0] + yg[0]);
        int dij;
        dij = int(RRR * std::acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);

        return dij;
    };

    int TsplibReader::distancePseudoEuclidean(vector<float> x, vector<float> y) {
        float xd = x[0] - y[0];
        float yd = x[1] - y[1];
        float rij = std::sqrt((xd * xd + yd * yd) / 10.0);
        int tij = nint(rij);
        if (tij < rij) {
            return tij + 1;
        };
        return tij;
    };

    void TsplibReader::readFullMatrix(std::vector<std::string> weights)
    {
        int id = 0;
        for (int i = 0; i < this->dimension; ++i)
        {
            for (int j = 0; j < this->dimension; ++j)
            {
                int weight = std::stoi(weights[id]);
                if (i < j & weight != 0)
                {
                    lemon::ListGraph::Edge edge = lemon::findEdge(this->graph, this->graph.nodeFromId(i), this->graph.nodeFromId(j));
                    this->edgeMap[edge] = std::stoi(weights[id]);
                };
                ++id;
            };
        };
    }

    void TsplibReader::readUpperRow(std::vector<std::string> weights)
    {
        int id = 0;
        for (int i = 0; i < this->dimension; ++i)
        {
            for (int j = i + 1; j < this->dimension; ++j)
            {
                lemon::ListGraph::Edge edge = lemon::findEdge(this->graph, this->graph.nodeFromId(i), this->graph.nodeFromId(j));
                this->edgeMap[edge] = std::stoi(weights[id]);
                ++id;
            };
        };
    };

    void TsplibReader::readLowerRow(std::vector<std::string> weights)
    {
        int id = 0;
        for (int i = 1; i < this->dimension; ++i)
        {
            for (int j = 0; j < i; ++j)
            {
                lemon::ListGraph::Edge edge = lemon::findEdge(this->graph, this->graph.nodeFromId(i), this->graph.nodeFromId(j));
                this->edgeMap[edge] = std::stoi(weights[id]);
                ++id;
            };
        };
    };
    
    void TsplibReader::readUpperDiagRow(std::vector<std::string> weights)
    {
        int id = 0;
        for (int i = 0; i < this->dimension; ++i)
        {
            for (int j = i; j < this->dimension; ++j)
            {
                if (i != j)
                {
                    lemon::ListGraph::Edge edge = lemon::findEdge(this->graph, this->graph.nodeFromId(i), this->graph.nodeFromId(j));
                    this->edgeMap[edge] = std::stoi(weights[id]);
                }
                ++id;
            };
        };
    };

    void TsplibReader::readLowerDiagRow(std::vector<std::string> weights)
    {
        int id = 0;
        for (int i = 0; i < this->dimension; ++i)
        {
            for (int j = 0; j < i + 1; ++j)
            {
                if (i != j)
                {
                    lemon::ListGraph::Edge edge = lemon::findEdge(this->graph, this->graph.nodeFromId(i), this->graph.nodeFromId(j));
                    this->edgeMap[edge] = std::stoi(weights[id]);
                }
                ++id;
            };
        };
    };

    void TsplibReader::getSparseGraph(int percent)
    {
        std::vector<lemon::ListGraph::Edge> solutionEdges;
        for (int i = 0; i < this->solution.size() - 1; ++i)
        {
            lemon::ListGraph::Edge edge = lemon::findEdge(this->graph, 
                                                          this->graph.nodeFromId(this->solution[i] - 1), 
                                                          this->graph.nodeFromId(this->solution[i+1] - 1));
            solutionEdges.push_back(edge);
        };

        int numEdges = (int) (lemon::countEdges(this->graph) * percent / 100); 
        std::cout << "number edges of sparse graph is " << numEdges << std::endl;
        std::random_shuffle(this->edges.begin(), this->edges.end());

        for (int i = 0; i < solutionEdges.size(); ++i)
        {
            lemon::ListGraph::Node n1 = this->graph.u(solutionEdges[i]);
            lemon::ListGraph::Node n2 = this->graph.v(solutionEdges[i]);
            lemon::ListGraph::Edge edge = this->sparseGraph.addEdge(n1, n2);
            this->edgeMapSparse[edge] = this->edgeMap[solutionEdges[i]];
        };

        while (numEdges > 0)
        {
            lemon::ListGraph::Edge e = this->edges[numEdges];
            if (std::find(solutionEdges.begin(), solutionEdges.end(), e) == solutionEdges.end())
            {
                lemon::ListGraph::Node n1 = this->graph.u(e);
                lemon::ListGraph::Node n2 = this->graph.v(e);
                lemon::ListGraph::Edge edge = this->sparseGraph.addEdge(n1, n2);
                this->edgeMapSparse[edge] = this->edgeMap[e];
            };
            numEdges--;
        };
    };

}; // namespace tsp
