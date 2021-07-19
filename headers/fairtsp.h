#ifndef FAIRTSP_H
#define FAIRTSP_H

#include "utils.h"

#include <vector>
#include <math.h>
#include <algorithm>
#include <map>
#include <string>

#include <lemon/list_graph.h>
#include <lemon/gomory_hu.h>

// Magic tricks to have CPLEX behave well:
#ifndef IL_STD
#define IL_STD
#endif
#include <cstring>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
// End magic tricks

namespace fairtsp {
    std::vector<std::vector<int>> getSubtours(const lemon::ListGraph &graph, std::vector<int> edgesUsed);
    void sortEdges(std::vector<int> indexEdgesUsed, std::vector<int> valueEdgesUsed, std::vector<int> &sortedEdges);

    template <class ET, class WT>
    void createILPfairtsp(IloModel model, IloIntVarArray x, IloNumVarArray v, const lemon::ListGraph& graph, const lemon::ListGraph::EdgeMap<ET>& edgeMap, std::vector<WT> weightFair)
    {
        IloEnv env = model.getEnv();
        int numNodes = lemon::countNodes(graph);
        int numEdges = lemon::countEdges(graph);
        char varName[100];

        IloExpr objExpr(env);
        IloExprArray degreeExpr(env, numNodes);
        for (int i = 0; i < numNodes; ++i)
        {
            v[i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
            std::sprintf(varName, "v.%d", i);
            v[i].setName(varName);

            objExpr += weightFair[i] * v[i];
            degreeExpr[i] = IloExpr(env);
        };
        model.add(IloMinimize(env, objExpr));

        for (int i = 0; i < numEdges; ++i)
        {            
            lemon::ListGraph::Edge edge = graph.edgeFromId(i);
            IloInt endpoint_1 = graph.id(graph.u(edge));
            IloInt endpoint_2 = graph.id(graph.v(edge));

            x[i] = IloIntVar(env, 0, 1);
            std::sprintf(varName, "x.%d.%d", endpoint_1 + 1, endpoint_2 + 1);
            x[i].setName(varName);

            degreeExpr[endpoint_1] += x[i];
            degreeExpr[endpoint_2] += x[i];
        };
        model.add(IloRangeArray(env, 2, degreeExpr, 2));
    };

    ILOUSERCUTCALLBACK2(SubtourUserCallback, IloIntVarArray, edgeVar, const lemon::ListGraph&, graph)
    {
        time_t start, end;
        start = clock();
        IloEnv masterEnv = getEnv();
        IloInt numEdges = lemon::countEdges(graph);
        IloNumArray xSol(masterEnv, numEdges);

        IloInt nNodes = getNnodes();
        if (nNodes % 10 != 0)
            return;

        if ( !isAfterCutLoop() )
            return;

        lemon::ListGraph::EdgeMap<double> edgeVarValue(graph);
        IloBool isIntegral = true;
        
        for (int i = 0; i < numEdges; ++i)
        {
            xSol[i] = getValue(edgeVar[i]);
            edgeVarValue[graph.edgeFromId(i)] = xSol[i] + 1;
            
            if (std::abs(xSol[i] - (int)xSol[i]) > TOLERANCE)
                isIntegral = false;
        };

        if (isIntegral)
            return;

        lemon::GomoryHu<lemon::ListGraph, lemon::ListGraph::EdgeMap<double>> tree(graph, edgeVarValue);
        tree.run();

        int count = 0;
        for (lemon::ListGraph::NodeIt it(graph); it != lemon::INVALID; ++it)
        {
            lemon::ListGraph::Node jt = tree.predNode(it);
            if (2 - tree.predValue(it) > TOLERANCE)
            {
                count++;
                IloExpr expr(masterEnv);
                for (lemon::GomoryHu<lemon::ListGraph, lemon::ListGraph::EdgeMap<double>>::MinCutEdgeIt e(tree, it, jt); e != lemon::INVALID; ++e)
                {
                    expr += edgeVar[graph.id((lemon::ListGraph::Edge)e)];
                };
                add(expr >= 2).end();
                expr.end();
            };
        };

        xSol.end();
        end = clock();
        std::cout << "Number of min cut " << count << std::endl;
        std::cout << "Time to add user cut " << (double)(end - start) / CLOCKS_PER_SEC << std::endl;

        return;
    };

    ILOLAZYCONSTRAINTCALLBACK4(FairLazyCallback, IloIntVarArray, edgesVar, IloNumVarArray, v, const lemon::ListGraph&, graph, const lemon::ListGraph::EdgeMap<int>&, edgeMap)
    {
        IloEnv masterEnv = getEnv();
        IloInt numEdges = lemon::countEdges(graph);
        IloInt numNodes = lemon::countNodes(graph);
        std::vector<int> edgesUsed;
        std::vector<int> valEdgesUsed;
        IloNumArray xSol(masterEnv, numEdges);
        IloNumArray vSol(masterEnv, numNodes);

        // get value of variables x
        for (int i = 0; i < numEdges; ++i) 
        {
            xSol[i] = getValue(edgesVar[i]);
            lemon::ListGraph::Edge edge = graph.edgeFromId(i);
            if (std::abs(xSol[i] - 1) < TOLERANCE)
            {
                edgesUsed.push_back(i);
                valEdgesUsed.push_back(edgeMap[edge]);
            };
        };

        // get value of variables v
        for (int i = 0; i < numNodes; ++i)
        {
            vSol[i] = getValue(v[i]);
        };

        // get subtours from optimal solution
        std::vector<std::vector<int>> listSubtour;
        listSubtour = getSubtours(graph, edgesUsed);

        if (listSubtour.size() != 1)
        {
            for (std::vector<int> subtour : listSubtour)
            {
                IloExpr expr(masterEnv);
                for (int i : subtour)
                {
                    expr += edgesVar[i];
                };
                add(expr <= (IloInt)subtour.size() - 1);
                expr.end();
            };
        };

        // add sort condition
        std::vector<int> orderIndex;
        sortEdges(edgesUsed, valEdgesUsed, orderIndex);

        IloExpr expr(masterEnv);
        float sumV, sumCX;
        for (int i = orderIndex.size() - 1; i != 0; --i)
        {
            int index = orderIndex.size() - 1 - i;
            sumV += vSol[index];
            lemon::ListGraph::Edge edge = graph.edgeFromId(orderIndex[index]);
            sumCX += edgeMap[edge];
            expr += v[index] - edgeMap[edge] * edgesVar[index];
            if (sumV < sumCX)
                add(expr >= 0);
        };

        xSol.end();
        return;
    };


    class FairTsp
    {
        const lemon::ListGraph &graph;
        const lemon::ListGraph::EdgeMap<int> &edgeMap;
        std::vector<float> weightFair;

    public:
        explicit FairTsp(const lemon::ListGraph& graph, const lemon::ListGraph::EdgeMap<int>& edgeMap) : graph{graph}, edgeMap{edgeMap} {};
        //explicit MinCut(const lemon::ListGraph& graph, const lemon::ListGraph::EdgeMap<WT>& edgeMap) : graph{graph}, edgeMap{edgeMap}  {};
        void readWeightFair(std::string fileName){
            std::ifstream file(fileName);
            std::string line;
            getline(file, line);
            std::vector<std::string> strWeight = splitSentence(line, ',');
            for (int i = 0; i < strWeight.size(); ++i)
            {
                std::cout << "weight fair is " << strWeight[i] << std::endl;
                this->weightFair.push_back(std::stof(strWeight[i]));
            };
        };

        void solve(){
            IloEnv masterEnv;
            IloModel masterMod(masterEnv, "fair-tsp");
            IloInt numEdges = lemon::countEdges(this->graph);
            IloInt numNodes = lemon::countNodes(this->graph);
            IloIntVarArray x(masterEnv, numEdges);
            IloNumVarArray v(masterEnv, numNodes);

            createILPfairtsp(masterMod, x, v, this->graph, this->edgeMap, this->weightFair);
            std::cout << "CREATED MASTER ILP" << std::endl;

            try
            {
                IloCplex masterCplex(masterMod);
                masterCplex.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);
                masterCplex.setParam(IloCplex::Param::Threads, 1);
                masterCplex.setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional);
                masterCplex.setParam(IloCplex::Param::Read::DataCheck, 2);
                masterCplex.setParam(IloCplex::Param::MIP::Display, 0);
                masterCplex.setParam(IloCplex::Param::MIP::Interval, 1);

                masterCplex.use(FairLazyCallback(masterEnv, x, v, this->graph, this->edgeMap));
                masterCplex.use(SubtourUserCallback(masterEnv, x, this->graph));

                masterCplex.exportModel("model.lp");

                if (masterCplex.solve())
                {
                    masterEnv.out() << endl << "Solution status: " << masterCplex.getStatus() << endl;

                    masterEnv.out() << "Objective value: "
                                    << masterCplex.getObjValue() << endl;

                    masterCplex.writeSolution("solution.sol");
                };
            }
            catch(const IloException& e)
            {
                std::cerr << "CPLEX found the following exception: " << e << std::endl;
            };
            
        };
    };

    std::vector<std::vector<int>> getSubtours(const lemon::ListGraph& graph, std::vector<int> edgesUsed)
    {
        std::vector<std::vector<int>> listSubtour;
        std::vector<int> subtour;
        lemon::ListGraph::Edge edge;
        int edgeId;
        lemon::ListGraph::Edge nextEdge;
        lemon::ListGraph::Node source, target;
        bool found = false;
        int searchEdgeId = 0;
        int nextEdgeId;

        while (edgesUsed.size() > 0) 
        {
            edgeId = edgesUsed[0];
            subtour.push_back(edgeId);
            edgesUsed.erase(edgesUsed.begin());
            
            edge = graph.edgeFromId(edgeId);
            source = graph.u(edge);
            target = graph.v(edge);

            while (source != target) 
            {
                found = false;

                for (lemon::ListGraph::IncEdgeIt e(graph, target); e != lemon::INVALID; ++e)
                {
                    if (edge != e)
                    {
                        searchEdgeId = graph.id(lemon::findEdge(graph, graph.u(e), graph.v(e)));
                        if (std::find(edgesUsed.begin(), edgesUsed.end(), searchEdgeId) != edgesUsed.end())
                        {
                            if (found) 
                            {
                                std::cout << "Invalid solution: There are 2 incident edges" << std::endl;
                            };
                            found = true;
                            nextEdge = e;
                            nextEdgeId = searchEdgeId;
                        };
                    };
                };

                if (found)
                {
                    edgesUsed.erase(std::remove(edgesUsed.begin(), edgesUsed.end(), nextEdgeId), edgesUsed.end());
                    subtour.push_back(nextEdgeId);
                    if (target != graph.u(nextEdge))
                    {
                        target = graph.u(nextEdge);
                    }
                    else
                    {
                        target = graph.v(nextEdge);
                    };                    
                }
                else
                {
                    break;
                }
                 
            };
            listSubtour.push_back(subtour);
            subtour.clear();
        };

        return listSubtour;
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
}


#endif