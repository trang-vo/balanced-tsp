#ifndef FAIROPT_H
#define FAIROPT_H

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

    class FairOpt {
        const lemon::ListGraph &graph;
        const lemon::ListGraph::EdgeMap<int> &edgeMap;
        std::string instance;
        const float M = 0.5;

    public:
        explicit FairOpt(std::string instance, const lemon::ListGraph& graph, const lemon::ListGraph::EdgeMap<int>& edgeMap) : instance{instance}, graph{graph}, edgeMap{edgeMap} {};
        void createLPformulationMaxminDirected(IloModel model, IloArray<IloIntVarArray> x, IloNumVar u, IloNumVar l, IloNumVar t);
        void createLPformulationBender(IloModel model, IloArray<IloIntVarArray> x, IloNumVar u, IloNumVar t);
        void createLPformulationAddVar(IloModel model, IloArray<IloIntVarArray> x, IloNumVar u, IloNumVar l, IloIntVarArray li);
        IloCplex solve(int version);
    };

};

#endif 