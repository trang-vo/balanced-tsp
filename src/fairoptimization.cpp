#include "fairoptimization.h"

namespace fairtsp {
    
    std::vector<std::vector<std::vector<int>>> getSubtours(std::vector<std::vector<int>> edgeUsed) {
        std::vector<std::vector<std::vector<int>>> ans;
        std::vector<std::vector<int>> subtour;

        while (edgeUsed.size() > 0) {
            std::vector<int> edge = edgeUsed[0];
            subtour.push_back(edge);
            edgeUsed.erase(edgeUsed.begin());

            int source = edge[0];
            int target = edge[1];
            int foundId;
            std::vector<int> nextEdge;
            int found;

            while (source != target) {
                found = 0;
                for (int j = 0; j < edgeUsed.size(); ++j) {
                    if (target == edgeUsed[j][0]) {
                        if (found > 2) {
                            std::cout << "Invalid solution: There are 2 incident edges" << std::endl;
                        }
                        found++;
                        nextEdge = edgeUsed[j];
                        foundId = j;
                    };
                };

                if (found > 0) {
                    edgeUsed.erase(edgeUsed.begin() + foundId);
                    subtour.push_back(nextEdge);
                    target = nextEdge[1];
                }
                else {
                    std::cout << "Can not find a tour" << std::endl;
                    break;
                };
            };
            ans.push_back(subtour);
            subtour.clear();
        };
        return ans;
    };

    ILOLAZYCONSTRAINTCALLBACK2(FairOptLazyCallback, IloArray<IloIntVarArray>, x, const lemon::ListGraph&, graph)
    {
        IloEnv masterEnv = getEnv();
        IloInt numEdges = lemon::countEdges(graph);
        IloInt numNodes = lemon::countNodes(graph);
        
        std::vector<std::vector<int>> edgesUsed;
        IloNumArray2 xSol(masterEnv, numNodes);

        for (int i = 0; i < numNodes; ++i) {
            xSol[i] = IloNumArray(masterEnv, numNodes);
            for (int j = 0; j < numNodes; ++j) {
                xSol[i][j] = getValue(x[i][j]);
                if (std::abs(xSol[i][j] - 1) < TOLERANCE) {
                    std::vector<int> tmp = {i, j};
                    edgesUsed.push_back(tmp);
                };
            };
        };

        std::vector<std::vector<std::vector<int>>> listSubtour;
        listSubtour = getSubtours(edgesUsed);

        if (listSubtour.size() != 1) {
            for (std::vector<std::vector<int>> subtour : listSubtour) {
                IloExpr expr(masterEnv);
                for (std::vector<int> edge : subtour)
                {
                    expr += x[edge[0]][edge[1]];
                };
                add(expr <= (IloInt)subtour.size() - 1);
                expr.end();
            };
        };
        xSol.end();
        return;
    };

    ILOLAZYCONSTRAINTCALLBACK5(BenderLazyCallback, IloNumVar, t, IloNumVar, u, IloArray<IloIntVarArray>, x, const lemon::ListGraph&, graph, const lemon::ListGraph::EdgeMap<int>&, edgeMap)
    {
        IloEnv masterEnv = getEnv();
        IloInt numEdges = lemon::countEdges(graph);
        IloInt numNodes = lemon::countNodes(graph);
        
        std::vector<std::vector<int>> edgesUsed;
        IloNumArray2 xSol(masterEnv, numNodes);
        std::vector<int> costs;

        for (int i = 0; i < numNodes; ++i) {
            xSol[i] = IloNumArray(masterEnv, numNodes);
            for (int j = 0; j < numNodes; ++j) {
                xSol[i][j] = getValue(x[i][j]);
                if (std::abs(xSol[i][j] - 1) < TOLERANCE) {
                    std::vector<int> tmp = {i, j};
                    edgesUsed.push_back(tmp);
                };
            };
        };

        std::vector<std::vector<std::vector<int>>> listSubtour;
        listSubtour = getSubtours(edgesUsed);

        if (listSubtour.size() > 1) {
            for (std::vector<std::vector<int>> subtour : listSubtour) {
                IloExpr expr(masterEnv);
                for (std::vector<int> edge : subtour)
                {
                    expr += x[edge[0]][edge[1]];
                };
                add(expr <= (IloInt)subtour.size() - 1);
                expr.end();
            };
        } else {
	    float maxSum = 0;
	    int maxInd = 0;
	    for (int i = 0; i < numNodes; ++i){
	      float sum_i = 0;
	      for(int j = 0; j < numNodes; ++j){
		lemon::ListGraph::Edge edge = lemon::findEdge(graph, graph.nodeFromId(i), graph.nodeFromId(j));
		sum_i += edgeMap[edge] * xSol[i][j];
	      };
	      if (sum_i > maxSum) {
		maxInd = i;
		maxSum = sum_i;
	      };
	    };
            IloExpr expr(masterEnv);
	    for (int j = 0; j < numNodes; ++j) {
	      lemon::ListGraph::Edge edge = lemon::findEdge(graph, graph.nodeFromId(maxInd), graph.nodeFromId(j));
	      expr += edgeMap[edge] * x[maxInd][j];
	    };
            add(t >= u - expr);
        };
        xSol.end();
        return;
    };

    ILOUSERCUTCALLBACK4(LocalCutCallback, IloArray<IloIntVarArray>, x, IloIntVarArray, li, const lemon::ListGraph&, graph, const lemon::ListGraph::EdgeMap<int>&, edgeMap) {
        if ( !isAfterCutLoop() )
            return;

        if (getNnodes() < 100)
            return;

        IloEnv masterEnv = getEnv();
        IloInt numEdges = lemon::countEdges(graph);
        IloInt numNodes = lemon::countNodes(graph);
        
        IloNumArray2 xSol(masterEnv, numNodes);
        IloNumArray m(masterEnv, numNodes);
        IloNumArray liSol(masterEnv, numNodes);

        std::cout << "CREATE LOCAL CUTS" << std::endl;
        IloExpr expr(masterEnv);
        for (int i = 0; i < numNodes; ++i)
        {
            xSol[i] = IloNumArray(masterEnv, numNodes);
            m[i] = 0;
            liSol[i] = getValue(li[i]);
            for (int j = 0; j < numNodes; ++j)
            {
                xSol[i][j] = getValue(x[i][j]);
                lemon::ListGraph::Edge edge = lemon::findEdge(graph, graph.nodeFromId(i), graph.nodeFromId(j));
                m[i] += edgeMap[edge] * xSol[i][j];
            };
            std::cout << "liSol[" << i << "] = " << liSol[i] << " mi = " << m[i] << std::endl;
            if (std::abs(liSol[i] - m[i]) < TOLERANCE) {
                addLocal(expr >= IloCeil(m[i]));
                std::cout << "Add local cut for vertex " << i << std::endl;
            };
        };

        return;
    };

    ILOMIPINFOCALLBACK1(InfoCallback, std::ofstream&, outputFile) {
        if (getNnodes() % 100 == 0) {
            // std::cout << "Objective = " << getBestObjValue() << std::endl;
            // std::cout << "Time = " << getCplexTime() - getStartTime() << std::endl;
            // std::cout << "GAP = " << getMIPRelativeGap() << std::endl;
            outputFile << getIncumbentObjValue() << "," << getCplexTime() - getStartTime() << "," << getMIPRelativeGap() << std::endl;
        };
    };

    void FairOpt::createLPformulationMaxminDirected(IloModel model, IloArray<IloIntVarArray> x, IloNumVar u, IloNumVar l, IloNumVar t) {
        typedef lemon::ListGraph Graph;

        IloEnv masterEnv = model.getEnv();
        int numEdges = lemon::countEdges(this->graph);
        int numNodes = lemon::countNodes(this->graph);

        IloRangeArray inbound_arcs(masterEnv, numNodes);
        IloRangeArray outbound_arcs(masterEnv, numNodes);
        std::stringstream name;

        // create variables x
        for(uint i = 0; i < numNodes; ++i) {
            x[i] = IloIntVarArray(masterEnv, numNodes);
            for(uint j = 0; j < numNodes; ++j) {
                name << "x." << i + 1 << "." << j + 1;
                x[i][j] = IloIntVar(masterEnv, 0, 1);
                x[i][j].setName(name.str().c_str());
                name.str(""); // Clean name
            };
        };

        // create inbound constraints
        IloExpr expr(masterEnv);
        for(auto i = 0u; i < numNodes; ++i) {
            for(auto j = 0u; j < numNodes; ++j) {
                expr += x[j][i];
            }

            name << "inbound_" << i;
            inbound_arcs[i] = IloRange(masterEnv, 1, expr, 1, name.str().c_str());
            name.str(""); // Clean name
            expr.clear(); // Clean expr
        };
        model.add(inbound_arcs);

        // create outbound constraints
        for(auto i = 0u; i < numNodes; ++i) {
            for(auto j = 0u; j < numNodes; ++j) {
                expr += x[i][j];
            }

            name << "outbound_" << i;
            outbound_arcs[i] = IloRange(masterEnv, 1, expr, 1, name.str().c_str());
            name.str(""); // Clean name
            expr.clear(); // Clean expr
        };
        model.add(outbound_arcs);

        IloExpr objExpr(masterEnv);
        for (auto i = 0u; i < numNodes; ++i)
        {
            for(auto j = 0u; j < numNodes; ++j) {
                Graph::Edge edge = lemon::findEdge(this->graph, graph.nodeFromId(i), graph.nodeFromId(j));
                model.add(u - edgeMap[edge] * x[i][j] >= 0);
                expr += edgeMap[edge] * x[i][j];
            };
            // model.add(expr - u <= 0);
            model.add(expr - l >= 0);
            expr.clear();
        };
        model.add(t >= u - l);
        IloObjective obj(masterEnv, t, IloObjective::Minimize);

        // Add the objective function to the model
        model.add(obj);

        // Free the memory used by expr
        expr.end();
    }

    void FairOpt::createLPformulationBender(IloModel model, IloArray<IloIntVarArray> x, IloNumVar u, IloNumVar t) {
        typedef lemon::ListGraph Graph;

        IloEnv masterEnv = model.getEnv();
        int numEdges = lemon::countEdges(this->graph);
        int numNodes = lemon::countNodes(this->graph);

        IloRangeArray inbound_arcs(masterEnv, numNodes);
        IloRangeArray outbound_arcs(masterEnv, numNodes);
        std::stringstream name;

        // create variables x
        for(uint i = 0; i < numNodes; ++i) {
            x[i] = IloIntVarArray(masterEnv, numNodes);
            for(uint j = 0; j < numNodes; ++j) {
                name << "x." << i + 1 << "." << j + 1;
                x[i][j] = IloIntVar(masterEnv, 0, 1);
                x[i][j].setName(name.str().c_str());
                name.str(""); // Clean name
            };
        };

        // create inbound constraints
        IloExpr expr(masterEnv);
        for(auto i = 0u; i < numNodes; ++i) {
            for(auto j = 0u; j < numNodes; ++j) {
                expr += x[j][i];
            }

            name << "inbound_" << i;
            inbound_arcs[i] = IloRange(masterEnv, 1, expr, 1, name.str().c_str());
            name.str(""); // Clean name
            expr.clear(); // Clean expr
        };
        model.add(inbound_arcs);

        // create outbound constraints
        for(auto i = 0u; i < numNodes; ++i) {
            for(auto j = 0u; j < numNodes; ++j) {
                expr += x[i][j];
            }

            name << "outbound_" << i;
            outbound_arcs[i] = IloRange(masterEnv, 1, expr, 1, name.str().c_str());
            name.str(""); // Clean name
            expr.clear(); // Clean expr
        };
        model.add(outbound_arcs);

        IloExpr objExpr(masterEnv);
        for (auto i = 0u; i < numNodes; ++i)
        {
            for(auto j = 0u; j < numNodes; ++j) {
                Graph::Edge edge = lemon::findEdge(this->graph, graph.nodeFromId(i), graph.nodeFromId(j));
                model.add(u - edgeMap[edge] * x[i][j] >= 0);
            };
        };
        IloObjective obj(masterEnv, t, IloObjective::Minimize);

        // Add the objective function to the model
        model.add(obj);

        // Free the memory used by expr
        expr.end();
    };

    void FairOpt::createLPformulationAddVar(IloModel model, IloArray<IloIntVarArray> x, IloNumVar u, IloNumVar l, IloIntVarArray li) {
        typedef lemon::ListGraph Graph;

        IloEnv masterEnv = model.getEnv();
        int numEdges = lemon::countEdges(this->graph);
        int numNodes = lemon::countNodes(this->graph);

        IloRangeArray inbound_arcs(masterEnv, numNodes);
        IloRangeArray outbound_arcs(masterEnv, numNodes);
        std::stringstream name;

        // create variables x
        for(uint i = 0; i < numNodes; ++i) {
            x[i] = IloIntVarArray(masterEnv, numNodes);
            for(uint j = 0; j < numNodes; ++j) {
                name << "x." << i + 1 << "." << j + 1;
                x[i][j] = IloIntVar(masterEnv, 0, 1);
                x[i][j].setName(name.str().c_str());
                name.str(""); // Clean name
            };
        };

        // create inbound constraints
        IloExpr expr(masterEnv);
        for(auto i = 0u; i < numNodes; ++i) {
            for(auto j = 0u; j < numNodes; ++j) {
                expr += x[j][i];
            }

            name << "inbound_" << i;
            inbound_arcs[i] = IloRange(masterEnv, 1, expr, 1, name.str().c_str());
            name.str(""); // Clean name
            expr.clear(); // Clean expr
        };
        model.add(inbound_arcs);

        // create outbound constraints
        for(auto i = 0u; i < numNodes; ++i) {
            for(auto j = 0u; j < numNodes; ++j) {
                expr += x[i][j];
            }

            name << "outbound_" << i;
            outbound_arcs[i] = IloRange(masterEnv, 1, expr, 1, name.str().c_str());
            name.str(""); // Clean name
            expr.clear(); // Clean expr
        };
        model.add(outbound_arcs);

        for (uint i = 0u; i < numNodes; ++i)
        {
            li[i] = IloIntVar(masterEnv, 0, IloInfinity);
            name << "l_" << i + 1  ;
            li[i].setName(name.str().c_str());
            name.str("");
            for (uint j = 0u; j < numNodes; ++j)
            {
                if (i != j) {
                    Graph::Edge edge = lemon::findEdge(this->graph, graph.nodeFromId(i), graph.nodeFromId(j));
                    model.add(u - edgeMap[edge] * x[i][j] >= 0);
                    expr += edgeMap[edge] * x[i][j];
                }
            };
            model.add(expr - li[i] >= 0);
            model.add(l <= li[i]);
            expr.clear();
        };
        IloObjective obj(masterEnv, u-l, IloObjective::Minimize);
        model.add(u >= l);

        // Add the objective function to the model
        model.add(obj);

        // Free the memory used by expr
        expr.end();
    };

    IloCplex FairOpt::solve(int useCutVersion) {
        IloEnv masterEnv;
        IloModel masterMod(masterEnv, "fair-tsp");
        IloInt numEdges = lemon::countEdges(this->graph);
        IloInt numNodes = lemon::countNodes(this->graph);

        IloArray<IloIntVarArray> x(masterEnv, numEdges * 2);
        IloNumVar u(masterEnv, 0, IloInfinity, ILOFLOAT);
        u.setName("u");
        IloNumVar l(masterEnv, 0, IloInfinity, ILOFLOAT);
        l.setName("l");
        IloNumVar t(masterEnv, 0, IloInfinity, ILOFLOAT);
        t.setName("t");
        IloIntVarArray li(masterEnv, numNodes);

        // createLPformulationMaxminDirected(masterMod, x, u, l, t);
        // createLPformulationAddVar(masterMod, x, u, l, li);
        createLPformulationBender(masterMod, x, u, t);
        std::cout << "CREATED MASTER ILP" << std::endl;

        IloCplex masterCplex(masterMod);
        // masterCplex.exportModel("model_opt.lp");
        masterCplex.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);
        masterCplex.setParam(IloCplex::Param::Threads, 1);
        masterCplex.setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional);
        
        // config for print log 
        masterCplex.setParam(IloCplex::Param::MIP::Display, 2);
        // masterCplex.setOut(masterEnv.getNullStream());
        // masterCplex.setParam(IloCplex::Param::MIP::Limits::CutsFactor, 0);
        // masterCplex.setParam(IloCplex::Param::MIP::Interval, 10);

        // useCutVersion = 0 (no cut), 1 (all cuts), 2 (selection cuts)
        if (useCutVersion == 0) {
            masterCplex.setParam(IloCplex::Param::MIP::Limits::CutsFactor, 0);
        } else if (useCutVersion == 2) {
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::BQP, -1);
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::Cliques, -1);
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::Covers, -1);
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::Disjunctive, -1);
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, -1);
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::Gomory, 1);
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::GUBCovers, -1);
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::LocalImplied, -1);
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::LiftProj, 1);
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::MIRCut, 1);
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::MCFCut, -1);
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::RLT, -1);
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::ZeroHalfCut, -1);
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::Implied, -1);
            masterCplex.setParam(IloCplex::Param::MIP::Cuts::PathCut, -1);
        };

        // masterCplex.use(FairOptLazyCallback(masterEnv, x, this->graph));
	masterCplex.use(BenderLazyCallback(masterEnv, t, u, x, this->graph, this->edgeMap));
        // masterCplex.use(LocalCutCallback(masterEnv, x, li, this->graph, this->edgeMap));

        //std::stringstream name;
        //name << "results/" << this->instance << "_" << useCutVersion << "_resultLogFile.txt";
        //std::ofstream outputFile(name.str().c_str());
        //masterCplex.use(InfoCallback(masterEnv, outputFile));
        //name.str("");

        try
        {
            bool check = masterCplex.solve();
            if (check)
            {
                masterEnv.out() << endl << "Solution status: " << masterCplex.getStatus() << endl;

                masterEnv.out() << "Objective value: "
                                << masterCplex.getObjValue() << endl;

                IloNumArray2 xSol(masterEnv, numNodes);

                for (int i = 0; i < numNodes; ++i) {
                    xSol[i] = IloNumArray(masterEnv, numNodes);
                    for (int j = 0; j < numNodes; ++j) {
                        xSol[i][j] = masterCplex.getValue(x[i][j]);
                        if (std::abs(xSol[i][j] - 1) < TOLERANCE) {
                            std::cout << x[i][j].getName() << '\t' << edgeMap[lemon::findEdge(graph, graph.nodeFromId(i), graph.nodeFromId(j))] << std::endl;
                        };
                    };
                };
                // name << "results/" << this->instance << "_" << useCutVersion << "_solution.sol";
                // masterCplex.writeSolution(name.str().c_str());
            } else {
                masterEnv.out() << endl << "Solution status: " << masterCplex.getStatus() << endl;
            };
        }
        catch(const IloException& e)
        {
            std::cerr << "CPLEX found the following exception: " << e << std::endl;
        }; 
        return masterCplex;
    };
};
