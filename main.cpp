/*
 * File:   main.cpp
 * Author: einstein
 *
 * Created on 18 de Maio de 2011, 15:43
 */

#include <cstdlib>
#include <cstdio>
#include <string>
#include "Timer.hpp"
#include "tabuSearch.hpp"
#include "CapacitadedGraph.hpp"


using namespace std;

/*
 *
 */

int main(int argc, char** argv) {
    CapacitadedGraph * g;
    //19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
    //83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149,
    //151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211,
    int tenure = 59;
    int inimax = 100;
    int maxIte = 200;
    string nomeDaInstancia = "";
    string caminhoArquivo = "";
    srand(7);
    if (argc == 1) {
    	nomeDaInstancia = "TA25.dat";
    	caminhoArquivo = "/home/hedley/workspace/CG/instances/";
        g = new CapacitadedGraph((caminhoArquivo+nomeDaInstancia).c_str());
        g->F = 5;
        g->P = 5;

    } else {
        g = new CapacitadedGraph(argv[1]);
        g->P = 0;
        if (argc >= 3) sscanf(argv[2], "%lf", &g->F);
        if (argc >= 4) sscanf(argv[3], "%d", &g->P);
        if (argc >= 5) sscanf(argv[4], "%d", &tenure);
        if (argc >= 6) sscanf(argv[5], "%d", &inimax);
        if (argc >= 7) sscanf(argv[6], "%d", &maxIte);
    }

    globalTimer.start();

    TabuSearch * tabu;

    if (g->P > 0){
        tabu = new TabuSearch_P(g);
    }else{
        tabu = new TabuSearch(g);
    }

    tabu->buscaTabu(tenure,inimax,maxIte);
    tabu->best->printToPs(g ,argv[1]);
    tabu->best->printSolutionToFile(nomeDaInstancia,globalTimer,
    		tabu->getRelacaoClustersVertices());

//    CG * cg = new CG(g,model);
//    cg->trivialInit();
//    cg->findIntSolution();
//    delete cg;

    delete tabu;
    delete g;
    cout << "Total time: " << globalTimer.getSec() << endl;

    //segundaInstancia(argc,argv);
    //terceiraInstancia(argc,argv);


    return 0;
}

