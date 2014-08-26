/*
 * File:   solution.hpp
 * Author: einstein
 *
 * Created on 13 de Fevereiro de 2012, 14:16
 */

#ifndef SOLUTION_HPP
#define	SOLUTION_HPP
#include <iostream>
#include "Box.hpp"

using namespace std;

struct Cluster {
	int id;
    /**current number of vertex in this cluster*/
    int size;
    /**vector of vertex of this cluster*/
    int * v;
    /**current cost of this cluster*/
    double cost;
    /**current anount free demand on this cluster*/
    int freeDemand;
    /**flag thar sinalize whethe this cluster was altered*/
    char flagNew;
    /**center of the cluster*/
    double cx, cy;
    /**caixa que envolve os pontos do cluster*/
    Box box;

    int getId(){
    	return id;
    }


};

/**Solução viavel*/
class Solution {
private:
    /**numero de clusters alocados na memoria*/
    int allocedCluster;
        /**Grafo de origem da solução*/
    CapacitadedGraph * graph;
public:
    /**numero de clusters*/
    int size;
    /**soma dos custos dos closter, sem custo de abertura*/
    double cost;
    /**clusters*/
    Cluster * cluster;
    /**cluster de cada vertice*/
    int * clusterOf;
    /*tempo até a melhor solução*/
    double timeToBest;

    Solution(CapacitadedGraph * g) {
        graph = g;
        clusterOf = new int[graph->size];

        if (graph->P > 0)
            allocedCluster = graph->P;
        else
            allocedCluster = graph->size;
        cluster = new Cluster[allocedCluster];

        for (int i = 0; i < allocedCluster; i++) {
            cluster[i].size = 0;
            cluster[i].v = new int[graph->size];
            cluster[i].cost = 0;
            cluster[i].freeDemand = graph->C;
            cluster[i].id = i;
            clusterOf[i] = -1;
        }
        size = 0;
    }

    virtual ~Solution() {
        delete [] clusterOf;
        for (int i = 0; i < allocedCluster; i++) free(cluster[i].v);
        delete [] cluster;
    }

    void printSolutionToFile(string nomeArquivo,MyTimer timer,
    		int* clustersDosVerticesDaMelhorSolucao){
    		int tempo = timer.getSec();
        	ofstream os;
        	os.open(nomeArquivo.c_str(),ios::out | ios::trunc);
        	os << "custo final: " << this->cost << endl;
        	os << "com tempo total de: " << tempo << "s" << endl;

        	for(int i = 0;i < graph->P; i++){
        		os << "Vértices do cluster "  << i << ": ";
        		for(int j = 0; j < cluster[i].size; j++){
        			os  << cluster[i].v[j] << " , ";
        		}
        		os << "\n";
        	}

        	//os << "time to best: " << timeToBest;
        	//os.flush();
        	os.close();
       }

    /**imprime a solução*/
    void printSolution(const CapacitadedGraph * g) {
        //    cout << "<#clusters>" << endl << "<cost>" << endl;
        //    cout << "<clusterCount>\t<#vertex>\t<demand>\t<cost>\t<cx>\t<cy>  :\t <V1>\t<V2> ..." << endl;
        //    cout << "..." << endl;

        cout << size << endl;
        if (g->P == 0)
            cout << cost << endl;
        else
            cout << cost + g->F * size << endl;
        for (int i = 0; i < size; i++) {
            cout << i << "\t" << cluster[i].size << "\t" << (g->C - cluster[i].freeDemand) << "\t" <<
                    cluster[i].cost << "\t" << cluster[i].cx <<
                    "\t" << cluster[i].cy << " :\t";
            for (int j = 0; j < cluster[i].size; j++) {
                cout << cluster[i].v[j] << "\t";
            }
            cout << endl;

        }

    }

    void printToPs(const CapacitadedGraph * g, char * name) {
        Box box;
        box.isNull = 0;
        box.maxX = box.minX = g->x[0];
        box.maxY = box.minY = g->y[0];
        for (int i = 1; i < g->size; i++) box.add(i, g);
        double scale = max((box.maxX - box.minX) / 570, (box.maxY - box.minY) / 810);

        ofstream psfile;
        char m[100];
        sprintf(m, "%s.ps", name);
        cout << m << endl;
        psfile.open(m);
        psfile << "\%!PS-Adobe-3.0\n\%\%Creator: flop\n\%\%Pages: 1\n\%\%Orientation: Portrait\n" <<
                "\%\%BoundingBox: 0 0 596 842\n\%\%HiResBoundingBox: 0 0 596 842\n\%\%DocumentMedia: plain 596 842 0 () ()\n" <<
                "\%\%EndComments\n\%\%Page: 1 1\n\%\%layOutinfor stands Dim 1600 stands Num 223\n";
        psfile << "0.2 setlinewidth" << endl;
        for (int i = 0; i < size; i++) {
            //psfile << "0." << i % 5 << " setgray\n" << endl;
            psfile << "0." << rand() % 7 << " 0." << rand() % 7 << " 0." << rand() % 7 << " setrgbcolor" << endl;
            for (int j = 0; j < cluster[i].size; j++) {
                psfile << 10 + (cluster[i].cx - box.minX) / scale << " " << 10 + (cluster[i].cy - box.minY) / scale << " moveto" << endl;
                psfile << 10 + (g->x[cluster[i].v[j]] - box.minX) / scale << " " << 10 + (g->y[cluster[i].v[j]] - box.minY) / scale << " lineto" << endl;
                psfile << "gsave\nstroke" << endl;
                psfile << 10 + (g->x[cluster[i].v[j]] - box.minX) / scale << " " << 10 + (g->y[cluster[i].v[j]] - box.minY) / scale << " 4 0 360 arc closepath " << endl;
                psfile << "fill" << endl;
            }
        }


        psfile.close();
    }

    int getClusterOfPorIndice(int indice){
       	return clusterOf[indice];
    }

};

int comparaCresSol(const void * a, const void * b) {
    if (((Solution*) a)->cost > ((Solution*) b)->cost)
        return 1;
    if (((Solution*) a)->cost < ((Solution*) b)->cost)
        return -1;
    return 0;
}


#endif	/* SOLUTION_HPP */

