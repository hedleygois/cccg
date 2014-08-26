// -*- C++ -*-
/*
 * File:   CapacitadedGraph.hpp
 * Author: einstein
 *
 * Created on 18 de Maio de 2011, 16:00
 */

#ifndef CAPACITADEDGRAPH_HPP
#define	CAPACITADEDGRAPH_HPP

//#define useOMP 10
#include <iostream>
#include <fstream>
#include <cctype>
#include <vector>
#include <string>
#include "utils.hpp"


#ifdef useOMP //ativa versão openmp
#include "omp.h"
#endif

using namespace std;

class CapacitadedGraph {
public:
    /**numero de vertices no problema*/
    int size;
    /** a capacidade máxima dos clusters*/
    int C;
    /** zero se não há numero de clusters fixado (g-CCCP) */
    int P;
    /**custo de abertura dos cluster na versão g-CCCP*/
    double F;
    /**peso de cada vertice*/
    int * w;
    /**coordenadas de cada vertice*/
    double * x, *y;
    /**custo dual usado na geração de coluna, (tirar daqui depois)*/
    double * dualCost;
    /**matrix de distancia do grafo*/
    double ** d;
    /**distância maxima*/
    double dmax;
    /**distância maxima*/
    double dmin;
    /**peso maximo*/
    double wmax;
    /**peso minimo*/
    double wmin;
    /**maxCluster*/
    int maxNCluster;

    CapacitadedGraph(const char * path) {
        int i;
        float f;
        int k;
        ifstream myfile;

#ifdef useOMP
        omp_set_num_threads(omp_get_num_procs());
#endif

        F = 0;
        P = 0;
        myfile.open(path);
        if (myfile.is_open()) {
            vector<int> X, Y, W;
            i = 0;
            myfile >> f; //
            size = (int) f;
            myfile >> f;
            while (i < size) {
                //                cout << str << " ";
                myfile >> f;
                X.push_back(f);
                //                cout << X[i] << " ";
                myfile >> f;
                Y.push_back(f);
                //                cout << Y[i] << " ";
                myfile >> f;
                C = (int) f;

                myfile >> f;
                k = (int) f;
                W.push_back(k);
                //                cout << W[i] << " ";
                //                cout << endl;
                i++;
            }
            size = i;
            x = new double[size];
            y = new double[size];
            w = new int[size];
            dualCost = new double[size];
            d = new double*[size];
            for (int i = 0; i < size; ++i) {
                x[i] = X[i];
                y[i] = Y[i];
                w[i] = W[i];
                d[i] = new double[size];
            }
            updateData();
            myfile.close();
            cout << "vertices : " << size << endl;
            cout << "capacity : " << C << endl;

        }

    }

    ~CapacitadedGraph() {
        if (x != NULL)
            delete [] x;
        if (y != NULL)
            delete [] y;
        if (dualCost != NULL)
            delete [] dualCost;
        if (w != NULL)
            delete [] w;

        if (d != NULL) {
            for (int i = 0; i < size; i++) {
                delete [] d[i];
            }
            delete [] d;
        }

    }

    /**preenche a matrix de custo*/
    void updateData() {
        dmax = 0;
        wmax = w[0];
        wmin = w[0];

        for (int i = 0; i < size; i++) {
            d[i][i] = 0;
            for (int j = i + 1; j < size; j++) {
                d[i][j] = d[j][i] = dist(i, j);
                if (dmax < d[i][j]) dmax = d[i][j];
            }
            if (wmax < w[i])
                wmax = w[i];
            else
                if (wmin > w[i])
                wmin = w[i];
        }
        cout << "dmax: " << dmax << endl;
        cout << "wmax: " << wmax << endl;
        cout << "wmin: " << wmin << endl;

        int temp[size];
        for (int i = 0; i < size; i++) temp[i] = i;
        qsort_r(temp,size,sizeof(int),comparaCresInt,w);



    }

    /**clacula a distancia entre o par de coordenadas e o vertice*/
    double dist(double cx, double cy, int index) {
        cx -= x[index];
        cx *= cx;
        cy -= y[index];
        cy *= cy;
        return sqrt(cx + cy);
    }

    /**distancia entre duas coordenadas */
    double dist(double x1, double y1, double x2, double y2) {
        x1 -= x2;
        x1 *= x1;
        y1 -= y2;
        y1 *= y1;
        return sqrt(x1 + y1);
    }

    /**distancia entre dois vertices */
    double dist(int i, int j) {
        double xx = x[i] - x[j];
        double yy = y[i] - y[j];
        xx *= xx;
        yy *= yy;
        return sqrt(xx + yy);
    }

    /**calcula a soma das distancia ao centroide do subconjunto*/
    double cost(int * subSet, int size) {
        double cx = x[*subSet], cy = y[*subSet], r = 0;
        subSet++;
        for (int i = 1; i < size; i++) {
            cx += x[*subSet];
            cy += y[*subSet];
            subSet++;
        }
        cx /= (double) size;
        cy /= (double) size;

        for (int i = 0; i < size; i++) {
            subSet--;
            r += dist(cx, cy, *subSet);
        }
        return r;
    }

    /**calcula a soma das distancia ao centroide do subconjunto e retorna
     as coordenadas do centroid*/
    double cost(int * subSet, int size, double &cx, double &cy) {
        double r = 0;
        cx = x[*subSet];
        cy = y[*subSet];
        subSet++;
        for (int i = 1; i < size; i++, subSet++) {
            cx += x[*subSet];
            cy += y[*subSet];
        }
        cx /= (double) size;
        cy /= (double) size;

        subSet--;
        for (int i = 0; i < size; i++, subSet--) {
            r += dist(cx, cy, *subSet);
        }

        return r;
    }

    /**calcula o custo do cluster apos a troca de dois vertice,
     retorna as coordenadas do novo centroid (cx,cy),
     vIn é o indice do vertice no GRAFO que entrará no cluster,
     out é o indice do vertice no CLUSTER que sairá*/
    double simulateCostSwap(int * subSet, const int size, double &cx, double &cy, const int vIn, const int out) {
        cx = (cx * size - x[subSet[out]] + x[vIn]) / size;
        cy = (cy * size - y[subSet[out]] + y[vIn]) / size;


#ifdef useOMP
        int i;
        double r = 0;
        int nt = omp_get_num_threads();

#pragma omp parallel for default(shared) private(i) reduction(+:r) schedule(static)
        for (i = 0; i < size; i++)
            r += dist(cx, cy, subSet[i]);


        r += dist(cx, cy, vIn) - dist(cx, cy, subSet[out]);
        //cout << r << endl;
#else
        double r = dist(cx, cy, vIn);
        for (int i = 0; i < out; i++, subSet++)
            r += dist(cx, cy, *subSet);
        subSet++;
        for (int i = out + 1; i < size; i++, subSet++)
            r += dist(cx, cy, *subSet);
        //cout << r << endl;
#endif

        return r;
    }

    /**calcula o custo do cluster apos a adição de um vertice,
     retorna as coordenadas do novo centroid (cx,cy),
     vIn é o indice do vertice no GRAFO que entrará no cluster */
    double simulateCostAdd(int * subSet, int size, double &cx, double &cy, int vIn) {
        cx = (cx * size + x[vIn]) / (size + 1);
        cy = (cy * size + y[vIn]) / (size + 1);

#ifdef useOMP
        double r = 0;
        int i;
#pragma omp parallel for default(shared) private(i) reduction(+:r) schedule(static)
        for (i = 0; i < size; i++) {
            r += dist(cx, cy, subSet[i]);
        }

        r += dist(cx, cy, vIn);
        //        cout << r << endl;
#else
        double r = dist(cx, cy, vIn);

        for (int i = 0; i < size; i++, subSet++) {
            r += dist(cx, cy, *subSet);
        }
        //cout << r << endl;
#endif
        return r;
    }

    /**calcula o custo do cluster apos a remoção de um vertice,
     retorna as coordenadas do novo centroid (cx,cy)
          out é o indice do vertice no CLUSTER que sairá*/
    double simulateCostRemove(int * subSet, int size, double &cx, double &cy, int out) {
        cx = (cx * size - x[subSet[out]]) / (size - 1);
        cy = (cy * size - y[subSet[out]]) / (size - 1);
        double r = 0;
        int i;
#ifdef useOMP
#pragma omp parallel for default(shared) private(i) reduction(+:r) schedule(static)
        for (i = 0; i < size; i++)
            r += dist(cx, cy, subSet[i]);

        r -= dist(cx, cy, subSet[out]);
        //        cout << r << endl;
#else
        for (i = 0; i < out; i++, subSet++)
            r += dist(cx, cy, *subSet);
        subSet++;
        for (i = out + 1; i < size; i++, subSet++)
            r += dist(cx, cy, *subSet);
        //        cout << r << endl;
#endif
        return r;
    }

    /**soma dos custos duais, (tirar daqui depois)*/
    double prizes(int * subSet, int size) {
        double r = dualCost[*subSet++];
        for (int i = 1; i < size; i++, subSet++) {
            r += dualCost[*subSet];
        }
        return r;
    }

}; //class CapacitadedGraph



#endif	/* CAPACITADEDGRAPH_HPP */
