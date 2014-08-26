/*
 * File:   facilityModel.hpp
 * Author: einstein
 *
 * Created on 30 de Setembro de 2011, 20:08
 */

#ifndef FACILITYMODEL_HPP
#define	FACILITYMODEL_HPP
#include <cstdio>
#include <iostream>
#include "tabuSearch.hpp"
#include "gurobi_c++.h"
#include "CapacitadedGraph.hpp"
#include "utils.hpp"


class FacilityModel {
private:
    CapacitadedGraph * graph;
    GRBEnv * env;
    GRBModel * grb;
    GRBVar * var;

    int x(int i, int j) {
        return i * jsize + j;
    }


public:
    int jsize;
    Cluster * localC;
    Tabu * tabu;

    FacilityModel(CapacitadedGraph * g) {
        graph = g;
        try {
            env = new GRBEnv();

            env->set(GRB_IntParam_Method, 4);
            //            env->set(GRB_IntParam_Threads,2);
            env->set(GRB_IntParam_OutputFlag, 0);
            // env->set(GRB_StringParam_LogFile, "log.txt");
            env->set(GRB_IntParam_Presolve, 1);
            //env->set(GRB_DoubleParam_MIPGap, 0.2 / 100);
            //env->set(GRB_IntParam_SolutionLimit, 2);
            //env->set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_OPTIMALITY); //focus on bound

//            env->set(GRB_IntParam_CoverCuts, 2);
            grb = NULL;
            tabu = NULL;

            //            printModel();
            // env->set(GRB_IntParam_MIPFocus,3);
            localC = new Cluster[graph->size];
            for (int i = 0; i < graph->size; i++) {
                localC[i].v = new int[graph->size];
                localC->size = 0;
            }

        } catch (GRBException e) {
            cout << "Erro on gurobi model init" << endl;
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
            exit(e.getErrorCode());
        }
    }

    ~FacilityModel() {
        for (int i = 0; i < graph->size; i++) {
            delete [] localC[i].v;
        }
        delete [] localC;

        if (grb != NULL) delete grb;
        delete env;
    }

    void printModel() {
        grb->update();
        grb->write("model.lp");
    }

    GRBVar varx(int i, int j) {
        return var[x(i, j)];
    }

    char search(double &inicost, Cluster * cluster, int &cSize) {
        double max = 0;
        for (int i = 0; i < cSize; i++) {
            for (int j = 0; j < cluster[i].size; j++) {
                if (max < graph->dist(cluster[i].cx, cluster[i].cy, cluster[i].v[j]))
                    max = graph->dist(cluster[i].cx, cluster[i].cy, cluster[i].v[j]);
            }

        }
        char flag = 0;
        while (runModel(inicost, cluster, cSize, max + EPSLON)) {
            flag = 1;
            cout << "Facility " << inicost << endl;
            max = 0;
            for (int i = 0; i < cSize; i++) {
                for (int j = 0; j < cluster[i].size; j++) {
                    if (max < graph->dist(cluster[i].cx, cluster[i].cy, cluster[i].v[j]))
                        max = graph->dist(cluster[i].cx, cluster[i].cy, cluster[i].v[j]);
                }

            }

        }

        return flag;
    }

    void init(double &inicost, Cluster * cluster, int &cSize) {
        double max = INFINITY;
        char flag = 0;
        inicost = INFINITY;
        while (runModel(inicost, cluster, cSize, max + EPSLON)) {

            flag = 1;
            //cout << "Facility " << inicost << endl;
            max = 0;
            for (int i = 0; i < cSize; i++) {
                for (int j = 0; j < cluster[i].size; j++) {
                    if (max < graph->dist(cluster[i].cx, cluster[i].cy, cluster[i].v[j]))
                        max = graph->dist(cluster[i].cx, cluster[i].cy, cluster[i].v[j]);
                }

            }

        }


    }

    char runModel(double &inicost, Cluster * cluster, int &cSize, double cut) {
        //tem q limpar a solução antes de usar esse metodo
        //        cout << "Creating model" << endl;
        char in[graph->size][cSize];
        char name[20];
        double d;
        try {
            jsize = cSize;
            if (grb) {
                delete grb;
            }
            grb = new GRBModel(*env);
            grb->set(GRB_IntAttr_ModelSense, 1); //minimizar + 1



            var = grb->addVars(graph->size * cSize, GRB_BINARY);
            grb->update();
            GRBLinExpr exp;

            //exp = graph->F * jsize;
            for (int j = 0; j < jsize; j++)
                if (cluster[j].size > 0 || inicost == INFINITY) {
                    for (int i = 0; i < graph->size; i++) {
                        d = graph->dist(cluster[j].cx, cluster[j].cy, i);

                        if (d < cut) {
                            sprintf(name, "x%d_%d", i, j);
                            varx(i, j).set(GRB_StringAttr_VarName, name);

                            in[i][j] = 1;
                            exp += d*d * varx(i, j);

                        } else {
                            in[i][j] = 0;
                            grb->remove(varx(i, j));
                        }

                    }
                } else {
                    cout << "empty cluster in model erro" << endl;
                    exit(1);
                }



            grb->setObjective(exp);
            if (tabu != NULL)
                for (int i = 0; i < MAXTENURE; i++) {
                    if (tabu[i].C > -1 && tabu[i].C < jsize && tabu[i].V > -1 && in[tabu[i].V][tabu[i].C]) {
                        grb->remove(varx(tabu[i].V, tabu[i].C));
                        in[tabu[i].V][tabu[i].C] = 0;
                    }

                }

            //cout << "\tadding cover constraints" << endl;
            for (int i = 0; i < graph->size; i++) {
                exp.clear();
                sprintf(name, "Cover%d", i);
                for (int j = 0; j < jsize; j++)
                    if (in[i][j]) {
                        exp += varx(i, j);
                    }
                grb->addConstr(exp, '=', 1, name);
            }



            //        cout << "\tadding pack constraints" << endl;

            for (int j = 0; j < jsize; j++) {
                sprintf(name, "Pack%d", j);
                exp.clear();
                exp = -graph->C;
                for (int i = 0; i < graph->size; i++)
                    if (in[i][j]) {
                        exp += graph->w[i] * varx(i, j);
                    }
                grb->addConstr(exp, '<', 0, name);
            }




            grb->update();
            if (inicost < INFINITY)
                for (int j = 0; j < cSize; j++) {

                    for (int i = 0; i < cluster[j].size; i++)
                        if (in[cluster[j].v[i]][j]) {
                            varx(cluster[j].v[i], j).set(GRB_DoubleAttr_Start, 1.0);
                        }
                }



            //            printModel();
            //grb->update();
            grb->optimize();

            //            if (grb->get(GRB_DoubleAttr_ObjVal) + EPSLON > inicost) {
            //                return 0;
            //            }

            double cost = 0;
            cSize = 0;
            for (int j = 0; j < jsize; j++) {

                localC[cSize].size = 0;
                localC[cSize].flagNew = 1;
                localC[cSize].freeDemand = graph->C;


                for (int i = 0; i < graph->size; i++)
                    if (in[i][j] && varx(i, j).get(GRB_DoubleAttr_X) + EPSLON > 1) {
                        localC[cSize].freeDemand -= graph->w[i];
                        localC[cSize].v[localC[cSize].size++] = i;
                    }
                localC[cSize].cost = graph->cost(localC[cSize].v, localC[cSize].size, localC[cSize].cx, localC[cSize].cy);
                cost += localC[cSize].cost;
                cSize++;

            }

            //cout << "facility " << cSize << "*" << graph->F << " + " << cost - (cSize * graph->F) << " = " << cost << endl;
            //if (inicost > cost + EPSLON) {
                inicost = cost;
                //cout << "*";
                for (int i = 0; i < cSize; i++) {
                    if (cluster[i].size != localC[i].size
                            || cluster[i].freeDemand != localC[i].freeDemand
                            || fabs(cluster[i].cx - localC[i].cx) > EPSLON
                            || fabs(cluster[i].cy - localC[i].cy) > EPSLON)
                        cluster[i].flagNew = 1;

                    cluster[i].cost = localC[i].cost;
                    cluster[i].cx = localC[i].cx;
                    cluster[i].cy = localC[i].cy;
                    cluster[i].freeDemand = localC[i].freeDemand;
                    cluster[i].size = localC[i].size;
                    for (int j = 0; j < localC[i].size; j++) {
                        cluster[i].v[j] = localC[i].v[j];
                    }
                }

               // return 1;
           // }
            return 0;
        } catch (GRBException e) {
            cout << "Erro on gurobi model " << endl;
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
            exit(e.getErrorCode());
        }
    }
};




#endif	/* FACILITYMODEL_HPP */

