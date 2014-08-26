/*
 * File:   utils.hpp
 * Author: einstein
 *
 * Created on 26 de Setembro de 2011, 09:27
 */

#ifndef UTILS_HPP
#define	UTILS_HPP

#define EPSLON 0.0001
//#define INFINITY 999999
//#define min(A,B) ((A<B)?A:B)
//#define max(A,B) ((A>B)?A:B)

//#define debug

#include <iostream>
#include <cctype>
#include <cmath>


using namespace std;


int comparaDec(const void * a, const void * b, void * cost) {
    double * score = (double*) cost;
    if (score[*(int*) a] > score[*(int*) b])
        return -1;
    if (score[*(int*) a] + EPSLON < score[*(int*) b])
        return 1;
    return 0;
}

int comparaCres(const void * a, const void * b, void * cost) {
    double * score = (double*) cost;
    if (score[*(int*) a] > score[*(int*) b])
        return 1;
    if (score[*(int*) a] < score[*(int*) b])
        return -1;
    return 0;
}


int comparaCresInt(const void* a, const void* b, void* cost) {
    int * score = (int*) cost;
    if (score[*(int*) a] > score[*(int*) b])
        return 1;
    if (score[*(int*) a] < score[*(int*) b])
        return -1;
    return 0;
}

int comparaDecInt(const void* a, const void* b, void* cost) {
    int * score = (int*) cost;
    if (score[*(int*) a] > score[*(int*) b])
        return -1;
    if (score[*(int*) a] < score[*(int*) b])
        return 1;
    return 0;
}


void swap(int &a, int &b) {
    int aux = a;
    a = b;
    b = aux;
}

void swap(double &a, double &b) {
    double aux = a;
    a = b;
    b = aux;
}

double dist(double x1, double y1, double x2, double y2) {
    x1 -= x2;
    y1 -= y2;
    x1 *= x1;
    y1 *= y1;
    return sqrt(x1 + y1);
}

double quadrado(double x) {
    return x*x;
}

#endif	/* UTILS_HPP */

