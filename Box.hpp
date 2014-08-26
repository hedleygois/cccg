// -*- C++ -*-
/*
 * File:   Box.hpp
 * Author: einstein
 *
 * Created on 9 de Junho de 2011, 18:47
 */

#ifndef BOX_HPP
#define	BOX_HPP
#include "CapacitadedGraph.hpp"
#include "utils.hpp"



//#include "CG.hpp"


class Box {
public:
    double maxX, maxY, minX, minY;
    char isNull;

    Box(double max_x, double max_y, double min_x, double min_y) {
        maxX = max_x;
        maxY = max_y;
        minX = min_x;
        minY = min_y;
        isNull = 0;
    }

    Box() {
        isNull = 1;
    }

    void inline copy(Box const &b) {
        maxX = b.maxX;
        maxY = b.maxY;
        minX = b.minX;
        minY = b.minY;
        isNull = b.isNull;
    }

    void inline setInfinity() {
        maxX = INFINITY;
        maxY = INFINITY;
        minX = -INFINITY;
        minY = -INFINITY;
        isNull = 0;
    }

    void inline inteception(Box const &b1, Box const &b2) {
        if (b1.isNull || b2.isNull) {
            isNull = 1;
            return;
        }
        maxX = min(b1.maxX, b2.maxX);
        maxY = min(b1.maxY, b2.maxY);
        minX = max(b1.minX, b2.minX);
        minY = max(b1.minY, b2.minY);
        if (maxX <= minX || maxY <= minY) isNull = 1;
        else isNull = 0;
    }

    void inline inteception(Box const &b2) {
        if (isNull || b2.isNull) {
            isNull = 1;
            return;
        }
        if (b2.maxX < maxX)
            maxX = b2.maxX;
        if (b2.maxY < maxY)
            maxY = b2.maxY;
        if (b2.minX > minX)
            minX = b2.minX;
        if (b2.minY > minY)
            minY = b2.minY;

        if (maxX <= minX || maxY <= minY) isNull = 1;
        else isNull = 0;
    }

    char inline intercept(Box const &b) {
        if (isNull || b.isNull) return 0;

        if (min(maxX, b.maxX) < max(minX, b.minX))
            return 0;

        if (min(maxY, b.maxY) < max(minY, b.minY))
            return 0;

        return 1;

    }

    void inline unionBox(Box const &b2) {
        if (isNull) {
            copy(b2);
            return;
        }
        if (b2.isNull) {
            return;
        }
        if (b2.maxX > maxX)
            maxX = b2.maxX;
        if (b2.maxY > maxY)
            maxY = b2.maxY;
        if (b2.minX < minX)
            minX = b2.minX;
        if (b2.minY < minY)
            minY = b2.minY;

        if (maxX <= minX || maxY <= minY) isNull = 1;
        else isNull = 0;
    }

    char has(const double x, const double y) {
        if (isNull) return 0;
        if (maxX >= x && minX <= x && maxY >= y && minY <= y) {
            return 1;
        }
        return 0;
    }


    void inline setBox(double x, double y, double r) {
        if (r < EPSLON) {
            isNull = 1;
            return;
        }
        isNull = 0;
        maxX = x + r;
        minX = x - r;
        maxY = y + r;
        minY = y - r;
    }

    void add(const int i, const CapacitadedGraph * graph) {
    if (maxX < graph->x[i])maxX = graph->x[i];
    if (maxY < graph->y[i])maxY = graph->y[i];
    if (minX > graph->x[i])minX = graph->x[i];
    if (minY > graph->y[i])minY = graph->y[i];
}
private:
};

#endif	/* BOX_HPP */

