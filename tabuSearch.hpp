// -*- C++ -*-
/*
 * File:   tabuSearch.hpp
 * Author: einstein
 *
 * Created on 16 de Junho de 2011, 09:48
 */

#ifndef TABUSEARCH_HPP
#define	TABUSEARCH_HPP

//#define TENURE 7 //TA instances
#define MAXTENURE 1000 //SJC instances
//#define TENURE 100
#define UM_DIA 84200


#include <iostream>
#include <iomanip> //setprecision
#include "Timer.hpp"
#include "utils.hpp"
#include "CapacitadedGraph.hpp"
#include "solution.hpp"

/**Regra tabu, vertice V não pode ser inserido no cluster C*/
struct Tabu {
	/**cluster*/
	int C;
	/**vertice*/
	int V;
};

/**Movimento*/
struct Move {
	double delta; // Variação no custo global causado no movimento
	int C_1; // cluster 1
	int C_2; // cluster 2
	int V_1; // vertice no cluster 1 a ser movido para o cluster 2
	int V_2; // vertice no cluster 2 a ser movido para o cluster 1, ou -1 se não se tratar de um swap.
};

class TabuSearch {
public:
	/**grafo*/
	CapacitadedGraph * graph;
	/**contador de cluster vazios na versão g-CCCP*/
	int emptyCt;
	/**lista de regras tabu*/
	Tabu * tabu;
	/**indice da lisa tabu*/
	int tabuHead;
	/**numero maximo de clusters */
	int numMaxCluster;
	/**Solução corrente*/
	Solution * curr;
	/**Melhor solução*/
	Solution * best;

	/**tamanho da lista tabu*/
	int tenure;

	TabuSearch() {}

	TabuSearch(CapacitadedGraph * g) {
		graph = g;
		numMaxCluster = graph->size;
		curr = new Solution(g);
		best = new Solution(g);
		tabu = new Tabu[MAXTENURE];
		tenure = MAXTENURE;
		cout << "Tabu data initialized " << globalTimer.getMillSec() << "s." << endl;
	}

	~TabuSearch() {
		delete best;
		delete curr;
		delete [] tabu;

	}

	void printCluster(int k) {
		cout << curr->cluster[k].cost << ": ";
		for (int i = 0; i < curr->cluster[k].size; i++) {
			cout << curr->cluster[k].v[i] << " ";
		}
		cout << endl;

	}

	/**seta todos os cluster como novos*/
	void renewCurr() {
		Cluster * c = curr->cluster;
		for (int i = 0; i < curr->size; i++, c++)
			c[i].flagNew = 1;
	}

	/**recalcula os custos de uma solução*/
	void updateCosts(Solution * sol) {
		double localCost = 0;
		for (int i = 0; i < sol->size; i++) {
			setBox(sol->cluster[i]);
			if (sol->cluster[i].size == 0)
				sol->cluster[i].cost = 0;
			else localCost += (sol->cluster[i].cost = graph->cost(sol->cluster[i].v, sol->cluster[i].size, sol->cluster[i].cx, sol->cluster[i].cy));
		}
		sol->cost = localCost;
	}

	/**debugar*/
	char check(int * v, int size, int free) {
		int w = 0;
		for (int i = 0; i < size; i++) {
			w += graph->w[v[i]];
		}
		if (free != graph->C - w) {
			cout << "erro in freedemand. haltin" << endl;
			return 0;
			//exit(1);
		}
		if (w > graph->C) {
			cout << "erro in w. haltin" << endl;
			return 0;
			//exit(1);
		}
		for (int i = 0; i < size; i++) {
			for (int j = i + 1; j < size; j++) {
				if (v[i] == v[j]) {
					cout << "erro douled vertex in a cluster. haltin" << endl;
					return 0;
					//exit(1);
				}
			}
		}

		return 1;
	}

	/**debugar */
	char check_ALLcurrente() {
		double cost, ac = 0;
		int cont = 0;
		for (int i = 0; i < curr->size; i++) {
			if (!check(curr->cluster[i].v, curr->cluster[i].size, curr->cluster[i].freeDemand))
				return 0;
			if (curr->cluster[i].size > 0)
				cost = graph->cost(curr->cluster[i].v, curr->cluster[i].size);
			else cost = 0;
			ac += cost;
			if (fabs(cost - curr->cluster[i].cost) > EPSLON) {
				cout << "wrong cost  erro" << endl;
				return 0;
			}
			cont += curr->cluster[i].size;
		}
		if (fabs(ac - curr->cost) > EPSLON) {
			cout << "wrong total cost  erro" << endl;
			return 0;
		}

		if (cont != graph->size) {
			cout << "wrong nuber of vertex in solution" << endl;
			return 0;
		}

		for (int i = 0; i < curr->size; i++) {
			for (int j = i + 1; j < curr->size; j++) {
				for (int vi = 0; vi < curr->cluster[i].size; vi++) {
					for (int vj = 0; vj < curr->cluster[j].size; vj++) {
						if (curr->cluster[i].v[vi] == curr->cluster[j].v[vj]) {
							cout << "vertice in two corrente.clusterSol erro" << endl;
							return 0;
						}
					}
				}
			}
		}
		return 1;
	}

	/** particiona o vetor entre menores e maiores que o pivot*/
	int partition(int * v, int size, double * w, double pivot) {
		int i = 0, j = size - 1;
		do {
			while (i < j && w[v[i]] < pivot)i++;
			while (i < j && w[v[j]] > pivot)j--;
			if (i < j)
				swap(v[i++], v[j--]);
		} while (i < j);

		return j;
	}

	/**add vertex v to cluster c, updating the freeDemand and size*/
	void inline addVertex(Cluster &c, int v) {
		if (c.size == 0)
			c.freeDemand = graph->C - graph->w[v];
		else
			c.freeDemand -= graph->w[v];

		c.v[c.size] = v;
		c.size++;

		c.cx = (c.cx * (c.size - 1) + graph->x[v]) / c.size;
		c.cy = (c.cy * (c.size - 1) + graph->y[v]) / c.size;
		//não é necessário calcular o custo, pois este método é utilizado somente na função
		//de construção
	}

	/**add vertex v to cluster c, updating the freeDemand, clusterof and size*/
	void inline addVertex(Solution &sol, int c, int v) {
		sol.clusterOf[v] = c;
		if (sol.cluster[c].size == 0)
			sol.cluster[c].freeDemand = graph->C - graph->w[v];
		else
			sol.cluster[c].freeDemand -= graph->w[v];

		sol.cluster[c].v[sol.cluster[c].size] = v;
		sol.cluster[c].size++;

		sol.cluster[c].cx = (sol.cluster[c].cx * (sol.cluster[c].size - 1) + graph->x[v]) / sol.cluster[c].size;
		sol.cluster[c].cy = (sol.cluster[c].cy * (sol.cluster[c].size - 1) + graph->y[v]) / sol.cluster[c].size;
		//não é necessário calcular o custo, pois este método é utilizado somente na função
		//de construção
	}

	/** realiza o movimento swap*/
	void commiteSwapMove(Move &move) {
		double iniCost;

		//cout << "swapTabu " << move.delta << endl;
		//#ifdef debug // incorreto para swap move
		//        if (corrente.cluster[move.C_1].size <= move.V_1 || corrente.cluster[move.C_2].freeDemand < graph->w[corrente.cluster[move.C_1].v[move.V_1]]) {
		//            cout << "Invalid trans move" << endl;
		//            exit(1);
		//        }
		//#endif
		iniCost = curr->cluster[move.C_1].cost + curr->cluster[move.C_2].cost;

		curr->cluster[move.C_1].flagNew = 1;
		curr->cluster[move.C_2].flagNew = 1;


		//if (move.delta > -EPSLON) {
		addTabu(move.C_1, curr->cluster[move.C_1].v[move.V_1]);
		addTabu(move.C_2, curr->cluster[move.C_2].v[move.V_2]);
		//}
		curr->cluster[move.C_1].freeDemand += graph->w[curr->cluster[move.C_1].v[move.V_1]] - graph->w[curr->cluster[move.C_2].v[move.V_2]];
		curr->cluster[move.C_2].freeDemand -= graph->w[curr->cluster[move.C_1].v[move.V_1]] - graph->w[curr->cluster[move.C_2].v[move.V_2]];

		swap(curr->cluster[move.C_1].v[move.V_1], curr->cluster[move.C_2].v[move.V_2]);
		curr->cluster[move.C_2].cost = graph->cost(curr->cluster[move.C_2].v, curr->cluster[move.C_2].size, curr->cluster[move.C_2].cx, curr->cluster[move.C_2].cy);
		curr->cluster[move.C_1].cost = graph->cost(curr->cluster[move.C_1].v, curr->cluster[move.C_1].size, curr->cluster[move.C_1].cx, curr->cluster[move.C_1].cy);
		curr->cost += curr->cluster[move.C_1].cost + curr->cluster[move.C_2].cost - iniCost;
		//        cout << corrente.cluster[move.C_1].cost + corrente.cluster[move.C_2].cost - iniCost << endl;
		setBox(curr->cluster[move.C_1]);
		setBox(curr->cluster[move.C_2]);
	}

	/**realiza o movimento transposição*/
	void commiteMove(Move &move) {
		double iniCost;

		//cout << "TRANS " << bestMove.C_1 << " " << bestMove.V_1 << " " << bestMove.C_2 << " " << endl;
#ifdef debug
		if (curr->cluster[move.C_1].size <= move.V_1 || curr->cluster[move.C_2].freeDemand < graph->w[curr->cluster[move.C_1].v[move.V_1]]) {
			cout << "Invalid trans move" << endl;
			exit(1);
		}
#endif
		iniCost = curr->cluster[move.C_1].cost + curr->cluster[move.C_2].cost;

		curr->cluster[move.C_1].flagNew = 1;
		curr->cluster[move.C_2].flagNew = 1;



		//if (move.delta > -EPSLON)
		addTabu(move.C_1, curr->cluster[move.C_1].v[move.V_1]);

		curr->cluster[move.C_1].freeDemand += graph->w[curr->cluster[move.C_1].v[move.V_1]];
		curr->cluster[move.C_2].freeDemand -= graph->w[curr->cluster[move.C_1].v[move.V_1]];
		curr->cluster[move.C_2].v[curr->cluster[move.C_2].size] = curr->cluster[move.C_1].v[move.V_1];
		curr->cluster[move.C_1].v[move.V_1] = curr->cluster[move.C_1].v[curr->cluster[move.C_1].size - 1];
		curr->cluster[move.C_1].size--;
		curr->cluster[move.C_2].size++;
		curr->cluster[move.C_2].cost = graph->cost(curr->cluster[move.C_2].v, curr->cluster[move.C_2].size, curr->cluster[move.C_2].cx, curr->cluster[move.C_2].cy);
		if (curr->cluster[move.C_1].size == 0) {
			curr->cluster[move.C_1].cost = 0;
			emptyCt++;
		} else curr->cluster[move.C_1].cost = graph->cost(curr->cluster[move.C_1].v, curr->cluster[move.C_1].size, curr->cluster[move.C_1].cx, curr->cluster[move.C_1].cy);
		if (curr->cluster[move.C_2].size == 1) emptyCt--;

		curr->cost += curr->cluster[move.C_1].cost + curr->cluster[move.C_2].cost - iniCost;

		setBox(curr->cluster[move.C_1]);
		setBox(curr->cluster[move.C_2]);
	}

	/**verifica se o vertice v pode ir para cluster c*/
	char isTabuTrans(int c, int v) {//atenção, este v é o indice global do vertice
		for (int i = 0; i < tenure; i++) {
			if (tabu[i].C == c && tabu[i].V == v) {
				return 1;
			}
		}

		return 0;
	}

	/**verifica se o vertice v1 pode irp ara o cluster c1 e o vertice v2 pode ir
     cluster c2*/
	char isTabuSwap(int c1, int c2, int v1, int v2) {
		for (int i = 0; i < tenure; i++) {
			if ((tabu[i].C == c2 && tabu[i].V == curr->cluster[c1].v[v1]) || (tabu[i].C == c1 && tabu[i].V == curr->cluster[c2].v[v2])) {
				return 1;
			}
		}
		return 0;
	}

	/**adiciona uma regra na lista tabu*/
	void addTabu(int c, int v) {
		tabu[tabuHead].C = c;
		tabu[tabuHead].V = v;
		if (tabuHead < tenure - 1) tabuHead++;
		else tabuHead = 0;
	}

	/**remove a regra mais antiga da lista*/
	void removeTabu() {
		tabu[tabuHead].C = -1;
		tabu[tabuHead].V = -1;
		if (tabuHead < tenure - 1) tabuHead++;
		else tabuHead = 0;
	}

	/**remove todas as regras da lista e marca os cluster relacionados
     como novos*/
	void resetTabu() {
		tabuHead = 0;
		for (int i = 0; i < tenure; i++)
			if (tabu[i].C != -1) {
				curr->cluster[tabu[i].C].flagNew = 1;
				tabu[i].C = tabu[i].V = -1;
			}
	}

	/**inicializa a lista tabu*/
	void initTabu() {
		tabuHead = 0;
		for (int i = 0; i < tenure; i++)
			tabu[i].C = tabu[i].V = -1;
	}

	/**busca um movimento de transposição entre os cluster a e b
     executando o movimento quando encontrado*/
	char tryTrans(int a, int b) {
		double cxA, cyA, cyB, cxB;
		double deltaA, deltaB;
		int vi;
		char flag = 0;
		if (curr->cluster[a].size == 1 && graph->P > 0) return 0;
		for (int i = 0; i < curr->cluster[a].size; i++)
			if (curr->cluster[b].freeDemand >= graph->w[curr->cluster[a].v[i]] && !isTabuTrans(b, curr->cluster[a].v[i])) {
				vi = curr->cluster[a].v[i];
				cxA = curr->cluster[a].cx;
				cyA = curr->cluster[a].cy;
				cxB = curr->cluster[b].cx;
				cyB = curr->cluster[b].cy;

				if (curr->cluster[a].size > 1)
					deltaA = graph->simulateCostRemove(curr->cluster[a].v, curr->cluster[a].size, cxA, cyA, i) - curr->cluster[a].cost;
				else deltaA = -graph->F;
				if (curr->cluster[a].size > 0)
					deltaB = graph->simulateCostAdd(curr->cluster[b].v, curr->cluster[b].size, cxB, cyB, vi) - curr->cluster[b].cost;
				else deltaB = graph->F;

				if (deltaA + deltaB < -EPSLON) {

					if (curr->cluster[a].size > 1 && curr->cluster[a].size > 0) curr->cost += deltaA + deltaB;
					else if (curr->cluster[a].size == 1) curr->cost += deltaB;
					else curr->cost += deltaA;

					curr->cluster[a].v[i] = curr->cluster[a].v[curr->cluster[a].size - 1];
					curr->cluster[a].flagNew = 1;
					curr->cluster[a].size--;
					curr->cluster[a].freeDemand += graph->w[vi];
					if (curr->cluster[a].size > 0) curr->cluster[a].cost += deltaA;
					else curr->cluster[a].cost = 0;
					if (curr->cluster[a].size == 0)emptyCt++;
					curr->cluster[a].cx = cxA;
					curr->cluster[a].cy = cyA;


					curr->cluster[b].v[curr->cluster[b].size] = vi;
					curr->cluster[b].flagNew = 1;
					curr->cluster[b].size++;
					curr->cluster[b].cost += deltaB;
					curr->cluster[b].freeDemand -= graph->w[vi];
					if (curr->cluster[b].size == 1)emptyCt--;
					curr->cluster[b].cx = cxB;
					curr->cluster[b].cy = cyB;

					boxRemove(curr->cluster[a], vi);
					boxAdd(curr->cluster[b], vi);
#ifdef debug
					check(curr->cluster[a].v, curr->cluster[a].size, curr->cluster[a].freeDemand);
					check(curr->cluster[b].v, curr->cluster[b].size, curr->cluster[b].freeDemand);
#endif
					//cout << "trans " << corrente.cost << endl;
					//return 1;
					flag = 1;


					//                    i--;
					//                    i = -1;
					//   continue;

				}
			}
		return flag;
	}

	/**calcula o valor aproximado de um movimento de swap*/
	double inline guessCost(Cluster &a, Cluster &b, int va, int vb) {
		double cax, cay, cbx, cby, d;
		//if (a.size * b.size < 10000) return -1;
		cax = (a.cx * a.size - graph->x[va] + graph->x[vb]) / a.size;
		cay = (a.cy * a.size - graph->y[va] + graph->y[vb]) / a.size;
		cbx = (b.cx * b.size - graph->x[vb] + graph->x[va]) / b.size;
		cby = (b.cy * b.size - graph->y[vb] + graph->y[va]) / b.size;
		d = graph->dist(cax, cay, vb) - graph->dist(b.cx, b.cy, vb) +
				graph->dist(cbx, cby, va) - graph->dist(a.cx, a.cy, va);
		//d -= graph->dist(graph->x[va], graph->y[va], graph->x[vb], graph->y[vb]) / (a.size + b.size);
		//d -= 0.2*(graph->dist(cax, cay, a.cx, a.cy)*(a.size - 1)+graph->dist(cbx, cby, b.cx, b.cy)*(b.size - 1));
		//cout<< graph->dist(cax, cay, a.cx, a.cy)*(a.size - 1)+graph->dist(cbx, cby, b.cx, b.cy)*(b.size - 1) <<" ";
		return d;
	}

	/**busca um movimento de swap entre os cluster a e b
     executando o movimento quando encontrado*/
	char trySwap(int a, int b) {
		double deltaB, deltaA;
		double cxA, cyA, cyB, cxB;
		char flag = 0;
		//double avg = corrente.cluster[a].cost / corrente.cluster[a].size;
		//        if (!curr->cluster[a].box.intercept(curr->cluster[b].box)) return 0;
		for (int i = 0; i < curr->cluster[a].size; i++)
			if (curr->cluster[b].box.has(graph->x[curr->cluster[a].v[i]], graph->y[curr->cluster[a].v[i]]))
				for (int j = 0; j < curr->cluster[b].size; j++)
					if (curr->cluster[a].box.has(graph->x[curr->cluster[b].v[j]], graph->y[curr->cluster[b].v[j]]) &&
							guessCost(curr->cluster[a], curr->cluster[b], curr->cluster[a].v[i], curr->cluster[b].v[j]) < 0 &&
							curr->cluster[a].freeDemand + graph->w[curr->cluster[a].v[i]] >= graph->w[curr->cluster[b].v[j]] &&
							curr->cluster[b].freeDemand + graph->w[curr->cluster[b].v[j]] >= graph->w[curr->cluster[a].v[i]] &&
							//                        graph->dist(corrente.cluster[a].cx, corrente.cluster[a].cy, corrente.cluster[b].v[j]) - graph->dist(corrente.cluster[b].cx, corrente.cluster[b].cy, corrente.cluster[b].v[j]) +
							//                        graph->dist(corrente.cluster[b].cx, corrente.cluster[b].cy, corrente.cluster[a].v[i]) - graph->dist(corrente.cluster[a].cx, corrente.cluster[a].cy, corrente.cluster[a].v[i]) < avg &&
							(graph->w[curr->cluster[a].v[i]] != graph->w[curr->cluster[b].v[j]] || graph->d[curr->cluster[a].v[i]][curr->cluster[b].v[j]] > EPSLON) &&
							!isTabuSwap(a, b, i, j)) {

						cxA = curr->cluster[a].cx;
						cyA = curr->cluster[a].cy;
						deltaA = graph->simulateCostSwap(curr->cluster[a].v, curr->cluster[a].size, cxA, cyA, curr->cluster[b].v[j], i) - curr->cluster[a].cost;

						cyB = curr->cluster[b].cy;
						cxB = curr->cluster[b].cx;
						deltaB = graph->simulateCostSwap(curr->cluster[b].v, curr->cluster[b].size, cxB, cyB, curr->cluster[a].v[i], j) - curr->cluster[b].cost;


						if (deltaA + deltaB < -EPSLON) {
							//                        cout << cA + cB - iniCost << "\t" << guessCost(corrente.cluster[a], corrente.cluster[b], corrente.cluster[a].v[i], corrente.cluster[b].v[j]) << "\t" << corrente.cluster[a].size << "\t" << corrente.cluster[b].size << endl;
							curr->cost += deltaA + deltaB;
							swap(curr->cluster[a].v[i], curr->cluster[b].v[j]);

							boxAdd(curr->cluster[a], curr->cluster[a].v[i]);
							boxRemove(curr->cluster[a], curr->cluster[b].v[j]);
							boxAdd(curr->cluster[b], curr->cluster[b].v[j]);
							boxRemove(curr->cluster[b], curr->cluster[a].v[i]);


							curr->cluster[a].cost += deltaA;
							curr->cluster[b].cost += deltaB;
							flag = curr->cluster[a].flagNew = curr->cluster[b].flagNew = 1;
							curr->cluster[a].freeDemand -= graph->w[curr->cluster[a].v[i]] - graph->w[curr->cluster[b].v[j]];
							curr->cluster[b].freeDemand -= graph->w[curr->cluster[b].v[j]] - graph->w[curr->cluster[a].v[i]];
							curr->cluster[a].cx = cxA;
							curr->cluster[a].cy = cyA;
							curr->cluster[b].cx = cxB;
							curr->cluster[b].cy = cyB;
#ifdef debug
							check(curr->cluster[a].v, curr->cluster[a].size, curr->cluster[a].freeDemand);
							check(curr->cluster[b].v, curr->cluster[b].size, curr->cluster[b].freeDemand);
#endif
							//cout << "swap " << corrente.cost << endl;
							//                        return 1;

							//                                                i = 0;
							//                                                j = -1;


							break;


						}//if
					}


		return flag;
	}

	/**look and commite tabu move*/
	void tabuMove() {
		double cxA, cyA, cyB, cxB;
		double deltaA, deltaB;
		int v;
		Move move;
		move.delta = INFINITY;
		move.C_1 = -1;
		for (int a = 0; a < curr->size; a++)
			if (curr->cluster[a].size > 1 || graph->P == 0) {
				for (int i = 0; i < curr->cluster[a].size; i++) {
					v = curr->cluster[a].v[i];
					cxA = curr->cluster[a].cx;
					cyA = curr->cluster[a].cy;
					if (curr->cluster[a].size > 1) deltaA = graph->simulateCostRemove(curr->cluster[a].v, curr->cluster[a].size, cxA, cyA, i) - curr->cluster[a].cost;
					else deltaA = -graph->F;

					for (int b = 0; b < curr->size; b++) if (a != b) {
						if (curr->cluster[b].freeDemand >= graph->w[curr->cluster[a].v[i]] && !isTabuTrans(b, curr->cluster[a].v[i])) {
							cxB = curr->cluster[b].cx;
							cyB = curr->cluster[b].cy;
							if (curr->cluster[b].size > 0) deltaB = graph->simulateCostAdd(curr->cluster[b].v, curr->cluster[b].size, cxB, cyB, v) - curr->cluster[b].cost;
							else deltaB = graph->F;

							if (deltaA + deltaB < move.delta - EPSLON) {
								move.delta = deltaA + deltaB;
								move.C_1 = a;
								move.C_2 = b;
								move.V_1 = i;
								move.V_2 = -1;
								if (move.delta < EPSLON) {
									commiteMove(move);
									//cout << "fogo" << endl;
									return;
								}

							}

						}
					}
				}
			}

		//swap
		if (move.C_1 == -1)
			for (int a = 0; a < curr->size; a++) if (curr->cluster[a].size > 0)
				for (int b = a + 1; b < curr->size; b++) if (curr->cluster[b].size > 0)
					for (int i = 0; i < curr->cluster[a].size; i++)//if (corrente.cluster[b].freeDemand < graph->w[corrente.cluster[a].v[i]])
						for (int j = 0; j < curr->cluster[b].size; j++)// if (corrente.cluster[a].freeDemand < graph->w[corrente.cluster[b].v[j]])
							if ((graph->w[curr->cluster[a].v[i]] != graph->w[curr->cluster[b].v[j]] || graph->d[curr->cluster[a].v[i]][curr->cluster[b].v[j]] > EPSLON) &&
									guessCost(curr->cluster[a], curr->cluster[b], curr->cluster[a].v[i], curr->cluster[b].v[j]) < move.delta - EPSLON &&
									curr->cluster[a].freeDemand + graph->w[curr->cluster[a].v[i]] >= graph->w[curr->cluster[b].v[j]] &&
									curr->cluster[b].freeDemand + graph->w[curr->cluster[b].v[j]] >= graph->w[curr->cluster[a].v[i]] &&
									!isTabuSwap(a, b, i, j)
							) {

								cxA = curr->cluster[a].cx;
								cyA = curr->cluster[a].cy;
								deltaA = graph->simulateCostSwap(curr->cluster[a].v, curr->cluster[a].size, cxA, cyA, curr->cluster[b].v[j], i) - curr->cluster[a].cost;

								cyB = curr->cluster[b].cy;
								cxB = curr->cluster[b].cx;
								deltaB = graph->simulateCostSwap(curr->cluster[b].v, curr->cluster[b].size, cxB, cyB, curr->cluster[a].v[i], j) - curr->cluster[b].cost;


								if (deltaA + deltaB < move.delta - EPSLON) {
									move.delta = deltaA + deltaB;
									move.C_1 = a;
									move.C_2 = b;
									move.V_1 = i;
									move.V_2 = j;
								}
							}
		if (move.delta == INFINITY) { // tratar isto depois
			cout << "no trans moviments found" << endl;
			exit(1);
		}

		if (move.V_2 == -1)
			commiteMove(move);
		else
			commiteSwapMove(move);
	}

	//    char trySplit(int c) {
	//        //if (corrente.clusterSol[c].size < 2) return 0;
	//        double cxA, cyA, cyB, cxB;
	//        char r;
	//        double cost1, cost2, max = 0, iniCost;
	//        int *v1 = new int[curr->cluster[c].size];
	//        int *v2 = new int[curr->cluster[c].size];
	//        int s1, s2, maxi = 0, maxj = 0, w1, w2;
	//
	//
	//        for (int i = 0; i < curr->cluster[c].size; i++) {
	//            for (int j = i + 1; j < curr->cluster[c].size; j++) {
	//                if (max < graph->d[curr->cluster[c].v[i]][curr->cluster[c].v[j]]) {
	//                    max = graph->d[curr->cluster[c].v[i]][curr->cluster[c].v[j]];
	//                    maxi = i;
	//                    maxj = j;
	//                }
	//            }
	//        }
	//        v1[0] = curr->cluster[c].v[maxi];
	//        v2[0] = curr->cluster[c].v[maxj];
	//        w1 = graph->w[curr->cluster[c].v[maxi]];
	//        w2 = graph->w[curr->cluster[c].v[maxj]];
	//        s1 = s2 = 1;
	//        for (int i = 0; i < curr->cluster[c].size; i++)
	//            if (i != maxi && i != maxj) {
	//                if (graph->d[v1[0]][curr->cluster[c].v[i]] < graph->d[v2[0]][curr->cluster[c].v[i]] && w1 + graph->w[curr->cluster[c].v[i]] <= graph->C) {
	//                    v1[s1] = curr->cluster[c].v[i];
	//                    s1++;
	//                    w1 += graph->w[curr->cluster[c].v[i]];
	//                } else {
	//                    v2[s2] = curr->cluster[c].v[i];
	//                    s2++;
	//                    w2 += graph->w[curr->cluster[c].v[i]];
	//                }
	//            }
	//        cost1 = graph->cost(v1, s1, cxA, cyA);
	//        cost2 = graph->cost(v2, s2, cxB, cyB);
	//        if (cost1 + cost2 + graph->F + EPSLON < curr->cluster[c].cost) {
	//
	//            iniCost = curr->cluster[c].cost;
	//            curr->cluster[c].size = s1;
	//            for (int i = 0; i < s1; i++)
	//                curr->cluster[c].v[i] = v1[i];
	//            curr->cluster[c].freeDemand = graph->C - w1;
	//            curr->cluster[c].cost = cost1;
	//            curr->cluster[c].flagNew = 1;
	//            curr->cluster[c].cx = cxA;
	//            curr->cluster[c].cy = cyA;
	//            for (maxj = 0; maxj < curr->size && curr->cluster[maxj].size != 0; maxj++); //do nothing!
	//            if (maxj == curr->size) curr->size++;
	//            //cout << maxj << " " << corrente.size << " - ";
	//            curr->cluster[maxj].size = s2;
	//
	//            for (int i = 0; i < s2; i++)
	//                curr->cluster[maxj].v[i] = v2[i];
	//            curr->cluster[maxj].freeDemand = graph->C - w2;
	//            curr->cluster[maxj].cost = cost2;
	//            curr->cluster[maxj].flagNew = 1;
	//            curr->cluster[maxj].cx = cxB;
	//            curr->cluster[maxj].cy = cyB;
	//
	//            curr->cost += cost1 + cost2 - iniCost;
	//            r = 1;
	//            //cout << "Split " << cost1 + cost2 - iniCost << " " << s1 << " " << s2 << endl;
	//        }
	//        delete [] v1;
	//        delete [] v2;
	//        return r;
	//    }
	//
	//    char tryMerge(int c1, int c2) {
	//        double iniCost, cost;
	//        char r = 0;
	//        double cx, cy;
	//        iniCost = curr->cluster[c1].cost + curr->cluster[c2].cost;
	//        if (curr->cluster[c1].freeDemand + curr->cluster[c2].freeDemand >= graph->C) {
	//            for (int i = 0; i < curr->cluster[c2].size; i++) {
	//                curr->cluster[c1].v[curr->cluster[c1].size + i] = curr->cluster[c2].v[i];
	//            }
	//            cost = graph->cost(curr->cluster[c1].v, curr->cluster[c1].size + curr->cluster[c2].size, cx, cy);
	//            if (cost + EPSLON < iniCost + graph->F) {
	//                //                    for (int i = 0; i < cluster[c1].size + cluster[c2].size; i++)
	//                //                        cout << cluster[c1].v[i] << " ";
	//                //                    cout << endl;
	//                curr->cluster[c1].size += curr->cluster[c2].size;
	//                curr->cluster[c1].cost = cost;
	//                curr->cluster[c1].flagNew = 1;
	//                curr->cluster[c1].freeDemand -= graph->C - curr->cluster[c2].freeDemand;
	//                curr->cluster[c1].cx = cx;
	//                curr->cluster[c1].cy = cy;
	//                curr->cluster[c2].size = 0;
	//                curr->cluster[c2].cost = 0;
	//                curr->cluster[c2].freeDemand = graph->C;
	//                curr->cluster[c2].flagNew = 1;
	//                emptyCt++;
	//                curr->cost += curr->cluster[c1].cost - iniCost;
	//                //cout << "Merge " << c1 << " " << c2 << " " << corrente.clusterSol[c1].cost - iniCost << endl;
	//            }
	//        }
	//        return r;
	//    }

	char localSearch() {
		char changed;
		char flag = 0;
		char flagNew[graph->size];

		do {
			changed = 0;
			for (int i = 0; i < curr->size; i++) {
				flagNew[i] = curr->cluster[i].flagNew;
				if (curr->cluster[i].flagNew) {
					sortCluster(curr->cluster[i]);
					curr->cluster[i].flagNew = 0;
				}

			}

			for (int i = 0; i < curr->size; i++) {
				for (int j = 0; j < curr->size; j++)
					if ((flagNew[j] || flagNew[i]) && i != j) {
						if (tryTrans(i, j)) {
							changed = 1;

						}
						if (j > i && curr->cluster[i].box.intercept(curr->cluster[j].box) && trySwap(i, j)) {
							changed = 1;

						}
						//  cout << "LS1: " << corrente.cost << endl;
					}
				//cout << "super" << endl;
				if (superMoveStart(i, 30)) {
					changed = 1;

					//                    cout <<i<< "super " << curr->cost << endl;
				}


				//                cout << "LS2: " << i << " " << curr->cost << "  " << globalTimer.getSec() << endl;


			}


			//            if (!changed) {
			//                for (int i = 0; i < corrente.size; i++) {
			//                    //cout << "super " << i << endl;
			//                    if (superMove(i, 0, 20)) {
			//                        changed = 1;
			//                        break;
			//                    }
			//                }
			//            }
			if (changed) flag = 1;
		} while (changed);
		//cout << "LS1: " << " " << corrente.cost << endl;
		return flag;
	}

	/**Ordena os vertices de um cluster pela distancia ao centroid*/
	void sortCluster(Cluster &c) {
		double d[graph->size];
		for (int i = 0; i < c.size; i++) {
			d[c.v[i]] = graph->dist(c.cx, c.cy, c.v[i]);
			//            cout << c.v[i] << endl;
		}
		//        cout << endl;
		qsort_r(c.v, c.size, sizeof (int), comparaDec, d);
	}

	void setBox(Cluster &c) {
		c.box.isNull = 0;
		c.box.maxY = c.box.maxX = -INFINITY;
		c.box.minY = c.box.minX = INFINITY;
		for (int i = 0; i < c.size; i++) {

			//            cout << c.v[i] << endl;
			if (c.box.maxX < graph->x[c.v[i]])c.box.maxX = graph->x[c.v[i]];
			if (c.box.maxY < graph->y[c.v[i]])c.box.maxY = graph->y[c.v[i]];

			if (c.box.minX > graph->x[c.v[i]])c.box.minX = graph->x[c.v[i]];
			if (c.box.minY > graph->y[c.v[i]])c.box.minY = graph->y[c.v[i]];

		}
		//        cout << endl;

	}

	void boxRemove(Cluster &c, int i) {
		if (c.box.maxX == graph->x[i] || c.box.minX == graph->x[i] || c.box.maxY == graph->y[i] || c.box.minY == graph->y[i])
			setBox(c);
	}

	void boxAdd(Cluster &c, int i) {
		if (c.box.maxX < graph->x[i])c.box.maxX = graph->x[i];
		if (c.box.maxY < graph->y[i])c.box.maxY = graph->y[i];
		if (c.box.minX > graph->x[i])c.box.minX = graph->x[i];
		if (c.box.minY > graph->y[i])c.box.minY = graph->y[i];
	}
	//    char preLS() {
	//        char changed, newI;
	//        do {
	//            changed = 0;
	//            for (int i = 0; i < curr->size; i++)
	//                if (curr->cluster[i].size > 0) {
	//                    newI = curr->cluster[i].flagNew;
	//                    curr->cluster[i].flagNew = 0;
	//                    if (newI && trySplit(i)) changed = 1;
	//                    for (int j = 0; j < i; j++)
	//                        if ((newI || curr->cluster[j].flagNew) && dist(curr->cluster[i].cx, curr->cluster[i].cy, curr->cluster[j].cx, curr->cluster[j].cy) < graph->F && curr->cluster[j].size > 0) {
	//                            if (tryMerge(i, j)) changed = 1;
	//                        }
	//                }
	//        } while (changed);
	//    }

	/**variaveis utilizada pelo super move*/
	double global_gain;
	int global_maxCall;

	/**movimento de propragação*/
	virtual int superMove(int c, int level, int maxLevel) {
		if (level >= maxLevel) return 0;
		double max = 0, costC, costI, cxC, cyC, cxI, cyI, delta;
		int maxV, maxArg = -1;
		char flagC;
		char flagI;
		global_maxCall--;
		if (global_maxCall < 0) {
			// cout << "super move maxcalls" << endl;
			return 0;
		}
		for (int i = 0; i < curr->cluster[c].size - 1; i++) {


			if (graph->w[curr->cluster[c].v[i]] >= -curr->cluster[c].freeDemand && max < graph->d[curr->cluster[c].v[curr->cluster[c].size - 1]][curr->cluster[c].v[i]]/**(graph->C - graph->w[corrente.cluster[c].v[i]])*/) {
				//                isolado = 1;
				//                for (int j = 0; j < curr->size; j++) {
				//                    if (c != j && curr->cluster[j].box.has(graph->x[curr->cluster[c].v[i]], graph->y[curr->cluster[c].v[i]])) {
				//                        isolado = 0;
				//                        break;
				//                    }
				//                }
				//                if (!isolado) {
				max = graph->d[curr->cluster[c].v[curr->cluster[c].size - 1]][curr->cluster[c].v[i]]/**(graph->C - graph->w[corrente.cluster[c].v[i]])*/;
				maxArg = i;
				//                }
				//                                else {
				//                                    cout << "isolado" << endl;
				//                                }
			}

		}
		if (maxArg < 0) {
			//acontece quando só quem poder tornar o cluster viavel é o último
			//vertice.
			return 0;
		}
		maxV = curr->cluster[c].v[maxArg];
		cxC = curr->cluster[c].cx;
		cyC = curr->cluster[c].cy;
		if (curr->cluster[c].size == 1) costC = 0;
		else costC = graph->simulateCostRemove(curr->cluster[c].v, curr->cluster[c].size, cxC, cyC, maxArg);

		curr->cluster[c].v[maxArg] = curr->cluster[c].v[curr->cluster[c].size - 1];


		double acuGain = curr->cluster[c].cost - costC - global_gain - EPSLON;
		if (curr->cluster[c].size == 1) acuGain += graph->F;
		for (int i = 0; i < curr->size && global_maxCall > 0; i++)
			if (i != c && curr->cluster[i].size > 0) {
				cxI = curr->cluster[i].cx;
				cyI = curr->cluster[i].cy;
				costI = graph->simulateCostAdd(curr->cluster[i].v, curr->cluster[i].size, cxI, cyI, maxV);
				if (acuGain > costI - curr->cluster[i].cost
						&& !isTabuTrans(i, maxV)) {
					curr->cluster[i].v[curr->cluster[i].size] = maxV;
					if (curr->cluster[i].freeDemand >= graph->w[maxV]) {
						curr->cost += costC + costI - (curr->cluster[i].cost + curr->cluster[c].cost);
						curr->cluster[i].size++;
						curr->cluster[i].flagNew = 1;
						curr->cluster[i].freeDemand -= graph->w[maxV];
						curr->cluster[i].cost = costI;
						curr->cluster[i].cx = cxI;
						curr->cluster[i].cy = cyI;

						curr->cluster[c].size--;
						curr->cluster[c].flagNew = 1;
						curr->cluster[c].freeDemand += graph->w[maxV];
						curr->cluster[c].cost = costC;
						curr->cluster[c].cx = cxC;
						curr->cluster[c].cy = cyC;
						if (curr->cluster[c].size == 0)
							emptyCt++;
#ifdef debug
						check(curr->cluster[c].v, curr->cluster[c].size, curr->cluster[c].freeDemand);
						check(curr->cluster[i].v, curr->cluster[i].size, curr->cluster[i].freeDemand);
#endif
						//                        cout << "SM " << level << endl;
						//cout << c << " " << i << " ";
						setBox(curr->cluster[c]);
						setBox(curr->cluster[i]);
						return 1;
					} else {
						//fake commite
						delta = costC + costI - (curr->cluster[i].cost + curr->cluster[c].cost);
						global_gain += delta;
						curr->cost += delta;

						curr->cluster[c].size--;
						curr->cluster[c].freeDemand += graph->w[maxV];
						swap(curr->cluster[c].cost, costC);
						swap(curr->cluster[c].cx, cxC);
						swap(curr->cluster[c].cy, cyC);
						flagC = curr->cluster[c].flagNew;
						curr->cluster[c].flagNew = 1;

						curr->cluster[i].size++;
						curr->cluster[i].freeDemand -= graph->w[maxV];
						swap(curr->cluster[i].cost, costI);
						swap(curr->cluster[i].cx, cxI);
						swap(curr->cluster[i].cy, cyI);
						flagI = curr->cluster[i].flagNew;
						curr->cluster[i].flagNew = 1;




						if (superMove(i, level + 1, maxLevel)) {
							curr->cluster[c].flagNew = curr->cluster[i].flagNew = 1;
							//cout << c << " ";
							setBox(curr->cluster[c]);
							return 1;
						}

						//undo the fake commite
						curr->cost -= delta;
						global_gain -= delta;

						curr->cluster[c].v[curr->cluster[c].size] = maxV;
						curr->cluster[c].size++;
						curr->cluster[c].freeDemand -= graph->w[maxV];
						swap(curr->cluster[c].cost, costC);
						swap(curr->cluster[c].cx, cxC);
						swap(curr->cluster[c].cy, cyC);
						curr->cluster[c].flagNew = flagC;

						curr->cluster[i].size--;
						curr->cluster[i].freeDemand += graph->w[maxV];
						swap(curr->cluster[i].cost, costI);
						swap(curr->cluster[i].cx, cxI);
						swap(curr->cluster[i].cy, cyI);
						curr->cluster[i].flagNew = flagI;


					}
				}
			}
		//add new cluster, just for gCCCP instances
		if (graph->P == 0 && -global_gain + (curr->cluster[c].cost - costC) > graph->F + EPSLON) {
			curr->cost += costC - curr->cluster[c].cost;
			curr->cluster[c].size--;
			curr->cluster[c].flagNew = 1;
			curr->cluster[c].freeDemand += graph->w[maxV];
			curr->cluster[c].cost = costC;
			curr->cluster[c].cx = cxC;
			curr->cluster[c].cy = cyC;

			curr->cluster[curr->size].cost = 0;
			curr->cluster[curr->size].flagNew = 1;
			curr->cluster[curr->size].cx = graph->x[maxV];
			curr->cluster[curr->size].cy = graph->y[maxV];
			curr->cluster[curr->size].freeDemand = graph->C - graph->w[maxV];
			curr->cluster[curr->size].size = 1;
			curr->cluster[curr->size].v[0] = maxV;
			curr->size++;

			setBox(curr->cluster[c]);
			setBox(curr->cluster[curr->size]);
#ifdef debug
			cout << "New cluster added " << curr->size << " " << curr->cost << endl;
			check_ALLcurrente();
			cout << "ok" << endl;
#endif
			return 1;
		}

		curr->cluster[c].v[curr->cluster[c].size - 1] = curr->cluster[c].v[maxArg];
		curr->cluster[c].v[maxArg] = maxV;

		return 0;
	}

	virtual int superMoveStart(int c, int maxLevel) {
		double costC, costI, cxC, cyC, cxI, cyI, delta;
		int maxV;
		char flagC;
		char flagI;

		global_maxCall = graph->size;
		for (int vi = 0; vi < curr->cluster[c].size; vi++) {
			if (graph->P > 0 && curr->cluster[c].size == 1) return 0;
			global_gain = 0;

			maxV = curr->cluster[c].v[vi];

			cxC = curr->cluster[c].cx;
			cyC = curr->cluster[c].cy;

			if (curr->cluster[c].size == 1) costC = 0;
			else costC = graph->simulateCostRemove(curr->cluster[c].v, curr->cluster[c].size, cxC, cyC, vi);

			curr->cluster[c].v[vi] = curr->cluster[c].v[curr->cluster[c].size - 1];

			double acuGain = curr->cluster[c].cost - costC - EPSLON;
			if (curr->cluster[c].size == 1) acuGain += graph->F;
			for (int i = 0; i < curr->size; i++)
				if (/*curr->cluster[i].box.has(graph->x[maxV],graph->y[maxV]) &&*/ curr->cluster[i].freeDemand < graph->w[maxV] && i != c && curr->cluster[i].size > 0 && !isTabuTrans(i, maxV)) {
					cxI = curr->cluster[i].cx;
					cyI = curr->cluster[i].cy;
					costI = graph->simulateCostAdd(curr->cluster[i].v, curr->cluster[i].size, cxI, cyI, maxV);
					if (acuGain > costI - curr->cluster[i].cost) {
						curr->cluster[i].v[curr->cluster[i].size] = maxV;
						//fake commite
						delta = costC + costI - (curr->cluster[i].cost + curr->cluster[c].cost);
						global_gain += delta;
						curr->cost += delta;

						curr->cluster[c].size--;
						curr->cluster[c].freeDemand += graph->w[maxV];
						swap(curr->cluster[c].cost, costC);
						swap(curr->cluster[c].cx, cxC);
						swap(curr->cluster[c].cy, cyC);
						flagC = curr->cluster[c].flagNew;
						curr->cluster[c].flagNew = 1;

						curr->cluster[i].size++;
						curr->cluster[i].freeDemand -= graph->w[maxV];
						swap(curr->cluster[i].cost, costI);
						swap(curr->cluster[i].cx, cxI);
						swap(curr->cluster[i].cy, cyI);
						flagI = curr->cluster[i].flagNew;
						curr->cluster[i].flagNew = 1;

						if (superMove(i, 1, maxLevel)) {
							curr->cluster[c].flagNew = curr->cluster[i].flagNew = 1;
							//cout << c << " ";
							setBox(curr->cluster[c]);
							return 1;

						}

						//undo the fake commite
						curr->cost -= delta;
						global_gain -= delta;

						curr->cluster[c].v[curr->cluster[c].size] = maxV;
						curr->cluster[c].size++;
						curr->cluster[c].freeDemand -= graph->w[maxV];
						swap(curr->cluster[c].cost, costC);
						swap(curr->cluster[c].cx, cxC);
						swap(curr->cluster[c].cy, cyC);
						curr->cluster[c].flagNew = flagC;

						curr->cluster[i].size--;
						curr->cluster[i].freeDemand += graph->w[maxV];
						swap(curr->cluster[i].cost, costI);
						swap(curr->cluster[i].cx, cxI);
						swap(curr->cluster[i].cy, cyI);
						curr->cluster[i].flagNew = flagI;


					}
				}
			//add new cluster, just for gCCCP instances
			if (graph->P == 0 && curr->cluster[c].cost - costC > graph->F + EPSLON) {
				curr->cost += costC - curr->cluster[c].cost;
				curr->cluster[c].size--;
				curr->cluster[c].flagNew = 1;
				curr->cluster[c].freeDemand += graph->w[maxV];
				curr->cluster[c].cost = costC;
				curr->cluster[c].cx = cxC;
				curr->cluster[c].cy = cyC;

				curr->cluster[curr->size].cost = 0;
				curr->cluster[curr->size].flagNew = 1;
				curr->cluster[curr->size].cx = graph->x[maxV];
				curr->cluster[curr->size].cy = graph->y[maxV];
				curr->cluster[curr->size].freeDemand = graph->C - graph->w[maxV];
				curr->cluster[curr->size].size = 1;
				curr->cluster[curr->size].v[0] = maxV;
				curr->size++;

				setBox(curr->cluster[c]);
				setBox(curr->cluster[curr->size]);
#ifdef debug
				cout << "New cluster added " << curr->size << " " << curr->cost << endl;
				check_ALLcurrente();
				cout << "ok" << endl;
#endif

				return 1;


			}

			curr->cluster[c].v[curr->cluster[c].size - 1] = curr->cluster[c].v[vi];
			curr->cluster[c].v[vi] = maxV;
		}
		return 0;
	}

	/**remove empty clustes*/
	void cleanBestSol() {
		int i = 0, j = best->size - 1;
		do {
			while (best->cluster[i].size > 0)i++;
			while (best->cluster[j].size == 0)j--;
			if (i < j) {
				best->cluster[i].cost = best->cluster[j].cost;
				best->cluster[i].cx = best->cluster[j].cx;
				best->cluster[i].cy = best->cluster[j].cy;
				best->cluster[i].freeDemand = best->cluster[j].freeDemand;
				best->cluster[i].size = best->cluster[j].size;
				for (int k = 0; k < best->cluster[i].size; k++) {
					best->cluster[i].v[k] = best->cluster[j].v[k];
				}
				best->cluster[j].size = best->cluster[j].cost = 0;
			}
		} while (i < j);
		best->size = i;
	}

	/**remove empty clustes*/
	void cleanCurrrent() {
		int i = 0, j = curr->size - 1;
		do {
			while (curr->cluster[i].size > 0)i++;
			while (curr->cluster[j].size == 0)j--;
			if (i < j) {
				curr->cluster[i].flagNew = 1;
				curr->cluster[i].cost = curr->cluster[j].cost;
				curr->cluster[i].cx = curr->cluster[j].cx;
				curr->cluster[i].cy = curr->cluster[j].cy;
				curr->cluster[i].freeDemand = curr->cluster[j].freeDemand;
				curr->cluster[i].size = curr->cluster[j].size;
				for (int k = 0; k < curr->cluster[i].size; k++) {
					curr->cluster[i].v[k] = curr->cluster[j].v[k];
				}
				curr->cluster[j].size = curr->cluster[j].cost = 0;
			}
		} while (i < j);
		curr->size = i;
		emptyCt = 0;
	}


	//    virtual void constructInitialSol() {
	//        timer.start();
	//        double min, max;
	//        int iArg = graph->far1, jArg = graph->far2, countV;
	//        int * clusterOfI;
	//        int * vizinhoArg;
	//        Cluster * clusterI;
	//        graph->calculateDensity();
	//
	//        clusterOfI = clusterOf;
	//        clusterI = corrente.cluster;
	//        for (int i = 0; i < graph->size; i++, clusterOfI++, clusterI++) {
	//            *clusterOfI = -1;
	//            clusterI->freeDemand = graph->C;
	//            clusterI->size = 0;
	//        }
	//
	//        countV = 2;
	//        corrente.size = 2;
	//        addVertex(0, iArg);
	//        addVertex(1, jArg);
	//        vizinhoArg = nearest[iArg];
	//
	//        for (int i = 1; i < graph->size && graph->d[iArg][*vizinhoArg] < graph->F; i++, vizinhoArg++)
	//            if (graph->w[*vizinhoArg] <= corrente.cluster[0].freeDemand && clusterOf[*vizinhoArg] == -1) {
	//                addVertex(0, *vizinhoArg);
	//                //                check(cluster[0].v, cluster[0].size, cluster[0].freeDemand);
	//                countV++;
	//            }
	//        vizinhoArg = nearest[jArg];
	//        for (int i = 1; i < graph->size && graph->d[jArg][*vizinhoArg] < graph->F; i++, vizinhoArg++)
	//            if (graph->w[*vizinhoArg] <= corrente.cluster[1].freeDemand && clusterOf[*vizinhoArg] == -1) {
	//                addVertex(1, *vizinhoArg);
	//                //                check(cluster[1].v, cluster[1].size, cluster[1].freeDemand);
	//                countV++;
	//            }
	//
	//        while (countV < graph->size) {
	//            max = 0;
	//            clusterOfI = clusterOf;
	//            for (int i = 0; i < graph->size; i++, clusterOfI++)
	//                if (*clusterOfI == -1) {
	//                    min = INFINITY;
	//                    vizinhoArg = nearest[i];
	//                    for (int j = 0; j < graph->size; j++, vizinhoArg++) {
	//                        if (clusterOf[*vizinhoArg] != -1) {
	//                            //if (min > graph->d[i][*vizinhoArg])
	//                            min = graph->d[i][*vizinhoArg];
	//                            break;
	//                        }
	//                    }
	//                    if (min > max) {
	//                        max = min;
	//                        iArg = i;
	//                    }
	//                }
	//            addVertex(corrente.size, iArg);
	//            countV++;
	//            vizinhoArg = nearest[iArg];
	//            for (int j = 1; j < graph->size && graph->d[iArg][*vizinhoArg] < graph->F; j++, vizinhoArg++)
	//                if (graph->w[*vizinhoArg] <= corrente.cluster[corrente.size].freeDemand && clusterOf[*vizinhoArg] == -1) {
	//                    addVertex(corrente.size, *vizinhoArg);
	//                    countV++;
	//                }
	//            //            check(cluster[corrente.size].v, cluster[corrente.size].size, cluster[corrente.size].freeDemand);
	//            corrente.size++;
	//        }
	//
	//        updateCosts();
	//        timer.stop();
	//#ifdef debug
	//        check_ALLcorrente.clusterSols();
	//#endif
	//        cout << "Initia solution created in " << timer.time << "s." << endl;
	//    }

	virtual void constructRandomInitialSol(int iniCount) {

		int * next = (int*) malloc(sizeof (int) *graph->size);
		Cluster * clusterI;

		best->cost = INFINITY;

		for (int k = 0; k < iniCount; k++) {
			for (int i = 0; i < graph->size; i++) {
				next[i] = i;
				curr->clusterOf[i] = -1;
			}
			for (int i = 0; i < graph->size; i++) swap(next[i], next[random() % graph->size]);

			clusterI = curr->cluster;
			curr->size = 0;
			for (int i = 0; i < graph->size; i++, clusterI++) {
				clusterI->freeDemand = graph->C;
				clusterI->size = 0;
				clusterI->cx = clusterI->cy = clusterI->cost = 0;
				clusterI->flagNew = 1;
			}


			double min;
			int minArg;
			for (int i = 0; i < graph->size; i++) {
				min = INFINITY;
				minArg = -1;
				for (int c = 0; c < curr->size; c++) {
					if (curr->cluster[c].freeDemand >= graph->w[next[i]] && min > graph->dist(curr->cluster[c].cx, curr->cluster[c].cy, next[i])) {
						min = graph->dist(curr->cluster[c].cx, curr->cluster[c].cy, next[i]);
						minArg = c;
					}
				}
				if (minArg == -1 || min + EPSLON > graph->F) {
					addVertex(curr->cluster[curr->size++], next[i]);
				} else {
					addVertex(curr->cluster[minArg], next[i]);
				}
			}

			updateCosts(curr);

#ifdef debug
			check_ALLcurrente();
#endif
			localSearch();
#ifdef debug
			check_ALLcurrente();
#endif

			cout << "random " << k << " " << curr->cost + curr->size * graph->F
					<< " " << best->cost + best->size * graph->F << "  "
					<< globalTimer.getSec() << "s" << endl;

			if (curr->cost + curr->size * graph->F < best->cost + best->size * graph->F) {
				best->cost = curr->cost;
				best->size = curr->size;
				cout << best->cost << endl;
				for (int i = 0; i < curr->size; i++) {
					best->cluster[i].cost = curr->cluster[i].cost;
					best->cluster[i].cx = curr->cluster[i].cx;
					best->cluster[i].cy = curr->cluster[i].cy;
					best->cluster[i].flagNew = 1;
					best->cluster[i].freeDemand = curr->cluster[i].freeDemand;
					best->cluster[i].size = curr->cluster[i].size;
					for (int j = 0; j < curr->cluster[i].size; j++) best->cluster[i].v[j] = curr->cluster[i].v[j];
				}
			}

		}

		curr->cost = best->cost;
		curr->size = best->size;
		for (int i = 0; i < best->size; i++) {
			curr->cluster[i].cost = best->cluster[i].cost;
			curr->cluster[i].cx = best->cluster[i].cx;
			curr->cluster[i].cy = best->cluster[i].cy;
			curr->cluster[i].flagNew = 0;
			curr->cluster[i].freeDemand = best->cluster[i].freeDemand;
			curr->cluster[i].size = best->cluster[i].size;
			for (int j = 0; j < best->cluster[i].size; j++) curr->cluster[i].v[j] = best->cluster[i].v[j];
		}
		emptyCt = 0;
		updateCosts(curr);
		free(next);

	}


	int* getRelacaoClustersVertices() {
			int* clustersDosVerticesDaMelhorSolucao = new int[graph->size];
			//mapear os clusters para os vertices
			for (int i = 0; i < graph->P; i++) {
				//pra cada cluster que eu tenho
				for (int j = 0; j < best->cluster[i].size; j++) {
					/** O vetor(clustersDosVertices) vai funcionar da seguinte forma:
					 * 		o vértice X do cluster Y ficará na posição X do vetor,
					 * 		e a posição X do vetor, conterá Y, que é o ID do cluster.
					 */
					clustersDosVerticesDaMelhorSolucao[best->cluster[i].v[j]] =
							best->cluster[i].id;
				}
			}

			//for(int i = 0; i < graph->size; i++){

			//}

			return clustersDosVerticesDaMelhorSolucao;
		}



	virtual void buscaTabu(int t, int iniCount, int maxIte) {
		//        FacilityModel * fm = new FacilityModel(graph);
		int count = 0;
		double bestTime;
		int trCount = 1;
		int nextTR = t;
		tenure = t;
		initTabu();


		//constructInitialSol();
		constructRandomInitialSol(iniCount);

		best->size = 0;
		best->cost = curr->cost;
		cout << "tenure " << tenure << endl;
		cout << "ini\t" << setprecision(10) << curr->cost + curr->size * graph->F << "\t" << curr->size << "\t" << globalTimer.getSec() << "s" << endl;

		renewCurr();
		//        preLS();
		//        renewcorrente.cluster();
		//        cout << "preLS\t" << corrente.cost << "\t" << corrente.size << "\t" << corrente.cost + (corrente.size - emptyCt) * graph->F << "\t" << globalTimer.getSec() << "s" << endl;

#ifdef debug
		check_ALLcurrente();
#endif


		best->cost = curr->cost;
		best->size = curr->size;
		for (int i = 0; i < curr->size; i++) {
			best->cluster[i].size = curr->cluster[i].size;
			best->cluster[i].freeDemand = curr->cluster[i].freeDemand;
			best->cluster[i].cost = curr->cluster[i].cost;
			for (int j = 0; j < curr->cluster[i].size; j++) {
				best->cluster[i].v[j] = curr->cluster[i].v[j];
			}
		}
		renewCurr();

		do {
			count++;

			localSearch();
			if (emptyCt > 0)
				cleanCurrrent();

#ifdef debug
			check_ALLcurrente();
#endif
			cout << count << "\t" << setprecision(10) << curr->cost + (curr->size - emptyCt) * graph->F << "\t" << curr->size << "\t" << setprecision(10) << best->cost + best->size * graph->F << "\t" << best->size << "\t" << globalTimer.getSec() << "s" << endl;
			if (curr->cost + curr->size * graph->F + EPSLON < best->cost + best->size * graph->F) {
				resetTabu();
				cleanCurrrent();
				best->cost = curr->cost;
				bestTime = globalTimer.getSec();
				count = 0;
				best->size = curr->size;
				for (int i = 0; i < curr->size; i++) {

					best->cluster[i].size = curr->cluster[i].size;
					best->cluster[i].freeDemand = curr->cluster[i].freeDemand;
					best->cluster[i].cost = curr->cluster[i].cost;
					best->cluster[i].cx = curr->cluster[i].cx;
					best->cluster[i].cy = curr->cluster[i].cy;

					for (int j = 0; j < curr->cluster[i].size; j++) {
						best->cluster[i].v[j] = curr->cluster[i].v[j];
					}
				}
				//                if (fm->search(corrente.cost, corrente.cluster, corrente.size))
				//                    continue;
			}
			//commiteMove(savedMove);
			if (count == nextTR) {
				resetTabu();
				trCount++;
				nextTR = tenure * (((1 + trCount) * trCount) / 2);


				cout << "resetTabu " << trCount << endl;

			} else
				tabuMove();
		} while (count < maxIte || (globalTimer.getSec() <= UM_DIA));

		best->printSolution(graph);

		cout << "time to best " << bestTime << endl;
		cout << "total time " << globalTimer.getSec() << endl;

		cleanBestSol();


	}


};

class TabuSearch_P : public TabuSearch {
public:

	TabuSearch_P(CapacitadedGraph * g) : TabuSearch(g) {
	}

	char bestFit(int * next) {
		Cluster * clusterI;
		clusterI = curr->cluster;
		//reset curr solution
		for (int i = 0; i < graph->P; i++, clusterI++) {
			clusterI->freeDemand = graph->C;
			clusterI->size = 0;
			clusterI->cx = clusterI->cy = clusterI->cost = 0;
			clusterI->flagNew = 1;
		}
		//distribui os pivot
		for (int c = 0; c < graph->P; c++) {
			addVertex(curr->cluster[c], next[c]);
		}
		double min;
		int minArg;

		for (int i = graph->P; i < graph->size; i++) {

			minArg = -1;
			//best fit
			min = INFINITY;
			for (int c = 0; c < graph->P; c++) {
				if (curr->cluster[c].freeDemand >= graph->w[next[i]] && min > graph->d[curr->cluster[c].v[0]][next[i]]) {
					min = graph->d[curr->cluster[c].v[0]][next[i]];
					minArg = c;
					break; //nextfit
				}
			}


			if (minArg == -1) {
				return 0;
			} else {
				addVertex(curr->cluster[minArg], next[i]);
			}
		}
		updateCosts(curr);
		return 1;
	}

	virtual void pathRelinkInitialSol(int starts) {
		cout << "chamou o PR" << endl;
		int * seed = (int*) ((malloc(sizeof(int) * graph->size)));
		int* bestseed = (int*) ((malloc(sizeof(int) * graph->size)));
		int contadorGeracoesSemMelhoria = 0;
		best->cost = INFINITY;
		best->size = curr->size = graph->P;
		for (int i = 0; i < graph->size; i++) {
			seed[i] = i;
		}
		for (int k = 0; k < starts; k++) {
			/*if(contadorGeracoesSemMelhoria == 20 || (globalTimer.getSec() >= UM_DIA)){
			 break;
			 }*/

			if ((globalTimer.getSec() >= UM_DIA)) {
				break;
			}

			for (int i = 0; i < graph->size; i++)
				swap(seed[i], seed[random() % graph->size]);

			if (!bestFit(seed))
				continue;

#ifdef debug
			check_ALLcurrente();
#endif
			localSearch();
#ifdef debug
			check_ALLcurrente();
#endif
			//            for (int i = 0; i < graph->P; i++)
			//                seed[i] = curr->cluster[i].v[curr->cluster[i].size - 1];
			//            for (int i = 0, u = graph->P; i < graph->P; i++) {
			//                for (int j = 0; j < curr->cluster[i].size - 1; j++) {
			//                    seed[u++] = curr->cluster[i].v[j];
			//                }
			//            }

			cout << "random " << k << " " << curr->cost << " " << best->cost
					<< "  " << globalTimer.getSec() << "s" << endl;

			//            exit(0);
			if (curr->cost + EPSLON < best->cost) {
				// k = 0;
				contadorGeracoesSemMelhoria = 0;
				best->cost = curr->cost;
				//                for (int i = 0; i < graph->P; i++)
				//                    bestseed[i] = curr->cluster[i].v[curr->cluster[i].size - 1];
				//                for (int i = 0, k = graph->P; i < graph->P; i++) {
				//                    for (int j = 0; j < curr->cluster[i].size - 1; j++) {
				//                        bestseed[k++] = curr->cluster[i].v[j];
				//                    }
				//                }
				for (int i = 0; i < graph->size; i++) {
					bestseed[i] = seed[i];
				}
				for (int i = 0; i < graph->P; i++) {
					best->cluster[i].cost = curr->cluster[i].cost;
					best->cluster[i].cx = curr->cluster[i].cx;
					best->cluster[i].cy = curr->cluster[i].cy;
					best->cluster[i].flagNew = 1;
					best->cluster[i].freeDemand = curr->cluster[i].freeDemand;
					best->cluster[i].size = curr->cluster[i].size;
					for (int j = 0; j < curr->cluster[i].size; j++)
						best->cluster[i].v[j] = curr->cluster[i].v[j];
				}
			} else { //path relink

				for (int i = 0; i < graph->P - 1; i++) {
					for (int m = i; m < graph->size; m += graph->P)
						if (seed[m] != bestseed[m]) {
							for (int j = 0; j < graph->size; j++) {
								if (seed[j] == bestseed[m]) {
									swap(seed[j], seed[m]);
									break;
								}
							}
						}

					if (!bestFit(seed))
						continue;

					localSearch();

					if (curr->cost + EPSLON < best->cost) {
						contadorGeracoesSemMelhoria = 0;
						cout << "pathRelink " << curr->cost << "  "
								<< globalTimer.getSec() << "s" << endl;

						//k = 0;
						best->cost = curr->cost;

						//                        for (int i = 0; i < graph->P; i++)
						//                            bestseed[i] = curr->cluster[i].v[curr->cluster[i].size - 1];
						//                        for (int i = 0, k = graph->P; i < graph->P; i++) {
						//                            for (int j = 0; j < curr->cluster[i].size - 1; j++) {
						//                                bestseed[k++] = curr->cluster[i].v[j];
						//                            }
						//                        }
						for (int i = 0; i < graph->size; i++) {
							bestseed[i] = seed[i];
						}
						for (int i = 0; i < graph->P; i++) {
							best->cluster[i].cost = curr->cluster[i].cost;
							best->cluster[i].cx = curr->cluster[i].cx;
							best->cluster[i].cy = curr->cluster[i].cy;
							best->cluster[i].flagNew = 1;
							best->cluster[i].freeDemand =
									curr->cluster[i].freeDemand;
							best->cluster[i].size = curr->cluster[i].size;
							for (int j = 0; j < curr->cluster[i].size; j++)
								best->cluster[i].v[j] = curr->cluster[i].v[j];
						}
						break;
					}
				}
			}
			if (curr->cost != best->cost) {
				++contadorGeracoesSemMelhoria;
			}

		}
		best->cost = curr->cost;
		for (int i = 0; i < graph->P; i++) {
			curr->cluster[i].cost = best->cluster[i].cost;
			curr->cluster[i].cx = best->cluster[i].cx;
			curr->cluster[i].cy = best->cluster[i].cy;
			curr->cluster[i].flagNew = 0;
			curr->cluster[i].freeDemand = best->cluster[i].freeDemand;
			curr->cluster[i].size = best->cluster[i].size;
			for (int j = 0; j < best->cluster[i].size; j++)
				curr->cluster[i].v[j] = best->cluster[i].v[j];
		}
		emptyCt = 0;
		graph->F = 0;
		updateCosts(curr);
		free(seed);
		free(bestseed);
	}

	virtual void constructRandomInitialSol(int iniCount) {
		int* next = (int*) ((malloc(sizeof(int) * graph->size)));
		Cluster* clusterI;
		//double cost;
		best->cost = INFINITY;
		curr->size = best->size = graph->P;
		for (int k = 0; k < iniCount; k++) {
			for (int i = 0; i < graph->size; i++) {
				next[i] = i;
				curr->clusterOf[i] = -1;
			}
			for (int i = 0; i < graph->size; i++)
				swap(next[i], next[random() % graph->size]);

			clusterI = curr->cluster;
			for (int i = 0; i < graph->P; i++, clusterI++) {
				clusterI->freeDemand = graph->C;
				clusterI->size = 0;
				clusterI->cx = clusterI->cy = clusterI->cost = 0;
				clusterI->flagNew = 1;
			}

			for (int c = 0; c < graph->P; c++) {
				addVertex(curr->cluster[c], next[c]);
			}
			double min;
			int minArg;
			char fail;
			fail = 0;
			for (int i = graph->P; i < graph->size; i++) {

				minArg = -1;
				//best fit
				min = INFINITY;
				for (int c = 0; c < graph->P; c++) {
					if (curr->cluster[c].freeDemand >= graph->w[next[i]]
							&& min > graph->d[curr->cluster[c].v[0]][next[i]]) {
						min = graph->d[curr->cluster[c].v[0]][next[i]];
						minArg = c;
					}
				}

				if (minArg == -1) {
					fail = 1;
					break;
				} else {
					addVertex(curr->cluster[minArg], next[i]);
				}
			}
			if (fail)
				continue;
			updateCosts(curr);

			localSearch();

			//cout << "random " << k << " " << curr->cost << " " << best->cost << "  " << globalTimer.getSec() << "s" << endl;

			if (curr->cost < best->cost) {
				best->cost = curr->cost;

				for (int i = 0; i < graph->P; i++) {
					best->cluster[i].cost = curr->cluster[i].cost;
					best->cluster[i].cx = curr->cluster[i].cx;
					best->cluster[i].cy = curr->cluster[i].cy;
					best->cluster[i].flagNew = 1;
					best->cluster[i].freeDemand = curr->cluster[i].freeDemand;
					best->cluster[i].size = curr->cluster[i].size;
					for (int j = 0; j < curr->cluster[i].size; j++)
						best->cluster[i].v[j] = curr->cluster[i].v[j];
				}
			} else { //path relink

			}

		}
		best->cost = curr->cost;
		for (int i = 0; i < graph->P; i++) {
			curr->cluster[i].cost = best->cluster[i].cost;
			curr->cluster[i].cx = best->cluster[i].cx;
			curr->cluster[i].cy = best->cluster[i].cy;
			curr->cluster[i].flagNew = 0;
			curr->cluster[i].freeDemand = best->cluster[i].freeDemand;
			curr->cluster[i].size = best->cluster[i].size;
			for (int j = 0; j < best->cluster[i].size; j++)
				curr->cluster[i].v[j] = best->cluster[i].v[j];
		}
		emptyCt = 0;
		graph->F = 0;
		updateCosts(curr);
		free(next);
	}


	void buscaTabu(int t, int iniCount, int maxIte) {
		int count = 0;
		double bestT, iniCost, iniT;
		int trCount = 1;
		int nextTR = t;
		tenure = t;
		initTabu();
		//constructRandomInitialSol(iniCount);
		pathRelinkInitialSol(iniCount);
		iniCost = best->cost = curr->cost;
		iniT = globalTimer.getSec();
		cout << "tenure " << tenure << endl;
		cout << "ini\t" << setprecision(10) << curr->cost << "\t"
				<< globalTimer.getMillSec() << "s" << endl;
		best->cost = curr->cost;
		best->size = curr->size;
		for (int i = 0; i < curr->size; i++) {
			best->cluster[i].size = curr->cluster[i].size;
			best->cluster[i].freeDemand = curr->cluster[i].freeDemand;
			best->cluster[i].cost = curr->cluster[i].cost;
			for (int j = 0; j < curr->cluster[i].size; j++) {
				best->cluster[i].v[j] = curr->cluster[i].v[j];
			}
		}

		//renewcorrente.cluster();
		//        //        //        fastlocalSearch();
		//                renewcorrente.cluster();
		//        updateF();
		//        updateCosts();
		//        cout << "preLS\t" << corrente.cost - graph->P * graph->F << "\t" << graph->F << "\t\t" << timer.getTime() << "s" << endl;
		do {
			count++;
			//            fm->search(corrente.cost, corrente.cluster, corrente.size);
			localSearch();
			//cout << count << "\t" << setprecision(10) << curr->cost << "\t" << setprecision(10) << best->cost << "\t" << globalTimer.getSec() << "s" << endl;
			//return ;
			if (curr->cost + EPSLON < best->cost) {
				trCount = 1;
				nextTR = tenure;
				//printToPs(curr,graph,"teste");

				bestT = globalTimer.getSec();
				resetTabu();
				trCount = 1;
				nextTR = tenure;
#ifdef debug
				check_ALLcurrente();
#endif
				best->cost = curr->cost;
				count = 0;
				best->size = curr->size;
				for (int i = 0; i < curr->size; i++) {

					best->cluster[i].size = curr->cluster[i].size;
					best->cluster[i].freeDemand = curr->cluster[i].freeDemand;
					best->cluster[i].cost = curr->cluster[i].cost;
					best->cluster[i].cx = curr->cluster[i].cx;
					best->cluster[i].cy = curr->cluster[i].cy;
					for (int j = 0; j < curr->cluster[i].size; j++) {
						best->cluster[i].v[j] = curr->cluster[i].v[j];
					}
				}
				//                if (fm->search(corrente.cost, corrente.cluster, corrente.size))
				//                    continue;
			}
			//commiteMove(savedMove);
			//updateF();
			//updateCosts();

			if (count == nextTR) { //|| (globalTimer.getSec() <= (UM_DIA / 2))) {
				resetTabu();
				trCount++;
				nextTR = tenure * (((1 + trCount) * trCount) / 2);

				//cout << "resetTabu " << trCount << endl;

			} else
				tabuMove();

		} while (count < maxIte);
		best->printSolution(graph);
		best->timeToBest = bestT;
		cout << "ini cost     " << iniCost << endl;
		cout << "time init    " << iniT << endl;
		cout << "cost         " << best->cost << endl;
		cout << "time to best " << bestT << endl;
		cout << "time total   " << globalTimer.getSec() << endl;
		//check_ALLcurrente();
		//getRelacaoClustersVertices();
	}
};


#endif	/* TABUSEARCH_HPP */


