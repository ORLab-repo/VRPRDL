#pragma once
#include "Param.h"

/*
* Based on paper "A unified solution framework for multi-attribute vehicle routing problems" 
*/

class SeqData
{
public:
	//attributes:
	int firstnode;//index location
	int lastnode;//index location
	int cost;// will set to oo if infeasible sub-sequence
	int load;
	int E;//earliest time completion
	int L;//latest starting time
	int T;//sum of travel and service time
	bool F;//check feasibility;
	Param* pr;
	
	//constructor
	SeqData() {};
	SeqData(Param* _pr) {
		pr = _pr;
	};
	~SeqData() {};

	//method:
	/*
	* Construct sequence containing only one node
	*/
	void init(int idxLoc) {
		firstnode = idxLoc;
		lastnode = idxLoc;
		T = 0;// changed to service time if have
		E = pr->listLoc[idxLoc].stTime;// plus service time if have.
		L = pr->listLoc[idxLoc].enTime;
		cost = 0;
		load = pr->listLoc[idxLoc].demand;
		F = true;
	}

	/**/
	void concatOneAfter(SeqData* seq, int idxLoc) {
		F = seq->F;
		
		if (seq->E + pr->times[seq->lastnode][idxLoc] > pr->listLoc[idxLoc].enTime) {
			F = false;
		}
		if (seq->load + pr->listLoc[idxLoc].demand > pr->Q) {
			F = false;
		}
		cost = seq->cost + pr->costs[seq->lastnode][idxLoc];		
		load = seq->load + pr->listLoc[idxLoc].demand;
		E = max(seq->E + pr->times[seq->lastnode][idxLoc], pr->listLoc[idxLoc].enTime);
		L = min(seq->L, pr->listLoc[idxLoc].enTime - pr->times[seq->lastnode][idxLoc] - seq->T);
		T = seq->T + pr->times[seq->lastnode][idxLoc];
		lastnode = idxLoc;
		firstnode = seq->firstnode;
	}
	/**/
	void concatOneBefore(SeqData* seq, int idxLoc) {
		F = seq->F;

		if (pr->listLoc[idxLoc].stTime + pr->times[idxLoc][seq->firstnode] > seq->L) {
			F = false;
		}
		if (seq->load + pr->listLoc[idxLoc].demand > pr->Q) {
			F = false;
		}
		cost = seq->cost + pr->costs[idxLoc][seq->lastnode];
		load = seq->load + pr->listLoc[idxLoc].demand;
		E = max(pr->listLoc[idxLoc].stTime + pr->times[idxLoc][seq->firstnode] + seq->T, seq->E);
		L = min(pr->listLoc[idxLoc].enTime, seq->L - pr->times[idxLoc][seq->firstnode]);
		T = seq->T + pr->times[idxLoc][seq->firstnode];
		firstnode = idxLoc;
		lastnode = seq->lastnode;
	}
	/**/
	int evaluation(vector<SeqData*> seqs) {
		int costR = seqs[0]->cost;
		int loadR = seqs[0]->load;
		bool totalF = seqs[0]->F;
		int totalE = seqs[0]->E;
		int totalL = seqs[0]->L;
		int totalT = seqs[0]->T;
		int u = -1, v = -1;
		for (int i = 0; i < seqs.size() - 1; ++i) {
			u = seqs[i]->lastnode;
			v = seqs[i + 1]->firstnode;
			costR += (pr->costs[u][v] + seqs[i + 1]->cost);
			loadR += seqs[i + 1]->load;
			totalF = (totalF ^ seqs[i]->F) ^ (totalE + pr->times[u][v] <= seqs[i + 1]->L);
			totalE = max(totalE + pr->times[u][v] + seqs[i + 1]->T, seqs[i + 1]->E);
			totalL = max(totalL, seqs[i + 1]->L - pr->times[u][v] - totalT);
			totalT += pr->times[u][v] + seqs[i + 1]->T;			
		}		
		if (!totalF)return oo;
		return costR;
	}
private:
};