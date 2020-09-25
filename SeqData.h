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
	int F;//check feasibility;
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
		if (seq == NULL) init(idxLoc);
		F = seq->F;
		
		int t_12 = pr->times[seq->lastnode][idxLoc];
		if (seq->E + t_12 > pr->listLoc[idxLoc].enTime) {
			F = false;
		}
		if (seq->load + pr->listLoc[idxLoc].demand > pr->Q) {
			F = false;
		}
		cost = seq->cost + pr->costs[seq->lastnode][idxLoc];		
		load = seq->load + pr->listLoc[idxLoc].demand;
		E = max(seq->E + t_12, pr->listLoc[idxLoc].stTime);
		L = min(seq->L, pr->listLoc[idxLoc].enTime - t_12 - seq->T);
		T = seq->T + t_12;
		lastnode = idxLoc;
		firstnode = seq->firstnode;
	}
	/**/
	void concatOneBefore(SeqData* seq, int idxLoc) {
		if (seq == NULL) init(idxLoc);
		F = seq->F;

		int t_12 = pr->times[idxLoc][seq->firstnode];
		if (pr->listLoc[idxLoc].stTime + t_12> seq->L) {
			F = false;
		}
		if (seq->load + pr->listLoc[idxLoc].demand > pr->Q) {
			F = false;
		}
		cost = seq->cost + pr->costs[idxLoc][seq->firstnode];
		load = seq->load + pr->listLoc[idxLoc].demand;
		E = max(pr->listLoc[idxLoc].stTime + t_12 + seq->T, seq->E);
		L = min(pr->listLoc[idxLoc].enTime, seq->L - t_12);
		T = seq->T + t_12;
		firstnode = idxLoc;
		lastnode = seq->lastnode;
	}
	/**/
	int evaluation(vector<SeqData*> seqs) {		
		if (seqs.front() == NULL)seqs.erase(seqs.begin());//remove null sequence
		if (seqs.back() == NULL)seqs.pop_back();//remove null sequence
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
			totalF = (totalF & seqs[i]->F) & (totalE + pr->times[u][v] <= seqs[i + 1]->L) & (loadR <= pr->Q) &(costR < oo);
			if (!totalF)return oo;// only use for checking feasible solution
			totalE = max(totalE + pr->times[u][v] + seqs[i + 1]->T, seqs[i + 1]->E);
			totalL = max(totalL, seqs[i + 1]->L - pr->times[u][v] - totalT);
			totalT += pr->times[u][v] + seqs[i + 1]->T;			
		}		
		return costR;
	}	

private:
};