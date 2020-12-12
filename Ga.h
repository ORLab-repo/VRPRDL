#pragma once
#include "Solution.h"
#include <ctime>

class GA
{
public:
	int n;
	int bestCost;
	const int nPop = 25;
	//const int delta = 100;
	const int delta = 83;
	vector<Solution*> pop;// [nPop + 200];
	Solution* valPop;
	vector<int> adapt;// [nPop + 200] ;
	int sumAdapt;
	vector<II> sortAdapt;// [nPop + 200] ;
	///SCP
	//map with hash value and index in route pool
	vector<map<II, int> > idWithLen; // [num of customers] 
	vector<int> lenR;
	vector<vector<int> > idRouBelong;
	vector<int*> routePool; 
	vector<int*> routePoolLoc;
	vector<int*> routePoolPrv;// previous location of node
	vector<int*> routePoolNxt;// next location of node
	vector<int> costR;// cost of rou
	int* valRou;
	int* valRouLoc;
	int valLength;
	int maxNumRou = 600;
	int curNumRou;
	/// 
	int nPop1;
	int numNotCha;
	int ddID[1000];
	int omega = 0;
	int threshold;	
	double pM = 0.2;
	Param* pr;
	GA();
	~GA();
	void init(Param* _pr) {
		pr = _pr;
		n = pr->numClient - 1;
		//init for SCP
		valRou = new int[n + 1];
		valRouLoc = new int[n + 1];
		idWithLen.resize(n + 1);
		idRouBelong.resize(n + 1);
		routePool.clear();
		routePoolLoc.clear();
		routePoolPrv.clear();
		routePoolNxt.clear();
		costR.clear();
		lenR.clear();
		for (int i = 1; i <= maxNumRou + n + 1; ++i) {
			routePool.push_back(new int[n + 1]);
			routePoolLoc.push_back(new int[n + 1]);
			routePoolPrv.push_back(new int[n + 1]);
			routePoolNxt.push_back(new int[n + 1]);
			costR.push_back(0);
			lenR.push_back(0);
		}		

		for (int i = 0; i <= nPop + 2*delta; ++i) {
			pop.push_back(new Solution(_pr));
			adapt.push_back(0);
			sortAdapt.push_back(II(0, 0));
		}
		Solution* val = new Solution(_pr);
		valPop = val;
	}
	///SCP FUNC
	void addRou(Solution* u);
	int solveSCP();
	///

    bool checkIdSol(Solution* u);
	int broken_pairs(Solution* u, Solution* v);
    bool CheckEqual(Solution* u, Solution* v);
    void equalSol(Solution* u, Solution* v);    
    void FindAdapt();    
    int getChild();    
    void choose(int& u, int& v);    
    void uni(Solution* u, Solution* v, Solution* u1, Solution* v1);
    void insertNew(Solution* u);
    void DelPopu();
	void InitPopu();
    void DiversifyPopu(Solution* bestSol);        
	void findGasSol(int maxNumGas = 5000);

private:
};

