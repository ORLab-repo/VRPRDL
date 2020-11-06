#pragma once
#include "Solution.h"
#include <ctime>

class GA
{
public:
	int n;
	const int nPop = 40;
	const int delta = 100;
	vector<Solution*> pop;// [nPop + 200];
	Solution* valPop;
	vector<int> adapt;// [nPop + 200] ;
	int sumAdapt;
	vector<II> sortAdapt;// [nPop + 200] ;
	int nPop1;
	int numNotCha;
	int ddID[1000];
	Param* pr;
	GA();
	~GA();
	void init(Param* _pr) {
		pr = _pr;
		n = pr->numClient - 1;
		for (int i = 0; i <= nPop + 200; ++i) {
			pop.push_back(new Solution(_pr));
			adapt.push_back(0);
			sortAdapt.push_back(II(0, 0));
		}
		Solution* val = new Solution(_pr);
		valPop = val;
	}
    bool checkIdSol(Solution* u);
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

