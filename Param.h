#pragma once
#include "Client.h"
#include "Location.h"

class Client;
class Location;

class Param {
public:
	std::string nameIns;
	ofstream fileOut;
	int Q;// capacity
	int T;// time horizon	
	int numVeh = int(1e9);// number of vehicle
	int numClient;//number of client (containing depot)	
	int numLoc;// number of node
	vector<Client> listCL;
	vector<Location> listLoc;
	vector<vector<int> > corDis;// correlation measure
	vector<vector<int> > costs;// distance	
	vector<vector<int> > new_costs;// distance	
	vector<vector<int> > times;// distance			
	bool isDebug = false;
	bool isTurnCkSol = false;
	int lambda = 2;	
	clock_t start, end;
	double total = 0;
	int ldTw = 5;// coff for TW
	//int maxNeibor = 40;// max size of neigbor vertex set for granular search (can be used for dynamically sertification).
	int maxNeibor = 120;// max size of neigbor vertex set for granular search (can be used for dynamically sertification).
	int nbR = 20;
	int nbIls = 100;
	int nbF = 20;
	int sizeSub = 10;// max size of subsequence for using concatenation
	int TL = -1;// time limit
	bool bi = true;//cheking best improvement
	int nI = 10;
	int nC = 10;
	int pMin = 1;
	int pMax = 2;
	int debugLS = false;	
	bool m2Opt = false;
	bool crosver = false;
	bool swap2 = false;
	bool moveor = false;
	bool stringex = false;
	string pathOut;
	Param() {
		cout << "clear param\n";
		listCL.clear();
		listLoc.clear();
		costs.clear();		
		times.clear();
	};
};
