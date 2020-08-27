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
	int numClient;//number of client	
	int numLoc;// number of node
	vector<Client> listCL;
	vector<Location> listLoc;
	vector<vector<int> > costs;// distance	
	vector<vector<int> > new_costs;// distance	
	vector<vector<int> > times;// distance			
	bool isDebug = false;
	int lambda = 2;
	clock_t start;
	int TL = -1;// time limit
	bool bi = true;//cheking best improvement
	int nI = 200;
	int nC = 200;
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
