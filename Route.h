#pragma once
#include "Node.h"

using namespace std;

class Node;

class Route{
public:
    Node* depot;
    //vector<double> C;// cumulative cost
    double* A = nullptr;// cumulative constant
    double* Arv;// cumulative constant
    double* E;// cumulative energy
    double* Erv;// cumulative energy
    double* B;// cumulative constant
    double* Brv;// cumulative constant
    int* Q;// cumulative load
    int* Qrv;// cumulative load
    int length = 0;// number of client plus depot on route
    //int posInSet;// position in set of route (using for update LS)
    Param* pr;
    Route(){}
    Route(Param* _pr) {
        pr = _pr;
        A = new double[_pr->numClient+3];
        Arv = new double[_pr->numClient + 3];
        B = new double[_pr->numClient+3];
        Brv = new double[_pr->numClient + 3];
        E = new double[_pr->numClient + 3];
        Erv = new double[_pr->numClient + 3];
        Q = new int[_pr->numClient + 3]; 
        Qrv = new int[_pr->numClient + 3];
    }

    void insertToRou(Node* u);
    void insertToRouPrev(Node* u);
    void clearRouteFrom(Node* u);
    void clearRouteFromRev(Node* u);
    double caculateDis();
    void updateRoute();
    double calAInRan(int u, int v);
    int calQInRan(int u, int v);
    double calBInRan(int u, int v);
    double calEInRan(int u, int v);
    void showR();
    void showR_rev();
    void reverse1Dep();
    void reverse();
    ~Route();
};