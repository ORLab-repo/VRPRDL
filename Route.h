#pragma once
#include "Param.h"
#include "Node.h"

using namespace std;

class Node;

class Route{
public:
    Node* depot;
    //vector<double> C;// cumulative cost       
    int length = 0;// number of client plus depot on route
    //int posInSet;// position in set of route (using for update LS)
    Param* pr;
    Route(){}
    Route(Param* _pr) {
        pr = _pr;        
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