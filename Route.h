#pragma once
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
    vector<bool> isNodeTested;// true if all move with this node are tested on this route
    Route(){}
    Route(Param* _pr) {
        pr = _pr;   
        for (int i = 1; i <= pr->numClient; ++i)isNodeTested.push_back(false);
    }
    
    void ckRoute();
    void insertToRou(Node* u);    
    void insertToRouPrev(Node* u);
    void clearRouteFrom(Node* u);
    void clearRouteFromRev(Node* u);
    int caculateDis();
    void updateRoute();    
    void showR();
    void showRLoc();
    void showR_rev();
    void reverse1Dep();
    void reverse();
    void clearNode();// remove all client in route
    ~Route();
};