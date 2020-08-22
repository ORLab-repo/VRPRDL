#pragma once
#include "Route.h"
#include <vector>

//class Route;
using namespace std;

class Node{
public:
    bool isDepot;
    int posInRoute;
    int posInSol;
    Route* rou;
    Node* pred;
    Node* suc;
    int stTime;//stating time
    int enTime;//ending time
    int demand;
    vector<int> moves; //init based on correlation measure
    int idxClient = -1;
    int idxLoc = -1;

    Node() {        
        isDepot = false;               
    }

    Node(int idx,int idl) {       
        idxClient = idx;
        idxLoc = idl;
        moves.clear();
    }

    bool ckNearDepot() {        
        return (pred->idxClient == 0 || suc->idxClient == 0);
    }  
    ~Node() {};    
};
