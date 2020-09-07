#pragma once
#include "Route.h"
#include "SeqData.h"
#include <vector>

class Route;
using namespace std;

class Node{
public:
    bool isDepot;
    int posInRoute;
    int posInSol;
    Route* rou;
    Node* pred;
    Node* suc;    
    SeqData* seq0_i;
    SeqData* seqi_0;
    SeqData* seqi_n;
    SeqData* seqn_i;
    vector <SeqData*> seqi_j; // data for (i,j) with j > i
    vector <SeqData*> seqj_i; // data for (j,i) (for the same subsequence as i_j, but reversed)    
    int demand;
    vector<int> movesClu; //init based on correlation measure for cluster
    vector<int> movesLoc; //init based on correlation measure for locations (still contain the index of customer)
    vector<bool> idxLocMoves; //index of locations that correspond to the index of customer (-1 if not exist).
    vector<bool> idxCluMoves; //index of locations that correspond to the index of customer (-1 if not exist).
    int idxClient = -1;
    int idxLoc = -1;

    Node() {        
        isDepot = false;               
    }    

    bool ckNearDepot() {        
        return (pred->idxClient == 0 || suc->idxClient == 0);
    }  
    ~Node() {
        movesClu.clear();
        movesLoc.clear();
        idxLocMoves.clear();
    };    
};
