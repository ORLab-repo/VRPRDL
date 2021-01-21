#pragma once
#include "Route.h"
#include "SeqData.h"
#include <vector>

using namespace std;

class Route;
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
    vector<bool> idxLocMoves; //true if index of customer is in the movesLoc (false if not exist).
    vector<bool> idxCluMoves; //true if index of customer is in the movesClu (false if not exist).
    int idxClient = -1;
    int idxLoc = -1;

    Node() {        
        isDepot = false;           
        rou = nullptr;
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
