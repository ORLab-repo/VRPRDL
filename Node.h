#pragma once
#include "Route.h"
#include "SeqData.h"
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
    SeqData* seq0_i;
    SeqData* seqi_0;
    SeqData* seqi_n;
    SeqData* seqn_i;
    vector <SeqData*> seqi_j; // data for (i,j) with j > i
    vector <SeqData*> seqj_i; // data for (j,i) (for the same subsequence as i_j, but reversed)
    int stTime;//stating time
    int enTime;//ending time
    int demand;
    vector<int> moves; //init based on correlation measure
    int idxClient = -1;
    int idxLoc = -1;

    Node() {        
        isDepot = false;               
    }    

    bool ckNearDepot() {        
        return (pred->idxClient == 0 || suc->idxClient == 0);
    }  
    ~Node() {
    };    
};
