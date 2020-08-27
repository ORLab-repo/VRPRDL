#include "lib.h"             
#include "Route.h"

double Route::caculateDis() {
    double res = 0.0;
    Node* val = this->depot;
    Node* valSuc = this->depot->suc;
    do
    {
        res += pr->costs[valSuc->idxClient][val->idxClient];
        val = valSuc;
        valSuc = valSuc->suc;
    } while (val ->idxClient);
    return res;
}

Route::~Route(void) {        
}               
void Route::updateRoute() {    
    length = 0;
    Node* val = depot;
    val->rou = this;
    int u, v;
    for (int i = 1;; ++i) {
        u = val->idxClient;
        val = val->suc;
        val->rou = this;
        v = val->idxClient;        
        length = i;
        //if(pr->isDebug)cout << v << " " << E[i] << endl;
        if (v == 0)break;
        val->posInRoute = i;
    }    
    val = depot->pred;
    for (int i = length-1; i >= 0; --i) {
        u = val->idxClient;
        val = val->pred;
        v = val->idxClient;        
        if (v == 0)break;
    }
}

void Route::showR() {
    Node* val = this->depot;
    do {
        cout << val->idxClient << " ";
        val = val->suc;
    } while (val != this->depot);
    cout << "\n";
}

void Route::showR_rev() {
    Node* val = this->depot->pred;
    do {
        cout << val->idxClient << " ";
        val = val->pred;
    } while (val != this->depot->pred);
    cout << "\n";
}
void Route::reverse() {
    Node* valU = this->depot->pred;
    Node* valV = this->depot->suc;
    do {
        Node* nxt = valV->suc;
        if (valV->idxClient == 0) {
            valV = this->depot;
        }
        valU->pred = valV;
        valV->suc = valU;
        if (valV->idxClient == 0)break;
        valU = valV;
        valV = nxt;
    } while (true);
    //updateRoute();
}

void Route::reverse1Dep() {
    Node* valU = this->depot;
    Node* valV = valU->suc;
    do {
        Node* nxt = valV->suc;
        valU->pred = valV;
        valV->suc = valU;
        if (valV->idxClient == 0)break;
        valU = valV;
        valV = nxt;
    } while (true);
}
void Route::clearRouteFrom(Node* val) {
    if (val->rou != this) {
        throw "Wrong here";
    }    
    val->suc = this->depot->pred;
    this->depot->pred->pred = val;
}

void Route::insertToRou(Node* u) {
    if (u->idxClient == 0)return;
    Node* depotEn = this->depot->pred;
    Node* nodeEn = depotEn->pred;
    nodeEn->suc = u;
    u->pred = nodeEn;
    u->suc = depotEn;
    depotEn->pred = u;
}

void Route::clearRouteFromRev(Node* u) {
    reverse();
    clearRouteFrom(u);
}
void Route::insertToRouPrev(Node* u) {
    if (u->idxClient == 0)return;
    Node* nxtDepot = this->depot->suc;
    this->depot->suc = u;
    u->pred = this->depot;
    u->suc = nxtDepot;
    nxtDepot->pred = u;
}

