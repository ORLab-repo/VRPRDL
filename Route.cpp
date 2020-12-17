#include "lib.h"             
#include "Route.h"

int Route::caculateDis() {
    int res = 0;
    Node* val = this->depot;
    Node* valSuc = this->depot->suc;
    do
    {
        res += pr->costs[val->idxLoc][valSuc->idxLoc];
        val = valSuc;
        valSuc = valSuc->suc;
    } while (val ->idxLoc);
    return res;
}

Route::~Route(void) {        
    isNodeTested.clear();
}               
void Route::updateRoute() {    
    length = 0;
    Node* val = depot;
    val->posInRoute = 0;
    int ckstg = 0;
    val->seq0_i->init(0);
    val->seqi_0->init(0);
    val->rou = this;
    int u, v;
    //update seqdata (0->i, i->0):
    for (int i = 1;; ++i) {
        u = val->idxLoc;
        val = val->suc;
        val->rou = this;
        v = val->idxLoc;        
        if(pr->isDebug)cout << v << " " << val->demand <<" "<<pr->listLoc[v].demand<<"\n";
        val->seq0_i->concatOneAfter(val->pred->seq0_i, v);
        val->seqi_0->concatOneBefore(val->pred->seqi_0, v);
        length = i;       
        //if(pr->isDebug)cout << v << " " << E[i] << endl;
        val->posInRoute = i;
        if (v == 0)break;        
    }    
    //update seqdata(i->n, n->i):
    val = depot->pred;
    if (!(val->idxLoc == 0)) {
        throw "bug update route";
    }
    val->seqi_n->init(0);
    val->seqn_i->init(0);
    while (true)
    {
        val = val->pred;
        v = val->idxLoc;
        val->seqn_i->concatOneAfter(val->suc->seqn_i, v);
        val->seqi_n->concatOneBefore(val->suc->seqi_n, v);
        if (v == 0)break;
    }

    //update seqdata(i->j)(j->i) 
    //|i-j| < pr->sizeSub (i and j is not a depot).    
    Node* curNode = depot->suc;
    val = curNode;
    while (true)
    {
        if (curNode->idxLoc == 0)break;
        curNode->seqi_j[0]->init(curNode->idxLoc);
        curNode->seqj_i[0]->init(curNode->idxLoc);
        val = curNode->suc;
        for (int i = 1; i < pr->sizeSub; ++i) {            
            if (val->idxLoc == 0)break;
            curNode->seqi_j[i]->concatOneAfter(curNode->seqi_j[i - 1], val->idxLoc);
            curNode->seqj_i[i]->concatOneBefore(curNode->seqj_i[i - 1], val->idxLoc);
            val = val->suc;
        }
        curNode = curNode->suc;
    }
    //reset nodes tested:
    for (int i = 1; i < pr->numClient; ++i)isNodeTested[i] = false;
}

int Route::getCliInRou(int* arr, int* arrLoc)
{
    Node* val = this->depot;
    int numCus = 0;
    do {
        if (val->idxClient) {
            arr[++numCus] = val->idxClient;            
            arrLoc[numCus] = val->idxLoc;
        }
        val = val->suc;
    } while (val != this->depot);
    return numCus;
}

void Route::showR() {
    Node* val = this->depot;
    do {
        cout << val->idxClient << " ";
        val = val->suc;
    } while (val != this->depot);
    cout << "\n";
}

void Route::showRInFile(ostream& os)
{
    Node* val = this->depot;
    do {
        os << val->idxClient << " ";
        val = val->suc;
    } while (val != this->depot);
    os << "\n";
}

void Route::showRLoc()
{
    Node* val = this->depot;
    do {
        cout << val->idxLoc << " ";
        val = val->suc;
    } while (val != this->depot);
    cout << "\n";
}

void Route::showRLocInFile(ostream& os)
{
    Node* val = this->depot;
    do {
        os << val->idxLoc << " ";
        val = val->suc;
    } while (val != this->depot);
    os << "\n";
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

void Route::clearNode()
{
    depot->suc = depot->pred;
    depot->pred->pred = depot;
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

void Route::ckRoute()
{
    //assert(depot->pred->seq0_i->F == true);    
    int cost = 0;
    int time = 0;
    int load = 0;
    Node* val = this->depot;
    Node* valSuc = this->depot->suc;
    int u, v;
    do
    {
        u = val->idxLoc;
        v = valSuc->idxLoc;
        cost += pr->costs[u][v];        
        time += pr->times[u][v];        
        time = max(time, pr->listLoc[v].stTime);        
        if (!(time <= pr->listLoc[v].enTime)) {
            throw "bug TW ckRoute";
        }
        val = valSuc;
        load += val->demand;
        if (!(load <= pr->Q)) {
            throw "bug capacity";
        }
        valSuc = valSuc->suc;
    } while (val->idxLoc);
    //if (pr->debugLS)cout << cost << "\n";
    if (!(depot->pred->seq0_i->cost == depot->seqi_n->cost)) {
        throw "bug cost compare 1";
    }
    if (!(cost == depot->pred->seq0_i->cost)) {
        throw "bug cost compare 2";
    }
    if (!(load == depot->pred->seq0_i->load)) {
        throw "bug load compare";
    }
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
