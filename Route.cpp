#include "lib.h"             
#include "Route.h"

double Route::caculateDis() {
    double res = 0.0;
    Node* val = this->depot;
    Node* valSuc = this->depot->suc;
    do
    {
        res += pr->Costs[valSuc->idxClient][val->idxClient];
        val = valSuc;
        valSuc = valSuc->suc;
    } while (val ->idxClient);
    return res;
}

Route::~Route(void) {    
    if (A != nullptr) {
        delete[] A;
        delete[] Q;
        delete[] B;
        delete[] E;
        delete[] Arv;
        delete[] Qrv;
        delete[] Brv;
        delete[] Erv;
    }
}               
void Route::updateRoute() {
    A[0] = 0;
    B[0] = 0;
    E[0] = 0;
    Q[0] = 0;
    length = 0;
    Node* val = depot;
    val->rou = this;
    int u, v;
    for (int i = 1;; ++i) {
        u = val->idxClient;
        val = val->suc;
        val->rou = this;
        v = val->idxClient;
        Q[i] = Q[i - 1] + pr->listCL[v].demand;
        B[i] = B[i - 1] + pr->listCL[u].calEnergy(pr->listCL[v], 0);
        A[i] = A[i - 1] + pr->listCL[u].calA(pr->listCL[v]);
        E[i] = E[i - 1] + Q[i - 1] * (A[i] - A[i - 1]) + (B[i] - B[i - 1]);
        length = i;
        //if(pr->isDebug)cout << v << " " << E[i] << endl;
        if (v == 0)break;
        val->posInRoute = i;
    }
    Arv[length] = 0;
    Brv[length] = 0;
    Erv[length] = 0;
    Qrv[length] = 0;
    val = depot->pred;
    for (int i = length-1; i >= 0; --i) {
        u = val->idxClient;
        val = val->pred;
        v = val->idxClient;
        Qrv[i] = Qrv[i + 1] + pr->listCL[v].demand;
        Brv[i] = Brv[i + 1] + pr->listCL[u].calEnergy(pr->listCL[v], 0);
        Arv[i] = Arv[i + 1] + pr->listCL[u].calA(pr->listCL[v]);
        Erv[i] = Erv[i + 1] + Qrv[i + 1] * (Arv[i] - Arv[i + 1]) + (Brv[i] - Brv[i + 1]);
        if (v == 0)break;
    }
}

double Route::calAInRan(int u, int v) {
    if (u == v)return 0;
    if (u > v)return Arv[v] - Arv[u];
    return A[v] - A[u];
}

int Route::calQInRan(int u, int v) {
    //if (u == v)return pr->listCL[u].demand;
    if (u > v)swap(u, v);
    return Q[v] - Q[max(0, u - 1)];
}

double Route::calBInRan(int u, int v) {
    if (u == v)return 0;
    if (u > v)return Brv[v] - Brv[u];
    return B[v] - B[u];
}

double Route::calEInRan(int u, int v) {
    if (u == v)return 0;
    if (u > v) {//reverse:
        return Erv[v] - Erv[u] - (Qrv[min(length, u + 1)] * calAInRan(u, v));
    }    
    return E[v] - E[u] - (Q[max(0, u - 1)] * calAInRan(u, v));
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

