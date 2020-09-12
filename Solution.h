#pragma once
#include "lib.h"
#include "Rng.h"
#include "Node.h"
#include "Route.h"

bool notDominatePop(II &label, II &labelA, II &labelB, II &labelC)
{
    if (label.first >= labelA.first && label.second >= labelA.second)
        if (label.first != labelA.first || label.second != labelA.second) return false;
    if (label.first >= labelB.first && label.second >= labelB.second)
        if (label.first != labelB.first || label.second != labelB.second) return false;
    if (label.first >= labelC.first && label.second >= labelC.second)
        if (label.first != labelC.first || label.second != labelC.second) return false;
    return true;
}
bool notDominatePush(II &label, II &labelA, II &labelB, II &labelC)
{
    if (label.first >= labelA.first && label.second >= labelA.second) return false;
    if (label.first >= labelB.first && label.second >= labelB.second) return false;
    if (label.first >= labelC.first && label.second >= labelC.second) return false;
    return true;
}
void updateLabel(II label, II& labelA, II& labelB, II& labelC)
{
    if (label.first <= labelA.first && label.second <= labelA.second) labelA = label;
    if (label.first < labelB.first) labelB = label;
    if (label.second < labelC.second) labelC = label;
}

class Solution {
public:
    vector<int> giantT;// giant tour
    vector<int> solT;// solution tour (contains index of clients)
    vector<Node*> nodes;// for nodes
    vector<Route*> setR;// for set Route      
    vector<SeqData*> myseqs;// for concatenation in LS
    SeqData* seqSet;// init sequences for each node.
    Param* pr;
    int cost;// objective
    int n;// number of customer
    int m;// number of vehicle (use when limit number of vehicle)
    /*
    * using when limiting the num of Veh
    double** F;
    int** pred;
    */
    int* F;
    int* pred;

    Solution(Param* _pr) {
        pr = _pr;
        n = pr->numClient - 1;// not contain depot
        if (pr->numVeh != oo)m = pr->numVeh;
        else m = n;
        m++;//do for initialization
        // create giant tour
        for (int i = 1; i <= n + 1; ++i)giantT.pb(0);//indexing from 1        
        //create node        
        /*
        * indexing from 1
        * each route contains 2 depot         
        */
        for (int i = 1; i <= n + 2 * m + 1; ++i)nodes.pb(new Node());
        //for customer        
        for (int i = 1; i <= n; ++i) {
            nodes[i]->idxClient = i;
            nodes[i]->demand = pr->listCL[i].demand;
        }
        //for depot
        for (int i = n + 1; i <= 2 * m + n; ++i) {
            nodes[i]->idxClient = 0;
            nodes[i]->demand = 0;
            nodes[i]->idxLoc = 0;
        }
        //init neibor list:

        solT.push_back(0);
        set<DI> sDis;
        for (int i = 1; i <= m + n; ++i) {
            solT.push_back(0);
            /*sDis.clear();
            for (int j = 1; j <= m + n; ++j)if (i != j) {
                sDis.insert(DI(pr->costs[nodes[i]->idxClient][nodes[j]->idxClient], j));
            }
            for (auto val : sDis) {
                nodes[i]->moves.push_back(val.second);
            }*/
        }
        solT.push_back(0);

        for (int i = n + 1; i <= m + n; ++i) {
            nodes[i]->pred = nodes[i + m];
            nodes[i]->suc = nodes[i + m];
            nodes[i + m]->suc = nodes[i];
            nodes[i + m]->pred = nodes[i];
        }

        //for route        
        for (int i = 1; i <= m + 1; ++i)setR.pb(new Route(pr));
        for (int i = 1; i <= m; ++i) {
            setR[i]->depot = nodes[i + n];            
        }       
        //for split:
        /* using when num of vehicles is limit.
        F = new double* [m + 1];
        pred = new int* [m + 1];
        for (int k = 0; k <= m; ++k) {
            F[k] = new double[n + 1];
            pred[k] = new int[n + 1];                       
        }
        */
        //for LS
        /*
        * init seqdata for each node
        */    
        int maxSizeSeq = (pr->sizeSub * 2 + 4) * n + 2 * m * 4;
        SeqData* myseqDatas = new SeqData[maxSizeSeq];
        for (int i = 0; i < maxSizeSeq; ++i)myseqDatas[i].pr = _pr;
        seqSet = myseqDatas;
        int posSeq = 0;

        //for client
        for (int i = 1; i <= n; ++i) {
            nodes[i]->seq0_i = &myseqDatas[posSeq];
            nodes[i]->seqi_0 = &myseqDatas[posSeq + 1];
            nodes[i]->seqi_n = &myseqDatas[posSeq + 2];
            nodes[i]->seqn_i = &myseqDatas[posSeq + 3];
            for (int j = 0; j < pr->sizeSub; ++j) {
                nodes[i]->seqi_j.push_back(&myseqDatas[posSeq + 4 + j]);
            }
            for (int j = 0; j < pr->sizeSub; ++j) {
                nodes[i]->seqj_i.push_back(&myseqDatas[posSeq + 4 + pr->sizeSub + j]);
            }
            posSeq += 4 + 2 * pr->sizeSub;
        }
        //for depot:
        for (int i = n + 1; i <= 2 * m + n; ++i) {
            nodes[i]->seq0_i = &myseqDatas[posSeq];
            nodes[i]->seqi_0 = &myseqDatas[posSeq + 1];
            nodes[i]->seqi_n = &myseqDatas[posSeq + 2];
            nodes[i]->seqn_i = &myseqDatas[posSeq + 3];
            posSeq += 4;
        }

        /*
        * init set of moves for each node
        * fixed and flexible version.
        */
        for (int i = 1; i <= n; ++i) {            
            nodes[i]->idxLocMoves.resize(n + 1, false);
            nodes[i]->movesLoc.resize(pr->maxNeibor);
            nodes[i]->idxCluMoves.resize(n + 1, false);            
            nodes[i]->movesClu.resize(pr->maxNeibor);                        
        }

        /*
        * init for split
        */
        F = new int[n + 1];
        pred = new int[n + 1];
    }

    void genGiantT() {
        for (int i = 1; i <= n; ++i)giantT[i] = i;
        random_shuffle(giantT.begin() + 1, giantT.end());
    }

    void cvGiantT() {
        int pos = 0;
        for (int i = 1; i <= m; ++i) {
            Node* val = setR[i]->depot;
            do
            {
                if (val->idxClient != 0)giantT[++pos] = val->idxClient;
                val = val->suc;
            } while (val != setR[i]->depot);
        }
        if (pr->isDebug) {
            int* dd = new int[n + 1];
            for (int i = 1; i <= n; ++i)dd[i] = 0;
            for (int i = 1; i <= n; ++i)dd[giantT[i]] = 1;
            for (int i = 1; i <= n; ++i)if (dd[i] == 0) {
                throw "Wrong here";
            }
            if (pos != n) {
                throw "Wrong here";
            }
            delete[] dd;
        }
    }

    
    bool ckSol() {
        vector<int> arrSol;
        int* dd = new int[n + 1];
        for (int i = 1; i <= n; ++i)dd[i] = 0;
        int totalCost = 0;
        for (int i = 1; i <= m; ++i) {
            setR[i]->showR();
            Node* val = setR[i]->depot;
            setR[i]->ckRoute();
            while (true)
            {
                val = val->suc;
                if (val->idxClient == 0)break;
                dd[val->idxClient] = 1;                
            }
            totalCost += setR[i]->depot->pred->seq0_i->cost;
        }
        for (int i = 1; i <= n; ++i)if (dd[i] == 0) {
            throw "Wrong here";
        }
        assert(totalCost == cost);
        delete [] dd;
    }    

    //void cvSolT() {
    //    if (pr->isDebug) {
    //        cout << "status local search: " << endl;
    //        if (pr->m2Opt)cout << "2OPT\n";
    //        if (pr->crosver)cout << "CrossOver\n";
    //        if (pr->swap2)cout << "Swap2Node\n";
    //        if (pr->stringex)cout << "StringExchange\n";
    //        cout << endl;
    //    }
    //    int pos = 0;     
    //    double testCost = 0.0;
    //    double testE = 0.0;
    //    if (pr->debugLS) {
    //        testCost = 0.0;
    //        testE = 0.0;
    //        for (int i = 1; i <= m; ++i) {
    //            testCost += setR[i]->caculateDis();
    //            testE = setR[i]->calEInRan(0, setR[i]->length);
    //            if (testE > pr->maxE) {
    //                cout << testE << endl;
    //                throw "Wrong here";
    //            }
    //        }            
    //        if (abs(testCost - cost) >= EP) {
    //            throw "Wrong here";
    //        }
    //    }
    //    if (pr->isDebug) {
    //        testCost = 0.0;
    //        testE = 0.0;            
    //        cout << "solution tour:\n";
    //        for (auto val : solT) cout << val << " ";
    //        cout << endl;
    //        cout << endl;
    //        cout << "detail route: \n";
    //        for (int i = 1; i <= m; ++i) {
    //            setR[i]->updateRoute();
    //            if (pr->isDebug) {
    //                setR[i]->showR();
    //                setR[i]->showR_rev();
    //            }
    //            testCost += setR[i]->caculateDis();
    //            testE = setR[i]->calEInRan(0, setR[i]->length);
    //            if (pr->isDebug) cout << testE << endl;
    //            if (testE > pr->maxE) {
    //                //cout << testE << endl;
    //                throw "Wrong here";
    //            }
    //            //setR[i]->showR_rev();
    //        }           
    //        
    //        if (testCost - cost >= EP) {
    //            if (pr->isDebug) {
    //                //cout << testE << endl;
    //                cout << testCost << " "<< cost<<endl;
    //            }
    //            throw "Wrong here";

    //        }
    //    }

    //    for (int i = 1; i <= m; ++i) {
    //        solT[++pos] = i + n;
    //        Node* val = nodes[i + n];
    //        nodes[i + n]->posInSol = pos;
    //        do {
    //            val = val->suc;
    //            if (val->idxClient == 0)break;
    //            solT[++pos] = val->idxClient;
    //            val->posInSol = pos;
    //        } while (1);
    //    }
    //    //debug:
    //    if (pr->debugLS) {
    //        if (!ckSolT()) {
    //            throw "Wrong here";
    //        }
    //    }
    //    if (pr->isDebug) {
    //        cout << "after transfer: \n";
    //        for (auto val : solT)cout << val << " ";
    //        cout << "\n\n";            
    //        if (!ckSolT()) {
    //            throw "Wrong here";
    //        }
    //    }
    //    if (pos != m + n) {
    //        /*for (auto val : giantT)cout << val << ", ";
    //        exit(0);*/
    //        throw "Wrong size";
    //    }        
    //    solT[0] = solT[m + n];
    //    solT[m + n + 1] = solT[1];        
    //}    

    //reset all moves in route r
    void reinitSingleMoveInRou(Route* r) {
        for (int i = 1; i <= n; ++i)r->isNodeTested[i] = false;
        for (Node* tempNode = r->depot->suc; tempNode->idxClient != 0; tempNode = tempNode->suc) {
            for (int i = 1; i <= min(m + 1, pr->numVeh); ++i)setR[i]->isNodeTested[tempNode->idxClient] = false;
        }
    }

    //reset all moves on all routes
    void reinitAllMovesInRou() {
        for (int i = 1; i <= min(m + 1, pr->numVeh); ++i) {
            for (int j = 1; j <= n; ++j)setR[i]->isNodeTested[j] = false;
        }
    }

    
    //reset neigborhood cluster set (run when changing location)
    void reinitCluSet(Node* nod) {
        for (auto val : nod->movesClu) {
            nod->idxCluMoves[val] = false;
        }        
        nod->movesClu = nod->rou->pr->listLoc[nod->idxLoc].moves;
    }

    //reset neigborhood set for each node
    void reinitNegiborSet() {
        set<II> sDis;
        int u, v;
        for (int i = 1; i <= n; ++i) {
            u = nodes[i]->idxLoc;
            //cluster set:
            reinitCluSet(nodes[i]);
            //location set: 
            /// (can update on the fly by using set but trading off when processing)
            //reset all existing index
            for (auto val : nodes[i]->movesLoc) {
                nodes[i]->idxLocMoves[pr->listLoc[val].idxClient] = false;
            }

            //checking based on the corelation measure:
            for (int j = 1; j <= n; ++j) if(j!=i){
                v = nodes[j]->idxLoc;
                sDis.insert(II(pr->corDis[u][v], v));
            }
            //adding new index:
            nodes[i]->movesLoc.clear();
            for (auto val : sDis) {
                nodes[i]->movesLoc.push_back(pr->listLoc[val.sc].idxClient);
                nodes[i]->idxLocMoves[nodes[i]->movesLoc.back()] = true;
                if (nodes[i]->movesLoc.size() == pr->maxNeibor)break;
            }
        }
        sDis.clear();
    }
    // split without limit number of vehicles
    void Split() {                
        /*
        * (cost, time)
        * idx of location
        */
        cout << m << "\n";
        if (m != pr->numVeh + 1 && m != n + 1)reinitAllMovesInRou();
        else m--;                    
        vector<III> lstLabel; //((cost, time), loc)
        vector<III> curLabel; //((cost, time), prv)
        vector<int> prvIdLb; //previous label in vector
        //the index that cannot exceed in a route due to the capacity.
        int* maxIdx = new int[n + 1];
        int curLoad = 0;
        maxIdx[0] = 0;       
        for (int i = 1; i <= n; ++i) {
            maxIdx[i] = max(maxIdx[i - 1], i);
            int prvIdx = giantT[i - 1];
            curLoad -= pr->listCL[prvIdx].demand;            
            if (maxIdx[i - 1] < i)curLoad = pr->listCL[giantT[i]].demand;            
            if (maxIdx[i - 1] == n){
                maxIdx[i] = n;
                continue;
            }
            for (int j = max(i + 1, maxIdx[i - 1] + 1); j <= n; ++j) {
                if (curLoad + pr->listCL[giantT[j]].demand <= pr->Q) {
                    curLoad += pr->listCL[giantT[j]].demand;
                    maxIdx[i] = j;
                }
                else {
                    break;
                }
            }            
        }        
        //init split phase:
        for (int i = 1; i <= n; ++i) {
            F[i] = oo;
            pred[i] = -1;
        }        
        F[0] = 0;
        pred[0] = -1;
        int stLb = -1, enLb = -1;
        III lbU, lbV;
        II label1, label2, label3;
        int idCus, costU, timeU, idLocU;
        int timeV, costV;
        for (int st = 1; st <= n; st++) {          
            //init first label:            
            lstLabel.pb(III(II(0, 0), 0));
            prvIdLb.pb(pred[st - 1]);
            stLb = lstLabel.size() - 1;
            enLb = lstLabel.size() - 1;
            for (int v = st; v <= maxIdx[st]; ++v) {
                int idCus = giantT[v];                
                for (auto idLocV : pr->listCL[idCus].listLoc) {
                    label1 = II(oo, oo);
                    label2 = II(oo, 0);
                    label3 = II(0, oo);
                    curLabel.clear();
                    if (pr->costs[0][idLocV] == oo) {//eliminated node:
                        continue;
                    }
                    //construct label of next layer
                    for (int i = stLb; i <= enLb; ++i) {
                        III lbU = lstLabel[i];
                        costU = lbU.ft.ft;
                        timeU = lbU.ft.sc;
                        idLocU = lbU.sc;
                        timeV = max(timeU + pr->times[idLocU][idLocV], pr->listLoc[idLocV].stTime);                        
                        if (timeV > pr->listLoc[idLocV].enTime
                            || timeV + pr->times[idLocV][0] > pr->T // only true when travel times satisfy triangle inequality 
                            ) {
                            continue;
                        }
                        costV = costU + pr->costs[idLocU][idLocV];
                        lbV = III(II(costV, timeV), i);
                        if (notDominatePush(lbV.first, label1, label2, label3)) {
                            updateLabel(lbV.first, label1, label2, label3);
                            curLabel.pb(lbV);
                        }                        
                    }                    
                    //dominate label
                    sort(curLabel.begin(), curLabel.end());                    
                    int minT = oo;
                    //for (auto val : curLabel) {
                    for(int i=0;i<curLabel.size();++i){
                        III val = curLabel[i];
                        if (val.ft.sc < minT)
                        {
                            minT = val.ft.sc;
                            lstLabel.pb(III(val.ft, idLocV));
                            prvIdLb.pb(val.sc);
                            //update F function
                            if (//val.ft.sc + pr->times[idLocV][0] < pr->T // uncomment when travel times does not satisfy triangle inequality
                                F[v] > val.ft.ft + pr->costs[idLocV][0] + F[st - 1]) {
                                F[v] = val.ft.ft + pr->costs[idLocV][0] + F[st - 1];
                                pred[v] = lstLabel.size() - 1;                                
                            }                            
                        }                        
                    }
                }
                stLb = enLb + 1;
                enLb = lstLabel.size() - 1;
                if (stLb > enLb)break;
            }
        }        
        //cout << F[n] << "\n";
        cost = F[n];
        if (cost == oo)return;
        int indexLb = pred[n];// index of last label.
        ///construct solution
        vector<int> tourLoc;        
        int numVeh = 0;
        while (indexLb != -1)
        {                        
            tourLoc.pb(lstLabel[indexLb].second);
            if (lstLabel[indexLb].second == 0)numVeh++;
            indexLb = prvIdLb[indexLb];
        }
        m = numVeh;
        assert(tourLoc.size() == m + n);
        reverse(tourLoc.begin(), tourLoc.end());        
        numVeh = 0;
        for (auto val : tourLoc) {
            if (val == 0) {
                numVeh++;                
            }
            else {                
                idCus = pr->listLoc[val].idxClient;
                nodes[idCus]->idxClient = idCus;
                nodes[idCus]->idxLoc = val;
                setR[numVeh]->insertToRou(nodes[idCus]);
            }
        }                
        for (int i = 1; i <= m; ++i) {            
            setR[i]->updateRoute();
            //setR[i]->showR();
        }
        reinitNegiborSet();
        //ckSol();
        //cvSolT(); // uncomment when need tracking the specific postion in solution (used in sequential search)
        lstLabel.clear();
        curLabel.clear();
        prvIdLb.clear();
        tourLoc.clear();
        delete[] maxIdx;
    }
   
//    void initSol() {
//        vector<Route*> setOfR;
//        setOfR.push_back(nullptr);
//        Node** setDepot = new Node * [2*n + 1];
//        for (int i = 1; i <= n; ++i) {
//            setDepot[i] = new Node(0);
//            setDepot[i+n] = new Node(0);            
//            setDepot[i]->suc = setDepot[i + n];
//            setDepot[i + n]->pred = setDepot[i];
//            setDepot[i + n]->suc = setDepot[i];
//            setDepot[i]->pred = setDepot[i + n];
//        }
//        for (int i = 1; i <= n; ++i) {
//            setOfR.push_back(new Route(pr));
//            //init depot:
//            setOfR[i]->depot = setDepot[i];
//            //set direct:                                    
//            setOfR[i]->insertToRou(nodes[i]);
//            setOfR[i]->updateRoute();
//        }
//
//        //symmetric ?? 
//        vector<DII> saving;
//        for (int i = 1; i <= n; ++i)
//            for (int j = i + 1; j <= n; ++j) {
//                saving.push_back(DII(pr->costs[i][0] + pr->costs[j][0] - pr->costs[i][j], II(i, j)));
//                //saving.push_back(DII(pr->costs[i][0] + pr->costs[j][0] - pr->costs[i][j], II(j, i)));
//            }
//        sort(saving.begin(), saving.end(), greater<>());
//        int u, v;
//        double eUV, eVU;
//        for (auto val : saving) {
//            u = val.sc.ft;
//            v = val.sc.sc;
//            if (nodes[u]->rou == nodes[v]->rou)continue;
//            if (!nodes[u]->ckNearDepot() || !nodes[v]->ckNearDepot())continue;
//            eUV = ckEnerCWH(u, v);
//            eVU = ckEnerCWH(v, u);
//            if(min(eUV, eVU) <= pr->maxE){                
//                if (eUV < eVU)ckEnerCWH(u, v, true);
//                else ckEnerCWH(v, u, true);
//            }
//        }
//
//        int pos = 0;
//        for (int i = 1; i <= n; ++i)
//        {
//            Route* val = setOfR[i];            
//            if (val->depot == nullptr)continue;            
//            val->updateRoute();
//            Node* test = val->depot;
//            do {
//                test = test->suc;
//                if (test->idxClient == 0)break;
//                giantT[++pos] = test->idxClient;// construc giant tour
//            } while (1);
//            //cout << val->E[val->length]<<endl;
//        }
//        delete[] setDepot;
//        for (int i = 1; i <= n; ++i)delete setOfR[i];
//        Split();        
//    }
//    
    ///Local Search:

    
    ///ELSALGO:
    void exchange() {
        int posU = Rng::getNumInRan(1, n);
        int posV = Rng::getNumInRan(1, n);
        while (posU == posV)
        {
            posV = Rng::getNumInRan(1, n);
        }
        swap(giantT[posU], giantT[posV]);
    }

    void interchange() {
        int posU = Rng::getNumInRan(1, n);
        int posV = Rng::getNumInRan(1, n);
        while (posU == posV)
        {
            posV = Rng::getNumInRan(1, n);
        }
        if (posU < posV) {
            for (int i = posU + 1; i <= posV; ++i) {
                swap(giantT[i], giantT[i - 1]);
            }
        }
        else {
            for (int i = posU - 1; i >= posV; --i) {
                swap(giantT[i], giantT[i + 1]);
            }
        }
    }

    void mutate(int numP) {
        for (int i = 1; i <= numP; ++i) {
            if (Rng::generator() % 2 == 0) {
                exchange();
            }
            else {
                interchange();
            }
        }
    }        
    void ELS() {
        //initSol();
        int* resGiantT = new int[n + 1];        
        int* curGiantT = new int[n + 1];
        for (int i = 1; i <= n; ++i) {
            curGiantT[i] = giantT[i];
            resGiantT[i] = giantT[i];
            //cout << giantT[i] << " " << endl;
        }        
        double bestObj=cost;
        cout << cost << endl;
        double curObj;
        int curP = pr->pMin;
        for (int i = 1; i <= pr->nI; ++i) {
            curObj = bestObj;
            for (int j = 1; j <= pr->nC; ++j) {                
                for (int i1 = 1; i1 <= n; ++i1)giantT[i1] = curGiantT[i1];
                mutate(curP);                
                Split();                
                if (cost == oo)continue;                
                try {
                    //updateObj();
                    continue;
                }
                catch (...) {
                    cout << "bug here\n";
                    for (auto val : giantT)cout << val << ", ";
                    cout << endl;
                    system("pause");
                    exit(0);
                }
                if (cost + EP < curObj) {
                    curObj = cost;
                    for (int i1 = 1; i1 <= n; ++i1)resGiantT[i1] = giantT[i1];
                    curP = pr->pMin;
                }
                else {
                    curP = min(pr->pMax, curP + 1);
                }

                /*if ((clock() - pr->start) / CLOCKS_PER_SEC >= pr->TL) {
                    if (curObj + EP < bestObj) {
                        bestObj = curObj;
                        goto thispro;
                    }
                }*/
                //cout << i << " " << j << "\n";
            } 
            //cout << i << endl;
            if (curObj + EP < bestObj) {
                bestObj = curObj;
                for (int i1 = 1; i1 <= n; ++i1)curGiantT[i1] = resGiantT[i1];
            }
        }        
        //thispro:
        cout << "best found obj: " << Util::round2num(bestObj) << endl;
        pr->fileOut<< "best found obj: " << Util::round2num(bestObj) << endl;
        //if (pr->debugLS) {            
            for (int i = 1; i <= n; ++i)giantT[i] = resGiantT[i];
            Split();            
            if (cost != oo && cost - bestObj > EP) {
                cout << "bug here\n";
                system("pause");
                exit(0);
            }
            cout << Util::round2num(cost) << endl;
            pr->fileOut << "Giant Tour:\n";
            for (int i = 1; i <= n; ++i)pr->fileOut << giantT[i] << ", ";
            pr->fileOut << "\n";
            pr->fileOut<<"cost split: " << Util::round2num(cost) << endl;
        //}
        delete[] resGiantT;
        delete[] curGiantT;
    }

    //deconstructor:
    ~Solution(){        
        /*for (int k = 0; k <= m; ++k) {
            delete[] F[k];
            delete[] pred[k];
        }*/
        delete[] seqSet;
        delete[] F;
        delete[] pred;
        giantT.clear();
        solT.clear();                
        for (int i = 0; i < n + 2 * pr->numVeh + 1; ++i)delete nodes[i];
        for (int i = 1; i <= m; ++i)delete setR[i];
        delete pr;
    }
};
