#pragma once
#include "lib.h"
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
        // create giant tour
        for (int i = 1; i <= n + 1; ++i)giantT.pb(0);//indexing from 1        
        //create node        
        /*
        * indexing from 1
        * each route contains 2 depot         
        */
        for (int i = 1; i <= n + 2 * m + 1; ++i)nodes.pb(new Node());
        //for customer        
        for (int i = 1; i <= n; ++i)nodes[i]->idxClient = i;
        //for depot
        for (int i = n + 1; i <= 2 * m + n; ++i)nodes[i]->idxClient = 0;        
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

        for(int i=n+1; i<= m+n; ++i){
            nodes[i]->pred = nodes[i+m];
            nodes[i+m]->suc = nodes[i];
        }

        //for route        
        for (int i = 1; i <= m; ++i)setR.pb(new Route(pr));
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

    /*
    bool ckSolT() {
        int* dd = new int[n + 1];
        for (int i = 1; i <= n; ++i)dd[i] = 0;
        for (auto val : solT) {
            if(val <= n)dd[pr->listLoc[val]] = 1;
        }
        for (int i = 1; i <= n; ++i)if (dd[i] == 0) {
            throw "Wrong here";
        }
        delete [] dd;
    }
    */

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

    // split without limit number of vehicles
    void Split() {                
        /*
        * (cost, time)
        * idx of location
        */        
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
        int indexLb = pred[n];
        ///construct soluton
        while (true)
        {
            
        }
        int st, en = n;
        //reset:
        for (int i = 1; i <= m; ++i) {
            setR[i]->depot = nodes[i + n];
            nodes[i + n]->pred = nodes[i + n + m];
            nodes[i + n + m]->suc = nodes[i + n];
            nodes[i + n]->suc = nodes[i + n + m];
            nodes[i + n + m]->pred = nodes[i + n];
        }

        for (int i = numVeh; i >= 1; --i) {
            st = pred[i][en] + 1;
            for (int j = st; j <= en; ++j) {
                setR[i]->insertToRou(nodes[giantT[j]]);
            }
            en = st - 1;
        }
        for (int i = 1; i <= m; ++i) {            
            setR[i]->updateRoute();
            //setR[i]->showR();
        }
        cvSolT();
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
//    int randPos() {
//        return rand() % n + 1;
//    }
//    void exchange() {
//        int posU = randPos();
//        int posV = randPos();
//        while (posU == posV)
//        {
//            posV = randPos();
//        }
//        swap(giantT[posU], giantT[posV]);
//    }
//
//    void interchange() {
//        int posU = randPos();
//        int posV = randPos();
//        while (posU == posV)
//        {
//            posV = randPos();
//        }
//        if (posU < posV) {
//            for (int i = posU + 1; i <= posV; ++i) {
//                swap(giantT[i], giantT[i - 1]);
//            }
//        }
//        else {
//            for (int i = posU - 1; i >= posV; --i) {
//                swap(giantT[i], giantT[i + 1]);
//            }
//        }
//    }
//
//    void mutate(int numP) {
//        for (int i = 1; i <= numP; ++i) {
//            if (rand() % 2 == 0) {
//                exchange();
//            }
//            else {
//                interchange();
//            }
//        }
//    }
//
//    //Local Search:
//    bool localSearch() {
//        bool flag = true;
//        while (flag) {
//            //if(pr->debugLS)cout << cost << endl;
//            flag = Move2Opt(pr->bi) 
//                || CrossOver(pr->bi) 
//                || Swap2Node(pr->bi) 
//                || OrOptMove(pr->lambda, pr->bi) 
//                || StringExchange(pr->lambda, pr->bi);
//        }
//        //if(!flag)cout << caculateE() << endl;
//        cvGiantT();
//        return flag;
//    }
//
//    void updateObj() {
//        /*double oldObj = oo;
//        double ckCost = oldObj + 1;
//        while (true)
//        {            
//            Split();
//            oldObj = cost;
//            if (oldObj > ckCost + EP) {                
//                cout << oldObj << " " << ckCost << endl;
//                throw "Wrong here";
//            }                       
//            localSearch();                        
//            if (abs(oldObj - cost) <= EP)break;
//            ckCost = cost;
//        }
//        cvGiantT();*/
//        localSearch();
//    }
//
//    double calC(Node* u, Node* v) {
//        return pr->costs[u->idxClient][v->idxClient];
//    }
//
//    double calA(Node* u, Node* v) {
//        return pr->listCL[u->idxClient].calA(pr->listCL[v->idxClient]);
//    }
//
//    double calB(Node* u, Node* v) {
//        return pr->listCL[u->idxClient].calEnergy(pr->listCL[v->idxClient], 0);
//    }
//
//    //2OPT:
//    bool Move2Opt(bool bi) {
//        if (pr->isDebug) {
//            cout << "In 2 Opt\n";
//            pr->m2Opt = true;
//        }
//        int* resPos = new int[5];// best move
//        int* pos = new int[5];// pos in solT
//        int* posR = new int[5];// pos in route
//        double bestG = 0.0, curGain;
//        vector<II> setTours(2, II(0, 0));
//        vector<II> bestTours(2, II(0, 0));
//        //double posi1, posi2, posi3, posi4;
//        int i1, i2, ap, t1, t2, t3, t4, idVec;
//        double c12, c23;
//        Route* r1;
//        Route* r2;
//        double sumE, sumL;
//        for (i1 = 1; i1 <= m + n; ++i1) 
//            for (ap = -1; ap <= 1; ap += 2) {
//                t1 = solT[i1];
//                t2 = solT[i1 + ap];
//                c12 = calC(nodes[t1], nodes[t2]);
//                for (idVec = 0; idVec< nodes[t2]->moves.size();++idVec) {
//                    t3 = nodes[t2]->moves[idVec];
//                    c23 = calC(nodes[t2], nodes[t3]);
//                    if (c23 >= c12 - bestG / 2 - EP)break;
//                    i2 = nodes[t3]->posInSol;
//                    t4 = solT[i2 - ap];
//                    curGain = c12 - c23
//                        + calC(nodes[t3], nodes[t4])
//                        -calC(nodes[t4], nodes[t1]);
//                    if (curGain - EP <= bestG)continue;
//                    //config order:
//                    pos[1] = i1;
//                    pos[2] = i1 + ap;
//                    pos[3] = i2 - ap;
//                    pos[4] = i2;                                        
//                               
//                    if (ap == -1) {
//                        swap(pos[1], pos[2]);
//                        swap(pos[3], pos[4]);
//                    }                    
//                    //convert to index nodes
//                    for (int i = 1; i <= 4; ++i) pos[i] = solT[pos[i]];                    
//                    if (pos[1] == pos[2]
//                        || pos[2] == pos[3]
//                        || pos[3] == pos[4])continue;                   
//                    r1 = nodes[pos[1]]->rou;
//                    r2 = nodes[pos[3]]->rou;
//                    //dealing with special case intra-route (ie. one of successor is first depot)
//                    if(r1 == r2 && nodes[pos[1]]->posInRoute > nodes[pos[3]]->posInRoute) {
//                        swap(pos[1], pos[3]);
//                        swap(pos[2], pos[4]);
//                    }
//                    //create pos in route:
//                    for (int i = 1; i <= 4; ++i)posR[i] = nodes[pos[i]]->posInRoute;
//                    if (pos[2] > n)posR[2] = r1->length;
//                    if (pos[4] > n)posR[4] = r2->length;                   
//                    //same route:
//                    if(r1 == r2){                        
//                        //capacity constrain still satisfy
//                        //so we only need to check energy constraint:                       
//                        
//                        //segment 1 : depot->pos1->pos3->pos2                        
//                        //double sumE1;
//                        sumE = r1->calEInRan(0, posR[1])
//                            + r1->calEInRan(posR[3], posR[2])
//                            + r1->calQInRan(0, posR[1]) * (
//                                r1->calAInRan(posR[3], posR[2])
//                                + calA(nodes[pos[1]], nodes[pos[3]])
//                                )
//                            + calB(nodes[pos[1]], nodes[pos[3]]);
//                        //sumE1 = sumE;
//                        sumL = r1->calQInRan(0, posR[3]);// sum load                                                              
//                        //segment2:
//                        sumE = sumE + r1->calEInRan(posR[4], r1->length)
//                            + sumL * (r1->calAInRan(posR[4], r1->length) + calA(nodes[pos[2]], nodes[pos[4]]) )
//                            +calB(nodes[pos[2]], nodes[pos[4]]);                        
//                        if (sumE > pr->maxE)continue;                                           
//                        bestG = curGain;
//                        for (int i = 1; i <= 4; ++i)resPos[i] = pos[i];  
//                        if (!bi)goto thispro;
//                    }
//                    //different route:
//                    else{                                       
//                        //cal load and energy of new route 1:
//                        sumL = r1->calQInRan(0, posR[1]) + r2->calQInRan(posR[3], 0);
//                        sumE = r1->calEInRan(0, posR[1]) + r2->calEInRan(posR[3], 0)
//                            + r1->calQInRan(0, posR[1]) * (r2->calAInRan(posR[3], 0))
//                            + r1->calQInRan(0, posR[1]) * calA(nodes[pos[1]], nodes[pos[3]])
//                            + calB(nodes[pos[1]], nodes[pos[3]]);                         
//                        if (sumE > pr->maxE || sumL > pr->cap)continue;                        
//                        //cal load and energy of new route 2:
//                        sumL = r1->calQInRan(r1->length, posR[2]) + r2->calQInRan(posR[4], r2->length);
//                        sumE = r1->calEInRan(r1->length, posR[2]) + r2->calEInRan(posR[4], r2->length)
//                            + r1->calQInRan(r1->length, posR[2]) * r2->calAInRan(posR[4], r2->length)
//                            + r1->calQInRan(r1->length, posR[2]) * calA(nodes[pos[2]], nodes[pos[4]])
//                            + calB(nodes[pos[2]], nodes[pos[4]]);                                                
//                        if (sumE > pr->maxE || sumL > pr->cap)continue;                         
//                        bestG = curGain;                        
//                        for (int i = 1; i <= 4; ++i)resPos[i] = pos[i]; 
//                        if (!bi)goto thispro;
//                    }
//                }
//        }
//        //update:        
//        thispro:
//        if (bestG != 0.0) {
//            if (pr->isDebug) {
//                cout << "improve 2 opt\n";
//                for (auto val : solT)cout << val << " ";
//                cout << endl;
//                for (int i = 1; i <= 4; ++i)cout << resPos[i] << " ";
//                cout << endl; cout << endl;
//            }
//            cost -= bestG;          
//            for (int i = 1; i <= 4; ++i) pos[i] = resPos[i];
//            Route* r1 = nodes[pos[1]]->rou;
//            Route* r2 = nodes[pos[3]]->rou;
//            //create pos in route:            
//            Node* n1 = nodes[pos[1]];
//            Node* n2 = nodes[pos[2]];
//            Node* n3 = nodes[pos[3]];
//            Node* n4 = nodes[pos[4]];
//            if (pos[2] > n)n2 = r1->depot->pred;
//            if (pos[4] > n)n4 = r2->depot->pred;
//            //same route:            
//            if (r1 == r2) {
//                n1->suc = n3;
//                n3->pred = n1;
//                //reverse:
//                Node* valU = n2;
//                Node* valV = n2->suc;
//                do {
//                    Node* nxt = valV->suc;
//                    valU->pred = valV;
//                    valV->suc = valU;
//                    if (valV == n3)break;
//                    valU = valV;
//                    valV = nxt;
//                } while (true);
//                //
//                n2->suc = n4;
//                n4->pred = n2;
//                r1->updateRoute();
//                cvSolT();
//            }
//            //different route:
//            else {               
//                //update route 1:
//                n1->suc = r1->depot->pred;
//                r1->depot->pred->pred = n1;
//                Node* val = n3;
//                Node* cur;
//                while (val->idxClient)
//                {
//                    cur = val->pred;
//                    r1->insertToRou(val);
//                    val = cur;
//                }                
//                //update route2:
//                n4->pred = r2->depot;
//                r2->depot->suc = n4;
//                val = n2;
//                while (val->idxClient)
//                {
//                    cur = val->suc;
//                    r2->insertToRouPrev(val);
//                    val = cur;
//                }
//                r1->updateRoute();
//                //r1->showR();
//                r2->updateRoute();
//                //r2->showR();
//                cvSolT();
//            }
//        }
//        delete[] resPos;
//        delete[] pos;
//        delete[] posR;        
//        pr->m2Opt = false;
//        return (bestG != 0.0);
//    }
//
//    //2OPT*:
//    bool CrossOver(bool bi) {
//        if (pr->isDebug) {
//            cout << "in crossover\n";
//            pr->crosver = true;
//        }
//        int* resPos = new int[5];// best move
//        int* pos = new int[5];// pos in solT
//        int* posR = new int[5];// pos in route
//        double bestG = 0.0;
//        int i1, i2, idVec;
//        int t1, t2, t3, t4;
//        double c12, c23;
//        double sumE;
//        int sumL;
//        double curGain;
//        //double posi1, posi2, posi3, posi4;
//        for (i1 = 1; i1 <= m + n; ++i1) {
//            t1 = solT[i1];
//            t2 = solT[i1 + 1];
//            c12 = calC(nodes[t1], nodes[t2]);
//            for (idVec = 0; idVec < nodes[t2]->moves.size(); ++idVec) {
//                t3 = nodes[t2]->moves[idVec];
//                c23 = calC(nodes[t2], nodes[t3]);
//                if (c23 >= c12 - bestG / 2 - EP)break;
//                i2 = nodes[t3]->posInSol;
//                t4 = solT[i2 + 1];
//                curGain = c12 - c23
//                    + calC(nodes[t3], nodes[t4])
//                    - calC(nodes[t4], nodes[t1]);
//                if (curGain - EP <= bestG)continue;
//                Route* r1 = nodes[t1]->rou;
//                Route* r2 = nodes[t3]->rou;
//                if (r1 == r2)continue;
//                pos[1] = i1;
//                pos[2] = i1 + 1;
//                pos[3] = i2;
//                pos[4] = i2 + 1;                
//                for (int i = 1; i <= 4; ++i) pos[i] = solT[pos[i]];
//                for (int i = 1; i <= 4; ++i) posR[i] = nodes[pos[i]]->posInRoute;
//                if (pos[2] > n)posR[2] = r1->length;
//                if (pos[4] > n)posR[4] = r2->length;
//                //check new route 1:
//                sumL = r1->calQInRan(0, posR[1]) + r2->calQInRan(posR[4], r2->length);
//                sumE = r1->calEInRan(0, posR[1]) + r2->calEInRan(posR[4], r2->length)
//                    + r1->calQInRan(0, posR[1]) * (r2->calAInRan(posR[4], r2->length) + calA(nodes[pos[1]], nodes[pos[4]]))
//                    + calB(nodes[pos[1]], nodes[pos[4]]);
//                if (sumE > pr->maxE || sumL > pr->cap)continue;
//                //check new route 2:
//                sumL = r1->calQInRan(posR[2], r1->length) + r2->calQInRan(0, posR[3]);
//                sumE = r2->calEInRan(0, posR[3]) + r1->calEInRan(posR[2], r1->length)
//                    + r2->calQInRan(0, posR[3]) * (r1->calAInRan(posR[2], r1->length) + calA(nodes[pos[3]], nodes[pos[2]]))
//                    + calB(nodes[pos[3]], nodes[pos[2]]);
//                if (sumE > pr->maxE || sumL > pr->cap)continue;                
//                bestG = curGain;                                
//                for (int i = 1; i <= 4; ++i)resPos[i] = pos[i];
//                if (!bi)goto thispro;
//            }
//        }
//        thispro:
//        if (bestG != 0.0) {
//            cost -= bestG;
//            for (int i = 1; i <= 4; ++i) pos[i] = resPos[i];
//            if (pr->isDebug) {
//                cout << "improve 2opt*\n";
//                for (auto val : solT)cout << val << " ";
//                cout << endl;
//                for (int i = 1; i <= 4; ++i)cout << nodes[pos[i]]->rou << endl;
//                cout << endl;
//                for (int i = 1; i <= 4; ++i)cout << pos[i] << " ";
//                cout << endl; cout << endl;
//            }
//            Route* r1 = nodes[pos[1]]->rou;
//            Route* r2 = nodes[pos[3]]->rou;           
//            //create pos in route:            
//            Node* n1 = nodes[pos[1]];
//            Node* n2 = nodes[pos[2]];
//            Node* n3 = nodes[pos[3]];
//            Node* n4 = nodes[pos[4]];
//            if (pos[2] > n)n2 = r1->depot->pred;
//            if (pos[4] > n)n4 = r2->depot->pred;
//            Node* endNodeR1 = r1->depot->pred->pred;
//            Node* endNodeR2 = r2->depot->pred->pred;
//            //update route 1:
//            //remove from n2:
//            n1->suc = r1->depot->pred;
//            r1->depot->pred->pred = n1;           
//            //add n4
//            if (n4->idxClient != 0) {                
//                endNodeR2->suc = r1->depot->pred;
//                r1->depot->pred->pred = endNodeR2;
//                n1->suc = n4;
//                n4->pred = n1;
//            }           
//            //update route 2:
//            //remove from n4:
//            n3->suc = r2->depot->pred;
//            r2->depot->pred->pred = n3;
//            //add n2:
//            if (n2->idxClient != 0) {
//                endNodeR1->suc = r2->depot->pred;
//                r2->depot->pred->pred = endNodeR1;
//                n3->suc = n2;
//                n2->pred = n3;
//            }           
//            r1->updateRoute();
//            r2->updateRoute();
//            cvSolT();            
//        }
//        delete[] resPos;
//        delete[] pos;
//        delete[] posR;        
//        pr->crosver = false;
//        return (bestG != 0.0);
//    }
//
//    //SWAP:
//    //applied when u and v are clients
//    bool Swap2Node(bool bi) {
//        if (pr->isDebug) {
//            cout << "In Swap 2 node\n";
//            pr->swap2 = true;
//        }
//        int posI1, posI2;
//        double bestG = 0.0;
//        int v1, t1, w1;
//        int v2, t2, w2;
//        int i1, i2;        
//        double curGain;
//        int pos, posR1, posR2;
//        Route* r1;
//        Route* r2;
//        double sumE;
//        int sumL;
//        for (i1 = 1; i1 <= m + n; ++i1)if (solT[i1] <= n) {
//            for (int ap = -1; ap <= 1; ap += 2) {               
//                for (pos = 0; pos < nodes[solT[i1]]->moves.size(); ++pos) {
//                    v1 = solT[i1 - 1];
//                    t1 = solT[i1];
//                    w1 = solT[i1 + 1];
//                    i2 = nodes[nodes[t1]->moves[pos]]->posInSol + ap;
//                    i2 = gID(i2);
//                    if (calC(nodes[t1], nodes[nodes[t1]->moves[pos]]) >= 
//                    (calC(nodes[v1], nodes[t1]) + calC(nodes[t1], nodes[w1])) / 2 - bestG / 2 - EP)break;
//                    if (solT[i2] > n || abs(i1 - i2) == 1)continue;
//                    v2 = solT[i2 - 1]; t2 = solT[i2]; w2 = solT[i2 + 1];                    
//                     curGain = calC(nodes[v1], nodes[t1]) + calC(nodes[t1], nodes[w1])
//                        + calC(nodes[v2], nodes[t2]) + calC(nodes[t2], nodes[w2])
//                        - calC(nodes[v1], nodes[t2]) - calC(nodes[t2], nodes[w1])
//                        - calC(nodes[v2], nodes[t1]) - calC(nodes[t1], nodes[w2]);
//                    if (curGain - EP <= bestG)continue;                   
//                    if (i1 > i2) {
//                        swap(t1, t2);
//                        swap(v1, v2);
//                        swap(w1, w2);
//                    }
//                    r1 = nodes[t1]->rou;
//                    r2 = nodes[t2]->rou;
//                    posR1 = nodes[t1]->posInRoute;
//                    posR2 = nodes[t2]->posInRoute;
//                    if (r1 == r2) {
//                        // depot -> v
//                        if (posR1 > posR2)swap(posR1, posR2);
//                        sumE = r1->calEInRan(0, posR1 - 1)
//                            + r1->calQInRan(0, posR1 - 1) * calA(nodes[v1], nodes[t2])
//                            + calB(nodes[v1], nodes[t2]);
//                        sumL = r1->calQInRan(0, posR1 - 1) + pr->listCL[nodes[t2]->idxClient].demand;
//                        //depot -> behind u:
//                        sumE = sumE + r1->calEInRan(posR1 + 1, posR2 - 1)
//                            + sumL * (r1->calAInRan(posR1 + 1, posR2 - 1) + calA(nodes[t2], nodes[w1]))
//                            + calB(nodes[t2], nodes[w1]);
//                        sumL += r1->calQInRan(posR1 + 1, posR2 - 1);
//                        //depot->u:
//                        sumE += sumL * calA(nodes[v2], nodes[t1]) + calB(nodes[v2], nodes[t1]);
//                        sumL += pr->listCL[nodes[t1]->idxClient].demand;
//                        //depot->end:
//                        sumE += r1->calEInRan(posR2 + 1, r1->length)
//                            + sumL * (r1->calAInRan(posR2 + 1, r1->length) + calA(nodes[t1], nodes[w2]))
//                            + calB(nodes[t1], nodes[w2]);
//                        if (sumE > pr->maxE)continue;                        
//                        bestG = curGain;
//                        posI1 = i1;
//                        posI2 = i2;
//                        if (!bi)goto thispro;
//                    }
//                    else {
//                        //route 1:
//                        sumE = r1->calEInRan(0, posR1 - 1)
//                            + r1->calQInRan(0, posR1 - 1) * calA(nodes[v1], nodes[t2])
//                            + calB(nodes[v1], nodes[t2])//end depot->v
//                            + r1->calEInRan(posR1 + 1, r1->length)
//                            + ((r1->calQInRan(0, posR1 - 1) + pr->listCL[t2].demand)
//                                * (r1->calAInRan(posR1 + 1, r1->length) + calA(nodes[t2], nodes[w1])))
//                            + calB(nodes[t2], nodes[w1]);//end depot->depot;
//                        sumL = r1->calQInRan(0, r1->length) - pr->listCL[t1].demand + pr->listCL[t2].demand;
//                        if (sumE > pr->maxE || sumL > pr->cap)continue;
//                        //route 2:
//                        sumE = r2->calEInRan(0, posR2 - 1)
//                            + r2->calQInRan(0, posR2 - 1) * calA(nodes[v2], nodes[t1])
//                            + calB(nodes[v2], nodes[t1])//end depot->v
//                            + r2->calEInRan(posR2 + 1, r2->length)
//                            + ((r2->calQInRan(0, posR2 - 1) + pr->listCL[t1].demand)
//                                * (r2->calAInRan(posR2 + 1, r2->length) + calA(nodes[t1], nodes[w2])))
//                            + calB(nodes[t1], nodes[w2]);//end depot->depot;
//                        sumL = r2->calQInRan(0, r2->length) - pr->listCL[t2].demand + pr->listCL[t1].demand;
//                        if (sumE > pr->maxE || sumL > pr->cap)continue;                        
//                        //cout << "type 2:" << sumE << endl;                        
//                        bestG = curGain;
//                        posI1 = i1;
//                        posI2 = i2;
//                        if (!bi)goto thispro;
//                    }
//                }
//            }
//        }
//        thispro:
//        if (bestG != 0.0) {                        
//            cost -= bestG;
//            if (pr->isDebug) {
//                cout << "improve swap:\n";
//                cout << solT[posI1] << " " << solT[posI2]<<endl;
//                cout << endl;
//            }
//            //update:                                   
//            Node* u = nodes[solT[posI1]];
//            Node* v = nodes[solT[posI2]];
//            Node* predU = u->pred;
//            Node* sucU = u->suc;
//            Node* predV = v->pred;
//            Node* sucV = v->suc;
//            Route* rU = u->rou;
//            Route* rV = v->rou;
//            //update u:
//            u->suc = sucV; 
//            sucV->pred = u;
//            u->pred = predV;
//            predV->suc = u;
//            //update v:
//            v->suc = sucU;
//            sucU->pred = v;
//            v->pred = predU;
//            predU->suc = v;
//            //update route:           
//            if (rU == rV)rU->updateRoute();
//            else {                
//                rU->updateRoute();
//                rV->updateRoute();
//            }           
//            //update solT:                            
//            cvSolT();
//        }        
//        pr->swap2 = false;
//        return (bestG != 0.0);
//    }
//
//    //OROPT:
//    //check and update using concatenation (note this only can use for oropt and string exchange)
//    //can use reverse by adding (u, v) with u < 0 and v < 0
//    //(|u|, |v|) need to be in range [1, m+n] 
//    bool ckBySet(vector<II> setTour, bool isUpdate = false) {        
//        int numT = setTour.size();        
//        int* ddT = new int[numT];        
//        for (int i = 0; i < numT; ++i)ddT[i] = 0;
//        int fiID = -1; // id of first tour containing depot
//        int u, v;
//        for (int i = 0; i < numT; ++i) {            
//            u = abs(setTour[i].first);
//            v = abs(setTour[i].second);                        
//            if (nodes[solT[u]]->idxClient == 0 || nodes[solT[v]]->idxClient == 0)ddT[i] = 1; else
//                if (nodes[solT[u]]->rou != nodes[solT[v]]->rou)ddT[i] = 1; else
//                    if (setTour[i].first > setTour[i].second)ddT[i] = 1;
//            if (ddT[i] == 1 && fiID == -1)fiID = i;            
//        }                        
//        u = solT[abs(setTour[fiID].first)];
//        v = solT[abs(setTour[fiID].second)];
//        int posRu = nodes[u]->posInRoute;
//        int posRv = nodes[v]->posInRoute;
//        bool isReverse = false;
//        if (setTour[fiID].first < 0 || setTour[fiID].second < 0) isReverse = true;
//        Route* curRoute = nodes[v]->rou;
//        Route* upRoute = nodes[v]->rou;
//        if (isUpdate) {
//            if (!isReverse)upRoute->clearRouteFrom(nodes[v]);
//            else {                
//                upRoute->clearRouteFromRev(nodes[v]);                
//            }                      
//        }        
//        Node* prevV = nodes[v];
//        double sumE, sumL;
//        int idT = fiID, st, en;
//        int endLoop;
//        if (!isUpdate) {
//            if (!isReverse) {
//                sumE = curRoute->calEInRan(0, posRv);
//                sumL = curRoute->calQInRan(0, posRv);
//            }
//            else {
//                sumE = curRoute->calEInRan(curRoute->length, posRv);
//                sumL = curRoute->calQInRan(curRoute->length, posRv);
//                if (posRv == 0)sumE = sumL = 0;
//            }           
//        }                
//        for (int i = 0;; ++i) {
//            idT = (idT + 1) % numT;           
//            u = solT[abs(setTour[idT].first)];
//            v = solT[abs(setTour[idT].second)];
//            isReverse = false;
//            if (setTour[idT].first < 0 || setTour[idT].second < 0) isReverse = true;
//            if (!isUpdate) {                
//                curRoute = nodes[u]->rou;
//                posRu = nodes[u]->posInRoute;
//                posRv = nodes[v]->posInRoute;                
//                st = posRu;                
//                if (ddT[idT] == 0)en = posRv;
//                else if (!isReverse) en = curRoute->length;
//                else en = 0;                                                  
//                //check energy and load
//                sumE += curRoute->calEInRan(st, en)
//                    + sumL * (curRoute->calAInRan(st, en) + calA(prevV, nodes[u]))
//                    + calB(prevV, nodes[u]);
//                sumL += curRoute->calQInRan(st, en);                      
//                if (sumE > pr->maxE || sumL > pr->cap) {
//                    delete[] ddT;
//                    return false;
//                }
//                prevV = nodes[v];
//            }
//            else {      
//                int idSolT;                
//                endLoop = solT[abs(setTour[idT].second)];   
//                if (setTour[idT].second < setTour[idT].first) endLoop = -1;
//                if (!isReverse) {
//                    for (int j = abs(setTour[idT].first);; ++j) {
//                        idSolT = solT[gID(j)];
//                        if (idSolT > n)break;
//                        upRoute->insertToRou(nodes[idSolT]);
//                        if (idSolT == endLoop)break;
//                    }
//                }
//                else {
//                    for (int j = abs(setTour[idT].first);; --j) {
//                        idSolT = solT[gID(j)];
//                        if (idSolT > n)break;
//                        upRoute->insertToRou(nodes[idSolT]);
//                        if (idSolT == endLoop)break;
//                    }
//                }
//            }  
//            //exit(0);
//            if (idT == fiID) {
//                /*if (!isUpdate) {
//                    cout <<"Energy in LS: " << sumE << endl;
//                }*/
//                if (isUpdate) {                                                            
//                    upRoute->updateRoute();                    
//                }
//                break;
//            }
//            if (ddT[idT]) {            
//                if (isUpdate) {
//                    upRoute->updateRoute();                    
//                }                
//                upRoute = nodes[v]->rou;                                                       
//                //reset when find a new depot: 
//                if (isUpdate) {
//                    if (!isReverse)upRoute->clearRouteFrom(nodes[v]);
//                    else upRoute->clearRouteFromRev(nodes[v]);
//                }
//                else {
//                    if (!isReverse) {
//                        sumE = upRoute->calEInRan(0, posRv);
//                        sumL = upRoute->calQInRan(0, posRv);
//                    }
//                    else {
//                        sumE = upRoute->calEInRan(curRoute->length, posRv);
//                        sumL = upRoute->calQInRan(curRoute->length, posRv);
//                        if (posRv == 0)sumE = sumL = 0;
//                    }
//                }
//            }
//        }           
//        if (isUpdate)cvSolT();        
//        delete[] ddT;        
//        return true;
//    }
//    
//    int gID(int val) {
//        if (val > m + n)return val % (m + n);
//        if (val <= 0)return (val % (m + n)) + m + n;
//        return val;
//    }
//    bool OrOptMove(int k, bool bi) {
//        if (pr->isDebug) {
//            cout << "In or opt\n";
//            pr->moveor = true;
//        }
//        vector<II> setTours(3, II(0,0));        
//        vector<II> bestTours(3, II(0, 0));
//        double bestG = 0.0;
//        int i1, i2, i3;
//        int t1,  t2,  t3,  t4,  t5,  t6;
//        int valAdd;
//        double c12, c23, c34, c45, c56, c61;
//        int idVec;
//        double curGain;
//        //LOOP 1:
//        for (i1 = 1; i1 <= m + n; ++i1) {
//            t1 = solT[i1];
//            t2 = solT[gID(i1 + 1)];
//            c12 = calC(nodes[t1], nodes[t2]);            
//            for (idVec = 0; idVec < nodes[t2]->moves.size(); ++idVec) {
//                t3 = nodes[t2]->moves[idVec];
//                c23 = calC(nodes[t2], nodes[t3]);
//                if (c23 >= c12 - bestG / 3 - EP)break;
//                i2 = nodes[t3]->posInSol; 
//                t4 = solT[gID(i2 + 1)];
//                for (valAdd = 1; valAdd <= k; ++valAdd) {                    
//                    i3 = gID(i1 + valAdd);                    
//                    if (!(
//                        (i1 < i3 && i3 < i2)
//                        || (i3 < i2 && i2 < i1)
//                        || (i2 < i1 && i1 < i3)))continue;                    
//                    t5 = solT[i3], t6 = solT[gID(i3 + 1)];
//                    curGain = c12 - c23
//                        + calC(nodes[t3], nodes[t4]) - calC(nodes[t4], nodes[t5])
//                        + calC(nodes[t5], nodes[t6]) - calC(nodes[t6], nodes[t1]);
//                    if (bestG >= curGain - EP)continue;                    
//                    setTours[0] = II(gID(i2 + 1), i1);
//                    setTours[1] = II(gID(i3 + 1), i2);
//                    setTours[2] = II(gID(i1 + 1), i3);                    
//                    if (ckBySet(setTours)) {                                                
//                        bestG = curGain;
//                        for (int i = 0; i <= 2; ++i)bestTours[i] = setTours[i];
//                        if (!bi)goto thispro;
//                    }                                        
//                }
//            }
//        }
//
//        //LOOP2:
//        for (i2 = 1; i2 <= m + n; ++i2) {
//            t3 = solT[i2];
//            t4 = solT[gID(i2 + 1)];
//            c34 = calC(nodes[t3], nodes[t4]);
//            for (idVec = 0; idVec < nodes[t4]->moves.size(); ++idVec) {
//                t5 = nodes[t4]->moves[idVec];
//                c45 = calC(nodes[t4], nodes[t5]);
//                if (c45 >= c34 - bestG / 3 - EP)break;
//                i3 = nodes[t5]->posInSol;
//                t6 = solT[gID(i3+1)];
//                for (valAdd = -k; valAdd <= -1; ++valAdd) {
//                    i1 = gID(i3 + valAdd);
//                    if (!(
//                        (i1 < i3 && i3 < i2)
//                        || (i3 < i2 && i2 < i1)
//                        || (i2 < i1 && i1 < i3)))continue;
//                    t1 = solT[i1], t2 = solT[gID(i1+1)];
//                    curGain = calC(nodes[t1], nodes[t2]) - calC(nodes[t2], nodes[t3])
//                        + c34 - c45
//                        + calC(nodes[t5], nodes[t6]) - calC(nodes[t6], nodes[t1]);
//                    if (bestG >= curGain - EP)continue;
//                    setTours[0] = II(gID(i2 + 1), i1);
//                    setTours[1] = II(gID(i3 + 1), i2);
//                    setTours[2] = II(gID(i1 + 1), i3);
//                    if (ckBySet(setTours)) {                                                
//                        bestG = curGain;
//                        for (int i = 0; i <= 2; ++i)bestTours[i] = setTours[i];
//                        if (!bi)goto thispro;
//                    }
//                }
//            }
//        }
//
//        //LOOP3:
//        double B2;
//        for (i3 = 1; i3 <= m + n; ++i3) {
//            t5 = solT[i3];
//            t6 = solT[gID(i3 + 1)];
//            c56 = calC(nodes[t5], nodes[t6]);
//            for (valAdd = -k; valAdd <= -1; ++valAdd) {
//                i1 = gID(i3 + valAdd);
//                t1 = solT[i1];
//                c61 = calC(nodes[t1], nodes[t6]);
//                if (c61 >= c56 - bestG / 3 - EP)continue;
//                t2 = solT[gID(i1 + 1)];
//                c12 = calC(nodes[t1], nodes[t2]);
//                B2 = c56 - c61 + c12 - (2 * bestG) / 3;
//                for (idVec = 0; idVec < nodes[t2]->moves.size(); ++idVec) {
//                    t3 = nodes[t2]->moves[idVec];
//                    c23 = calC(nodes[t2], nodes[t3]);
//                    if (c23 >= B2 - EP)break;
//                    i2 = nodes[t3]->posInSol;
//                    t4 = solT[gID(i2 + 1)];
//                    if (!(
//                        (i1 < i3 && i3 < i2)
//                        || (i3 < i2 && i2 < i1)
//                        || (i2 < i1 && i1 < i3)))continue;
//                    curGain = c12 - c23
//                        + calC(nodes[t3], nodes[t4]) - calC(nodes[t4], nodes[t5])
//                        + c56 - c61;
//                    if (bestG >= curGain - EP)continue;
//                    setTours[0] = II(gID(i2 + 1), i1);
//                    setTours[1] = II(gID(i3 + 1), i2);
//                    setTours[2] = II(gID(i1 + 1), i3);
//                    if (ckBySet(setTours)) {                                                
//                        bestG = curGain;
//                        for (int i = 0; i <= 2; ++i)bestTours[i] = setTours[i];
//                        if (!bi)goto thispro;
//                    }
//                }
//            }
//        }
//
//        //UPDATE:        
//        thispro:
//        if (bestG != 0.0) {
//            if (pr->isDebug) {
//                cout << "pos in solT\n";
//                for (auto val : bestTours)cout << val.first << " " << val.second << endl;
//                cout << "index nodes\n";
//                for (auto val : bestTours)cout << solT[abs(val.first)] << " " << solT[abs(val.second)] << endl;
//            }
//            cost -= bestG;
//            ckBySet(bestTours, true);            
//        }
//        pr->moveor = false;
//        return (bestG != 0.0);
//    }
//
//    //String-exchange:
//    bool StringExchange(int k, bool bi) {
//        if (pr->isDebug) {
//            cout << "in string exchange\n";
//            pr->stringex = true;
//        }
//        double bestG = 0.0;
//        vector<II> setTours(4, II(0, 0));
//        vector<II> bestTours(4, II(0, 0));
//        for (int i1 = 1; i1 <= m + n; ++i1)
//            for (int ap = -1; ap <= 1; ap += 2) {
//                int t1 = solT[i1], t2 = solT[i1 + ap];
//                double c12 = calC(nodes[t1], nodes[t2]);
//                for (auto t3 : nodes[t2]->moves) {
//                    double c23 = calC(nodes[t2], nodes[t3]);
//                    if (c23 >= c12 - bestG / 4 - EP)break;
//                    int i2 = nodes[t3]->posInSol, t4 = solT[gID(i2-ap)];
//                    double c34 = calC(nodes[t3], nodes[t4]);
//                    double c41 = calC(nodes[t1], nodes[t4]);
//                    if (c12 - c23 + c34 - c41 - bestG / 2 + EP <= 0.0) continue;
//                    for (int j1 = 2; j1 <= k + 1; ++j1) {
//                        int i3 = gID(i1 + j1 * ap);
//                        int t5 = solT[i3], t6 = solT[gID(i3 - ap)];
//                        for (int j2 = 2; j2 <= k + 1; ++j2) {
//                            int i4 = gID(i2 - j2 * ap);
//                            int t7 = solT[i4], t8 = solT[gID(i4 + ap)];
//                            double c56 = calC(nodes[t5], nodes[t6]);
//                            double c67 = calC(nodes[t7], nodes[t6]);
//                            double c78 = calC(nodes[t7], nodes[t8]);
//                            double c85 = calC(nodes[t5], nodes[t8]);
//                            double curGain = c12 - c23 + c34 - c41 + c56 - c67 + c78 - c85;
//                            if (bestG >= curGain - EP)continue;
//                            if (ap == 1) {                                
//                                setTours[0] = II(i2, i1);
//                                setTours[1] = II(-gID(i2 - 1), -gID(i4 + 1));
//                                setTours[2] = II(i3, i4);
//                                setTours[3] = II(-gID(i3-1), -gID(i1+1));
//                            }
//                            else {
//                                setTours[0] = II(i1, i2);
//                                setTours[1] = II(-gID(i1 - 1), -gID(i3 + 1));
//                                setTours[2] = II(i4, i3);
//                                setTours[3] = II(-gID(i4 - 1), -gID(i2 + 1));
//                            }   
//                            int stU = abs(setTours[1].second), enU = abs(setTours[1].first);
//                            int stV = abs(setTours[3].second), enV = abs(setTours[3].first);
//                            if (gID(enU + 1) == stV || gID(enV + 1) == stU)continue;
//                            if (stU >= stV && stU <= enV)continue;
//                            if (stV >= stU && stV <= enU)continue;
//                            bool ckLegal = true;
//                            for (int i = 1; i <= 3; i += 2) {
//                                int u = abs(setTours[i].first);
//                                int v = abs(setTours[i].second);
//                                if (nodes[solT[u]]->idxClient == 0 || nodes[solT[v]]->idxClient == 0) ckLegal = false; else
//                                    if (nodes[solT[u]]->rou != nodes[solT[v]]->rou)ckLegal = false; else
//                                        if (setTours[i].first > setTours[i].second)ckLegal = false;
//                            }
//                            if (!ckLegal)continue;                            
//                            if (ckBySet(setTours)) {                                  
//                                bestG = curGain;                                    
//                                for (int i = 0; i < 4; ++i)bestTours[i] = setTours[i];
//                                if (!bi)goto thispro;
//                            }                                                        
//                        }
//                    }
//                }
//            }
//        //UPDATE:        
//        thispro:
//        if (bestG != 0.0) {
//            if (pr->isDebug) {
//                cout << "improve string exchange\n";
//                cout << "pos in solT: \n";
//                for (auto val : bestTours)cout << val.first << " " << val.second << endl;
//                cout << "\n";
//                cout << "index node\n";
//                for (auto val : bestTours)cout << solT[abs(val.first)] << " " << solT[abs(val.second)] << endl;
//                cout << "\n\n";
//            }
//            //for (int i = 0; i < 4; ++i)cout << bestTours[i].first << " " << bestTours[i].second << endl;
//            cost -= bestG;
//            ckBySet(bestTours, true);
//        }
//        pr->stringex = false;
//        return (bestG != 0.0);        
//    }
//
//    //ELSALGO:
//    void ELS() {
//        initSol();
//        int* resGiantT = new int[n + 1];        
//        int* curGiantT = new int[n + 1];
//        for (int i = 1; i <= n; ++i) {
//            curGiantT[i] = giantT[i];
//            resGiantT[i] = giantT[i];
//            //cout << giantT[i] << " " << endl;
//        }        
//        double bestObj=cost;
//        cout << cost << endl;
//        double curObj;
//        int curP = pr->pMin;
//        for (int i = 1; i <= pr->nI; ++i) {
//            curObj = bestObj;
//            for (int j = 1; j <= pr->nC; ++j) {                
//                for (int i1 = 1; i1 <= n; ++i1)giantT[i1] = curGiantT[i1];
//                mutate(curP);                
//                Split();                
//                if (cost == oo)continue;                
//                try {
//                    updateObj();
//                }
//                catch (...) {
//                    cout << "bug here\n";
//                    for (auto val : giantT)cout << val << ", ";
//                    cout << endl;
//                    system("pause");
//                    exit(0);
//                }
//                if (cost + EP < curObj) {
//                    curObj = cost;
//                    for (int i1 = 1; i1 <= n; ++i1)resGiantT[i1] = giantT[i1];
//                    curP = pr->pMin;
//                }
//                else {
//                    curP = min(pr->pMax, curP + 1);
//                }
//
//                /*if ((clock() - pr->start) / CLOCKS_PER_SEC >= pr->TL) {
//                    if (curObj + EP < bestObj) {
//                        bestObj = curObj;
//                        goto thispro;
//                    }
//                }*/
//                //cout << i << " " << j << "\n";
//            } 
//            //cout << i << endl;
//            if (curObj + EP < bestObj) {
//                bestObj = curObj;
//                for (int i1 = 1; i1 <= n; ++i1)curGiantT[i1] = resGiantT[i1];
//            }
//        }        
//        //thispro:
//        cout << "best found obj: " << Util::round2num(bestObj) << endl;
//        pr->fileOut<< "best found obj: " << Util::round2num(bestObj) << endl;
//        //if (pr->debugLS) {            
//            for (int i = 1; i <= n; ++i)giantT[i] = resGiantT[i];
//            Split();            
//            if (cost != oo && cost - bestObj > EP) {
//                cout << "bug here\n";
//                system("pause");
//                exit(0);
//            }
//            cout << Util::round2num(cost) << endl;
//            pr->fileOut << "Giant Tour:\n";
//            for (int i = 1; i <= n; ++i)pr->fileOut << giantT[i] << ", ";
//            pr->fileOut << "\n";
//            pr->fileOut<<"cost split: " << Util::round2num(cost) << endl;
//        //}
//        delete[] resGiantT;
//        delete[] curGiantT;
//    }
//
//    //deconstructor:
    ~Solution(){        
        /*for (int k = 0; k <= m; ++k) {
            delete[] F[k];
            delete[] pred[k];
        }*/
        delete[] F;
        delete[] pred;
        giantT.clear();
        solT.clear();                
        for (int i = 1; i <= n + m; ++i)nodes[i]->moves.clear();        
        for (int i = 1; i <= m; ++i)delete setR[i];
        delete pr;
    }
};
