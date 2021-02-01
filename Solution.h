#pragma once
#include "lib.h"
#include "Node.h"
#include "Route.h"

class Solution {
public:
    vector<int> giantT;// giant tour
    vector<int> solT;// solution tour (contains index of clients)
    vector<int> ordNodeLs; // using for LS
    vector<Node*> nodes;// for nodes
    vector<Route*> setR;// for set Route      
    vector<SeqData*> myseqs1;// for concatenation in LS    
    vector<SeqData*> myseqs2;// for concatenation in LS        
    vector<int>* lstNeig;// getting current list of neiborhood;
    SeqData* seqSet;// init sequences for each node.    
    SeqData* seqDep;// only contain depot
    //GA parts:
    vector<int> predecessors; // store prev node of each client  
    vector<int> successors; // store next node of each client  
    multiset<pair<int, Solution*>> indivsPerProximity;
    double biasedFitness; // Biased fitness of the solution
    /*
    * Variables for dealing with LS
    */
    long long count[4][4];
    int resGenInMoves[4][4];//using for general insert
    III resSubFlex1[4][4];//using for flexible version
    III resSubFlex2[4][4];//using for flexible version
    vector<int> traceLoc[4][4];
    vector<int> traceLoc1[4][4];
    vector<int> traceLoc2[4][4];
    vector<III> lstLabel; //((cost, time), loc)
    vector<III> curLabel; //((cost, time), prv)
    vector<int> prvIdLb; //previous label in vector
    vector<int> tourLoc;
    Node* nodeU;
    Node* uPred;
    Node* uSuc;
    Node* nodeV;
    Node* vPred;
    Node* vSuc;
    Route* routeU;
    Route* routeV;
    bool isFixed = false;// for identifying the type of LS;

    Param* pr;
    int cost;// objective
    int n;// number of customer
    int m;// number of vehicle (use when limit number of vehicle)
    int curVeh; // used for transfering solution 
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
        if (pr->maxNeibor > n)pr->maxNeibor = n - 1;// deal with special case (exploring all neiborhood)
        if (pr->numVeh != oo)m = pr->numVeh;
        else m = n;
        m++;//do for initialization
        // create giant tour        
        for (int i = 0; i <= n; ++i) {
            giantT.pb(0);//indexing from 1        
            ordNodeLs.push_back(i + 1);
        }
        ordNodeLs.pop_back();
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
        for (int i = 1; i <= m; ++i) {
            setR[i]->updateRoute();
        }
        //only contain depot:
        seqDep = new SeqData(_pr);
        seqDep->init(0);
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

        indivsPerProximity.clear();
        predecessors.resize(n + 1);
        successors.resize(n + 1);
        /*
        * init for split
        */               
        F = new int[n + 1];
        pred = new int[n + 1];
    }

    bool notDominatePop(II& label, II& labelA, II& labelB, II& labelC)
    {
        if (label.first >= labelA.first && label.second >= labelA.second)
            if (label.first != labelA.first || label.second != labelA.second) return false;
        if (label.first >= labelB.first && label.second >= labelB.second)
            if (label.first != labelB.first || label.second != labelB.second) return false;
        if (label.first >= labelC.first && label.second >= labelC.second)
            if (label.first != labelC.first || label.second != labelC.second) return false;
        return true;
    }
    bool notDominatePush(II& label, II& labelA, II& labelB, II& labelC)
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

    void genGiantT() {
        for (int i = 1; i <= n; ++i)giantT[i] = i;
        shuffle(giantT.begin() + 1, giantT.end(), pr->Rng.generator);
        //random_shuffle(giantT.begin() + 1, giantT.end());
    }

    /*diversity contribution in GA*/
    void removeProximity(Solution* indiv)
    {
        auto it = indivsPerProximity.begin();
        while (it->second != indiv) ++it;
        indivsPerProximity.erase(it);
    }


    double brokenPairsDistance(Solution* valSol) {
        int differences = 0;
        for (int j = 1; j <= n; j++)
        {
            if (successors[j] != valSol->successors[j] && successors[j] != valSol->predecessors[j]) differences++;
            if (predecessors[j] == 0 && valSol ->predecessors[j] != 0 && valSol->successors[j] != 0) differences++;
        }
        return (double)differences / (double)n;
    }

    double averageBrokenPairsDistanceClosest(int nbClosest) {
        double result = 0;
        int maxSize = min<int>(nbClosest, indivsPerProximity.size());
        auto it = indivsPerProximity.begin();
        for (int i = 0; i < maxSize; i++)
        {
            result += it->first;
            ++it;
        }
        return result / (double)maxSize;
    }
    /**/

    void cvGiantT() {
        int pos = 0;
        for (int i = 1; i <= min(m + 1, n); ++i) {
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
                throw "Wrong here 1";
            }
            if (pos != n) {
                throw "Wrong here 2";
            }
            delete[] dd;
        }
    }


    bool ckSol() {
        if (pr->isDebug)cout << "Current Check\n";
        vector<int> arrSol;
        int* dd = new int[n + 1];
        for (int i = 1; i <= n; ++i)dd[i] = 0;
        int totalCost = 0;
        for (int i = 1; i <= min(m + 1, n); ++i) {
            if (pr->isDebug)setR[i]->showR();
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
        if (!(totalCost == cost)) {
            throw "error total cost";
        }
        delete[] dd;
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
            for (int i = 1; i <= min(m + 1, n); ++i)setR[i]->isNodeTested[tempNode->idxClient] = false;
        }
    }

    //reset all moves on all routes
    void reinitAllMovesInRou() {
        for (int i = 1; i <= min(m + 1, n); ++i) {
            for (int j = 1; j <= n; ++j)setR[i]->isNodeTested[j] = false;
        }
    }

    //marked that node had tested on all routes:
    void markNodeTested(Node* nod) {
        for (int i = 1; i <= min(m + 1, n); ++i)setR[i]->isNodeTested[nod->idxClient] = true;
    }

    //reset neigborhood cluster set (run when changing location)
    void reinitCluSet(Node* nod) {
        for (auto val : nod->movesClu) {
            nod->idxCluMoves[val] = false;
        }
        nod->movesClu = nod->rou->pr->listLoc[nod->idxLoc].moves;
        for (auto val : nod->movesClu) {
            nod->idxCluMoves[val] = true;
        }
    }

    //reset neigborhood set for each node
    void reinitNegiborSet(bool isConClu = true) {
        set<II> sDis;
        int u, v;
        for (int i = 1; i <= n; ++i) {
            u = nodes[i]->idxLoc;
            //cluster set:
            if (isConClu)reinitCluSet(nodes[i]);
            //location set: 
            /// (can update on the fly by using data structure set but trading off when processing)
            //reset all existing index
            for (auto val : nodes[i]->movesLoc) {
                nodes[i]->idxLocMoves[pr->listLoc[val].idxClient] = false;
            }

            //checking based on the corelation measure:
            sDis.clear();
            for (int j = 1; j <= n; ++j) if (j != i) {
                v = nodes[j]->idxLoc;
                sDis.insert(II(pr->corDis[v][u], v));
            }
            //adding new index:
            nodes[i]->movesLoc.clear();
            for (auto val : sDis) {
                if (val.first >= oo)continue;
                nodes[i]->movesLoc.push_back(pr->listLoc[val.sc].idxClient);
                nodes[i]->idxLocMoves[nodes[i]->movesLoc.back()] = true;
                if (nodes[i]->movesLoc.size() == pr->maxNeibor)break;
            }
        }
    }
    // split without limit number of vehicles
    void Split() {
        /*
        * (cost, time)
        * idx of location
        */
        if (pr->isDebug)cout << m << "\n";
        if (m != pr->numVeh + 1 && m != n + 1) {
            //rerun split
            //clear all route:
            for (int i = 1; i <= m; ++i) {
                setR[i]->clearNode();
                setR[i]->updateRoute();
            }
            reinitAllMovesInRou();
        }
        else m--;
        /*vector<III> lstLabel; //((cost, time), loc)
        vector<III> curLabel; //((cost, time), prv)
        vector<int> prvIdLb; //previous label in vector*/        
        lstLabel.clear();
        curLabel.clear();
        prvIdLb.clear();
        //the index that cannot exceed in a route due to the capacity.
        int* maxIdx = new int[n + 1];
        int curLoad = 0;
        maxIdx[0] = 0;
        for (int i = 1; i <= n; ++i) {
            maxIdx[i] = max(maxIdx[i - 1], i);
            int prvIdx = giantT[i - 1];
            curLoad -= pr->listCL[prvIdx].demand;
            if (maxIdx[i - 1] < i)curLoad = pr->listCL[giantT[i]].demand;
            if (maxIdx[i - 1] == n) {
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
                            //|| timeV + pr->times[idLocV][0] > pr->T // only true when travel times satisfy triangle inequality 
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
                    for (int i = 0; i < curLabel.size(); ++i) {
                        III val = curLabel[i];
                        if (val.ft.sc < minT)
                        {
                            minT = val.ft.sc;
                            lstLabel.pb(III(val.ft, idLocV));
                            prvIdLb.pb(val.sc);
                            //update F function
                            if (val.ft.sc + pr->times[idLocV][0] < pr->T && // uncomment when travel times does not satisfy triangle inequality
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
        cost = F[n];
        if (cost == oo) {                
            delete[] maxIdx;
            return;
        }
        int indexLb = pred[n];// index of last label.
        ///construct solution
        //vector<int> tourLoc;
        tourLoc.clear();
        int numVeh = 0;
        while (indexLb != -1)
        {
            tourLoc.pb(lstLabel[indexLb].second);
            if (lstLabel[indexLb].second == 0)numVeh++;
            indexLb = prvIdLb[indexLb];
        }
        m = numVeh;
        if (!(tourLoc.size() == m + n)) {
            throw "split bug loc size";
        }
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
        curVeh = numVeh;
        for (int i = 1; i <= min(m + 1, n); ++i) {               
            setR[i]->updateRoute();
            //setR[i]->showR();
        }
        reinitNegiborSet();        
        for (int i = 1; i <= n; ++i) {
            predecessors[i] = nodes[i]->pred->idxClient;
            successors[i] = nodes[i]->suc->idxClient;
        }
        if (pr->isTurnCkSol) {
            try {
                ckSol();
            }
            catch (const char* msg) {
                cerr << msg << endl;
                for (int i = 1; i <= n; ++i)pr->fileOut << giantT[i] << ", ";
                pr->fileOut << "\n";
                for (auto val : ordNodeLs)pr->fileOut << val << ", ";
                pr->fileOut.close();
                system("pause");
                exit(0);                
            }
        }
        //cvSolT(); // uncomment when need tracking the specific postion in solution (used in sequential search)        
        delete[] maxIdx;
    }

    void initSol() {
        vector<DD> points(n + 1);
        set<DI> setAng;
        setAng.clear();
        points[0] = DD(pr->listLoc[0].x, pr->listLoc[0].y);
        double avgX, avgY;
        for (int i = 1; i <= n; ++i) {
            avgX = 0.0;
            avgY = 0.0;
            for (auto val : pr->listCL[i].listLoc) {
                avgX += pr->listLoc[val].x;
                avgY += pr->listLoc[val].y;
            }
            avgX /= pr->listCL[i].listLoc.size();
            avgY /= pr->listCL[i].listLoc.size();
            //caculate angle with depot:
            setAng.insert(DI(atan2(avgY - points[0].sc, avgX - points[0].ft), i));
        }
        int idPos = 0;
        for (auto val : setAng) {
            giantT[++idPos] = val.sc;
        }
        if (idPos != n) {
            throw "error init";
        }
        Split();       
        cout << "initial cost: " << cost << "\n";
        setAng.clear();
        points.clear();
    }
    ///Local Search:
    /// isFixed =  true if running fixed version LS
    /// Otherwise if running flexible version LS
    /// 
    /// 
    // insert node u after node v
    void insertNode(Node* u, Node* v)
    {
        if (u->pred != v && u != v)
        {
            u->pred->suc = u->suc;
            u->suc->pred = u->pred;
            v->suc->pred = u;
            u->pred = v;
            u->suc = v->suc;
            v->suc = u;
            u->rou = v->rou;
        }
    }

    bool isCorrelated(Node* u, Node* v) {
        return true;
        if (u->idxClient == 0 || v->idxClient == 0)return true;
        if (isFixed)return v->idxLocMoves[u->idxClient];
        else return v->idxCluMoves[u->idxClient];
    }

    void updateObjInter() {
        //shuffle the moves:  
        shuffle(ordNodeLs.begin(), ordNodeLs.end(), pr->Rng.generator);
        /*for (int i = 1; i <= n; ++i) {            
            shuffle(nodes[i]->movesClu.begin(), nodes[i]->movesClu.end(), pr->Rng.generator);
            shuffle(nodes[i]->movesLoc.begin(), nodes[i]->movesLoc.end(), pr->Rng.generator);
        }*/
        //reset all moves on route:
        reinitAllMovesInRou();
        //checking move:
        int isMoved = 0;
        bool isMetEmpRou = false;
        bool isEndSearch = false;
        Node* tempNode;
        while (!isEndSearch)
        {
            isEndSearch = true;
            isMoved = 0;                        
            for (int posU = 0; posU < ordNodeLs.size(); ++posU) {
                //cout << "posU: " << posU << " " << isMoved << "\n";
                posU -= isMoved;
                isMoved = 0;
                isMetEmpRou = false;
                nodeU = nodes[ordNodeLs[posU]];
                routeU = nodeU->rou;
                uSuc = nodeU->suc;
                uPred = nodeU->pred;
                if (isFixed)lstNeig = &nodeU->movesLoc;
                else lstNeig = &nodeU->movesClu;
                for (int posV = 0; posV < lstNeig->size() && isMoved == 0; ++posV) {
                    nodeV = nodes[lstNeig->at(posV)];
                    routeV = nodeV->rou;
                    if (nodeV->rou->isNodeTested[nodeU->idxClient])continue;
                    vPred = nodeV->pred;
                    vSuc = nodeV->suc;
                    if (routeU != routeV) {
                        //apply inter general insert:
                        if (isMoved != 1) {
                            tempNode = nodeV;
                            nodeV = nodeV->suc;
                            vSuc = nodeV->suc;
                            vPred = nodeV->pred;
                            isMoved = interRouteGeneralInsert();
                            nodeV = tempNode;
                            vSuc = nodeV->suc;
                            vPred = nodeV->pred;
                        }
                        //2-Opt*
                        if (isMoved != 1) {
                            isMoved = interRoute2Opt();
                        }
                        //2-Opt* Inverse
                        if (isMoved != 1) {
                            isMoved = interRoute2OptInv();
                        }
                    }
                    else /*if (isFixed)*/ {
                        //intra route general insert (only consider with fixed version):                        
                        //isFixed = true;
                        tempNode = nodeV;
                        nodeV = nodeV->suc;
                        vSuc = nodeV->suc;
                        vPred = nodeV->pred;
                        int old_cost = cost;
                        if (isMoved != 1) {
                            isMoved = intraRouteGeneralInsert();
                        }
                        //isFixed = false;
                        nodeV = tempNode;
                        routeV = nodeV->rou;
                        vSuc = nodeV->suc;
                        vPred = nodeV->pred;
                    }
                }

                //intra route 2Opt:
                //cout << "status " << isMoved << "\n";
                nodeV = nodeU->suc;
                routeV = nodeV->rou;
                if (!routeV->isNodeTested[nodeU->idxClient]) {
                    //isFixed = true;
                    while (isMoved != 1 && nodeV->idxClient)
                    {
                        vPred = nodeV->pred;
                        vSuc = nodeV->suc;
                        if (isCorrelated(uPred, nodeV) || isCorrelated(nodeU, vSuc))
                            isMoved = intraRoute2Opt();
                        nodeV = nodeV->suc;
                    }
                    //isFixed = false;
                }

                // Special cases : testing the insertions behind the depot, and the empty routes
                //cout << "status " << isMoved << "\n";
                for (int idRou = 1; idRou <= min(m + 1, n) && isMoved == 0; ++idRou) {
                    nodeV = setR[idRou]->depot;
                    routeV = nodeV->rou;
                    vSuc = nodeV->suc;
                    vPred = nodeV->pred;
                    if (routeV->isNodeTested[nodeU->idxClient])continue;
                    if (isMetEmpRou && vSuc->idxClient == 0)continue;
                    if (vSuc->idxClient == 0)isMetEmpRou = true;
                    if (routeU != routeV) {
                        //apply inter general insert:
                        if (isMoved != 1) {
                            tempNode = nodeV;
                            nodeV = nodeV->suc;
                            vSuc = nodeV->suc;
                            vPred = nodeV->pred;
                            isMoved = interRouteGeneralInsert();
                            nodeV = tempNode;
                            vSuc = nodeV->suc;
                            vPred = nodeV->pred;
                        }
                        //2-Opt*
                        if (isMoved != 1) {
                            isMoved = interRoute2Opt();
                        }
                        //2-Opt* Inverse
                        if (isMoved != 1) {
                            isMoved = interRoute2OptInv();
                        }
                    }
                    else /*if (isFixed)*/ {
                        //intra route general insert (only consider with fixed version):                        
                        //isFixed = true;
                        tempNode = nodeV;
                        nodeV = nodeV->suc;
                        vSuc = nodeV->suc;
                        vPred = nodeV->pred;
                        if (isMoved != 1) {
                            isMoved = intraRouteGeneralInsert();
                        }
                        //isFixed = false;
                        nodeV = tempNode;
                        routeV = nodeV->rou;
                        vSuc = nodeV->suc;
                        vPred = nodeV->pred;
                    }
                    if (isMoved != 0 && idRou == m + 1) {
                        m++;
                    }
                }
                //cout << "status " << isMoved << "\n";

                if (isMoved == 0) {
                    markNodeTested(nodeU);
                }
                else {                    
                    isEndSearch = false;                    
                }
            }
        }        
    }

    //flexSeq must have larger than 2 sequences
    //with flexible version, it will change locations of the last node of first sequence and first node of last sequence .
    III evalFlex(vector<SeqData*> flexSeq) {        
        III res = III(II(-1, -1), oo);// (location of 2 changed node and result)        
        SeqData* firstE = flexSeq.front();
        SeqData* lastE = flexSeq.back();
        int cliFirst = pr->listLoc[firstE->lastnode].idxClient;
        int cliBe;
        if (firstE->beforeLaNode != -1)cliBe = pr->listLoc[firstE->beforeLaNode].idxClient;
        int cliLast = pr->listLoc[lastE->firstnode].idxClient;
        int cliAf;
        if (lastE->afterFiNode != -1)cliAf = pr->listLoc[lastE->afterFiNode].idxClient;
        int valCli;
        SeqData* seqLocU = new SeqData(pr);
        SeqData* seqLocV = new SeqData(pr);
        SeqData* seqTempU = seqLocU;
        SeqData* seqTempV = seqLocV;
        int ckRes;
        for (auto locU : pr->listCL[cliFirst].listLoc) {
            if (locU == 0) {
                seqLocU = seqDep;
            }
            else {
                if (cliBe == 0)seqLocU->concatOneAfter(seqDep, locU);
                else if (nodes[cliFirst]->pred == nodes[cliBe]) seqLocU->concatOneAfter(nodes[cliBe]->seq0_i, locU);
                else seqLocU->concatOneAfter(nodes[cliBe]->seqn_i, locU);
            }
            if (!seqLocU->F)continue;
            //flexSeq[0] = seqLocU;
            for (auto locV : pr->listCL[cliLast].listLoc) {
                if (locV == 0) {
                    seqLocV = seqDep;
                }
                else {
                    if (cliAf == 0)seqLocV->concatOneBefore(seqDep, locV);
                    else if (nodes[cliLast]->suc == nodes[cliAf])seqLocV->concatOneBefore(nodes[cliAf]->seqi_n, locV);
                    else seqLocV->concatOneBefore(nodes[cliAf]->seqi_0, locV);
                }
                if (!seqLocV->F)continue;
                flexSeq.front() = seqLocU;
                flexSeq.back() = seqLocV;
                ckRes = seqDep->evaluation(flexSeq);
                if (ckRes < res.sc) {
                    res.sc = ckRes;
                    res.ft = II(locU, locV);
                }
            }
        }        
        delete seqTempU;
        delete seqTempV;        
        return res;
    }

    int evalSpecFlex(vector<SeqData*>& flexSeq, vector<int>& trace) {        
        int* virGiantT = new int[6];
        int curNum = 0;
        int totalLoad = 0;
        assert(flexSeq.size() == 2 || flexSeq.size() == 3);
        SeqData* firstE = flexSeq.front();
        SeqData* lastE = flexSeq.back();
        int cliFirst = pr->listLoc[firstE->lastnode].idxClient;
        int cliBe = -1;
        if (firstE->beforeLaNode != -1)cliBe = pr->listLoc[firstE->beforeLaNode].idxClient;
        int cliLast = pr->listLoc[lastE->firstnode].idxClient;
        int cliAf = -1;
        if (lastE->afterFiNode != -1)cliAf = pr->listLoc[lastE->afterFiNode].idxClient;
        //init giant tour:
        if (cliBe != -1)virGiantT[++curNum] = cliFirst;
        else cliBe = cliFirst;
        if (flexSeq.size() == 3) {
            virGiantT[++curNum] = pr->listLoc[flexSeq[1]->firstnode].idxClient;
            if (flexSeq[1]->firstnode != flexSeq[1]->lastnode) {
                virGiantT[++curNum] = pr->listLoc[flexSeq[1]->lastnode].idxClient;
            }            
        }        
        if (cliAf != -1)virGiantT[++curNum] = cliLast;
        else cliAf = cliLast;
        if (curNum == 0) {
            delete[] virGiantT;
            return 0;
        }        
        //find min:
        int stT = 0, enT = pr->T;
        int resCost = 0;
        if (cliBe) {
            if (nodes[cliFirst]->pred == nodes[cliBe]) {
                resCost += nodes[cliBe]->seq0_i->cost;
                stT = nodes[cliBe]->seq0_i->E;
                totalLoad += nodes[cliBe]->seq0_i->load;
            }
            else {
                resCost += nodes[cliBe]->seqn_i->cost;
                stT = nodes[cliBe]->seqn_i->E;
                totalLoad += nodes[cliBe]->seqn_i->load;
            }
            cliBe = nodes[cliBe]->idxLoc;
        }        
        if (cliAf) {
            if (nodes[cliLast]->suc == nodes[cliAf]) {
                resCost += nodes[cliAf]->seqi_n->cost;
                enT = nodes[cliAf]->seqi_n->L;
                totalLoad += nodes[cliAf]->seqi_n->load;
            }
            else {
                resCost += nodes[cliAf]->seqi_0->cost;
                enT = nodes[cliAf]->seqi_0->L;
                totalLoad += nodes[cliAf]->seqi_0->load;
            }
            cliAf = nodes[cliAf]->idxLoc;
        }

        lstLabel.clear();
        curLabel.clear();
        prvIdLb.clear();
        for (int i = 1; i <= curNum; ++i) {
            F[i] = oo;
            pred[i] = -1;
            totalLoad += pr->listCL[virGiantT[i]].demand;
        }        
        if (totalLoad > pr->Q) {            
            delete[] virGiantT;
            return oo;
        }
        F[0] = oo;
        pred[0] = -1;
        int stLb = -1, enLb = -1;
        III lbU, lbV;
        II label1, label2, label3;
        int idCus, costU, timeU, idLocU;
        int timeV, costV;
        int st = 1;
        //init first label:            
        lstLabel.pb(III(II(0, stT), cliBe));
        prvIdLb.pb(pred[st - 1]);
        stLb = lstLabel.size() - 1;
        enLb = lstLabel.size() - 1;
        for (int v = st; v <= curNum; ++v) {
            int idCus = virGiantT[v];
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
                        //|| timeV + pr->times[idLocV][0] > pr->T // only true when travel times satisfy triangle inequality 
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
                for (int i = 0; i < curLabel.size(); ++i) {
                    III val = curLabel[i];
                    if (val.ft.sc < minT)
                    {
                        minT = val.ft.sc;
                        lstLabel.pb(III(val.ft, idLocV));
                        prvIdLb.pb(val.sc);
                        //update F function (only consider last client)                        
                        if (v == curNum) {
                            if (val.ft.sc + pr->times[idLocV][cliAf] <= enT // uncomment when travel times does not satisfy triangle inequality
                                && F[v] > val.ft.ft + pr->costs[idLocV][cliAf]) {
                                F[v] = val.ft.ft + pr->costs[idLocV][cliAf];
                                pred[v] = lstLabel.size() - 1;
                            }
                        }
                    }
                }
            }
            stLb = enLb + 1;
            enLb = lstLabel.size() - 1;
            if (stLb > enLb)break;
        }                
        if (F[curNum] >= oo) {
            delete[] virGiantT;
            return oo;
        }        
        resCost += F[curNum];
        int indexLb = pred[curNum];// index of last label.
       ///construct solution
       //vector<int> tourLoc;       
        tourLoc.clear();        
        while (indexLb != 0)
        {            
            tourLoc.pb(lstLabel[indexLb].second);
            indexLb = prvIdLb[indexLb];            
        }
        reverse(tourLoc.begin(), tourLoc.end());        
        //change location
        /*for (auto val : tourLoc) {
            changeLocCli(val);
        }*/
        trace = tourLoc;
        delete[] virGiantT;               
        return resCost;
    }
    //split seq:
    bool evalFlexRou(vector<SeqData*> &flexSeq, vector<int>& trace) {                  
        int* virGiantT = new int[routeU->length + 3];
        int curNum = 0;
        //get order cus in route.            
        for (auto curSeq : flexSeq) {            
            if (curSeq == NULL) continue;
            if (curSeq->firstnode == 0 && curSeq->lastnode == 0)continue;
            int curCli = pr->listLoc[curSeq->firstnode].idxClient;
            bool isRev = false;
            if (curSeq->firstnode == curSeq->lastnode) {
                virGiantT[++curNum] = curCli;
                continue;
            }            
            if (curSeq->firstnode) {
                if (nodes[curCli]->pred == nodes[pr->listLoc[curSeq->afterFiNode].idxClient])isRev = true;
            }
            else {
                if (nodes[pr->listLoc[curSeq->lastnode].idxClient]->suc == nodes[pr->listLoc[curSeq->beforeLaNode].idxClient])isRev = true;
            }
            int tst = 1;            
            while (true)
            {
                if(curCli)virGiantT[++curNum] = curCli;                
                if (curCli == pr->listLoc[curSeq->lastnode].idxClient)break;
                if (curCli == 0) {
                    curCli = pr->listLoc[curSeq->afterFiNode].idxClient;
                    continue;
                }
                if(!isRev)curCli = nodes[curCli]->suc->idxClient;
                else curCli = nodes[curCli]->pred->idxClient;
            }
            //for (auto val : curSeq->idxCliNode)virGiantT[++curNum] = val;
        }                
        assert(routeU->length + 2 >= curNum);
        //optimizing route to optimal:
        /*vector<III> lstLabel; //((cost, time), loc)
        vector<III> curLabel; //((cost, time), prv)
        vector<int> prvIdLb; //previous label in vector*/
        lstLabel.clear();
        curLabel.clear();
        prvIdLb.clear();
        //cout << "giant in rou:\n";
        for (int i = 1; i <= curNum; ++i) {
            //cout << virGiantT[i] << " ";
            F[i] = oo;
            pred[i] = -1;
        }
        //cout << "\n";
        F[0] = oo;
        pred[0] = -1;
        int stLb = -1, enLb = -1;
        III lbU, lbV;
        II label1, label2, label3;
        int idCus, costU, timeU, idLocU;
        int timeV, costV;
        int st = 1;
        //init first label:            
        lstLabel.pb(III(II(0, 0), 0));
        prvIdLb.pb(pred[st - 1]);
        stLb = lstLabel.size() - 1;
        enLb = lstLabel.size() - 1;
        for (int v = st; v <= curNum; ++v) {
            int idCus = virGiantT[v];
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
                        //|| timeV + pr->times[idLocV][0] > pr->T // only true when travel times satisfy triangle inequality 
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
                for (int i = 0; i < curLabel.size(); ++i) {
                    III val = curLabel[i];
                    if (val.ft.sc < minT)
                    {
                        minT = val.ft.sc;
                        lstLabel.pb(III(val.ft, idLocV));
                        prvIdLb.pb(val.sc);
                        //update F function (only consider last client)                        
                        if (v == curNum) {
                            if (val.ft.sc + pr->times[idLocV][0] < pr->T&& // uncomment when travel times does not satisfy triangle inequality
                                F[v] > val.ft.ft + pr->costs[idLocV][0]) {
                                F[v] = val.ft.ft + pr->costs[idLocV][0];
                                pred[v] = lstLabel.size() - 1;
                            }
                        }
                    }
                }
            }
            stLb = enLb + 1;
            enLb = lstLabel.size() - 1;
            if (stLb > enLb)break;
        }           
        int oldCost = routeU->depot->seqi_n->cost;              
        if (F[curNum] >= oldCost) {                                       
            delete[] virGiantT;            
            return false;
        }        
        int indexLb = pred[curNum];// index of last label.
        ///construct solution
        //vector<int> tourLoc;       
        tourLoc.clear();
        while (indexLb != 0)
        {
            tourLoc.pb(lstLabel[indexLb].second);            
            indexLb = prvIdLb[indexLb];
        }        
        reverse(tourLoc.begin(), tourLoc.end());        
        //change location
        /*for (auto val : tourLoc) {                           
            changeLocCli(val);
        }*/
        trace = tourLoc;                
        delete[] virGiantT;
        return true;
    }
    //change location of client 
    //input: location used to change
    void changeLocCli(int idLoc) {        
        if (idLoc == 0)return;
        if (nodes[pr->listLoc[idLoc].idxClient]->idxLoc == idLoc)return;
        nodes[pr->listLoc[idLoc].idxClient]->idxLoc = idLoc;
        reinitCluSet(nodes[pr->listLoc[idLoc].idxClient]);        
    }

    int interRouteGeneralInsert() {        
        int iBest = 0, jBest = 0;
        int moveMin;
        // 0 -> send nothing
        // 1 -> send U
        // 2 -> send U,Unext
        // 3 -> send U,U_prev for first index         
        // 3 -> send Unext,U for second index                 
        /*
        * the first index can't be 0 due to the property of granular search
        */
        SeqData* seq = nodeU->seq0_i;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)resGenInMoves[i][j] = oo;
        if (!isFixed) {
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) {                    
                    resSubFlex1[i][j].sc = oo;
                    resSubFlex2[i][j].sc = oo;
                    traceLoc1[i][j].clear();
                    traceLoc2[i][j].clear();
                }
        }
        resGenInMoves[0][0] = routeU->depot->seqi_n->cost + routeV->depot->seqi_n->cost;
        //condition for remain cases
        if (nodeU->idxClient == 0)return 0;
        if (pr->isDebug)cout << "ck inter Insert\n";        
        if (pr->isDebug) {
            cout << nodeU->idxClient << " " << nodeV->idxClient << "\n";
            routeU->showR();
            routeV->showR();
        }
        //1-0                
        myseqs1.clear();
        myseqs1.push_back(uPred->seq0_i);
        myseqs1.push_back(uSuc->seqi_n);
        if (isFixed)resGenInMoves[1][0] = seq->evaluation(myseqs1);
        else {            
            //resSubFlex1[1][0] = evalFlex(myseqs1);
            int testCost1 = evalSpecFlex(myseqs1, traceLoc1[1][0]);                       
            assert(resSubFlex1[1][0].sc >= testCost1);
            resSubFlex1[1][0].sc = testCost1;            
        }                
        myseqs2.clear();
        myseqs2.push_back(vPred->seq0_i);
        myseqs2.push_back(nodeU->seqi_j[0]);
        myseqs2.push_back(nodeV->seqi_n);
        if (isFixed)resGenInMoves[1][0] += seq->evaluation(myseqs2);
        else {
            //resSubFlex2[1][0] = evalFlex(myseqs2);
            int testCost2 = evalSpecFlex(myseqs2, traceLoc2[1][0]);            
            assert(resSubFlex2[1][0].sc >= testCost2);
            resSubFlex2[1][0].sc = testCost2;            
        }        
        //2-0        
        myseqs2.clear();
        myseqs2.push_back(vPred->seq0_i);
        myseqs2.push_back(nodeU->seqi_j[0]);
        myseqs2.push_back(nodeV->seqi_n);
        if (uSuc->idxClient) {
            myseqs1.clear();
            myseqs1.push_back(uPred->seq0_i);
            myseqs1.push_back(uSuc->suc->seqi_n);
            if (isFixed)resGenInMoves[2][0] = seq->evaluation(myseqs1);
            else {
                //resSubFlex1[2][0] = evalFlex(myseqs1);
                int testCost1 = evalSpecFlex(myseqs1, traceLoc1[2][0]);
                assert(resSubFlex1[2][0].sc >= testCost1);
                resSubFlex1[2][0].sc = testCost1;                
            }
            myseqs2[1] = nodeU->seqi_j[1];
            if (isFixed)resGenInMoves[2][0] += seq->evaluation(myseqs2);
            else {
                //resSubFlex2[2][0] = evalFlex(myseqs2);
                int testCost2 = evalSpecFlex(myseqs2, traceLoc2[2][0]);
                assert(resSubFlex2[2][0].sc >= testCost2);
                resSubFlex2[2][0].sc = testCost2;                
            }
        }        
        //3-0        
        if (uPred->idxClient) {
            myseqs1.clear();
            myseqs1.push_back(uPred->pred->seq0_i);
            myseqs1.push_back(uSuc->seqi_n);
            if (isFixed)resGenInMoves[3][0] = seq->evaluation(myseqs1);
            else {
                //resSubFlex1[3][0] = evalFlex(myseqs1);
                int testCost1 = evalSpecFlex(myseqs1, traceLoc1[3][0]);
                assert(resSubFlex1[3][0].sc >= testCost1);
                resSubFlex1[3][0].sc = testCost1;                
            }
            myseqs2[1] = uPred->seqj_i[1];
            if (isFixed)resGenInMoves[3][0] += seq->evaluation(myseqs2);
            else {
                //resSubFlex2[3][0] = evalFlex(myseqs2);
                int testCost2 = evalSpecFlex(myseqs2, traceLoc2[3][0]);
                assert(resSubFlex2[3][0].sc >= testCost2);
                resSubFlex2[3][0].sc = testCost2;                
            }
        }        
        //condition for remain cases
        if (nodeV->idxClient == 0)goto line1;
        //1-1                    
        myseqs1.clear();
        myseqs1.push_back(uPred->seq0_i);
        myseqs1.push_back(nodeV->seqi_j[0]);
        myseqs1.push_back(uSuc->seqi_n);
        if (isFixed)resGenInMoves[1][1] = seq->evaluation(myseqs1);
        else {
            //resSubFlex1[1][1] = evalFlex(myseqs1);
            int testCost1 = evalSpecFlex(myseqs1, traceLoc1[1][1]);
            assert(resSubFlex1[1][1].sc >= testCost1);
            resSubFlex1[1][1].sc = testCost1;            
        }
        myseqs2.clear();
        myseqs2.push_back(vPred->seq0_i);
        myseqs2.push_back(nodeU->seqi_j[0]);
        myseqs2.push_back(vSuc->seqi_n);
        if (isFixed)resGenInMoves[1][1] += seq->evaluation(myseqs2);
        else {
            //resSubFlex2[1][1] = evalFlex(myseqs2);
            int testCost2 = evalSpecFlex(myseqs2, traceLoc2[1][1]);
            assert(resSubFlex2[1][1].sc >= testCost2);
            resSubFlex2[1][1].sc = testCost2;            
        }
        //2-1        
        myseqs2.clear();
        myseqs2.push_back(vPred->seq0_i);
        myseqs2.push_back(nodeU->seqi_j[0]);
        myseqs2.push_back(vSuc->seqi_n);
        if (uSuc->idxClient) {
            myseqs1.clear();
            myseqs1.push_back(uPred->seq0_i);
            myseqs1.push_back(nodeV->seqi_j[0]);
            myseqs1.push_back(uSuc->suc->seqi_n);
            if (isFixed)resGenInMoves[2][1] = seq->evaluation(myseqs1);
            else {
                //resSubFlex1[2][1] = evalFlex(myseqs1);
                int testCost1 = evalSpecFlex(myseqs1, traceLoc1[2][1]);
                assert(resSubFlex1[2][1].sc >= testCost1);
                resSubFlex1[2][1].sc = testCost1;                
            }
            myseqs2[1] = nodeU->seqi_j[1];
            if (isFixed)resGenInMoves[2][1] += seq->evaluation(myseqs2);
            else {
                //resSubFlex2[2][1] = evalFlex(myseqs2);
                int testCost2 = evalSpecFlex(myseqs2, traceLoc2[2][1]);
                assert(resSubFlex2[2][1].sc >= testCost2);
                resSubFlex2[2][1].sc = testCost2;                
            }
        }
        //3-1        
        if (uPred->idxClient) {
            myseqs1.clear();
            myseqs1.push_back(uPred->pred->seq0_i);
            myseqs1.push_back(nodeV->seqi_j[0]);
            myseqs1.push_back(uSuc->seqi_n);
            if (isFixed)resGenInMoves[3][1] = seq->evaluation(myseqs1);
            else {
                //resSubFlex1[3][1] = evalFlex(myseqs1);
                int testCost1 = evalSpecFlex(myseqs1, traceLoc1[3][1]);
                assert(resSubFlex1[3][1].sc >= testCost1);
                resSubFlex1[3][1].sc = testCost1;                
            }
            myseqs2[1] = uPred->seqj_i[1];
            if (isFixed)resGenInMoves[3][1] += seq->evaluation(myseqs2);
            else {
                //resSubFlex2[3][1] = evalFlex(myseqs2);
                int testCost2 = evalSpecFlex(myseqs2, traceLoc2[3][1]);
                assert(resSubFlex2[3][1].sc >= testCost2);
                resSubFlex2[3][1].sc = testCost2;                
            }
        }
        //condition for remain cases
        if (vSuc->idxClient == 0)goto line1;
        //with second indexes are 2 and 3, some concat operations may do at the same time
        //1-2                    
        myseqs1.clear();
        myseqs1.push_back(uPred->seq0_i);
        myseqs1.push_back(nodeV->seqi_j[1]);
        myseqs1.push_back(uSuc->seqi_n);
        if (isFixed)resGenInMoves[1][2] = seq->evaluation(myseqs1);
        else {
            //resSubFlex1[1][2] = evalFlex(myseqs1);
            int testCost1 = evalSpecFlex(myseqs1, traceLoc1[1][2]);
            assert(resSubFlex1[1][2].sc >= testCost1);
            resSubFlex1[1][2].sc = testCost1;            
        }
        myseqs2.clear();
        myseqs2.push_back(vPred->seq0_i);
        myseqs2.push_back(nodeU->seqi_j[0]);
        myseqs2.push_back(vSuc->suc->seqi_n);
        if (isFixed) {
            resGenInMoves[1][3] = seq->evaluation(myseqs2);//same route as 1-3
            resGenInMoves[1][2] += resGenInMoves[1][3];
        }
        else {
            //resSubFlex2[1][3] = evalFlex(myseqs2);
            int testCost2 = evalSpecFlex(myseqs2, traceLoc2[1][3]);
            assert(resSubFlex2[1][3].sc >= testCost2);
            resSubFlex2[1][3].sc = testCost2;            
            resSubFlex2[1][2] = resSubFlex2[1][3];
            traceLoc2[1][2] = traceLoc2[1][3];
        }
        //1-3:            
        myseqs1[1] = nodeV->seqj_i[1];
        if (isFixed)resGenInMoves[1][3] += seq->evaluation(myseqs1);
        else {
            //resSubFlex1[1][3] = evalFlex(myseqs1);
            int testCost1 = evalSpecFlex(myseqs1, traceLoc1[1][3]);
            assert(resSubFlex1[1][3].sc >= testCost1);
            resSubFlex1[1][3].sc = testCost1;            
        }
        //2-2        
        myseqs2.clear();
        myseqs2.push_back(vPred->seq0_i);
        myseqs2.push_back(nodeU->seqi_j[0]);
        myseqs2.push_back(vSuc->suc->seqi_n);
        if (uSuc->idxClient) {
            myseqs1.clear();
            myseqs1.push_back(uPred->seq0_i);
            myseqs1.push_back(nodeV->seqi_j[1]);
            myseqs1.push_back(uSuc->suc->seqi_n);
            if (isFixed)resGenInMoves[2][2] = seq->evaluation(myseqs1);
            else {
                //resSubFlex1[2][2] = evalFlex(myseqs1);
                int testCost1 = evalSpecFlex(myseqs1, traceLoc1[2][2]);
                assert(resSubFlex1[2][2].sc >= testCost1);
                resSubFlex1[2][2].sc = testCost1;                
            }
            myseqs2[1] = nodeU->seqi_j[1];
            if (isFixed) {
                resGenInMoves[2][3] = seq->evaluation(myseqs2);//same route as 2-3
                resGenInMoves[2][2] += resGenInMoves[2][3];
            }
            else {
                //resSubFlex2[2][3] = evalFlex(myseqs2);
                int testCost2 = evalSpecFlex(myseqs2, traceLoc2[2][3]);
                assert(resSubFlex2[2][3].sc >= testCost2);
                resSubFlex2[2][3].sc = testCost2;                
                resSubFlex2[2][2] = resSubFlex2[2][3];
                traceLoc2[2][2] = traceLoc2[2][3];
            }
            //2-3:
            myseqs1[1] = nodeV->seqj_i[1];
            if (isFixed)resGenInMoves[2][3] += seq->evaluation(myseqs1);
            else {
                //resSubFlex1[2][3] = evalFlex(myseqs1);
                int testCost1 = evalSpecFlex(myseqs1, traceLoc1[2][3]);
                assert(resSubFlex1[2][3].sc >= testCost1);
                resSubFlex1[2][3].sc = testCost1;
            }
        }
        //3-2
        if (uPred->idxClient) {
            myseqs1.clear();
            myseqs1.push_back(uPred->pred->seq0_i);
            myseqs1.push_back(nodeV->seqi_j[1]);
            myseqs1.push_back(uSuc->seqi_n);
            if (isFixed)resGenInMoves[3][2] = seq->evaluation(myseqs1);
            else {
                //resSubFlex1[3][2] = evalFlex(myseqs1);
                int testCost1 = evalSpecFlex(myseqs1, traceLoc1[3][2]);
                assert(resSubFlex1[3][2].sc >= testCost1);
                resSubFlex1[3][2].sc = testCost1;
            }
            myseqs2[1] = uPred->seqj_i[1];
            if (isFixed) {
                resGenInMoves[3][3] = seq->evaluation(myseqs2);//same route as 3-3
                resGenInMoves[3][2] += resGenInMoves[3][3];
            }
            else {
                //resSubFlex2[3][3] = evalFlex(myseqs2);
                int testCost2 = evalSpecFlex(myseqs2, traceLoc2[3][3]);
                assert(resSubFlex2[3][3].sc >= testCost2);
                resSubFlex2[3][3].sc = testCost2;
                resSubFlex2[3][2] = resSubFlex2[3][3];
                traceLoc2[3][2] = traceLoc2[3][3];
            }
            //3-3:
            myseqs1[1] = nodeV->seqj_i[1];
            if (isFixed)resGenInMoves[3][3] += seq->evaluation(myseqs1);
            else {
                //resSubFlex1[3][3] = evalFlex(myseqs1);
                int testCost1 = evalSpecFlex(myseqs1, traceLoc1[3][3]);
                assert(resSubFlex1[3][3].sc >= testCost1);
                resSubFlex1[3][3].sc = testCost1;
            }
        }
    line1:
        moveMin = oo;
        iBest = 0; jBest = 0;
        //deal with flexible version:
        if (!isFixed) {
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) {
                    if (i == 0 && j == 0)continue;
                    resGenInMoves[i][j] = resSubFlex1[i][j].sc + resSubFlex2[i][j].sc;
                }
        }
        //find best move:
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (resGenInMoves[i][j] < moveMin) {
                    moveMin = resGenInMoves[i][j];
                    iBest = i;
                    jBest = j;
                }
            }
        }        
        if (pr->isDebug) {
            cout << iBest << " " << jBest << "\n";
        }
        if (iBest == 0 && jBest == 0)
            return 0;
        //reinitSingleMoveInRou(routeU);
        //reinitSingleMoveInRou(routeV);                
        //Moving nodes in each changed route
        if (!isFixed) {
            //change location when using flexible verion:
            //change for routeU:
            /*changeLocCli(resSubFlex1[iBest][jBest].ft.ft);
            changeLocCli(resSubFlex1[iBest][jBest].ft.sc);*/
            for (auto val : traceLoc1[iBest][jBest])changeLocCli(val);
            //change for routeV:
            /*changeLocCli(resSubFlex2[iBest][jBest].ft.ft);
            changeLocCli(resSubFlex2[iBest][jBest].ft.sc);*/
            for (auto val : traceLoc2[iBest][jBest])changeLocCli(val);
        }
        Node* placeU = nodeU->pred;
        if (iBest == 3)placeU = uPred->pred;

        if (iBest)insertNode(nodeU, vPred);
        if (iBest == 2)insertNode(uSuc, nodeU);
        if (iBest == 3)insertNode(uPred, nodeU);

        if (jBest == 1 || jBest == 2)insertNode(nodeV, placeU);
        if (jBest == 2)insertNode(vSuc, nodeV);
        if (jBest == 3) {
            insertNode(vSuc, placeU);
            insertNode(nodeV, vSuc);
        }
        if (pr->isDebug) {
            cout << iBest << " " << jBest << "\n";//2 0
            cout << nodeU->idxClient << " " << nodeV->idxClient << "\n";//4 111
            cout << resSubFlex1[iBest][jBest].sc << " " << resSubFlex2[iBest][jBest].sc << "\n";
            cout << resGenInMoves[iBest][jBest] << "\n";//289
        }
        //update route data
        routeU->updateRoute();
        routeV->updateRoute();
        if (pr->isDebug) {
            routeU->showR();
            routeU->showRLoc();
            routeV->showR();
            routeV->showRLoc();
            cout << routeU->caculateDis() << " " << routeV->caculateDis() << "\n";
        }

        //if (nodeU->idxClient == 80 && nodeV->idxClient == 108)goto line2;
        cost += resGenInMoves[iBest][jBest] - resGenInMoves[0][0];
        //count[iBest][jBest]++;
        if (pr->isTurnCkSol) {
            try {
                ckSol();
            }
            catch (const char* msg) {
                cerr << msg << endl;
                for (int i = 1; i <= n; ++i)pr->fileOut << giantT[i] << ", ";
                pr->fileOut << "\n";
                for (auto val : ordNodeLs)pr->fileOut << val << ", ";
                pr->fileOut.close();
                system("pause");
                exit(0);                
            }
        }
        return 1;
        //line2:
        //return 0;
    }

    
    
    void addSeqInPieces(Node* st, Node* en, vector<SeqData*>& myseqs) {
        //st and en are in the same route.        
        //st and en can't be depot
        assert(st->idxClient && en->idxClient);
        if (st->posInRoute > en->posInRoute)return;
        int disInR = -1;
        Node* val = st;
        while (true)
        {
            disInR = en->posInRoute - val->posInRoute;
            if (disInR + 1 <= pr->sizeSub) {
                myseqs.push_back(val->seqi_j[disInR]);
                break;
            }
            SeqData* curSeq = val->seqi_j[pr->sizeSub - 1];
            myseqs.push_back(curSeq);
            val = nodes[pr->listLoc[curSeq->lastnode].idxClient]->suc;
        }
    }

    void addRevSeqInPieces(Node* st, Node* en, vector<SeqData*>& myseqs) {
        //st and en are in the same route.        
        //st and en can't be depot
        int preBe = myseqs.size();
        assert(st->idxClient && en->idxClient);
        if (st->posInRoute > en->posInRoute)return;
        int disInR = -1;
        Node* val = st;
        while (true)
        {
            disInR = en->posInRoute - val->posInRoute;
            if (disInR + 1 <= pr->sizeSub) {
                myseqs.push_back(val->seqj_i[disInR]);
                break;
            }
            SeqData* curSeq = val->seqj_i[pr->sizeSub - 1];
            myseqs.push_back(curSeq);
            val = nodes[pr->listLoc[curSeq->firstnode].idxClient]->suc;
        }
        reverse(myseqs.begin() + preBe, myseqs.end());
    }

    int intraRouteGeneralInsert() {      
        //nodeV in this case can be depot (arrival depot)
        int difDir = nodeV->posInRoute - nodeU->posInRoute;
        if (difDir >= 0 && difDir <= 1) {// move in this case is useless
            return 0;
        }
        //initialization:
        int iBest = 0, jBest = 0;
        int moveMin;
        // 0 -> send nothing
        // 1 -> send U
        // 2 -> send U,Unext
        // 3 -> send U,U_prev for first index         
        // 3 -> send Unext,U for second index         
        // the case with index 3 may reverse when isTurn is true
        /*
        * the first index can't be 0 due to the property of granular search
        */
        SeqData* seq = nodeU->seq0_i;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)resGenInMoves[i][j] = oo;
        if (!isFixed) {
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)traceLoc[i][j].clear();
        }
        resGenInMoves[0][0] = routeU->depot->seqi_n->cost;
        if (difDir == -1) {
            //V->U
            //only two case need to consider: 
            //1-1: intra2opt will deal with this case
            //2-0: 
            if (uSuc->idxClient && nodeV->idxClient) {
                myseqs1.clear();
                myseqs1.push_back(vPred->seq0_i);
                myseqs1.push_back(nodeU->seqi_j[1]);
                myseqs1.push_back(nodeV->seqi_j[0]);
                myseqs1.push_back(uSuc->suc->seqi_n);
                if(isFixed)resGenInMoves[0][2] = seq->evaluation(myseqs1);
                else {
                    if (evalFlexRou(myseqs1, traceLoc[0][2]))resGenInMoves[0][2] = F[traceLoc[0][2].size()];
                }
            }                                    
        }

        // If difDir < 0, this means that V is on the left of U in the route
        // We don't want to deal with this case, so we simply inverse the notations of U and V 
        // And we deal in the following with only the case where V is on the right
        Node* tempU = nodeU;
        Node* tempV = nodeV;
        bool isTurn = false;
        if (difDir < 0)
        {
            nodeU = tempV;
            nodeV = tempU;
            uSuc = nodeU->suc;
            vSuc = nodeV->suc;
            uPred = nodeU->pred;
            vPred = nodeV->pred;
            difDir = -difDir;
            isTurn = true;
        }        

        /*check in detail*/        
        if (pr->isDebug)cout << "ck intra Insert\n";                
        if (difDir == 1)goto line1;
        if (difDir == 2) {
            //U X V Y
            //only consider 1-2, 3-2. the remain is useless or same as 2-opt.
            //useless in here is based on granular search
            if (!isTurn) {
                //1-2:
                if (nodeV->idxClient == 0 || vSuc->idxClient == 0)goto line1;
                if(nodeU->idxClient == 0)goto line1;      
                myseqs1.clear();
                myseqs1.push_back(uPred->seq0_i);
                myseqs1.push_back(nodeV->seqi_j[1]);                
                myseqs1.push_back(nodeU->seqj_i[1]);
                myseqs1.push_back(vSuc->suc->seqi_n);
                if(isFixed)resGenInMoves[1][2] = seq->evaluation(myseqs1);
                else {
                    if (evalFlexRou(myseqs1, traceLoc[1][2]))resGenInMoves[1][2] = F[traceLoc[1][2].size()];
                }
                //3-2:
                if (uPred->idxClient == 0)goto line1;
                myseqs1[0] = uPred->pred->seq0_i;
                myseqs1[2] = uPred->seqj_i[2];
                if(isFixed)resGenInMoves[3][2] = seq->evaluation(myseqs1);
                else {
                    if (evalFlexRou(myseqs1, traceLoc[3][2]))resGenInMoves[3][2] = F[traceLoc[3][2].size()];
                }
            }                        
            //V Y U X -> U X V Y
            //3 for first is U_next->U and second is U->U_pred
            //the second index is always > 0
            else {   
                // in this case U and V can't be depot based on foor loop feature.
                if (!(nodeU->idxClient != 0 && nodeV->idxClient != 0)) {
                    throw "bug intra 1";
                }
                //0-1 UVYX -> VUXY                
                myseqs1.clear();
                myseqs1.push_back(uPred->seq0_i);
                myseqs1.push_back(nodeV->seqi_j[0]);
                myseqs1.push_back(nodeU->seqi_j[1]);                
                myseqs1.push_back(vSuc->seqi_n);
                if(isFixed)resGenInMoves[0][1] = seq->evaluation(myseqs1);
                else {
                    if (evalFlexRou(myseqs1, traceLoc[0][1]))resGenInMoves[0][1] = F[traceLoc[0][1].size()];
                }
                //1-1 UYVX -> VXUY
                myseqs1[2] = nodeU->seqj_i[1];
                if(isFixed)resGenInMoves[1][1] = seq->evaluation(myseqs1);                
                else {
                    if (evalFlexRou(myseqs1, traceLoc[1][1]))resGenInMoves[1][1] = F[traceLoc[1][1].size()];
                }
                //0-2 UXVY -> VYUX
                if(vSuc->idxClient == 0)goto line1;
                myseqs1[1] = nodeV->seqi_j[1];
                myseqs1[2] = nodeU->seqi_j[1];
                myseqs1[3] = vSuc->suc->seqi_n;
                if(isFixed)resGenInMoves[0][2] = seq->evaluation(myseqs1);
                else {
                    if (evalFlexRou(myseqs1, traceLoc[0][2]))resGenInMoves[0][2] = F[traceLoc[0][2].size()];
                }
                /*
                    0-3// 2opt                
                    1-2// 2opt
                    1-3 same with 1-1
                    2-1 same with 0-1
                    2-2 same with 0-2
                    2-3 // fault with this case
                */
            }
        }
        else {            
            //UXVY
            if (!isTurn) {
                if (!(nodeU->idxClient != 0)) {
                    throw "bug intra 2";
                }
                //1-0             
                myseqs1.clear();
                myseqs1.push_back(uPred->seq0_i);
                addSeqInPieces(uSuc, vPred, myseqs1);
                myseqs1.push_back(nodeU->seqi_j[0]);
                myseqs1.push_back(nodeV->seqi_n);
                if(isFixed)resGenInMoves[1][0] = seq->evaluation(myseqs1);
                else {
                    if (evalFlexRou(myseqs1, traceLoc[1][0]))resGenInMoves[1][0] = F[traceLoc[1][0].size()];
                }
                //2-0
                if (uSuc->idxClient) {
                    myseqs1.clear();
                    myseqs1.push_back(uPred->seq0_i);
                    addSeqInPieces(uSuc->suc, vPred, myseqs1);
                    myseqs1.push_back(nodeU->seqi_j[1]);
                    myseqs1.push_back(nodeV->seqi_n);
                    if(isFixed)resGenInMoves[2][0] = seq->evaluation(myseqs1);
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[2][0]))resGenInMoves[2][0] = F[traceLoc[2][0].size()];
                    }
                }
                //3-0
                if (uPred->idxClient) {
                    myseqs1.clear();
                    myseqs1.push_back(uPred->pred->seq0_i);
                    addSeqInPieces(uSuc, vPred, myseqs1);
                    myseqs1.push_back(uPred->seqj_i[1]);
                    myseqs1.push_back(nodeV->seqi_n);
                    if(isFixed)resGenInMoves[3][0] = seq->evaluation(myseqs1);
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[3][0]))resGenInMoves[3][0] = F[traceLoc[3][0].size()];
                    }
                }
                if (nodeV->idxClient == 0) goto line1;
                //1-1
                myseqs1.clear();                
                myseqs1.push_back(uPred->seq0_i);
                myseqs1.push_back(nodeV->seqi_j[0]);
                addSeqInPieces(uSuc, vPred, myseqs1);
                myseqs1.push_back(nodeU->seqi_j[0]);
                myseqs1.push_back(vSuc->seqi_n);
                if(isFixed)resGenInMoves[1][1] = seq->evaluation(myseqs1);
                else {
                    if (evalFlexRou(myseqs1, traceLoc[1][1]))resGenInMoves[1][1] = F[traceLoc[1][1].size()];
                }
                //2-1
                if (uSuc->idxClient) {
                    myseqs1.clear();
                    myseqs1.push_back(uPred->seq0_i);
                    myseqs1.push_back(nodeV->seqi_j[0]);
                    addSeqInPieces(uSuc->suc, vPred, myseqs1);
                    myseqs1.push_back(nodeU->seqi_j[1]);
                    myseqs1.push_back(vSuc->seqi_n);
                    if(isFixed)resGenInMoves[2][1] = seq->evaluation(myseqs1);
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[2][1]))resGenInMoves[2][1] = F[traceLoc[2][1].size()];
                    }
                }
                //3-1
                if (uPred->idxClient) {
                    myseqs1.clear();
                    myseqs1.push_back(uPred->pred->seq0_i);
                    myseqs1.push_back(nodeV->seqi_j[0]);
                    addSeqInPieces(uSuc, vPred, myseqs1);
                    myseqs1.push_back(uPred->seqj_i[1]);
                    myseqs1.push_back(vSuc->seqi_n);
                    if(isFixed)resGenInMoves[3][1] = seq->evaluation(myseqs1);
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[3][1]))resGenInMoves[3][1] = F[traceLoc[3][1].size()];
                    }
                }    
                if (vSuc->idxClient == 0)goto line1;
                //1-2
                myseqs1.clear();
                myseqs1.push_back(uPred->seq0_i);
                myseqs1.push_back(nodeV->seqi_j[1]);
                addSeqInPieces(uSuc, vPred, myseqs1);
                myseqs1.push_back(nodeU->seqi_j[0]);
                myseqs1.push_back(vSuc->suc->seqi_n);
                if(isFixed)resGenInMoves[1][2] = seq->evaluation(myseqs1);
                else {
                    if (evalFlexRou(myseqs1, traceLoc[1][2]))resGenInMoves[1][2] = F[traceLoc[1][2].size()];
                }
                //1-3 is almost same
                myseqs1[1] = nodeV->seqj_i[1];                
                if(isFixed)resGenInMoves[1][3] = seq->evaluation(myseqs1);
                else {
                    if (evalFlexRou(myseqs1, traceLoc[1][3]))resGenInMoves[1][3] = F[traceLoc[1][3].size()];
                }
                //2-2
                if (uSuc->idxClient) {
                    myseqs1.clear();
                    myseqs1.push_back(uPred->seq0_i);
                    myseqs1.push_back(nodeV->seqi_j[1]);
                    addSeqInPieces(uSuc->suc, vPred, myseqs1);
                    myseqs1.push_back(nodeU->seqi_j[1]);
                    myseqs1.push_back(vSuc->suc->seqi_n);
                    if(isFixed)resGenInMoves[2][2] = seq->evaluation(myseqs1);                    
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[2][2]))resGenInMoves[2][2] = F[traceLoc[2][2].size()];
                    }
                    //2-3 is almost same
                    myseqs1[1] = nodeV->seqj_i[1];
                    if(isFixed)resGenInMoves[2][3] = seq->evaluation(myseqs1);
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[2][3]))resGenInMoves[2][3] = F[traceLoc[2][3].size()];
                    }
                }
                //3-2   
                if (uPred->idxClient) {
                    myseqs1.clear();
                    myseqs1.push_back(uPred->pred->seq0_i);
                    myseqs1.push_back(nodeV->seqi_j[1]);
                    addSeqInPieces(uSuc, vPred, myseqs1);
                    myseqs1.push_back(uPred->seqj_i[1]);
                    myseqs1.push_back(vSuc->suc->seqi_n);
                    if(isFixed)resGenInMoves[3][2] = seq->evaluation(myseqs1);
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[3][2]))resGenInMoves[3][2] = F[traceLoc[3][2].size()];
                    }
                    //3-3 is almost same
                    myseqs1[1] = nodeV->seqj_i[1];
                    if(isFixed)resGenInMoves[3][3] = seq->evaluation(myseqs1);
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[3][3]))resGenInMoves[3][3] = F[traceLoc[3][3].size()];
                    }
                }
            }
            else
            {
                //VYUX->UXVY
                //second index is always >0
                //0-1 
                if (!(nodeV->idxClient != 0)) {
                    throw "bug intra 3";
                }
                myseqs1.clear();
                myseqs1.push_back(uPred->seq0_i);
                myseqs1.push_back(nodeV->seqi_j[0]);
                addSeqInPieces(nodeU, vPred, myseqs1);
                myseqs1.push_back(vSuc->seqi_n);
                if(isFixed)resGenInMoves[0][1] = seq->evaluation(myseqs1);                
                else {
                    if (evalFlexRou(myseqs1, traceLoc[0][1]))resGenInMoves[0][1] = F[traceLoc[0][1].size()];
                }
                //0-2
                if (vSuc->idxClient) {
                    myseqs1.clear();
                    myseqs1.push_back(uPred->seq0_i);
                    myseqs1.push_back(nodeV->seqi_j[1]);
                    addSeqInPieces(nodeU, vPred, myseqs1);
                    myseqs1.push_back(vSuc->suc->seqi_n);
                    if(isFixed)resGenInMoves[0][2] = seq->evaluation(myseqs1);
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[0][2]))resGenInMoves[0][2] = F[traceLoc[0][2].size()];
                    }
                }
                //0-3:
                if (vPred->idxClient) {
                    myseqs1.clear();
                    myseqs1.push_back(uPred->seq0_i);
                    myseqs1.push_back(vPred->seqj_i[1]);
                    addSeqInPieces(nodeU, vPred->pred, myseqs1);
                    myseqs1.push_back(vSuc->seqi_n);
                    if(isFixed)resGenInMoves[0][3] = seq->evaluation(myseqs1);
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[0][3]))resGenInMoves[0][3] = F[traceLoc[0][3].size()];
                    }
                }
                //1-1:
                if (nodeU->idxClient == 0)goto line1;
                myseqs1.clear();
                myseqs1.push_back(uPred->seq0_i);
                myseqs1.push_back(nodeV->seqi_j[0]);
                addSeqInPieces(uSuc, vPred, myseqs1);
                myseqs1.push_back(nodeU->seqi_j[0]);
                myseqs1.push_back(vSuc->seqi_n);
                if(isFixed)resGenInMoves[1][1] = seq->evaluation(myseqs1);
                else {
                    if (evalFlexRou(myseqs1, traceLoc[1][1]))resGenInMoves[1][1] = F[traceLoc[1][1].size()];
                }
                //1-2:
                if (vSuc->idxClient) {
                    myseqs1.clear();
                    myseqs1.push_back(uPred->seq0_i);
                    myseqs1.push_back(nodeV->seqi_j[1]);
                    addSeqInPieces(uSuc, vPred, myseqs1);
                    myseqs1.push_back(nodeU->seqi_j[0]);
                    myseqs1.push_back(vSuc->suc->seqi_n);
                    if(isFixed)resGenInMoves[1][2] = seq->evaluation(myseqs1);
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[1][2]))resGenInMoves[1][2] = F[traceLoc[1][2].size()];
                    }
                }
                //1-3:
                if (vPred->idxClient) {
                    myseqs1.clear();
                    myseqs1.push_back(uPred->seq0_i);
                    myseqs1.push_back(vPred->seqj_i[1]);
                    addSeqInPieces(uSuc, vPred->pred, myseqs1);
                    myseqs1.push_back(nodeU->seqi_j[0]);
                    myseqs1.push_back(vSuc->seqi_n);
                    if(isFixed)resGenInMoves[1][3] = seq->evaluation(myseqs1);
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[1][3]))resGenInMoves[1][3] = F[traceLoc[1][3].size()];
                    }
                }
                if (uSuc->idxClient == 0)goto line1;
                //2-1
                myseqs1.clear();
                myseqs1.clear();
                myseqs1.push_back(uPred->seq0_i);
                myseqs1.push_back(nodeV->seqi_j[0]);
                addSeqInPieces(uSuc->suc, vPred, myseqs1);
                myseqs1.push_back(nodeU->seqi_j[1]);
                myseqs1.push_back(vSuc->seqi_n);
                if(isFixed)resGenInMoves[2][1] = seq->evaluation(myseqs1);
                else {
                    if (evalFlexRou(myseqs1, traceLoc[2][1]))resGenInMoves[2][1] = F[traceLoc[2][1].size()];
                }
                //3-1 is almost same
                myseqs1[myseqs1.size() - 2] = nodeU->seqj_i[1];
                if(isFixed)resGenInMoves[3][1] = seq->evaluation(myseqs1);
                else {
                    if (evalFlexRou(myseqs1, traceLoc[3][1]))resGenInMoves[3][1] = F[traceLoc[3][1].size()];
                }
                //2-2
                if (vSuc->idxClient) {
                    myseqs1.clear();
                    myseqs1.push_back(uPred->seq0_i);
                    myseqs1.push_back(nodeV->seqi_j[1]);
                    addSeqInPieces(uSuc->suc, vPred, myseqs1);
                    myseqs1.push_back(nodeU->seqi_j[1]);
                    myseqs1.push_back(vSuc->suc->seqi_n);
                    if(isFixed)resGenInMoves[2][2] = seq->evaluation(myseqs1);
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[2][2]))resGenInMoves[2][2] = F[traceLoc[2][2].size()];
                    }
                    //3-2 is almost same:
                    myseqs1[myseqs1.size() - 2] = nodeU->seqj_i[1];
                    if(isFixed)resGenInMoves[3][2] = seq->evaluation(myseqs1);
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[3][2]))resGenInMoves[3][2] = F[traceLoc[3][2].size()];
                    }
                }
                //2-3
                if (vPred->idxClient) {
                    myseqs1.clear();
                    myseqs1.push_back(uPred->seq0_i);
                    myseqs1.push_back(vPred->seqj_i[1]);
                    addSeqInPieces(uSuc->suc, vPred->pred, myseqs1);
                    myseqs1.push_back(nodeU->seqi_j[1]);
                    myseqs1.push_back(vSuc->seqi_n);
                    if(isFixed)resGenInMoves[2][3] = seq->evaluation(myseqs1);                    
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[2][3]))resGenInMoves[2][3] = F[traceLoc[2][3].size()];
                    }
                    //3-3 is almost same:
                    myseqs1[myseqs1.size() - 2] = nodeU->seqj_i[1];
                    if(isFixed)resGenInMoves[3][3] = seq->evaluation(myseqs1);
                    else {
                        if (evalFlexRou(myseqs1, traceLoc[3][3]))resGenInMoves[3][3] = F[traceLoc[3][3].size()];
                    }
              }
            }
        }

        line1:
        moveMin = oo;
        iBest = 0; jBest = 0;        
        //find best move:
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (resGenInMoves[i][j] < moveMin) {
                    moveMin = resGenInMoves[i][j];
                    iBest = i;
                    jBest = j;
                }
            }
        }

        if (iBest == 0 && jBest == 0) {
            nodeU = tempU;
            nodeV = tempV;
            uPred = nodeU->pred;
            vPred = nodeV->pred;
            uSuc = nodeU->suc;
            vSuc = nodeV->suc;
            return 0;
        }
        if (pr->isDebug) {
            routeU->showR();
            routeU->showRLoc();
            cout << resGenInMoves[0][0] << "\n";
            cout << iBest << " " << jBest << " " << difDir << "\n";
            cout << nodeU->idxClient << " " << nodeV->idxClient << "\n";
            cout << resGenInMoves[iBest][jBest] << "\n";
            cout << boolalpha << isTurn << "\n";
        }
        if (!isFixed) {
            //cout << traceLoc[iBest][jBest].size() << "\n";
            for (auto val : traceLoc[iBest][jBest])changeLocCli(val);
        }
        reinitSingleMoveInRou(routeU);
        Node* placeU = uPred;
        Node* placeV = vPred;
        if (isTurn && jBest == 3)placeV = vPred->pred;
        if (!isTurn && iBest == 3)placeU = uPred->pred;       
        if (difDir >= 4 || (difDir == 3 && !isTurn)) {               
            if (iBest == 1 || iBest == 2)insertNode(nodeU, placeV);
            if (iBest == 2)insertNode(uSuc, nodeU);
            if (iBest == 3) {
                if (isTurn) {
                    insertNode(uSuc, placeV);
                    insertNode(nodeU, uSuc);
                }
                else {
                    insertNode(nodeU, placeV);
                    insertNode(uPred, nodeU);
                }
            }

            if (jBest == 1 || jBest == 2)insertNode(nodeV, placeU);
            if (jBest == 2)insertNode(vSuc, nodeV);
            if (jBest == 3) {
                if (isTurn) {
                    insertNode(nodeV, placeU);
                    insertNode(vPred, nodeV);
                }
                else {
                    insertNode(vSuc, placeU);
                    insertNode(nodeV, vSuc);
                }
            }
        }
        else if (difDir == 3) {            
              if (jBest == 1 || jBest == 2)insertNode(nodeV, placeU);                
              if (jBest == 2)insertNode(vSuc, nodeV);

              if (jBest == 3) {
                  insertNode(nodeV, placeU);
                  insertNode(vPred, nodeV);
                  if (iBest == 1 || iBest == 3)insertNode(nodeU, uSuc);
              }
              else {
                  if (iBest == 1 || iBest == 2)insertNode(nodeU, placeV);
                  if (iBest == 2)insertNode(uSuc, nodeU);
                  if (iBest == 3) {
                      insertNode(uSuc, placeV);
                      insertNode(nodeU, uSuc);
                  }
              }
        }
        else if (difDir == 2) {
            if (isTurn) {
                if (jBest)insertNode(nodeV, placeU);
                if (jBest == 2)insertNode(vSuc, nodeV);
                if (iBest == 1)insertNode(nodeU, placeV);
            }
            else
            {
                insertNode(nodeV, placeU);
                insertNode(vSuc, nodeV);
                insertNode(nodeU, placeV);
                if (iBest == 3)insertNode(uPred, nodeU);
            }
        }
        else {
            //difDir == 1:
            //only one case 0-2 is considered:
            insertNode(nodeV, placeU);
            insertNode(vSuc, nodeV);
        }

        //count[iBest][jBest]++;
        routeU->updateRoute();        
        nodeU = tempU;
        nodeV = tempV;
        uPred = nodeU->pred;
        vPred = nodeV->pred;
        uSuc = nodeU->suc;
        vSuc = nodeV->suc;
        cost += resGenInMoves[iBest][jBest] - resGenInMoves[0][0];        
        if (pr->isTurnCkSol) {
            try {
                ckSol();
            }
            catch (const char* msg) {
                cerr << msg << endl;
                for (int i = 1; i <= n; ++i)pr->fileOut << giantT[i] << ", ";
                pr->fileOut << "\n";
                for (auto val : ordNodeLs)pr->fileOut << val << ", ";
                pr->fileOut.close();
                system("pause");
                exit(0);                
            }
        }
        return 1;       
    }        

    int interRoute2Opt() {                
        SeqData* seq = nodeU->seq0_i;
        if (uPred->idxClient == 0 && nodeV->idxClient == 0)return 0;
        int oldCost = routeU->depot->seqi_n->cost + routeV->depot->seqi_n->cost;
        int newCost;
        III flexCost1;
        III flexCost2;
        //route 1:            
        myseqs1.clear();
        myseqs1.push_back(uPred->seq0_i);
        myseqs1.push_back(vSuc->seqi_n);        
        //route 2:            
        myseqs2.clear();
        myseqs2.push_back(nodeV->seq0_i);
        myseqs2.push_back(nodeU->seqi_n);
        if (isFixed) {            
            newCost = seq->evaluation(myseqs1) + seq->evaluation(myseqs2);
        }
        else {
            flexCost1 = evalFlex(myseqs1);
            flexCost2 = evalFlex(myseqs2);
            newCost = flexCost1.sc + flexCost2.sc;
        }   
        if (pr->isDebug)cout << "ck inter route 2Opt\n";
        //check changed cost:
        if (oldCost <= newCost)return 0;
        if (pr->isDebug) {
            cout << "improved: " << newCost << "\n";            
        }
        //improvement change
        reinitSingleMoveInRou(routeU);
        reinitSingleMoveInRou(routeV);
        if (!isFixed) {
            //change location for routeU:
            changeLocCli(flexCost1.ft.ft);
            changeLocCli(flexCost1.ft.sc);
            //change location for routeV:
            changeLocCli(flexCost2.ft.ft);
            changeLocCli(flexCost2.ft.sc);            
        }
        Node* lastNodeU = routeU->depot->pred->pred;// last client in routeU
        Node* lastNodeV = routeV->depot->pred->pred;// last client in routeV
        //changed for route 1: uPred->vSuc        

        if (vSuc->idxClient) {            
            uPred->suc = vSuc;
            vSuc->pred = uPred;
            lastNodeV->suc = routeU->depot->pred;
            routeU->depot->pred->pred = lastNodeV;
        }
        else {
            uPred->suc = routeU->depot->pred;
            routeU->depot->pred->pred = uPred;
        }
        //changed for route 2: nodeV->nodeU        
        if (nodeU->idxClient) {
            nodeV->suc = nodeU;
            nodeU->pred = nodeV;
            lastNodeU->suc = routeV->depot->pred;
            routeV->depot->pred->pred = lastNodeU;
        }
        else {
            nodeV->suc = routeV->depot->pred;
            routeV->depot->pred->pred = nodeV;
        }        
        //update route data:
        routeU->updateRoute();
        routeV->updateRoute();        
        //reset infor for futher search
        uPred = nodeU->pred;
        uSuc = nodeU->suc;
        routeU = nodeU->rou;
        cost += newCost - oldCost;
        //count[1][0]++;        
        if (pr->isTurnCkSol) {
            try {
                ckSol();
            }
            catch (const char* msg) {
                cerr << msg << endl;
                for (int i = 1; i <= n; ++i)pr->fileOut << giantT[i] << ", ";
                pr->fileOut << "\n";
                for (auto val : ordNodeLs)pr->fileOut << val << ", ";
                pr->fileOut.close();
                system("pause");
                exit(0);                
            }
        }
        return 1;        
    }

    //2opt*-Inv
    int interRoute2OptInv() {        
        SeqData* seq = nodeU->seq0_i;
        int oldCost = routeU->depot->seqi_n->cost + routeV->depot->seqi_n->cost;
        int newCost = 0;
        int valCost;
        int valCostRev;
        III flexValCost;
        III flexCost1;
        III flexCost2;        
        bool reverseRouteU = false, reverseRouteV = false;                
        ///routeU        
        myseqs1.clear();
        myseqs1.push_back(vSuc->seqn_i);
        myseqs1.push_back(uSuc->seqi_n);
        if (isFixed)valCost = seq->evaluation(myseqs1);
        else {
            flexCost1 = evalFlex(myseqs1);
            valCost = flexCost1.sc;
        }        
        //reverse it:
        myseqs1[0] = uSuc->seqn_i;
        myseqs1[1] = vSuc->seqi_n;
        if (isFixed)valCostRev = seq->evaluation(myseqs1);
        else {
            flexValCost = evalFlex(myseqs1);
            valCostRev = flexValCost.sc;
        }        
        //check what is better:
        if (valCost < valCostRev)newCost += valCost;
        else {
            newCost += valCostRev; 
            flexCost1 = flexValCost;
            reverseRouteU = true;
        }        
        
        ///routeV
        myseqs1.clear();
        myseqs1.push_back(nodeV->seq0_i);
        myseqs1.push_back(nodeU->seqi_0);
        valCost = seq->evaluation(myseqs1);
        if (isFixed)valCost = seq->evaluation(myseqs1);
        else {
            flexCost2 = evalFlex(myseqs1);
            valCost = flexCost2.sc;
        }        
        //reverse it:
        myseqs1[0] = nodeU->seq0_i;
        myseqs1[1] = nodeV->seqi_0;
        if (isFixed)valCostRev = seq->evaluation(myseqs1);
        else {
            flexValCost = evalFlex(myseqs1);
            valCostRev = flexValCost.sc;
        }
        //check what is better:
        if (valCost < valCostRev)newCost += valCost;
        else {
            newCost += valCostRev;
            flexCost2 = flexValCost;
            reverseRouteV = true;
        }        
        ///check whether has improvement:
        if (oldCost <= newCost)return 0;
        if(pr->isDebug) {
            cout << "check inter route 2Opt Inv\n";
            cout << nodeU->idxClient << " " << nodeV->idxClient << "\n";
            routeU->showR();
            routeU->showRLoc();
            routeV->showR();
            routeV->showRLoc();
            cout << boolalpha << reverseRouteU << " " << reverseRouteV << "\n";
            cout << "improved: " << newCost << "\n";            
        }        
        reinitSingleMoveInRou(routeU);
        reinitSingleMoveInRou(routeV);
        if (!isFixed) {
            //change location for routeU:
            changeLocCli(flexCost1.ft.ft);
            changeLocCli(flexCost1.ft.sc);
            //change location for routeV:
            changeLocCli(flexCost2.ft.ft);
            changeLocCli(flexCost2.ft.sc);
        }
        //routeU:
        Node* endNode = uSuc;
        Node* startNode = vSuc;
        while(startNode->idxClient){
            Node* tempNode = startNode->suc;
            startNode->suc = endNode;
            endNode->pred = startNode;
            endNode = startNode;
            startNode = tempNode;
        }
        endNode->pred = routeU->depot;
        routeU->depot->suc = endNode;
        if (reverseRouteU)routeU->reverse();        
        //routeV:
        startNode = nodeV;
        endNode = nodeU;
        while (endNode->idxClient)
        {
            Node* tempNode = endNode->pred;
            startNode->suc = endNode;
            endNode->pred = startNode;
            startNode = endNode;
            endNode = tempNode;
        }
        startNode->suc = routeV->depot->pred;
        routeV->depot->pred->pred = startNode;
        if (reverseRouteV)routeV->reverse();

        //update route data:
        routeU->updateRoute();
        routeV->updateRoute();        
        //reset infor for futher search
        uPred = nodeU->pred;
        uSuc = nodeU->suc;
        routeU = nodeU->rou;
        cost += newCost - oldCost;
        //count[1][1]++;
        if (pr->isTurnCkSol) {
            try {
                ckSol();
            }
            catch (const char* msg) {
                cerr << msg << endl;
                for (int i = 1; i <= n; ++i)pr->fileOut << giantT[i] << ", ";
                pr->fileOut << "\n";
                for (auto val : ordNodeLs)pr->fileOut << val << ", ";
                pr->fileOut.close();
                system("pause");
                exit(0);                
            }
        }
        return 1;    
    }

    //intra route 2Opt:
    int intraRoute2Opt() {
        //nodeV is always on the right side of nodeU
        int oldCost = routeV->depot->seqi_n->cost;
        SeqData* seq = nodeU->seq0_i;
        myseqs1.clear();
        myseqs1.push_back(uPred->seq0_i);        
        addRevSeqInPieces(nodeU, nodeV, myseqs1);
        myseqs1.push_back(vSuc->seqi_n);
        int newCost = oo; 
        int newCostFixed = oo;
        III flexRes;
        if(isFixed)newCost = seq->evaluation(myseqs1);
        else {
            /*flexRes = evalFlex(myseqs1);
            newCost = flexRes.sc;*/           
            traceLoc[0][1].clear();
            if (evalFlexRou(myseqs1, traceLoc[0][1])) {                                
                newCost = F[traceLoc[0][1].size()];
            }
        }                
        if (oldCost <= newCost)return 0;                     
        if (pr->isDebug) {
            cout << "ck intra route 2Opt\n"; 
            cout << boolalpha << isFixed << "\n";
            routeU->showR();
            cout << nodeU->idxClient << " " << nodeV->idxClient << "\n";
        }
        reinitSingleMoveInRou(routeV);
        if (!isFixed) {
            /*changeLocCli(flexRes.ft.ft);
            changeLocCli(flexRes.ft.sc);*/
            for (auto val : traceLoc[0][1])changeLocCli(val);
        }
        //update route structure
        Node* startNode = uPred;
        Node* endNode = nodeV;
        while (startNode != nodeU)
        {
            Node* tempNode = endNode->pred;
            startNode->suc = endNode;
            endNode->pred = startNode;
            startNode = endNode;
            endNode = tempNode;
        }
        nodeU->suc = vSuc;
        vSuc->pred = nodeU;
        //update route data:        
        routeV->updateRoute();
        cost += newCost - oldCost;
        uSuc = nodeU->suc;
        uPred = nodeU->pred;        
        //count[1][1]++;
        //check sol
        if (pr->isTurnCkSol) {
            try {
                ckSol();                
            }
            catch (const char* msg) {
                cerr << msg << endl;                
                for (int i = 1; i <= n; ++i)pr->fileOut << giantT[i] << ", ";
                pr->fileOut << "\n";
                for (auto val : ordNodeLs)pr->fileOut << val << ", ";
                pr->fileOut.close();
                system("pause");
                exit(0);                
            }
        }
        return 1;
    }
    ///ELSALGO:
    void exchange() {
        int posU = pr->Rng.getNumInRan(1, n);
        int posV = pr->Rng.getNumInRan(1, n);
        while (posU == posV)
        {
            posV = pr->Rng.getNumInRan(1, n);
        }
        swap(giantT[posU], giantT[posV]);
    }

    void interchange() {
        int posU = pr->Rng.getNumInRan(1, n);
        int posV = pr->Rng.getNumInRan(1, n);
        while (posU == posV)
        {
            posV = pr->Rng.getNumInRan(1, n);
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
            exchange();
        }
    }        
    
    /*void LocalSearch() {
        for (int i = 1; i <= n; ++i) {
            shuffle(ordNodeLs.begin(), ordNodeLs.end(), pr->Rng.generator);
            shuffle(nodes[i]->movesClu.begin(), nodes[i]->movesClu.end(), pr->Rng.generator);
            shuffle(nodes[i]->movesLoc.begin(), nodes[i]->movesLoc.end(), pr->Rng.generator);
        }
        int isImproved = 1;
        while (isImproved)
        {
            isImproved = 0;
            for (int i = 1; i <= 6; ++i) {                
                switch (i) {
                case 1:
                    isImproved = updateObjInter(1, 1);                    
                    break;
                case 2:
                    isImproved = updateObjInter(1, 0);
                    break;
                case 3:
                    isImproved = updateObjInter(2, 1);
                    break;
                case 4:
                    isImproved = updateObjInter(1, 2);
                    break;
                case 5:
                    isImproved = updateObjInter(2, 2);
                    break;
                case 6:
                    isImproved = updateObjInter(2, 0);
                    break;                    
                default:
                    break;
                }
            }
        }
    }*/

    void R_ILS() {
        int bestCost = oo;
        int budget = pr->nbR * pr->nbIls;
        int countB = 0;
        int* initGiantT = new int[n];
        initSol();
        for (int i = 1; i <= n; ++i)initGiantT[i - 1] = giantT[i];
        int stIdx = -1;
        int* resGiantT = new int[n + 1];
        int* localGiantT = new int[n + 1];
        int* curGiantT = new int[n + 1];        
        while (countB < budget)
        {                        
            cout << countB << "\n";
            stIdx = (stIdx + 1) % n;
            int curIdx = stIdx;
            for (int i = 1; i <= n; ++i) {
                giantT[i] = initGiantT[curIdx];                
                curIdx = (curIdx + 1) % n;
            }
            Split();
            //updateTotal();
            for (int i = 1; i <= n; ++i) {                
                curGiantT[i] = giantT[i];                
            }
            int oldCost = cost;
            int countILS = 0, countF = 0;
            while (countILS < pr->nbIls && countF < pr->nbF && countB < budget) {
                countILS++;
                for (int i = 1; i <= n; ++i)giantT[i] = curGiantT[i];
                exchange();
                Split();
                //updateTotal();
                countB++;
                if (cost < oldCost) {
                    oldCost = cost;
                    for (int i = 1; i <= n; ++i) {
                        curGiantT[i] = giantT[i];
                        localGiantT[i] = giantT[i];
                    }
                    countF = 0;
                }
                else
                {
                    countF++;
                }
            }
            if (oldCost < bestCost) {
                for (int i = 1; i <= n; ++i) {
                    resGiantT[i] = localGiantT[i];
                }
                bestCost = oldCost;
            }
        }
        cout << "best cost: " << bestCost << "\n";
    }

    void updateTotal() {
        updateObjInter();
        cvGiantT();
        Split();
        /*isFixed = true;
        for (int i = 0; i <= 1; ++i) {
            updateObjInter();
            cvGiantT();
            Split();
            isFixed = !isFixed;
        }*/
    }

    void ELS() {
        //genGiantT();
        //Split();
        initSol();//        
        updateTotal();
        cout << "improved: " << cost << "\n";
        int* resGiantT = new int[n + 1];
        int* curGiantT = new int[n + 1];
        for (int i = 1; i <= n; ++i) {
            curGiantT[i] = giantT[i];
            resGiantT[i] = giantT[i];
            //cout << giantT[i] << " " << endl;
        }
        int bestObj = cost;
        cout << cost << endl;
        int curObj;
        int curP = pr->pMin;
        for (int i = 1; i <= pr->nI; ++i) {
            cout << "ILS" << i << " " << bestObj << "\n";
            curObj = bestObj;
            for (int j = 1; j <= pr->nC; ++j) {
                for (int i1 = 1; i1 <= n; ++i1)giantT[i1] = curGiantT[i1];
                mutate(curP);
                Split();
                if (cost == oo)continue;
                try {
                    updateTotal();
                }
                catch (...) {
                    cout << "bug here\n";
                    for (auto val : giantT)cout << val << ", ";
                    cout << endl;
                    system("pause");
                    exit(0);
                }
                if (cost < curObj) {
                    curObj = cost;
                    for (int i1 = 1; i1 <= n; ++i1)resGiantT[i1] = giantT[i1];
                }
                /*if ((clock() - pr->start) / CLOCKS_PER_SEC >= pr->TL) {
                    if (curObj + EP < bestObj) {
                        bestObj = curObj;
                        goto thispro;
                    }
                }*/
                //cout << i << " " << j << "\n";
            }
            if (curObj < bestObj) {
                bestObj = curObj;
                for (int i1 = 1; i1 <= n; ++i1)curGiantT[i1] = resGiantT[i1];
                curP = pr->pMin;
                i = 1;
            }
            else {
                curP = min(pr->pMax, curP + 1);
            }
            //cout << i << endl;            
        }
        //thispro:
        cout << "best found obj: " << bestObj << endl;
        //pr->fileOut<< "best found obj: " << Util::round2num(bestObj) << endl;
        //if (pr->debugLS) {            
        for (int i = 1; i <= n; ++i)giantT[i] = resGiantT[i];
        Split();
        if (cost != oo && cost - bestObj > 0) {
            cout << "bug here\n";
            for (int i = 1; i <= n; ++i)pr->fileOut << giantT[i] << ", ";
            pr->fileOut << "\n";
            for (auto val : ordNodeLs)pr->fileOut << val << ", ";
            pr->fileOut.close();
            system("pause");
            exit(0);
        }
        cout << "final cost \n" << cost << endl;
        /*pr->fileOut << "Giant Tour:\n";
        for (int i = 1; i <= n; ++i)pr->fileOut << giantT[i] << ", ";
        pr->fileOut << "\n";
        pr->fileOut<<"cost split: " << Util::round2num(cost) << endl;*/
        //}
        delete[] resGiantT;
        delete[] curGiantT;
    }

    void printSol(ostream& os) {
        os << "Cost: " << cost << "\n\n";
        os << "Giant Tour\n";
        for (int i = 1; i <= n; ++i)os << giantT[i] << ", ";
        os << "\n\nClients in Route:\n";
        for (int i = 1; i <= m; ++i) {
            setR[i]->showRInFile(os);
        }
        os << "\nLocations in Route:\n";
        for (int i = 1; i <= m; ++i) {
            setR[i]->showRLocInFile(os);
        }
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
        ordNodeLs.clear();
        //for (int i = 0; i < n + 2 * pr->numVeh + 1; ++i)delete nodes[i];
        //for (int i = 1; i <= m; ++i)delete setR[i];
        //delete pr;
    }
};
