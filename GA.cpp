#include "Ga.h"

GA::GA()
{
}

GA::~GA()
{
}

void GA::addRou(Solution* u)
{   
    for (int idR = 1; idR <= u->m; ++idR) {
    //for (auto rou : u->setR) {
        Route* rou = u->setR[idR];
        valLength = rou->getCliInRou(valRou, valRouLoc);        
        if (valLength == 0)continue;        
        II hashVal = Util::getHash(valRou, valLength);    
        int existID = idWithLen[valLength][hashVal];
        if (existID == 0 || costR[existID] > rou->depot->seqi_n->cost){            
            if (existID == 0) {
                existID = ++curNumRou;
                idWithLen[valLength][hashVal] = curNumRou;
            }            
            routePoolPrv[existID][1] = -1;
            routePoolNxt[existID][valLength] = -1;
            for (int i = 1; i <= valLength; ++i) {               
                routePool[existID][i] = valRou[i];
                routePoolLoc[existID][i] = valRouLoc[i];
                if (i != 1)routePoolPrv[existID][i] = i - 1;
                if (i != valLength)routePoolNxt[existID][i] = i + 1;
            }            
            costR[existID] = rou->depot->seqi_n->cost;
            if(rou->caculateDis() != rou->depot->seqi_n->cost){
                throw "error cost route";
            }
            lenR[existID] = valLength;
        }        
    }
}

int GA::solveSCP()
{    
    IloEnv env;
    int res = -1;
    try {
        IloModel model(env);
        IloRangeArray cons(env);
        IloNumVarArray x(env, curNumRou + 1);
        //IloExpr sumX(env);
        for (int i = 1; i <= curNumRou; ++i) {
            x[i] = IloNumVar(env, 0, 1, ILOINT);
            //sumX += x[i];
        }
        IloExprArray sum(env, n + 1);
        for (int i = 1; i <= n; ++i) {
            sum[i] = IloExpr(env);
        }
        //construct sum
        int idCli;
        for (int i = 1; i <= curNumRou; ++i) {
            for (int j = 1; j <= lenR[i]; ++j) {
                idCli = routePool[i][j];
                sum[idCli] += x[i];
            }
        }
        //add constraint:
        for (int i = 1; i <= n; ++i)cons.add(sum[i] == 1);// set to "=1" for using SPP
        //add objective
        IloExpr sumObj(env);
        for (int i = 1; i <= curNumRou; ++i)sumObj += x[i] * costR[i];
        model.add(IloMinimize(env, sumObj));
        model.add(cons);
        IloCplex cplex(model);
        if (!cplex.solve()) {
            //res = oo;
            throw "can't sol";
        }
        //cplex.setOut(pr->fileOut);
        //cplex.setParam(IloCplex::TiLim, 10800);
        cplex.out() << "Solution value = " << Util::round2num(cplex.getObjValue()) << endl;                
        //cost = Util::round2num(cplex.getObjValue());
        cplex.out() << "Solution status = " << cplex.getStatus() << endl;        
        //construct new solution:                      
        for (int i = 1; i <= curNumRou; ++i)if (cplex.getValue(x[i]) >= 0.9) {
            //cout << i << "\n";
            /// get route containing duplicate customer.    
            for (int j = 1; j <= lenR[i]; ++j) {
                idRouBelong[routePool[i][j]].push_back(II(i,j));
                //cout << routePool[i][j] << " ";
            }
            //cout << "\n";
        }   
        /// remove duplicate customer
        int minIdRou = -1;
        int minImproved = oo;
        int idValRou, idValPos, valReduced;
        int idPrv, idNxt;
        int nodeU, prvU, nxtU;        
        for (int i = 1; i <= n; ++i)if (idRouBelong[i].size() > 1) {
            //cout << "client: " << i << "\n";
            minIdRou = -1;
            minImproved = oo;
            //find min route
            for (auto val : idRouBelong[i]) {
                idValRou = val.ft;
                idValPos = val.sc;
                nodeU = routePoolLoc[idValRou][idValPos];
                idPrv = routePoolPrv[idValRou][idValPos];
                idNxt = routePoolNxt[idValRou][idValPos];
                prvU = (idPrv == -1) ? 0 : routePoolLoc[idValRou][idPrv];
                nxtU = (idNxt == -1) ? 0 : routePoolLoc[idValRou][idNxt];
                /*cout << idValPos << " " << idValRou << " ";
                cout << prvU << " " << nodeU << " " << nxtU << "\n";*/
                valReduced = pr->costs[prvU][nxtU] - pr->costs[prvU][nodeU] - pr->costs[nodeU][nxtU];
                if (minImproved > valReduced) {
                    minImproved = valReduced;
                    minIdRou = idValRou;
                }
            }            
            //keep this route and remove the remains            
            for (auto val : idRouBelong[i]) {
                idValRou = val.ft;
                idValPos = val.sc;
                if (idValRou == minIdRou)continue;
                /*valPop->setR[1]->clearNode();
                valPop->setR[1]->updateRoute();
                cout << valPop->setR[1]->depot->pred->seq0_i->load << "\n";
                for (int i1 = 1; i1 <= lenR[idValRou]; ++i1) if(routePool[idValRou][i1] != -1){
                    valPop->nodes[routePool[idValRou][i1]]->idxClient = routePool[idValRou][i1];
                    valPop->nodes[routePool[idValRou][i1]]->idxLoc = routePoolLoc[idValRou][i1];
                    valPop->setR[1]->insertToRou(valPop->nodes[routePool[idValRou][i1]]);                                        
                }                
                valPop->setR[1]->showR();
                valPop->setR[1]->showRLoc();
                valPop->setR[1]->updateRoute();                
                valPop->setR[1]->ckRoute();                */
                idPrv = routePoolPrv[idValRou][idValPos];
                idNxt = routePoolNxt[idValRou][idValPos];
                routePoolNxt[idValRou][idPrv] = idNxt;
                routePoolPrv[idValRou][idNxt] = idPrv;
                routePool[idValRou][idValPos] = -1;//remove it
            }            
        }
        int valCount = 0;
        int ckCost = 0;
        for (int i = 1; i <= curNumRou; ++i)if (cplex.getValue(x[i]) >= 0.9) {
            ckCost += costR[i];
            for (int j = 1; j <= lenR[i]; ++j) if (routePool[i][j] != -1) {
                valPop->giantT[++valCount] = routePool[i][j];
            }
        }        
        if (valCount != n) {
            throw "Error in solve SCP";
        }
        else {
            for (int i = 1; i <= n; ++i)pr->fl << valPop->giantT[i] << ", ";
            pr->fl.close();
        }
        valPop->Split();        
        res = valPop->cost;
        if (res > (int)ceil(cplex.getObjValue())) {
            throw "error cost in SCP";
        }
    }
    catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
        throw "Error 7";
    }
    catch (const char* msg) {
        cerr << msg << endl;
        exit(0);
        system("pause");
    }
    catch (...)
    {
        cerr << "Unknown exception caught" << endl;
        throw "Error 7";
    }
    env.end();    
    for (int i = 1; i <= n; ++i) {
        idWithLen[i].clear();
        idRouBelong[i].clear();
    }   
    curNumRou = 0;
    return res;
}

bool GA::checkIdSol(Solution* u)
{
    for (int i = 1; i <= n; ++i)ddID[i] = 0;
    for (int i = 1; i <= n; ++i)ddID[u->giantT[i]] = 1;
    for (int i = 1; i <= n; ++i)if (ddID[i] == 0)return false;
    return true;
}

int GA::broken_pairs(Solution* u, Solution* v)
{
    //computing broken pair dis of u with v
    int* nxt = new int[n + 1];
    for (int i = 0; i < n - 1; ++i)nxt[u->giantT[i]] = u->giantT[i + 1];
    nxt[u->giantT[n - 1]] = -1;
    int numEqualEdge = 0;
    for (int i = 0; i < n - 1; ++i)if (v->giantT[i + 1] == nxt[v->giantT[i]])numEqualEdge++;
    delete[] nxt;
    return n - 1 - numEqualEdge;
}

bool GA::CheckEqual(Solution* u, Solution* v)
{    
    if (abs(u->cost - v->cost) <= omega) return true;
    /*for (int i = 1; i <= n; ++i)
        if(u->giantT[i] != v->giantT[i])return false;*/
    return false;
    /*return(abs(u->cost - v->cost) <= omega)
        || (broken_pairs(u, v) < threshold);*/
}

void GA::equalSol(Solution* u, Solution* v)
{
    u->cost = v->cost;
    for (int i = 1; i <= n; ++i) {
        u->giantT[i] = v->giantT[i];
    }
}

void GA::FindAdapt()
{
    sumAdapt = 0;
    for (int i = 1; i <= nPop1; ++i)
        sumAdapt += pop[i]->cost;
    for (int i = 1; i <= nPop1; ++i)
        adapt[i] = sumAdapt - pop[i]->cost;
    sumAdapt = 0;
    for (int i = 1; i <= nPop1; ++i)
        sumAdapt += adapt[i];
}

int GA::getChild()
{
    int u = pr->Rng.getNumInRan(1, nPop1);
    int v = pr->Rng.getNumInRan(1, nPop1);
    return (pop[u]->cost < pop[v]->cost) ? u : v;
    /*double percent = pr->Rng.genRealInRang01();
    int sum1 = 0;
    for (int i = 1; i <= nPop1; ++i) {
        sum1 += adapt[i];
        if (sumAdapt * percent <= (double)sum1)return i;
    }
    return 1;*/
}

void GA::choose(int& u, int& v)
{
    /// new selection
    u = getChild();
    do {
        v = getChild();
    } while (u == v);
    /*double avg=sumAdapt/2;
    percent=randomdN(0,1);
    double sum1=percent*avg;
    u=1;sumAdapt=0;*/   
}

void GA::uni(Solution* u, Solution* v, Solution* u1, Solution* v1)
{
    deque<int> q;
    q.clear();
    int* idu = new int[n + 1];
    int* idv = new int[n + 1];
    for (int i = 1; i <= n; ++i) {
        idu[i] = u->giantT[i];
        idv[i] = v->giantT[i];
    }
    /*for(int i=0;i<n;++i)out<<u.id[i]<<" ";
    out<<endl;
    for(int i=0;i<n;++i)out<<v.id[i]<<" ";
    out<<endl;*/
    //int rn=randomR(2*n/5,3*n/5);
    //int be=randomR(0,n-rn);
    /// test:
    //int en=be+rn-1;
    int be = pr->Rng.getNumInRan(1, n);
    //int en = pr->Rng.getNumInRan(be, n);
    int en = pr->Rng.getNumInRan(1, n);    
    while (be == en)
    {
        en = pr->Rng.getNumInRan(1, n);
    }
    int* dd = new int[n + 1];
    int vt;
    //born u:
    for (int i = 0; i <= n; ++i)
        dd[i] = 0;   
    int idI = be;
    /*for (int i = be; i <= en; ++i) {        
        u1->giantT[i] = idu[i];
        dd[idu[i]] = 1;
    }*/
    while (true)
    {        
        if ((en + 1) % n == idI % n)break;
        u1->giantT[idI] = idu[idI];
        dd[idu[idI]] = 1;
        idI++;        
        if (idI > n) idI = 1;
    }    
    /*for (int i = 1; i <= n; ++i)
        if (dd[idv[i]] == 0)
            q.push_back(idv[i]);
    vt = en + 1;
    while (!q.empty())
    {
        if (vt > n)
            vt = 1;
        u1->giantT[vt] = q.front();
        q.pop_front();
        vt++;
    }*/    
    int curId = en;
    for (int i = 1; i <= n; ++i) {                
        curId++;
        if (curId > n)curId = 1;
        if (dd[idv[curId]] == 0) {
            u1->giantT[idI] = idv[curId];
            idI = (idI + 1) % n;
            if (idI == 0)idI = n;
        }
    }
    //born v:
    /*for (int i = 0; i <= n; ++i)
        dd[i] = 0;
    for (int i = be; i <= en; ++i) {
        v1->giantT[i] = idv[i];
        dd[idv[i]] = 1;
    }
    for (int i = 1; i <= n; ++i)
        if (dd[idu[i]] == 0)
            q.push_back(idu[i]);
    vt = en + 1;
    while (!q.empty())
    {
        if (vt > n)
            vt = 1;
        v1->giantT[vt] = q.front();
        q.pop_front();
        vt++;
    }*/
    u1->Split();
    if (pr->Rng.genRealInRang01_muta() > pM) {//u1->interchange();
        u1->updateTotal();
    }
    //v1->Split(); v1->updateTotal();
    /*if(rand()%2==0){
        u1.updateObj();v1.updateObj();
    }else{
        u1.newSplit();v1.newSplit();
    }*/
    delete[] dd;
    delete[] idu;
    delete[] idv;
}

void GA::insertNew(Solution* u)
{
    int maxObj = oo;
    int posIns = nPop1 + 1;

    if (nPop1 == 0) {
        equalSol(pop[++nPop1], u);        
        return;
    }    
    if (u->cost < pop[1]->cost)posIns = 1;
    else {
        for (int i = 1; i <= nPop1; ++i)
        {
            if (CheckEqual(u, pop[i]))
            {                                
                return;
            }
        }

        for (int i = 1; i <= nPop1; ++i)
        {
            if (pop[i]->cost > u->cost) {
                posIns = i;
                break;
            }
        }
    }        
    //pop[nPop1 + 1] = u;
    equalSol(pop[nPop1 + 1], u);
    for (int i = nPop1; i >= posIns; --i) {
        swap(pop[i], pop[i + 1]);
    }
    nPop1++;
}

void GA::DelPopu()
{
    nPop1 = nPop;
}

void GA::InitPopu()
{        
    while (nPop1 != nPop) {                
        valPop->genGiantT();
        valPop->Split();
        //valPop->updateTotal();
        //cout<<valPop->cost<<endl;
        insertNew(valPop);
    }
}

void GA::DiversifyPopu(Solution* bestSol)
{
    int newPop = nPop / 3;    
    nPop1 = newPop;
    InitPopu();
    if (bestSol->cost > pop[1]->cost) {
        //pop[i].printSol();
        equalSol(bestSol, pop[1]);
    }
}

void GA::findGasSol(int maxNumGas)
{    
    int scpSol = oo;
    curNumRou = 0;
    for (int i = 1; i <= n; ++i) {
        idWithLen[i].clear();
        idRouBelong[i].clear();
    }
    int idFa, idMo;
    Solution* bestSol = new Solution(pr);
    clock_t be = clock();
    // generate population with 50 sols
    nPop1 = 0;
    InitPopu();
    /*for(int i=1;i<=25;++i){
        for(int j=0;j<n;++j)out<<pop[i].id[j]<<" ";
        out<<endl;
    }*/
    // hybrid 50 times:
    bestSol->cost = oo;
    equalSol(bestSol, pop[1]);
    //maxNumGas=500;
    numNotCha = 0;
    threshold = (n - 1) / 2;
    Solution* child1 = new Solution(pr);
    Solution* child2 = new Solution(pr);    
    for (int numga = 1;; ++numga)
    {        
        numNotCha++;
        //cout<<numga<<":"<<endl;
        cout << "name ins: " << pr->nameIns<<"\n";
        cout << numNotCha << " " << numga <<"{"<< endl;
        //out<<numga<<"{\n";
        
        //exit(0);
        // calc adaptation
        FindAdapt();
        // selection
        choose(idFa, idMo);
        /*cout<<"parent"<<endl;
        cout<<pop[idFa].obj<<" "<<pop[idMo].obj<<endl;*/
        /*cout<<idFa<<" "<<idMo<<endl;
        for(int i=0;i<n;++i)cout<<pop[idFa].id[i]<<" ";
        cout<<endl;
        for(int i=0;i<n;++i)cout<<pop[idMo].id[i]<<" ";
        cout<<endl;*/
        // hybrid
        uni(pop[idFa], pop[idMo], child1, child2);
        /*cout<<"children"<<endl;
        cout<<child1.obj<<" "<<child2.obj<<endl;*/
        /*for(int i=0;i<n;++i)cout<<child1.id[i]<<" ";
        cout<<endl;
        for(int i=0;i<n;++i)cout<<child2.id[i]<<" ";
        cout<<endl;*/
        //insert child:
        insertNew(child1);       
        addRou(child1);       
        if (nPop1 == nPop + delta)
            DelPopu();
        /*insertNew(child2);
        if (nPop1 == nPop + delta)
            DelPopu();*/
        if (bestSol->cost > pop[1]->cost) {
            numNotCha = 0;
            equalSol(bestSol, pop[1]);
        }
        // if bestSol don't change 100 times change 25 worst sols by 25 new sols
        //if(numNotCha>=200){numNotCha=0;addPopu();}
        if (curNumRou >= maxNumRou) {
            scpSol = min(scpSol, solveSCP());
            insertNew(valPop);
            if (nPop1 == nPop + delta)
                DelPopu();
            /*pop[1]->Split();
            addRou(pop[1]);*/
        }
        if (numNotCha == 3000)
        {
            threshold = min(threshold - 1, 1);
            int oldBestObj = bestSol->cost;
            DiversifyPopu(bestSol);
            if (bestSol->cost != oldBestObj)numNotCha = 0;
        }
        //cout<<"best obj:\n";cout<<bestSol.obj<<endl<<endl;
        if (numNotCha == 5000) {
            bestSol->Split();
            /*bestSol.printSol();
            cout << id_test << " " << bestSol.obj << " " << (double)(clock() - be) / CLOCKS_PER_SEC << endl;
            fl << bestSol.obj << " " << (double)(clock() - be) / CLOCKS_PER_SEC << "\n";*/
            bestCost = bestSol->cost;
            cout << bestSol->cost;
            pr->fileOut<< bestSol->cost<<"\n";            
            for (int i = 1; i <= n; ++i)pr->fileOut << bestSol->giantT[i] << ", ";
            pr->fileOut << "\n";
            pr->fileOut <<(double)(clock() - be) / CLOCKS_PER_SEC << "\n";
            //break;
        }
        //if((double)(clock()-be)/CLOCKS_PER_SEC>=600)break;
        cout << curNumRou << "\n";
        cout << "SCP sol: " << scpSol << "\n";
        cout<<bestSol->cost<<"}"<<endl;
    }
}
