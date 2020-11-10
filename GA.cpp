#include "Ga.h"

GA::GA()
{
}

GA::~GA()
{
}

bool GA::checkIdSol(Solution* u)
{
    for (int i = 1; i <= n; ++i)ddID[i] = 0;
    for (int i = 1; i <= n; ++i)ddID[u->giantT[i]] = 1;
    for (int i = 1; i <= n; ++i)if (ddID[i] == 0)return false;
    return true;
}

bool GA::CheckEqual(Solution* u, Solution* v)
{    
    for (int i = 1; i <= n; ++i)
        if(u->giantT[i] != v->giantT[i])return false;
    return true;    
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
    int en = pr->Rng.getNumInRan(be, n);
    int* dd = new int[n + 1];
    int vt;
    //born u:
    for (int i = 0; i <= n; ++i)
        dd[i] = 0;
    for (int i = be; i <= en; ++i) {
        u1->giantT[i] = idu[i];
        dd[idu[i]] = 1;
    }
    for (int i = 1; i <= n; ++i)
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
    }
    //born v:
    for (int i = 0; i <= n; ++i)
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
    }
    u1->Split(); u1->updateTotal();
    v1->Split(); v1->updateTotal();
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
    //threshold = (n - 1) / 2;
    Solution* child1 = new Solution(pr);
    Solution* child2 = new Solution(pr);    
    for (int numga = 1;; ++numga)
    {        
        numNotCha++;
        //cout<<numga<<":"<<endl;
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
        if (nPop1 == nPop + delta)
            DelPopu();
        insertNew(child2);
        if (nPop1 == nPop + delta)
            DelPopu();
        if (bestSol->cost > pop[1]->cost) {
            numNotCha = 0;
            equalSol(bestSol, pop[1]);
        }
        // if bestSol don't change 100 times change 25 worst sols by 25 new sols
        //if(numNotCha>=200){numNotCha=0;addPopu();}
        if (numNotCha == 4000)
        {
            //threshold = min(threshold - 1, 1);
            int oldBestObj = bestSol->cost;
            DiversifyPopu(bestSol);
            if (bestSol->cost != oldBestObj)numNotCha = 0;
        }
        //cout<<"best obj:\n";cout<<bestSol.obj<<endl<<endl;
        if (numNotCha == 15000) {
            bestSol->Split();
            /*bestSol.printSol();
            cout << id_test << " " << bestSol.obj << " " << (double)(clock() - be) / CLOCKS_PER_SEC << endl;
            fl << bestSol.obj << " " << (double)(clock() - be) / CLOCKS_PER_SEC << "\n";*/
            cout << bestSol->cost;
            pr->fileOut<< bestSol->cost<<"\n";
            break;
        }
        //if((double)(clock()-be)/CLOCKS_PER_SEC>=600)break;
        cout<<bestSol->cost<<"}"<<endl;
    }
}
