#pragma once
#include "lib.h"
#include "Param.h"

using namespace std;

Param* read_Ins(string path) {
    Param* pr = new Param();
    ifstream filein(path);
    int countLine = 0;
    vector<string> lineCont;// content of Line
    pr->listCL.clear();
    pr->listLoc.clear();
    while (filein) {
        string line;
        getline(filein, line);
        lineCont = Util::splitString(line, "\t");
        if (lineCont[0] == "NAME:") {
            //countLine++;
            pr->nameIns = lineCont[1];
            continue;
        }
        /*
        if (lineCont[0] == "TYPE:") {
            countLine++;
            continue;
        }*/
        if (lineCont[0] == "DIMENSION:") {
            pr->numLoc = stoi(lineCont[1]);
            continue;
        }
        if (lineCont[0] == "CAPACITY:") {
            pr->Q = stoi(lineCont[1]);
            continue;
        }
        if (lineCont[0] == "TIME_HORIZON:") {
            pr->T = stoi(lineCont[1]);
            cout << pr->T << endl;
            continue;
        }
        if (lineCont[0] == "NUM_CUSTOMERS:") {
            pr->numClient = stoi(lineCont[1]);
            pr->listCL.resize(pr->numClient);
            continue;
        }
        //1 NODE
        if (lineCont[0] == "NODE_COORD_SECTION") {
            countLine++;
            continue;
        }
        //2 WEIGHT
        if (lineCont[0] == "EDGE_WEIGHT_SECTION") {
            countLine++;
            continue;
        }
        //3 TRAVEL TIME
        if (lineCont[0] == "EDGE_TRAVEL_TIME_SECTION") {
            countLine++;
            continue;
        }
        //4 TIME WINDOW
        if (lineCont[0] == "TIME_WINDOW_SECTION") {
            countLine++;
            continue;
        }
        //5 CLUSTER
        if (lineCont[0] == "CLUSTER_SECTION") {
            countLine++;
            continue;
        }
        //6 DEMAND
        if (lineCont[0] == "DEMAND_SECTION") {
            countLine++;
            continue;
        }
        //7 DEPOT
        if (lineCont[0] == "DEPOT_SECTION") {
            countLine++;
            continue;
        }

        //get coordinate:
        if (countLine == 1) {
            pr->listLoc.push_back(Location(stod(lineCont[1]),
                stod(lineCont[2])
            )
            );
        }

        //get distance:
        if (countLine == 2) {
            assert(size(lineCont) == pr->numLoc);
            pr->costs.push_back(vector<int>(pr->numLoc));
            for (int i = 0; i < pr->numLoc; ++i) {
                pr->costs.back()[i] = stoi(lineCont[i]);
            }
        }

        //get travel time        
        if (countLine == 3) {
            assert(size(lineCont) == pr->numLoc);
            pr->times.push_back(vector<int>(pr->numLoc));
            for (int i = 0; i < pr->numLoc; ++i) {
                pr->times.back()[i] = stoi(lineCont[i]);
            }
        }

        //get TW:
        if (countLine == 4) {
            pr->listLoc[stoi(lineCont[0]) - 1].stTime = stoi(lineCont[1]);
            pr->listLoc[stoi(lineCont[0]) - 1].enTime = stoi(lineCont[2]);
        }

        //get cluster:
        if (countLine == 5) {
            int idx = stoi(lineCont[0]) - 1;
            for (int i = 1; i < lineCont.size(); i++) {
                pr->listCL[idx].listLoc.pb(stoi(lineCont[i]) - 1);
                pr->listLoc[stoi(lineCont[i]) - 1].idxClient = idx;

            }
        }

        //get demand:
        if (countLine == 6) {
            pr->listCL[stoi(lineCont[0]) - 1].demand = stoi(lineCont[1]);
        }

        ///uncomment if index of depot is not one
        /*if(countLine == 7){
            swap(pr->listCL[0], pr->listCL[stoi(lineCont[0])-1]);
        }*/

        if (lineCont[0] == "EOF")break;
    }
    filein.close();
    return pr;
}

Param* read_Ins_GVRPTW(string path) {
    Param* pr = new Param();
    ifstream filein(path);
    int countLine = 0;
    vector<string> lineCont;// content of Line
    pr->listCL.clear();
    pr->listLoc.clear();
    while (filein) {
        string line;
        getline(filein, line);
        lineCont = Util::splitString(line, "\t");
        if (lineCont[0] == "NAME:") {
            //countLine++;
            pr->nameIns = lineCont[1];
            continue;
        }
        /*
        if (lineCont[0] == "TYPE:") {
            countLine++;
            continue;
        }*/
        if (lineCont[0] == "DIMENSION:") {
            pr->numLoc = stoi(lineCont[1]);
            continue;
        }
        if (lineCont[0] == "VEHICLES:") {
            pr->maxVeh = stoi(lineCont[1]);
            continue;
        }
        if (lineCont[0] == "CAPACITY:") {
            pr->Q = stoi(lineCont[1]);
            continue;
        }
        if (lineCont[0] == "TIME_HORIZON:") {
            pr->T = stoi(lineCont[1]);
            cout << pr->T << endl;
            continue;
        }
        if (lineCont[0] == "NUM_CUSTOMERS:") {
            pr->numClient = stoi(lineCont[1]);
            pr->listCL.resize(pr->numClient);
            continue;
        }
        //1 NODE
        if (lineCont[0] == "NODE_COORD_SECTION") {
            countLine++;
            continue;
        }
        //2 WEIGHT
        if (lineCont[0] == "EDGE_WEIGHT_SECTION") {
            countLine++;
            continue;
        }
        //3 TRAVEL TIME
        if (lineCont[0] == "EDGE_TRAVEL_TIME_SECTION") {
            countLine++;
            continue;
        }
        //4 TIME WINDOW
        if (lineCont[0] == "TIME_WINDOW_SECTION") {
            countLine++;
            continue;
        }
        //5 CLUSTER
        if (lineCont[0] == "CLUSTER_SECTION") {
            countLine++;
            continue;
        }
        //6 DEMAND
        if (lineCont[0] == "DEMAND_SECTION") {
            countLine++;
            continue;
        }
        //7 DEPOT
        if (lineCont[0] == "DEPOT_SECTION") {
            countLine++;
            continue;
        }

        //get coordinate:
        if (countLine == 1) {
            pr->listLoc.push_back(Location(stod(lineCont[1]),
                stod(lineCont[2])
            )
            );
        }

        //get distance:
        if (countLine == 2) {
            assert(size(lineCont) == pr->numLoc);
            pr->costs.push_back(vector<int>(pr->numLoc));
            for (int i = 0; i < pr->numLoc; ++i) {
                pr->costs.back()[i] = stoi(lineCont[i]);
            }
        }

        //get travel time        
        if (countLine == 3) {
            assert(size(lineCont) == pr->numLoc);
            pr->times.push_back(vector<int>(pr->numLoc));
            for (int i = 0; i < pr->numLoc; ++i) {
                pr->times.back()[i] = stoi(lineCont[i]);
            }
        }

        //get TW:
        if (countLine == 4) {
            pr->listLoc[stoi(lineCont[0]) - 1].stTime = stoi(lineCont[1]);
            pr->listLoc[stoi(lineCont[0]) - 1].enTime = stoi(lineCont[2]);
        }

        //get cluster:
        if (countLine == 5) {
            int idx = stoi(lineCont[0]) - 1;
            for (int i = 1; i < lineCont.size(); i++) {
                pr->listCL[idx].listLoc.pb(stoi(lineCont[i]) - 1);
                pr->listLoc[stoi(lineCont[i]) - 1].idxClient = idx;

            }
        }

        //get demand:
        if (countLine == 6) {
            pr->listCL[stoi(lineCont[0]) - 1].demand = stoi(lineCont[1]);
        }

        ///uncomment if index of depot is not one
        /*if(countLine == 7){
            swap(pr->listCL[0], pr->listCL[stoi(lineCont[0])-1]);
        }*/

        if (lineCont[0] == "EOF")break;
    }
    filein.close();
    return pr;
}

void ckData(Param* pr) {
    cout << "Number of clients: " << pr->numClient << "\n";
    cout << "Number of locations: " << pr->numLoc << "\n";
    cout << "Capacity: " << pr->Q << "\n";
    cout << "Time horizon" << pr->T << "\n";
    //locations:
    cout << "List of locations:\n";
    for (auto val : pr->listLoc) {
        cout << fixed << "Coordinate: " << val.x << " " << val.y << endl;
        cout << "TW: " << val.stTime << " " << val.enTime << endl;
    }
    cout << "List of clients:\n";
    for (auto val : pr->listCL) {
        cout << "demand: " << val.demand << "\n";
        cout << "locations: ";
        for (auto loc : val.listLoc)cout << loc << " ";
        cout << "\n";
        assert(val.demand != -1);
    }
    cout << "costs:\n";
    for (int i = 0; i < pr->numLoc; ++i) {
        for (int j = 0; j < pr->numLoc; ++j)cout << pr->costs[i][j] << " ";
        cout << "\n";
    }
    cout << "times: \n";
    for (int i = 0; i < pr->numLoc; ++i) {
        for (int j = 0; j < pr->numLoc; ++j)cout << pr->times[i][j] << " ";
        cout << "\n";
    }
}
void init(Param* pr) {
    //init
    for (int i = 0; i < pr->listCL.size(); ++i) {
        for (auto idxLoc : pr->listCL[i].listLoc) {
            pr->listLoc[idxLoc].demand = pr->listCL[i].demand;
        }
    }
    set<II> sDis;
    pr->corDis.resize(pr->numLoc);
    for (int i = 0; i < pr->numLoc; ++i) {
        pr->corDis[i].resize(pr->numLoc);
        for (int j = 0; j < pr->numLoc; ++j) {
            pr->corDis[i].push_back(pr->costs[i][j]
                + pr->ldTw * max(pr->listLoc[i].stTime + pr->times[i][j] - pr->listLoc[j].enTime, 0));
        }        
    }
    for (int i = 1; i < pr->numLoc; ++i) {
        // calculating correlation measure for each cluster
        int valMin = oo;
        sDis.clear();
        for (int j = 1; j < pr->numClient; ++j)if (j != pr->listLoc[i].idxClient) {
            valMin = oo;
            for (auto val : pr->listCL[j].listLoc) {               
                valMin = min(valMin, pr->corDis[val][i]);
            }
            sDis.insert(II(valMin, j));
        }
        // getting maxNeibor nearest cluster
        pr->listLoc[i].moves.clear();
        for (auto val : sDis) {
            if (val.first >= oo)continue;
            pr->listLoc[i].moves.pb(val.second);
            if (pr->listLoc[i].moves.size() == pr->maxNeibor)break;
        }
    }

    ///preprocessing:  (only used for VRPHRDL or VRPRDL)  
    ///comment when travel time does not satisfy the triangle inequality
    /*for (int i = 0; i < pr->numLoc; ++i) {
        pr->new_costs.push_back(vector<int>(pr->numLoc));
    }
    for (int i = 0; i < pr->numLoc; ++i)
        for (int j = 0; j < pr->numLoc; ++j) {
            pr->new_costs[i][j] = pr->listLoc[i].calDis(pr->listLoc[j]);
        }
    for (int k = 0; k < pr->numLoc; ++k)
        for (int i = 0; i < pr->numLoc; ++i)
            for (int j = 0; j < pr->numLoc; ++j)pr->new_costs[i][j] = min(pr->new_costs[i][j], pr->new_costs[i][k] + pr->new_costs[k][j]);
    */

    //<- Compute cost
    //eliminate nodes:          
    int num_reduced = 0;
    for (int i = 1; i < pr->numLoc; ++i) {
        int val = pr->listLoc[0].stTime + pr->times[0][i];//e0+t0i
        if (val > pr->listLoc[i].enTime //li
            || max(val, pr->listLoc[i].stTime) + pr->times[i][0] > pr->T
            ) 
        {            
            for (int j = 0; j < pr->numLoc; ++j) {                
                pr->costs[i][j] = oo;
                pr->costs[j][i] = oo;
            }
            num_reduced++;
        }//eliminate arcs
        else {
            for (int j = 1; j < pr->numLoc; ++j)if (j != i) {
                int eStime = max(val, pr->listLoc[i].stTime);
                if (pr->listLoc[i].stTime + pr->times[i][j] > pr->listLoc[j].enTime //ei+tij>lj
                    || eStime + pr->times[i][j] > pr->listLoc[j].enTime //0-i-j-0
                    || max(eStime + pr->times[i][j], pr->listLoc[j].stTime) + pr->times[j][0] > pr->T
                    ) {                    
                    pr->costs[i][j] = oo; 
                    num_reduced++;
                }
            }
        }
    }
    cout << num_reduced << endl;
    //tighten TW:
    num_reduced = 0;
    while (true)
    {
        bool flag = true;// can't reduce
        for (int i = 1; i < pr->numLoc; ++i)if (pr->costs[i][0] != oo) {//not a eliminated node
            int valMinIJ = oo,valMinJI = oo;
            int valMaxIJ = -oo, valMaxJI = -oo;
            for (int j = 0; j < pr->numLoc; ++j)if (pr->costs[i][j] != oo && j!=i) {//not a eliminated arc                
                valMinIJ = min(valMinIJ, pr->listLoc[j].stTime - pr->times[i][j]);
                valMinJI = min(valMinJI, pr->listLoc[j].stTime + pr->times[j][i]);
                valMaxJI = max(valMaxIJ, pr->listLoc[j].enTime + pr->times[j][i]);
                valMaxIJ = max(valMaxIJ, pr->listLoc[j].enTime - pr->times[i][j]);
            }
            valMinIJ = min(valMinIJ, valMinJI);
            if (pr->listLoc[i].stTime < min(valMinIJ, pr->listLoc[i].enTime)) {
                pr->listLoc[i].stTime = min(valMinIJ, pr->listLoc[i].enTime);
                flag = false;
                num_reduced++;
            }
            valMaxIJ = max(valMaxIJ, valMaxJI);
            if (pr->listLoc[i].enTime > max(valMaxIJ,pr->listLoc[i].stTime)) {
                pr->listLoc[i].enTime = max(valMaxIJ, pr->listLoc[i].stTime);
                flag = false;
                num_reduced++;
            }
        }
        /*cout << num_reduced << endl;
        cout << flag << endl;*/
        if (flag)break;
    }    
}

void init_GVRPTW(Param* pr) {
    //init
    for (int i = 0; i < pr->listCL.size(); ++i) {
        for (auto idxLoc : pr->listCL[i].listLoc) {
            pr->listLoc[idxLoc].demand = pr->listCL[i].demand;
        }
    }
    set<II> sDis;
    pr->corDis.resize(pr->numLoc);
    for (int i = 0; i < pr->numLoc; ++i) {
        pr->corDis[i].resize(pr->numLoc);
        for (int j = 0; j < pr->numLoc; ++j) {
            pr->corDis[i].push_back(pr->costs[i][j]
                + pr->ldTw * max(pr->listLoc[i].stTime + pr->times[i][j] - pr->listLoc[j].enTime, 0));
        }
    }
    for (int i = 1; i < pr->numLoc; ++i) {
        // calculating correlation measure for each cluster
        int valMin = oo;
        sDis.clear();
        for (int j = 1; j < pr->numClient; ++j)if (j != pr->listLoc[i].idxClient) {
            valMin = oo;
            for (auto val : pr->listCL[j].listLoc) {
                valMin = min(valMin, pr->corDis[val][i]);
            }
            sDis.insert(II(valMin, j));
        }
        // getting maxNeibor nearest cluster
        pr->listLoc[i].moves.clear();
        for (auto val : sDis) {
            if (val.first >= oo)continue;
            pr->listLoc[i].moves.pb(val.second);
            if (pr->listLoc[i].moves.size() == pr->maxNeibor)break;
        }
    }

    ///preprocessing:  (only used for VRPHRDL or VRPRDL)  
    ///comment when travel time does not satisfy the triangle inequality
    /*for (int i = 0; i < pr->numLoc; ++i) {
        pr->new_costs.push_back(vector<int>(pr->numLoc));
    }
    for (int i = 0; i < pr->numLoc; ++i)
        for (int j = 0; j < pr->numLoc; ++j) {
            pr->new_costs[i][j] = pr->listLoc[i].calDis(pr->listLoc[j]);
        }
    for (int k = 0; k < pr->numLoc; ++k)
        for (int i = 0; i < pr->numLoc; ++i)
            for (int j = 0; j < pr->numLoc; ++j)pr->new_costs[i][j] = min(pr->new_costs[i][j], pr->new_costs[i][k] + pr->new_costs[k][j]);
    */

    //<- Compute cost
    //eliminate nodes:          
    int num_reduced = 0;
    for (int i = 0; i < pr->numLoc; ++i) {
        for (int j = 0; j < pr->numLoc; ++j)if (j != i) {            
            if (pr->listLoc[i].stTime + pr->times[i][j] > pr->listLoc[j].enTime //ei+tij>lj 
                ) {
                pr->costs[i][j] = oo;
                num_reduced++;
            }
        }
    }
    cout << num_reduced << endl;
    //tighten TW:
    num_reduced = 0;
    while (true)
    {
        bool flag = true;// can't reduce
        for (int i = 1; i < pr->numLoc; ++i) {
            int valMinIJ = oo, valMinJI = oo;
            int valMaxIJ = -oo, valMaxJI = -oo;
            for (int j = 0; j < pr->numLoc; ++j)if (pr->costs[i][j] != oo && j != i) {//not a eliminated arc                
                valMinIJ = min(valMinIJ, pr->listLoc[j].stTime - pr->times[i][j]);
                valMinJI = min(valMinJI, pr->listLoc[j].stTime + pr->times[j][i]);
                valMaxJI = max(valMaxIJ, pr->listLoc[j].enTime + pr->times[j][i]);
                valMaxIJ = max(valMaxIJ, pr->listLoc[j].enTime - pr->times[i][j]);
            }
            valMinIJ = min(valMinIJ, valMinJI);
            if (pr->listLoc[i].stTime < min(valMinIJ, pr->listLoc[i].enTime)) {
                pr->listLoc[i].stTime = min(valMinIJ, pr->listLoc[i].enTime);
                flag = false;
                num_reduced++;
            }
            valMaxIJ = max(valMaxIJ, valMaxJI);
            if (pr->listLoc[i].enTime > max(valMaxIJ, pr->listLoc[i].stTime)) {
                pr->listLoc[i].enTime = max(valMaxIJ, pr->listLoc[i].stTime);
                flag = false;
                num_reduced++;
            }
        }
        /*cout << num_reduced << endl;
        cout << flag << endl;*/
        if (flag)break;
    }
}