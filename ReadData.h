#pragma once
#include "lib.h"
#include "Param.h"

using namespace std;

Param* read_Ins(string path){
    Param* pr = new Param();
    ifstream filein(path);
    int countLine = 0;
    vector<string> lineCont;// content of Line
    pr->listCL.clear();
    pr->listLoc.clear();
    while(filein){
        string line;
        getline(filein,line);
        lineCont = Util::splitString(line, "\t");        
        /*if (lineCont[0] == "NAME:") {
            countLine++;
            continue;
        }
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
            continue;
        }
        if (lineCont[0] == "NUM_CUSTOMERS:") {            
            pr->numClient = stoi(lineCont[1]);
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
       /* if (lineCont[0] == "DEPOT_SECTION") {
            countLine++;
            continue;
        }*/
        //get coordinate:
        if(countLine == 1){                        
            pr->listLoc.push_back(Location(stod(lineCont[1]),
                              stod(lineCont[2])                                                            
                              )
                       );
         }
        //get distance:
        if(countLine == 2){
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

        }
        //if(ck_depot){
        //    swap(pr->listCL[0], pr->listCL[Util::convertStringToNum(lineCont[0])-1]);
        //}
        if (lineCont[0] == "EOF")break;      
    }
    return pr;
}

void init(Param* pr) {
    //init
    /*pr->Costs.resize(pr->numClient);
    pr->Slopes.resize(pr->numClient);
    for (int i = 0; i < pr->numClient; ++i) {
        pr->Costs[i] = vector<double>(pr->numClient);
        pr->Slopes[i] = vector<double>(pr->numClient);
    }

    for (int i = 0; i < pr->numClient; ++i)
        for (int j = 0; j < pr->numClient; ++j) {
            pr->Costs[i][j] = pr->listCL[i].countDis(pr->listCL[j]);
            if (i != j)pr->Slopes[i][j] = pr->listCL[i].calSlope(pr->listCL[j]);
            else pr->Slopes[i][j] = 0;
        }*/
}