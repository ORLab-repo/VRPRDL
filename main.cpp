#include "lib.h"
#include "ReadData.h"
#include "Solution.h"

using namespace std;

int arr[] = {
    0, 59, 9, 51, 39, 66, 35, 108, 96, 52, 49, 74, 79, 27, 95, 53, 57, 43, 29, 94, 46, 101, 42, 8, 73, 22, 16, 44, 87, 31, 88, 37, 25, 13, 3, 115, 70, 5, 107, 81, 92, 77, 83, 86, 45, 34, 54, 91, 17, 14, 15, 119, 71, 55, 100, 41, 80, 93, 36, 103, 110, 82, 72, 63, 67, 98, 24, 50, 89, 114, 68, 1, 7, 105, 60, 65, 85, 12, 97, 40, 69, 109, 84, 21, 75, 48, 58, 62, 118, 19, 120, 2, 78, 56, 32, 38, 6, 111, 61, 20, 4, 76, 26, 18, 23, 102, 117, 33, 116, 113, 28, 10, 11, 104, 64, 90, 112, 47, 99, 106, 30
};

int seed[] = {
    18319894,
    23390422,
    36197069,
    45945346,
    54500951,
    63196367,
    71110057,
    89578146,
    96527670,
    10415237,
};

//exe -ni -nc -pmin -pmax -ld - bi -TL -method
string typeIns[] = { "C", "R", "RC" };
int _numI = 200, _numC = 200, _pMin = 1, _pMax = 2, _ld = 2, timeLimit = oo;
string method = "ELS";
bool _bi = true;
int main(int argc, char* argv[]) {    
    /*for (int i = 1; i < argc; ++i) {        
        if (string(argv[i]) == "-method") {            
            method = argv[i + 1];
        }
        if (string(argv[i]) == "-ni") {
            _numI = atoi(argv[i + 1]);
        }
        if (string(argv[i]) == "-nc") {
            _numC = atoi(argv[i + 1]);
        }
        if (string(argv[i]) == "-pmin") {
            _pMin = atoi(argv[i + 1]);
        }
        if (string(argv[i]) == "-pmax") {
            _pMax = atoi(argv[i + 1]);
        }
        if (string(argv[i]) == "-ld") {
            _ld = atoi(argv[i+1]);
        }
        if (string(argv[i]) == "-bi") {
            _bi = (atoi(argv[i + 1])) ? true : false;
        }
        if (string(argv[i]) == "-TL") {
            timeLimit = atoi(argv[i + 1]);
        }
    }*/
    cin.tie(0); cout.tie(0);
    ios::sync_with_stdio(0);
    string pathIn, pathOut;    
    pathIn = "instances\\instance_30-triangle.vrp";
    //pathOut = "solution_"+ to_string(idSed)+"\\" + typeIns[idType] + "\\" + "sol_" + to_string(idx) + ".txt";
    //cin >> pathIn;
    Param* pr = read_Ins(pathIn);    
    //set up param for algo:
    pr->nI = _numI;
    pr->nC = _numC;
    pr->pMin = _pMin;
    pr->pMax = _pMax;
    pr->lambda = _ld;
    pr->bi = _bi;                
    //pr->maxE *= 10;             
    //init(pr);
    //pr->fileOut.open(pathOut);
    cout.precision(6);                   
    //ckData(pr);
    init(pr);
    Solution bestSol(pr);     
    for (int i = 1; i <= bestSol.n; ++i)bestSol.giantT[i] = arr[i];   
    bestSol.Split();
    /*bestSol.solT.clear();
    for (int i = 0; i < sizeof(arr) / sizeof(int); ++i)bestSol.solT.push_back(arr[i]);    
    int ckCosts = 0;    
    for (int i = 1; i < bestSol.solT.size(); ++i) {
        cout<<bestSol.solT[i - 1]<<" "<<bestSol.solT[i]<<" "<< pr->costs[bestSol.solT[i - 1]][bestSol.solT[i]]<<"\n";
        ckCosts += pr->costs[bestSol.solT[i - 1]][bestSol.solT[i]];
    }
    cout << ckCosts << endl;*/
    //for (int i = 0; i <= bestSol.n; i++)bestSol.giantT[i] = arr[i];
    //bestSol.Split();
    //bestSol.ELS();
    ////cout << Util::round2num(bestSol.cost) << endl;    
    //cout << (clock()-start) / CLOCKS_PER_SEC << endl;
    //pr->fileOut << (clock() - start) / CLOCKS_PER_SEC << endl;
    //pr->fileOut.close();           
    system("pause");
    exit(0);    
    return 0;
}
