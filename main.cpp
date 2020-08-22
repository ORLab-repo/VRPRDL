#include "lib.h"
#include "ReadData.h"
//#include "Solution.h"

using namespace std;

int arr[] = {
    0, 8, 10, 2, 4, 6, 3, 5, 1, 9, 7,
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
    for (int i = 1; i < argc; ++i) {        
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
    }
    cin.tie(0); cout.tie(0);
    ios::sync_with_stdio(0);
    string pathIn, pathOut;    
    pathIn = "instances\\instance_0-triangle.vrp";
    //pathOut = "solution_"+ to_string(idSed)+"\\" + typeIns[idType] + "\\" + "sol_" + to_string(idx) + ".txt";
    //cin >> pathIn;
    Param* pr = read_Ins(pathIn);
    //pr->cap = 10000000;
    //pr->start = start;
    //set up param for algo:
    pr->nI = _numI;
    pr->nC = _numC;
    pr->pMin = _pMin;
    pr->pMax = _pMax;
    pr->lambda = _ld;
    pr->bi = _bi;                
    //pr->maxE *= 10;             
    init(pr);
    //pr->fileOut.open(pathOut);
    cout.precision(6);                   
    //Solution bestSol(pr);                       
    //bestSol.ELS();
    ////cout << Util::round2num(bestSol.cost) << endl;    
    //cout << (clock()-start) / CLOCKS_PER_SEC << endl;
    //pr->fileOut << (clock() - start) / CLOCKS_PER_SEC << endl;
    //pr->fileOut.close();           
    system("pause");
    exit(0);    
    return 0;
}
