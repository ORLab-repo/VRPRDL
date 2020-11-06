#include "lib.h"
#include "ReadData.h"
//#include "Solution.h"
#include "Ga.h"

using namespace std;

int arr[] = {
    0, 54, 113, 61, 37, 44, 35, 16, 73, 43, 17, 38, 52, 110, 69, 55, 116, 50, 23, 6, 83, 68, 98, 112, 41, 75, 26, 14, 65, 115, 97, 120, 86, 11, 29, 42, 3, 18, 47, 108, 27, 15, 109, 46, 2, 106, 119, 28, 78, 84, 88, 79, 21, 51, 89, 103, 72, 48, 59, 91, 1, 5, 92, 111, 99, 63, 9, 104, 76, 94, 66, 77, 24, 40, 95, 10, 39, 90, 102, 33, 105, 82, 71, 101, 19, 67, 8, 20, 96, 34, 57, 93, 107, 31, 100, 53, 64, 58, 32, 56, 12, 7, 25, 87, 80, 4, 114, 30, 49, 22, 74, 13, 117, 60, 118, 70, 81, 36, 45, 85, 62
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

vector<vector<int>> samArr;
void getSamples() {
    ifstream fileIn("samples.txt");
    vector<string> lineCont;// content of Line
    while (fileIn)
    {
        string line;
        getline(fileIn, line);        
        lineCont = Util::splitString(line, ",");
        if (lineCont.size() != 120)break;
        samArr.push_back(vector<int>(120));
        for (int i = 0; i < 120; ++i){
            samArr.back()[i] = stoi(lineCont[i]);
        }
    }
    fileIn.close();
}
//exe -ni -nc -pmin -pmax -ld - bi -TL -method
string typeIns[] = { "C", "R", "RC" };
int _numI = 100, _numC = 20, _pMin = 1, _pMax = 2, _ld = 2, timeLimit = oo;
string method = "ELS";
bool _bi = true;
void ckChanged(vector<int> arr) {
    arr.push_back(0);
}
int main(int argc, char* argv[]) {
    //freopen("log.txt", "w", stdout);
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
    //cin.tie(0); cout.tie(0);    
    int numRun = 1;
    //cout << "Num of runs: ";
    //cin >> numRun;    
    int isFlex = 1;
    //cout << "isFlex: ";
    //cin >> isFlex;
    ios::sync_with_stdio(0);    
    //srand(seed[0]);    
    string pathIn, pathOut;
    pathIn = "instances\\instance_33-triangle.vrp";
    //getSamples();    
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
    pr->Rng.config(seed[0]);
    //pr->maxE *= 10;             
    //init(pr);
    //pr->fileOut.open(pathOut);
    cout.precision(6);
    //ckData(pr);    
    init(pr);    
    Solution bestSol(pr);
    ///check solution    
    ofstream fl("ckSol.txt");    
    /*std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    for (int i = 1; i <= 100000; ++i) {
        for (int j = 1; j <= bestSol.n; ++j)bestSol.giantT[j] = samArr[i - 1][j - 1];
        bestSol.Split();
        fl << bestSol.cost << ", ";
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
    fl.close();
    */
    //for (int i = 1; i <= bestSol.n; ++i)bestSol.giantT[i] = i;                 
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    bestSol.isFixed = false;        
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)bestSol.count[i][j] = 0LL;
    pr->isDebug = false;
    pr->debugLS = false;
    pr->isTurnCkSol = true;
    int minCost = oo;
    int oldCost;   
    //bestSol.R_ILS();        
    //for (int j = 1; j <= bestSol.n; ++j)bestSol.giantT[j] = arr[j];           
    //for (int i = 1; i <= bestSol.n; ++i)fl << bestSol.giantT[i] << ", ";
    //fl.close();
    //bestSol.ELS();
    GA Algo;
    Algo.init(pr);
    Algo.findGasSol();
    cout << "strategy: " << boolalpha << bestSol.isFixed << "\n";
    cout << "times: " << pr->total << "\n";
    cout << "cost: " << bestSol.cost << "\n";
    for (int i = 1; i < 4; ++i)
        for (int j = 0; j < 4; ++j)cout << i << " " << j << " " << bestSol.count[i][j] << "\n";
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
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
