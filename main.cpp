#include "lib.h"
#include "ReadData.h"
#include "Solution.h"

using namespace std;

int arr[] = {
    0, 82, 72, 63, 67, 24, 75, 48, 62, 118, 81, 102, 14, 15, 101, 42, 73, 46, 30, 106, 99, 112, 47, 87, 16, 19, 66, 59, 9, 51, 39, 111, 2, 20, 4, 44, 31, 88, 17, 95, 53, 57, 43, 29, 94, 22, 52, 74, 49, 79, 27, 40, 26, 23, 97, 41, 80, 93, 36, 103, 110, 76, 18, 1, 7, 60, 65, 85, 12, 92, 77, 83, 86, 45, 34, 54, 91, 10, 11, 104, 64, 90, 108, 96, 38, 55, 35, 78, 61, 120, 107, 68, 69, 109, 84, 21, 98, 50, 89, 114, 37, 25, 13, 3, 115, 70, 5, 28, 58, 56, 6, 32, 119, 71, 8, 117, 100, 105, 33, 116, 113
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
int _numI = 30, _numC = 100, _pMin = 1, _pMax = 2, _ld = 2, timeLimit = oo;
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
    Rng::config(seed[0]);
    //srand(seed[0]);    
    string pathIn, pathOut;
    pathIn = "instances\\instance_30-triangle.vrp";
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
    //pr->maxE *= 10;             
    //init(pr);
    //pr->fileOut.open(pathOut);
    cout.precision(6);
    //ckData(pr);
    init(pr);
    Solution bestSol(pr);
    ///check solution
    
    ofstream fl("ckSol.txt");
    fl << "abc";
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
    bestSol.genGiantT();
    bestSol.Split();
    cout << "initial: " << bestSol.cost << "\n";
    //bestSol.updateTotal();   
    bestSol.LocalSearch();
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
