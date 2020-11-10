#include "lib.h"
#include "ReadData.h"
//#include "Solution.h"
#include "Ga.h"

using namespace std;

int arr[] = {
    0, 83, 2, 34, 86, 11, 45, 120, 6, 26, 105, 79, 72, 57, 1, 110, 107, 51, 60, 78, 63, 97, 23, 91, 20, 24, 30, 116, 22, 16, 89, 85, 111, 38, 114, 7, 74, 106, 76, 54, 68, 59, 94, 35, 67, 9, 73, 77, 101, 18, 25, 29, 56, 13, 62, 82, 5, 61, 69, 103, 3, 96, 14, 39, 99, 64, 21, 98, 102, 58, 117, 90, 84, 15, 95, 28, 53, 92, 71, 112, 40, 75, 4, 80, 66, 41, 100, 37, 87, 10, 12, 46, 27, 47, 93, 49, 32, 17, 108, 44, 55, 88, 119, 19, 115, 36, 104, 8, 52, 109, 81, 65, 113, 31, 70, 43, 118, 42, 33, 48, 50
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
    cin.tie(0); cout.tie(0);    
    int numRun = 1;
    int idxIns;
    cin >> idxIns;
    //cout << "Num of runs: ";
    //cin >> numRun;    
    int isFlex = 1;
    //cout << "isFlex: ";
    //cin >> isFlex;
    ios::sync_with_stdio(0);    
    //srand(seed[0]);    
    string pathIn, pathOut;
    pathIn = "instances\\instance_" + to_string(idxIns) + "-triangle.vrp";
    //getSamples();    
    pathOut = "solution\\sol_" + to_string(idxIns) + ".txt";
    //cin >> pathIn;
    Param* pr = read_Ins(pathIn);
    //set up param for algo:
    pr->nI = _numI;
    pr->nC = _numC;
    pr->pMin = _pMin;
    pr->pMax = _pMax;
    pr->lambda = _ld;
    pr->bi = _bi;
    //pr->Rng.config(seed[0]);   
    //pr->maxE *= 10;             
    //init(pr);
    pr->fileOut.open(pathOut);
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
    GA Algo;
    Algo.init(pr);
    for (int numRun = 0; numRun < 10; ++numRun) {
        pr->Rng.config(seed[numRun]);
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
        int sumCost = 0;
        //bestSol.R_ILS();  
        /*for (int i = 1; i <= 100; ++i) {
            for (int j = 1; j <= bestSol.n; ++j)bestSol.giantT[j] = arr[j];
            bestSol.Split();
            cout << bestSol.cost << "\n";
            bestSol.updateTotal();
            cout << bestSol.cost << "\n\n";
            sumCost += bestSol.cost;
            minCost = min(minCost, bestSol.cost);
        }
        cout << "min: " << minCost << "\n";
        cout << "avg: " << (double)sumCost / 100<<"\n";*/
        //for (int i = 1; i <= bestSol.n; ++i)fl << bestSol.giantT[i] << ", ";
        //fl.close();
        //bestSol.ELS();        
        Algo.findGasSol();

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);

        std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
        pr->fileOut << "elapsed time: " << elapsed_seconds.count() << "s\n\n";
    }
    pr->fileOut.close();
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
