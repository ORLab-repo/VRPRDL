#include "lib.h"
#include "ReadData.h"
//#include "Solution.h"
#include "Ga.h"

using namespace std;

int arr[] = {
    0, 25, 59, 91, 1, 38, 52, 5, 92, 111, 23, 10, 39, 101, 41, 11, 64, 53, 80, 100, 89, 51, 103, 61, 113, 56, 88, 8, 76, 85, 36, 117, 94, 104, 66, 77, 14, 29, 13, 9, 24, 49, 47, 87, 30, 4, 114, 95, 60, 118, 70, 81, 45, 50, 19, 67, 108, 20, 96, 69, 63, 99, 58, 109, 7, 12, 54, 32, 105, 78, 79, 21, 28, 62, 98, 112, 110, 115, 44, 74, 35, 107, 31, 16, 22, 97, 37, 75, 57, 26, 65, 33, 82, 6, 71, 83, 93, 34, 120, 86, 42, 3, 18, 40, 119, 106, 2, 15, 43, 90, 17, 48, 102, 73, 55, 116, 68, 27, 46, 72, 84
};

int arrLS[] = {
    120, 68, 89, 33, 102, 23, 12, 116, 62, 4, 85, 108, 45, 56, 64, 21, 87, 5, 22, 65, 95, 52, 106, 59, 67, 25, 36, 41, 88, 77, 18, 3, 42, 80, 39, 44, 31, 53, 78, 13, 15, 81, 75, 29, 43, 74, 115, 61, 10, 111, 98, 112, 28, 17, 110, 9, 93, 20, 79, 24, 71, 8, 113, 50, 32, 92, 103, 69, 105, 34, 70, 46, 83, 49, 72, 2, 76, 86, 47, 101, 57, 26, 118, 11, 16, 48, 27, 73, 60, 82, 91, 117, 38, 30, 54, 66, 109, 84, 100, 119, 90, 37, 107, 35, 99, 51, 7, 19, 97, 96, 6, 58, 114, 55, 14, 104, 1, 63, 94, 40
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
        for (int i = 0; i < 120; ++i) {
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
    pathIn = "instances_VRPHRDL\\instance_" + to_string(idxIns) + "-triangle.vrp";
    //pathIn = "instances\\instance_" + to_string(idxIns) + "-triangle.vrp";
    //pathIn = "instances\\" + to_string(idxIns) + "-v2.vrp";
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
    int minCost = oo;
    for (int numRun = 0; numRun < 5; ++numRun) {
        pr->Rng.config(seed[numRun]);
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        bestSol.isFixed = false;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)bestSol.count[i][j] = 0LL;
        pr->isDebug = false;
        pr->debugLS = false;
        pr->isTurnCkSol = true;
        //int minCost = oo;
        int oldCost;
        int sumCost = 0;
        //bestSol.R_ILS();         
        /*for (int i = 1; i <= 1; ++i) {
            for (int j = 1; j <= bestSol.n; ++j) {
                bestSol.giantT[j] = arr[j];
                bestSol.ordNodeLs[j - 1] = arrLS[j - 1];
            }
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
        try {
            Algo.findGasSol();
        }
        catch (const char* msg) {
            cerr << msg << endl;                       
            exit(0);
            system("pause");
        }
        catch (...) {
            cout << "error\n";
            exit(0);
            system("pause");
        }

        minCost = min(minCost, Algo.bestCost);
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
        pr->fileOut << "elapsed time: " << elapsed_seconds.count() << "s\n\n";
    }
    pr->fileOut << "best run: " << minCost;
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
    return 0;
}