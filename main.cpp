#include "lib.h"
#include "ReadData.h"
#include "Solution.h"

using namespace std;

int arr[] = {
    0, 29, 73, 87, 35, 86, 84, 114, 53, 88, 100, 28, 66, 102, 106, 36, 94, 5, 12, 45, 104, 39, 77, 117, 67, 26, 81, 95, 1, 7, 97, 103, 16, 11, 69, 80, 57, 33, 89, 43, 79, 50, 24, 109, 85, 99, 40, 115, 25, 30, 2, 47, 14, 13, 51, 56, 112, 71, 91, 83, 70, 60, 37, 116, 96, 46, 68, 17, 108, 31, 18, 54, 101, 8, 21, 9, 32, 75, 119, 110, 27, 55, 52, 93, 49, 62, 41, 92, 78, 64, 58, 6, 98, 61, 59, 34, 65, 107, 120, 4, 74, 90, 10, 63, 76, 105, 19, 118, 72, 48, 3, 113, 23, 38, 22, 15, 111, 44, 20, 82, 42
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
int _numI = 200, _numC = 200, _pMin = 1, _pMax = 2, _ld = 2, timeLimit = oo;
string method = "ELS";
bool _bi = true;
int main(int argc, char* argv[]) {        
    srand(18319894);
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
    ios::sync_with_stdio(0);
    string pathIn, pathOut;    
    pathIn = "instances\\instance_30-triangle.vrp";
    getSamples();    
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
    /*
    ofstream fl("ckSol.txt");  
    std::chrono::time_point<std::chrono::system_clock> start, end;
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
    /*for (int i = 1; i <= bestSol.n; ++i)bestSol.giantT[i] = arr[i];   
    bestSol.Split();
    cout << bestSol.cost << "\n";*/
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
