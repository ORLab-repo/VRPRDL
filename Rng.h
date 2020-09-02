#pragma once
#include <cstdlib>
#include<random>

namespace Rng{
    unsigned seed;
    std::mt19937 generator;// using for gen num in range [0,1]    

    void config(unsigned int sd){
        seed = sd;
        generator.seed(seed);        
    }

    //Gen real number in range [0,1]
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double genRealInRang01(){
        return distribution(generator);
    }

    //using for mutation_reversal
    std::uniform_real_distribution<double> distribution1(0.0,1.0);
    double genRealInRang01_muta(){
        return distribution1(generator);
    }

    int getIntRand(){
        return generator();
    }

    //gen number in range [u,v]
    int getNumInRan(int u,int v){
        return (generator()%(v-u+1))+u;
    }
}
