#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>
#include <gmp.h>

using namespace std;

int main(){
    ifstream MODULUS;
    MODULUS.open("/Users/williamhill/Documents/Coursework/Fall Semester 2013/CSC 562/Modulus Generation/Modulus Generation/build/Debug/Modulus Sets.txt");
    
    FILE * productfile;
    productfile = fopen("Modulus Product Set WORKING.txt", "w");
    
    uint64_t modulus;
    string modulus_collection;
    while(getline(MODULUS, modulus_collection)){
        vector<uint64_t> modulusSET;
        istringstream mstream (modulus_collection);
        while(mstream >> modulus){
            modulusSET.push_back(modulus);
        }
        
        mpf_t M;
        mpf_init(M);
        mpf_set_ui(M, 1L);
        for(int i=0; i < modulusSET.size(); i++){
            mpf_mul_ui(M, M, modulusSET[i]);
        }
        
        mpf_out_str(productfile, 10, 0, M);
        fwrite("\n", sizeof(char), sizeof(char), productfile);
    }
    MODULUS.close();
    fclose(productfile);
    return 0;
}

