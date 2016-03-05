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
    // Modulus Set/Product Input Files
    ifstream MODULUS, PRODUCT;
    MODULUS.open("/Users/williamhill/Documents/Coursework/Fall Semester 2013/CSC 562/Modulus Generation/Modulus Generation/build/Debug/Modulus Sets.txt");
    PRODUCT.open("/Users/williamhill/Documents/Coursework/Fall Semester 2013/CSC 562/Modulus Generation/Modulus Product Set/build/Debug/Modulus Product Set.txt");
    
    // Output File
    FILE * excludedfile;
    excludedfile = fopen("Excluded Modulus Product Set.txt", "w");
    
    
    uint64_t modulus;
    const char * prod;
    string modulus_collection, modulus_product;

    while(getline(MODULUS, modulus_collection) && getline(PRODUCT, modulus_product)){
        // Modulus Set Construction
        vector<uint64_t> modulusSET;
        istringstream mstream (modulus_collection);
        while(mstream >> modulus){
            modulusSET.push_back(modulus);
        }
        
        // Modulus Product Declaration
        mpz_t M;
        mpz_init(M);
        prod = modulus_product.c_str();
        mpz_set_str(M, prod, 10);
        
        // Excluded Modulus Product Construction
        for(int i = 0; i < modulusSET.size(); i++){
            mpz_t excluded;
            mpz_init(excluded);
            mpz_div_ui(excluded, M, modulusSET[i]);
            mpz_out_str(excludedfile, 10, excluded);
            mpz_clear(excluded);
            
            fwrite("\n", sizeof(char), sizeof(char), excludedfile);
        }
        
    }
    MODULUS.close();
    PRODUCT.close();
    fclose(excludedfile);
    return 0;
}
