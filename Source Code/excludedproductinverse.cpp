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
    // Modulus Set/Excluded Modulus Product Input Files
    ifstream MODULUS, EXCLUDED;
    MODULUS.open("/Users/williamhill/Documents/Coursework/Fall Semester 2013/CSC 562/Modulus Generation/Modulus Generation/build/Debug/Modulus Sets.txt");
    EXCLUDED.open("/Users/williamhill/Documents/Coursework/Fall Semester 2013/CSC 562/Modulus Generation/Excluded Modulus Product Set/build/Debug/Excluded Modulus Product Set.txt");
    
    // Output File
    FILE * inversefile;
    inversefile = fopen("Excluded Modulus Product Inverse Set.txt", "w");
    
    uint64_t modulus;
    string modulus_collection, exclude_collection;
    while(getline(MODULUS, modulus_collection)){
        vector<string>   excludeSET;
        vector<uint64_t> modulusSET;
        istringstream mstream (modulus_collection);
        while(mstream >> modulus && getline(EXCLUDED, exclude_collection)){
            modulusSET.push_back(modulus);
            excludeSET.push_back(exclude_collection);
        }
        
        mpz_t * modset = new mpz_t[modulusSET.size()];
        mpz_t * exprod = new mpz_t[excludeSET.size()];
        const char * exclu;
        for(int i=0; i < modulusSET.size(); i++){
            mpz_init(exprod[i]);
            mpz_init(modset[i]);
            exclu = excludeSET[i].c_str();
            mpz_set_str(exprod[i], exclu, 10);
            mpz_set_ui(modset[i], modulusSET[i]);
        }
        
        // This portion of code is written to be an application of the Extended Euclid Algorithm to GMP integers.
        for(int i=0; i < modulusSET.size(); i++){
            mpz_t x, y, u, v, a, b, tempx, tempy, quot, remn;
            
            mpz_init(x);
            mpz_init(y);
            mpz_init(u);
            mpz_init(v);
            mpz_init(a);
            mpz_init(b);
            mpz_init(quot);
            mpz_init(remn);
            mpz_init(tempx);
            mpz_init(tempy);
            
            mpz_set_ui(x, 1L);
            mpz_set_ui(v, 1L);
            mpz_set_ui(y, 0L);
            mpz_set_ui(u, 0L);
            mpz_set(a, modset[i]);
            mpz_set(b, exprod[i]);
            
            while(mpz_cmp_ui(b, 0L) > 0){
                mpz_div(quot, a, b);
                mpz_mod(remn, a, b);
                
                mpz_set(a, b);
                mpz_set(b, remn);
                mpz_set(tempx, x);
                mpz_set(tempy, y);
                mpz_set(x, u);
                mpz_set(y, v);
                
                mpz_mul(u, quot, u);
                mpz_sub(u, tempx, u);
                mpz_mul(v, quot, v);
                mpz_sub(v, tempy, v);
            }
            
            while(mpz_cmp_ui(y, 0L) < 0){
                mpz_add(y, y, modset[i]);
            }

            mpz_out_str(inversefile, 10, y);
            fwrite(" ", sizeof(char), sizeof(char), inversefile);
            
            mpz_clear(x);
            mpz_clear(y);
            mpz_clear(u);
            mpz_clear(v);
            mpz_clear(a);
            mpz_clear(b);
            mpz_clear(quot);
            mpz_clear(remn);
            mpz_clear(tempx);
            mpz_clear(tempy);
        }
        fwrite("\n", sizeof(char), sizeof(char), inversefile);
    }
    
    
    
    MODULUS.close();
    EXCLUDED.close();
    fclose(inversefile);
    return 0;
}

