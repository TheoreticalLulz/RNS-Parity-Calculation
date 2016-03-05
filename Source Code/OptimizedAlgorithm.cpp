#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>
#include <gmp.h>

using namespace std;

int main(){
    // Input Files
    ifstream MODULUS, RESIDUE, PRODUCT, EXCLUDE, INVERSE;

    MODULUS.open("/Users/williamhill/Documents/Coursework/Fall Semester 2013/CSC 562/Modulus Generation/Modulus Generation/build/Debug/Modulus Sets.txt");
    RESIDUE.open("/Users/williamhill/Documents/Coursework/Fall Semester 2013/CSC 562/Modulus Generation/Modulus Generation/build/Debug/Residue Sets.txt");
    PRODUCT.open("/Users/williamhill/Documents/Coursework/Fall Semester 2013/CSC 562/Modulus Generation/Modulus Product Set/build/Debug/Modulus Product Set.txt");
    EXCLUDE.open("/Users/williamhill/Documents/Coursework/Fall Semester 2013/CSC 562/Modulus Generation/Excluded Modulus Product Set/build/Debug/Excluded Modulus Product Set.txt");
    INVERSE.open("/Users/williamhill/Documents/Coursework/Fall Semester 2013/CSC 562/Modulus Generation/Excluded Modulus Product Inverse Set/build/Debug/Excluded Modulus Product Inverse Set.txt");
    
    // Output Files
    FILE * modifiedfile;
    modifiedfile = fopen("Modified Algorithm Output.txt", "w");
    
    
    // Main Program Execution Loop
    string modulus_collection, residue_collection, modulus_product, exclude_collection, inverse_collection;
    while(getline(MODULUS, modulus_collection) && getline(RESIDUE, residue_collection) && getline(PRODUCT, modulus_product) && getline(INVERSE, inverse_collection)){
        mpz_t modprod;
        mpz_init(modprod);
        const char * product = modulus_product.c_str();
        mpz_set_str(modprod, product, 10);
        
        uint64_t         modulus, residue, inverse;
        vector<string>   excludeSET;
        vector<uint64_t> modulusSET, residueSET, inverseSET;
        istringstream mstream(modulus_collection), rstream(residue_collection), istream(inverse_collection);
        while(mstream >> modulus && rstream >> residue && istream >> inverse && getline(EXCLUDE, exclude_collection)){
            modulusSET.push_back(modulus);
            residueSET.push_back(residue);
            excludeSET.push_back(exclude_collection);
            inverseSET.push_back(inverse);
        }
        
        uint64_t * modset = new uint64_t[modulusSET.size()];
        uint64_t * resset = new uint64_t[residueSET.size()];
        mpz_t    * excset = new mpz_t[excludeSET.size()];
        const char * temp;
        for(int i=0; i < modulusSET.size(); i++){
            mpz_init(excset[i]);
            temp = excludeSET[i].c_str();
            mpz_set_str(excset[i], temp, 10);
            
            modset[i] = modulusSET[i];
            resset[i] = residueSET[i];
        }
        
        // Construct 't' Value for Rounding
        double   bitspace = log(modulusSET.size());
        uint64_t maximums = modset[0];
        for(int i=1; i < modulusSET.size(); i++){
            if(modset[i] > maximums)
                maximums = modset[i];
        }
        mpz_t rounding;
        mpz_init(rounding);
        bitspace += log(maximums);
        uint64_t t = ceil(bitspace);
        mpz_ui_pow_ui(rounding, 2L, t);
        
        // Construct Inverse Set
        uint64_t invers[modulusSET.size()][modulusSET.size()];
        for(int i=0; i < modulusSET.size(); i++){
            invers[i][0] = inverseSET[i];
            for(int j = 1; j < modulusSET.size(); j++){
                mpz_t inv;
                mpz_init_set_ui(inv, invers[i][j-1]);
                mpz_mul_ui(inv, inv, modset[modulusSET.size()- j]);
                mpz_mod_ui(inv, inv, modset[i]);
                invers[i][j] = mpz_get_ui(inv);
            }
        }
        
        
        // Begin Modified Algorithm
        struct start, end;
        double  elapsed = 0;
        mpz_t   parity;
        mpf_t   elwrite;
        
        mpz_t * S = new mpz_t[modulusSET.size()];
        mpz_t * U = new mpz_t[modulusSET.size()];
        mpz_t * V = new mpz_t[modulusSET.size()];
        mpz_t ur, vr;
        mpz_init_set_ui(ur, 0L);
        mpz_init_set_ui(vr, 0L);
        mpz_init(parity);
        
        for(int i = (int) modulusSET.size()-1; i >= 0; i--){
            for(int j=0; j <= i; j++){
                mpz_init(S[j]);
                mpz_init(U[i]);
                mpz_init(V[i]);
                
                mpz_init(S[j]);
                mpz_init(U[i]);
                mpz_init(V[i]);
                
                mpz_set_ui(S[j], resset[j]);
                mpz_mul_ui(S[j], S[j], invers[j][i]);
                mpz_mod_ui(S[j], S[j], modset[j]);
                mpz_mul(U[i], S[j], rounding);
                mpz_tdiv_q_ui(U[i], U[i], modset[i]);
                if(mpz_cmp_ui(U[i], 0L) != 0)
                    mpz_add_ui(V[i], U[i], 1L);
                else
                    mpz_set(V[i], U[i]);
            }
            
            //<------------------------------- Start Time
            clock_gettime(CLOCK_REALTIME, &start);
            for(int j=0; j <= i; j++){
                mpz_add(ur, ur, U[j]);
                mpz_add(vr, vr, V[j]);
            }
            
            mpz_tdiv_q(ur, ur, rounding);
            mpz_tdiv_q(vr, vr, rounding);
            clock_gettime(CLOCK_REALTIME, &end);
            elapsed += (end.tv_nsec - start.tv_nsec) + ((double)(end.tv_sec - start.tv_sec) * 1e9);
            //<------------------------------- End Time
            if(mpz_cmp(ur, vr) == 0)
                break;
        }
            //<------------------------------- Final Start Time
        clock_gettime(CLOCK_REALTIME, &start);
        mpz_mod_ui(ur, ur, 2L);
        mpz_set(parity, ur);
        for(int i=0; i < modulusSET.size(); i++){
            mpz_mod_ui(S[i], S[i], 2L);
            mpz_add(parity, parity, S[i]);
        }
        mpz_mod_ui(parity, parity, 2L);
        clock_gettime(CLOCK_REALTIME, &end);
        elapsed += (end.tv_nsec - start.tv_nsec) + ((double)(end.tv_sec - start.tv_sec) * 1e9);
            //<------------------------------- Final End Time
        
        mpf_init(elwrite);
        mpf_set_d(elwrite, elapsed);
        mpf_out_str(modifiedfile, 10, 0, elwrite);
        fwrite("\n", sizeof(char), sizeof(char), modifiedfile);

        mpz_clear(ur);
        mpz_clear(vr);
        mpz_clear(parity);
        mpz_clear(rounding);
    }
    
    MODULUS.close();
    RESIDUE.close();
    PRODUCT.close();
    EXCLUDE.close();
    INVERSE.close();
    fclose(modifiedfile);
    
    return 0;
}

