#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <stdint.h>
#include <vector>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <sys/time.h>

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
    FILE * controlfile;
    controlfile = fopen("Mi Lu & Jen-Shiun Chiang Algorithm Output.txt", "w");


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
        uint64_t * invset = new uint64_t[inverseSET.size()];
        mpz_t    * excset = new mpz_t[excludeSET.size()];
        const char * temp;
        for(int i=0; i < modulusSET.size(); i++){
            mpz_init(excset[i]);
            temp = excludeSET[i].c_str();
            mpz_set_str(excset[i], temp, 10);

            modset[i] = modulusSET[i];
            resset[i] = residueSET[i];
            invset[i] = inverseSET[i];
        }

        // Construct Additional Algorithm Components
        double   bitspace = log(modulusSET.size());
        for(int i=0; i < modulusSET.size(); i++){
            bitspace += log(modset[i]);
        }
        mpz_t rounding;
        mpz_init(rounding);
        uint64_t t = ceil(bitspace);
        mpz_ui_pow_ui(rounding, 2L, t);

        // Begin Mi Lu & Jen-Shiun Chiang Algorithm
        struct timespec start, end;
        double  elapsed;
        mpz_t   parity;
        mpf_t   elwrite;
        


        mpz_t * U = new mpz_t[modulusSET.size()];
        mpz_t * S = new mpz_t[modulusSET.size()];
        mpz_t r;
        mpz_init(r);
        mpz_init(parity);
        for(int i=0; i < modulusSET.size(); i++){
            mpz_init(S[i]);
            mpz_init(U[i]);
            mpz_set_ui(S[i], resset[i]);
            mpz_mul_ui(S[i], S[i], invset[i]);
            mpz_mod_ui(S[i], S[i], modset[i]);
            mpz_mul(U[i], S[i], rounding);
            mpz_tdiv_q_ui(U[i], U[i], modset[i]);
        }
        clock_gettime(CLOCK_REALTIME, &start);  //<----- Start Time
        for(int i=0; i < modulusSET.size(); i++){
            mpz_add(r, r, U[i]);
        }
        mpz_tdiv_q(r, r, rounding);
        mpz_mod_ui(r, r, 2L);
        mpz_set(parity, r);
        for(int i=0; i < modulusSET.size(); i++){
            mpz_mod_ui(S[i], S[i], 2L);
            mpz_add(parity, parity, S[i]);
        }
        mpz_mod_ui(parity, parity, 2L);
        clock_gettime(CLOCK_REALTIME, &end); //<-------- End Time

        elapsed = (end.tv_nsec - start.tv_nsec) + ((double)(end.tv_sec - start.tv_sec) * 1e9);
        mpf_init(elwrite);
        mpf_set_d(elwrite, elapsed);
        mpf_out_str(controlfile, 10, 0, elwrite);
        fwrite("\n", sizeof(char), sizeof(char), controlfile);


    }

    MODULUS.close();
    RESIDUE.close();
    PRODUCT.close();
    EXCLUDE.close();
    INVERSE.close();
    fclose(controlfile);

    return 0;
}

