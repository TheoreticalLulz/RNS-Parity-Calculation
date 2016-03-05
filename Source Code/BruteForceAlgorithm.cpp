#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>
#include <stdint.h>
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
    FILE * bruteforcefile, * residuefile;
    bruteforcefile = fopen("Brute Force Algorithm Output.txt", "w");
    residuefile    = fopen("Integer Decimal Representation.txt", "w");


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

        // Begin Brute Force Algorithm
        //clock_t start, end;
        struct timespec start, end;
        double  elapsed;
        mpf_t   elwrite;
        mpz_t   X, parity, integer;

        //start = clock(); //<--------------------- Start Time
        clock_gettime(CLOCK_REALTIME, &start);
        mpz_init(X);
        mpz_init(parity);
        for(int i=0; i < residueSET.size(); i++){
            mpz_init(integer);
            mpz_set_ui(integer, resset[i]);
            mpz_mul_ui(integer, integer, invset[i]);
            mpz_mod_ui(integer, integer, modset[i]);
            mpz_mul(integer, integer, excset[i]);

            mpz_add(X, X, integer);
        }
        mpz_mod(X, X, modprod);
        mpz_mod_ui(parity, X, 2L);
        //end = clock();   //<--------------------- End Time
        clock_gettime(CLOCK_REALTIME, &end);


        //elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
        elapsed = (end.tv_nsec - start.tv_nsec) + ((double)(end.tv_sec - start.tv_sec) * 1e9);
        mpf_init(elwrite);
        mpf_set_d(elwrite, elapsed);
        mpf_out_str(bruteforcefile, 10, 0, elwrite);
        mpz_out_str(residuefile, 10, X);
        fwrite("\n", sizeof(char), sizeof(char), bruteforcefile);
        fwrite("\n", sizeof(char), sizeof(char), residuefile);

    }

    MODULUS.close();
    RESIDUE.close();
    PRODUCT.close();
    EXCLUDE.close();
    INVERSE.close();
    fclose(bruteforcefile);
    fclose(residuefile);
    return 0;
}
