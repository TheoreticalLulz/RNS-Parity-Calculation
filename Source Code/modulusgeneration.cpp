#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <time.h>
using namespace std;

/*********************************************************************************************************************************************************************************************************/
class MODULUS{
public:
    MODULUS();
    void getPrimes(int x);
    void generateModuli(int y);
    
    vector<uint64_t> primes;
};
/*********************************************************************************************************************************************************************************************************/
int main(int argc, const char * argv[])
{
    srand(time(NULL));
    
    MODULUS set;
    set.getPrimes(1000);
    set.generateModuli(1000);
    
    return 0;
}
/*********************************************************************************************************************************************************************************************************/
MODULUS::MODULUS(){
    
}
/*********************************************************************************************************************************************************************************************************/
void MODULUS::getPrimes(int x){
    primes.push_back(2);
    for(int i = 3; primes.size() < x; i++){
        bool prime = true;
        for(int j = 0; j < primes.size(); j++){
            if(i % primes[j] == 0)
                prime = false;
        }
        if(prime == true)
            primes.push_back(i);
    }
    primes.erase(primes.begin());
}
/*********************************************************************************************************************************************************************************************************/
void MODULUS::generateModuli(int y){
    ofstream modulusfile;
    ofstream residuefile;
    modulusfile.open("Modulus Sets.txt");
    residuefile.open("Residue Sets.txt");
    for(int i = 0; i < y; i++){
        vector<uint64_t> temp = primes;
        vector<uint64_t> moduli;
        for(int mod = 0; mod < (rand()%100 + 1); mod++){
            moduli.push_back(1);
            for (int num = 0; num < (rand()%5 + 1); num++) {
                int loc = rand()%temp.size();
                moduli[mod] *= temp[loc];
                temp.erase(temp.begin() + loc);
            }
        }
        
        for(int i=0; i < moduli.size(); i++){
            modulusfile << moduli[i] << " ";
            residuefile << rand()%moduli[i] << " ";
        }
        
        //End Modulus Set
        modulusfile << "\n";
        residuefile << "\n";
    }
    modulusfile.close();
    residuefile.close();
}
/*********************************************************************************************************************************************************************************************************/
