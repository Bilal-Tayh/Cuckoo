#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <sstream>
#include <string>
#include <cstring>
#include <unordered_map>
// #include <windows.h>
// #include "Cuckoo.h"
//#include "Cuckoo_water_level.h"
#include <stdint.h>
#include <xmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <sys/timeb.h>
#include "Utils.hpp"
#include "HeavyHitters_Cuckoo.h"
#include "MurmurHash3.h"

#define CLK_PER_SEC CLOCKS_PER_SEC
#define CAIDA16_SIZE 152197439
#define CAIDA18_SIZE 175880896
#define UNIV1_SIZE 17323447

typedef unsigned key;
typedef unsigned short val;

using namespace std;


// void getKeysFromFile(string filename, vector<key*> &keys, int size) {
//   ifstream stream;
//   stream.open(filename, fstream::in | fstream::out | fstream::app);
//   if (!stream) {
//     throw invalid_argument("Could not open " + filename + " for reading.");
//   }
// 
//   key* data = (key*) malloc(sizeof(key) * size);
//   string line;
//   for (int i = 0; i< size; ++i){
//     getline(stream, line);
//     try {
//       data[i] = stoull(line);
//     } catch (const invalid_argument& ia) {
//       cerr << "Invalid argument: " << ia.what() << " at line " << i << endl;
//       cerr << line << endl;
//       --i;
//     }
//   }
// 
//   keys.push_back(data);
// 
//   stream.close();
// }




void getKeysAndWeightsFromFile(string filename, vector<key*> &keys, vector<val*> &value, int size) {
  ifstream stream;
  stream.open(filename, fstream::in | fstream::out | fstream::app);
  if (!stream) {
    throw invalid_argument("Could not open " + filename + " for reading.");
  }

  key* file_keys = (key*) malloc(sizeof(key) * size);
  val* file_ws = (val*) malloc(sizeof(val) * size);

  string line;
  string len;
  string id;
  for (int i = 0; i < size; ++i){
    getline(stream, line);
    std::istringstream iss(line);
    iss >> len;
    iss >> id;
    try {
      file_keys[i] = stoul(id);
      file_ws[i] = stoi(len);
    } catch (const std::invalid_argument& ia) {
      cerr << "Invalid argument: " << ia.what() << " at line " << i << endl;
      cerr << len << " " << id << endl;;
      --i;
      exit(1);
    }
  }
  keys.push_back(file_keys);
  value.push_back(file_ws);

  stream.close();
}




double run(int traceSize, key* kk, val* vv, int runs, double maxl, double gamma){
    double time = 0;
    for (int run = 0; run < runs; run++) {
            
            
            struct timeb begintb, endtb;
            clock_t begint, endt;
            Cuckoo_waterLevel_HH_no_FP_SIMD_256<uint32_t, uint32_t> hh(1, 2, 3, maxl,gamma,traceSize);
            begint = clock();
            for (int i = 0; i < traceSize; ++i) {
                hh.insert(kk[i], vv[i]);
//                     printf("%p\n",index);
            }
            endt = clock();
            cout << ((double)(endt-begint))/CLK_PER_SEC << endl;
            time += ((double)(endt-begint))/CLK_PER_SEC;
            
        
    }
    return time/runs;
}




float runNrmse(int traceSize, key* kk, val* vv, int runs, double maxl, double gamma){
    float sum=0.0;
    for (int run = 0; run < runs; run++) {
        uint64_t c =0; 
        uint64_t vol =0;
        unordered_map<key, val> map;
        Cuckoo_waterLevel_HH_no_FP_SIMD_256<uint32_t, uint32_t> hh(1, 2, 3, maxl,gamma,traceSize);
        for (int i = 0; i < traceSize; ++i) {
            map[kk[i]] +=vv[i];
            hh.insert(kk[i], vv[i]);
                uint64_t err = map[kk[i]] - (hh.query(kk[i]));
                c+= err*err;
                vol+=vv[i];
        }
        
        float mse = c/traceSize;
        float rmse = sqrt(mse);
        float nrmse = rmse/vol;
        sum+=nrmse;
        
    }

    return sum/runs;
}



double run1(double maxl, double gamma){
    fprintf(stderr, "run1\n");
    double time = 0;
    int traceSize = 30;
            
    struct timeb begintb, endtb;
    clock_t begint, endt;
    Cuckoo_waterLevel_HH_no_FP_SIMD_256<uint32_t, uint32_t> hh(1, 2, 3, maxl,gamma,traceSize);
    begint = clock();
    
    
    
    for (int i = 0; i < traceSize; ++i) {
        val v = std::rand();
        printf("%d  %d\n",i+1,v);
        hh.insert(i+1,  v);
    }
    
    
    
    for (int i = 0; i < traceSize; ++i) {
        val v = std::rand();
        key k = (std::rand()%19) + 1;
        printf("%d  %d\n",k,v);
        hh.insert(k,  v);
    }
    
    
    
    
    endt = clock();
    cout << ((double)(endt-begint))/CLK_PER_SEC << endl;
    time += ((double)(endt-begint))/CLK_PER_SEC;

    return time;
}

int main(int argc, char* argv[])
{
//     vector<key*> keys;
//     vector<val*> values;
//     getKeysAndWeightsFromFile("../../datasets/CAIDA16/mergedAggregatedPktlen_Srcip", keys, values, CAIDA16_SIZE);
// //     sizes.push_back(CAIDA16_SIZE);
// //     datasets.push_back("caida");
// 
// 
//     
//     
//     vector<key*>::iterator k_it1 = keys.begin();
//     key* k1 = *k_it1;
// 
//     int size2 =sizeof(key);
//     int traceSize = CAIDA16_SIZE;
//     ofstream f1;
// 
//     f1.open("ca16_k.dat",ios::binary|ios::out); 
// 
//     for(int j=0;j<traceSize;j++){
//         f1.write((char *)&k1[j],size2);
//     }
// 
//     f1.close();
//     
//     
//     
//     
//     vector<val*>::iterator v_it1 = values.begin();
//     val* v1 = *v_it1;
// 
// 
//     ofstream ff1;
//     int size1 =sizeof(val);
//     ff1.open("ca16_v.dat",ios::binary|ios::out); 
// 
//     for(int j=0;j<traceSize;j++){
//         ff1.write((char *)&v1[j],size1);
//     }
// 
//     ff1.close();
//     
    
    //##################
    

    
    

    
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    ifstream ff1_k;
    ifstream ff1_v;
    key kk[CAIDA16_SIZE]= {};
    val vv[CAIDA16_SIZE] = {};
    ff1_k.open("ca16_k.dat",ios::binary|ios::in);
    ff1_v.open("ca16_v.dat",ios::binary|ios::in);

    for(int i=0;i<CAIDA16_SIZE;i++){
        ff1_k.read((char *)&kk[i],sizeof(key));
        ff1_v.read((char *)&vv[i],sizeof(val));

    }
    
    ff1_k.close();
    ff1_v.close();
    
    
    
    
    
    
    
/*    
    ifstream ff2_k;
    ifstream ff2_v;
    key ca18_k[CAIDA18_SIZE]= {};
    val ca18_v[CAIDA18_SIZE] = {};
    ff2_k.open("ca18_k.dat",ios::binary|ios::in);
    ff2_v.open("ca18_v.dat",ios::binary|ios::in);

    for(int i=0;i<CAIDA18_SIZE;i++){

        ff2_k.read((char *)&ca18_k[i],sizeof(key*));
        ff2_v.read((char *)&ca18_v[i],sizeof(key*));

    }
    
    ff2_k.close();
    ff2_v.close();
    
    
    
    vector<key*> keys = {univ1_k, ca16_k, ca18_k};
    vector<val*> values = {univ1_v, ca16_v, ca18_v};*/

    



    
    
    
    int k = 1;
    double maxl = 0.6;
    double gamma = 1.0;
    
    ofstream ostream;
    setupOutputFile("time.raw_res", ostream, false);
    ofstream ostream1;
    setupOutputFile("nrmse.raw_res", ostream1, false);
    
//      vector<key*>::iterator k_it = keys.begin();
//     vector<val*>::iterator v_it = values.begin();
//     vector<int>::iterator s_it = sizes.begin();
//     vector<string>::iterator d_it = datasets.begin();
    
    
//     key* kk = *k_it;
//     val* vv = *v_it;
//     int size = *s_it;   
//     string dataset = *d_it;

    /*
    for (int run = 0; run < k; run++) {
        vector<key*>::iterator k_it = keys.begin();
        vector<val*>::iterator v_it = values.begin();
        vector<int>::iterator s_it = sizes.begin();
        vector<string>::iterator d_it = datasets.begin();
        
        for (int trc = 0; trc < 1; trc++) {
            key* kk = *k_it;
            val* vv = *v_it;
            int size = *s_it;
            string dataset = *d_it;
            */
//             cout<< "****  "<<dataset<<"  ****" << endl ;






int size = CAIDA16_SIZE;
double t = run(size, kk , vv,3, maxl, gamma);
float nrmse = runNrmse(size, kk , vv,3, maxl, gamma);

            
            
            
            
            
    int hh_nr_bins = 1 << 18;
    int q = hh_nr_bins * 12;
    int byteSize = 8*q;

    cout << "nrmse " << nrmse << endl;
    printf("time %f\n", t);
    cout<< byteSize << endl;
            
//             ostream1 << "Cuckoo_MC" << "      " << maxl <<  "         " << 0.83  << "        "  << byteSize <<  "        " << nrmse << endl;
//             
//             ostream << "Cuckoo_MC" << "      " << maxl <<  "         " << 0.83  << "        "  << byteSize <<  "        " << t << endl;
//             
            
       ostream1 << byteSize <<  "        " << nrmse << endl;
            
       ostream  << byteSize <<  "        " << t << endl;
            
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
            
//             ++k_it;
//             ++v_it;
//             ++s_it;
//             ++d_it;
            
//         }
//         
//     }
    ostream.close();
    ostream1.close();
    return 0;   
    
}





