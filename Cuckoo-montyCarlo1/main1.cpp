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

#define CLK_PER_SEC CLOCKS_PER_SEC
#define CAIDA16_SIZE 152197439
#define CAIDA18_SIZE 175880896
#define UNIV1_SIZE 17323447

typedef unsigned long long key;
typedef double val;

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
      file_keys[i] = stoull(id);
      file_ws[i] = stod(len);
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




int main(int argc, char* argv[])
{
//     vector<key*> keys;
//     vector<val*> values;
//     vector<int> sizes;
//     vector<string> datasets;
    

//     getKeysAndWeightsFromFile("../../datasets/UNIV1/mergedAggregatedPktlen_Srcip", keys, values, UNIV1_SIZE);
//     sizes.push_back(UNIV1_SIZE);
//     datasets.push_back("univ1");
//     


//     getKeysAndWeightsFromFile("../../datasets/CAIDA16/mergedAggregatedPktlen_Srcip", keys, values, CAIDA16_SIZE);
//     sizes.push_back(CAIDA16_SIZE);
//     datasets.push_back("caida");



//     getKeysAndWeightsFromFile("../../datasets/CAIDA18/mergedAggregatedPktlen_Srcip", keys, values, CAIDA18_SIZE);
//     sizes.push_back(CAIDA18_SIZE);
//     datasets.push_back("caida18");


    
    
    
    
    
    
    
/*    
    getKeysAndWeightsFromFile("../../datasets/UNIV1/mergedAggregatedPktlen_Srcip", keys, values, UNIV1_SIZE);
    sizes.push_back(UNIV1_SIZE);
    datasets.push_back("univ1");

    
    
    vector<key*>::iterator k_it = keys.begin();
    key* k = *k_it;

    int size =sizeof(key*);
    int traceSize = UNIV1_SIZE;
    ofstream f;

    f.open("univ1_k.dat",ios::binary|ios::out); 

    for(int j=0;j<traceSize;j++){
        f.write((char *)&k[j],size);
    }

    f.close();
    
    
    
    
    vector<val*>::iterator v_it = values.begin();
    val* v = *v_it;

    int size1 =sizeof(val*);
    ofstream ff;

    ff.open("univ1_v.dat",ios::binary|ios::out); 

    for(int j=0;j<traceSize;j++){
        ff.write((char *)&v[j],size1);
    }

    ff.close();
    */
    
    
    
    
    
//##################
    
    
    
    
    
    
//     getKeysAndWeightsFromFile("../../datasets/CAIDA16/mergedAggregatedPktlen_Srcip", keys, values, CAIDA16_SIZE);
//     sizes.push_back(CAIDA16_SIZE);
//     datasets.push_back("caida");
// 
// 
//     
//     
//     vector<key*>::iterator k_it = keys.begin();
//     key* k = *k_it;
// 
//     int size =sizeof(key*);
//     int traceSize = CAIDA16_SIZE;
//     ofstream f1;
// 
//     f1.open("ca16_k.dat",ios::binary|ios::out); 
// 
//     for(int j=0;j<traceSize;j++){
//         f1.write((char *)&k[j],size);
//     }
// 
//     f1.close();
//     
//     
//     
//     
//     vector<val*>::iterator v_it = values.begin();
//     val* v = *v_it;
// 
// 
//     ofstream ff1;
//     int size1 =sizeof(val*);
//     ff1.open("ca16_v.dat",ios::binary|ios::out); 
// 
//     for(int j=0;j<traceSize;j++){
//         ff1.write((char *)&v[j],size1);
//     }
// 
//     ff1.close();
    
    
    //##################
    
    
    
    
//     getKeysAndWeightsFromFile("../../datasets/CAIDA18/mergedAggregatedPktlen_Srcip", keys, values, CAIDA18_SIZE);
//     sizes.push_back(CAIDA18_SIZE);
//     datasets.push_back("caida18");
//     
//     
//     vector<key*>::iterator k_it = keys.begin();
//     key* k = *k_it;
// 
//     int size =sizeof(key*);
//     int traceSize = CAIDA18_SIZE;
//     ofstream f1;
// 
//     f1.open("ca18_k.dat",ios::binary|ios::out); 
// 
//     for(int j=0;j<traceSize;j++){
//         f1.write((char *)&k[j],size);
//     }
// 
//     f1.close();
//     
//     
//     
//     
//     vector<val*>::iterator v_it = values.begin();
//     val* v = *v_it;
// 
// 
//     ofstream ff1;
//     int size1 =sizeof(val*);
//     ff1.open("ca18_v.dat",ios::binary|ios::out); 
// 
//     for(int j=0;j<traceSize;j++){
//         ff1.write((char *)&v[j],size1);
//     }
// 
//     ff1.close();
    
    
    
    
    
    
    
    
    
    
    
    
    vector<string> datasets;
    datasets.push_back("univ1");
    datasets.push_back("caida16");
    datasets.push_back("caida18");
    vector<int> sizes;
    sizes.push_back(UNIV1_SIZE);
    sizes.push_back(CAIDA16_SIZE);
    sizes.push_back(CAIDA16_SIZE);
        
    
    
    
    
    
    
    ifstream ff_k;
    ifstream ff_v;
    key univ1_k[UNIV1_SIZE]= {};
    val univ1_v[UNIV1_SIZE] = {};
    ff_k.open("univ1_k.dat",ios::binary|ios::in);
    ff_v.open("univ1_v.dat",ios::binary|ios::in);

    for(int i=0;i<UNIV1_SIZE;i++){

        ff_k.read((char *)&univ1_k[i],sizeof(key*));
        ff_v.read((char *)&univ1_v[i],sizeof(key*));

    }
    
    ff_k.close();
    ff_v.close();
    
    
    
    
    
    
    
    
    ifstream ff1_k;
    ifstream ff1_v;
    key ca16_k[CAIDA16_SIZE]= {};
    val ca16_v[CAIDA16_SIZE] = {};
    ff1_k.open("ca16_k.dat",ios::binary|ios::in);
    ff1_v.open("ca16_v.dat",ios::binary|ios::in);

    for(int i=0;i<CAIDA16_SIZE;i++){
        ff1_k.read((char *)&ca16_k[i],sizeof(key*));
        ff1_v.read((char *)&ca16_v[i],sizeof(key*));

    }
    
    ff1_k.close();
    ff1_v.close();
    
    
    
    
    

    
    
    
    
    
    
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
    vector<val*> values = {univ1_v, ca16_v, ca18_v};

    



    
    
    
    int k = 1;
    double maxl = 0.5;
    double gamma = 1.0;
    
    ofstream ostream;
    setupOutputFile("time.raw_res", ostream, false);
    ofstream ostream1;
    setupOutputFile("nrmse.raw_res", ostream1, false);

    
    for (int run = 0; run < k; run++) {
        vector<key*>::iterator k_it = keys.begin();
        vector<val*>::iterator v_it = values.begin();
        vector<int>::iterator s_it = sizes.begin();
        vector<string>::iterator d_it = datasets.begin();
        
        for (int trc = 0; trc < 3; trc++) {
            key* kk = *k_it;
            val* vv = *v_it;
            int size = *s_it;
            string dataset = *d_it;
            
//             cout<< "****  "<<dataset<<"  ****" << endl ;
            
            
            struct timeb begintb, endtb;
            clock_t begint, endt;
            Cuckoo_waterLevel_HH_no_FP_SIMD_256<uint32_t, uint32_t> hh(1, 2, 3, maxl,gamma,size);
            begint = clock();
//             ftime(&begintb);
            for (int i = 0; i < size; ++i) {
                hh.insert(kk[i], vv[i]);
            }
            endt = clock();
//             ftime(&endtb);
            double time = ((double)(endt-begint))/CLK_PER_SEC;

            
                
//             double c =0.0; 
//             uint64_t vol =0;
//             unordered_map<key, val> map;
//             Cuckoo_waterLevel_HH_no_FP_SIMD_256<uint32_t, uint32_t> hh1(1, 2, 3, maxl,gamma,size);
//             for (int i = 0; i < size; ++i) {
//                     map[kk[i]] +=vv[i];
//                     hh1.insert(kk[i], vv[i]);
//                     uint32_t err = map[kk[1]] - hh1.query(kk[i]);
//                     c+= err*err;
//                     vol+=vv[i];
//             }
//             double mse = c/size;
//             double rmse = sqrt(mse);
//             double nrmse = rmse/vol;
//             
//             
//             
            int hh_nr_bins = 1 << 14;
            int hh_nr_bins_over_two = hh_nr_bins >> 1;
            int q = hh_nr_bins_over_two * 8;
            
            
            ostream << "Cuckoo_MC" << dataset<<   "     " << q <<  "    " << time << endl;
//             ostream1 << "Cuckoo_MC" << dataset<<  "     " << q << "     " << nrmse << endl;
            
            
            ++k_it;
            ++v_it;
            ++s_it;
            ++d_it;
            
        }
        
    }
    ostream.close();
    ostream1.close();
    return 0;   
    
}





