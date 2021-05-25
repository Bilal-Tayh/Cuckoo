#pragma once
#include <cstdint>
#ifndef CUCKOO_HH_H
#define CUCKOO_HH_H

#include <emmintrin.h>
#include "xxhash.h"
#include <unordered_map>
#include <iostream>
#include <chrono>
#include <algorithm>
#include "RngFast.hpp"
#include <queue>

using namespace std;

__m256i m1;                                            // multiplier used in fast division
__m256i d1;
uint32_t shf1;
uint32_t shf2;

__m256i m11;                                            // multiplier used in fast division
__m256i d11;
uint32_t shf11;
uint32_t shf12;



union U256 {
    int32_t i[8];
    __m256i ymm;
};



void initFastMod(uint32_t d) {                                 // Set or change divisor, calculate parameters
    uint32_t L, L2, sh1, sh2, m;
    switch (d) {
    case 0:
        m = sh1 = sh2 = 1 / d;                         // provoke error for d = 0
        break;
    case 1:
        m = 1; sh1 = sh2 = 0;                          // parameters for d = 1
        break;
    case 2:
        m = 1; sh1 = 1; sh2 = 0;                       // parameters for d = 2
        break;
    default:                                           // general case for d > 2
        L = ceil(log2(d));//bit_scan_reverse(d - 1) + 1;                  // ceil(log2(d))
        L2 = L < 32 ? 1 << L : 0;                      // 2^L, overflow to 0 if L = 32
        m = 1 + uint32_t((uint64_t(L2 - d) << 32) / d); // multiplier
        sh1 = 1;  sh2 = L - 1;                         // shift counts
    }
    m1 = _mm256_set1_epi32(m);
    d1 = _mm256_set1_epi32(d);
    shf1 = sh1;
    shf2 = sh2;
}


void initFastMod1(uint32_t d) {                                 // Set or change divisor, calculate parameters
    uint32_t L, L2, sh1, sh2, m;
    switch (d) {
    case 0:
        m = sh1 = sh2 = 1 / d;                         // provoke error for d = 0
        break;
    case 1:
        m = 1; sh1 = sh2 = 0;                          // parameters for d = 1
        break;
    case 2:
        m = 1; sh1 = 1; sh2 = 0;                       // parameters for d = 2
        break;
    default:                                           // general case for d > 2
        L = ceil(log2(d));//bit_scan_reverse(d - 1) + 1;                  // ceil(log2(d))
        L2 = L < 32 ? 1 << L : 0;                      // 2^L, overflow to 0 if L = 32
        m = 1 + uint32_t((uint64_t(L2 - d) << 32) / d); // multiplier
        sh1 = 1;  sh2 = L - 1;                         // shift counts
    }
    m11 = _mm256_set1_epi32(m);
    d11 = _mm256_set1_epi32(d);
    shf11 = sh1;
    shf12 = sh2;
}

constexpr int hh_nr_bins = 1 << 21;
constexpr int hh_nr_bins_over_two = hh_nr_bins >> 1;

template<class K, class V>
class Cuckoo_waterLevel_HH_no_FP_SIMD_256 {
private:
	K _keys[hh_nr_bins_over_two][8] = {};
	V _values[hh_nr_bins_over_two][8] = {};
	V _water_level_plus_1;
	__m256i _wl_item;

	int _insertions_till_maintenance;
	int _insertions_till_maintenance_restart;
	int _max_load;
    int counter;
	int _seed1;
	int _seed2;
	int _seed3;

	int _mask;
    
    int _Z;
    int _k;
    double _alpha;
    long double _delta;
    float _gamma;
    U256 eight_samples;
    U256 eight_samples1;
    rng::rng128 gen_arr;
    
    
    void updateZK();
    void get8samples();
    void get8samples1();
	void move_to_bin(K key, V value, int bin);
	void set_water_level(V wl);
	void maintenance();
public:
	V _water_level;
	void test_correctness();
	void test_speed();
	Cuckoo_waterLevel_HH_no_FP_SIMD_256(int seed1, int seed2, int seed3, float max_load,float gamma);
	~Cuckoo_waterLevel_HH_no_FP_SIMD_256();
	void insert(K key, V value);
	V query(K key);
};


template<class K, class V>
Cuckoo_waterLevel_HH_no_FP_SIMD_256<K, V>::Cuckoo_waterLevel_HH_no_FP_SIMD_256(int seed1, int seed2, int seed3, float max_load, float gamma) {

    max_load = max_load/(1+gamma);
    
	_seed1 = seed1;
	_seed2 = seed2;
	_seed3 = seed3;
    counter=0;
	_mask = hh_nr_bins_over_two - 1;
	_water_level = 0;
	_water_level_plus_1 = 1;
	_wl_item = _mm256_set1_epi32((int)_water_level_plus_1);
	_max_load = hh_nr_bins_over_two * 8 * max_load;
	_insertions_till_maintenance_restart = hh_nr_bins_over_two * 8 * (max_load);
	_insertions_till_maintenance = _max_load;
	srand(42);
    
        
    
    _delta = 0.99;
    _alpha = 0.83;
    _gamma = gamma;
    
    
    
    
    
    _insertions_till_maintenance_restart *= (2*_alpha) - 1;
    
    _k = ceil(((_alpha * _gamma * (2 + _gamma - _alpha * _gamma)) / (pow(_gamma - _alpha * _gamma, 2))) * log(1 / _delta));
    _Z = (int)((_k * (1 + _gamma)) / (_alpha * _gamma));
    if (_Z & 0b111) // integer multiple of 8
        _Z = _Z - (_Z & 0b111) + 8;
    _k = _Z * (_alpha * _gamma) / (1 + _gamma)+0.5;
    
    
    initFastMod(hh_nr_bins_over_two);
    initFastMod1(8);
}

template<class K, class V>
Cuckoo_waterLevel_HH_no_FP_SIMD_256<K, V>::~Cuckoo_waterLevel_HH_no_FP_SIMD_256() {

}

template<class K, class V>
void Cuckoo_waterLevel_HH_no_FP_SIMD_256<K, V>::insert(K key, V value)
{
	uint64_t hash = XXH64((void *) &key, sizeof(key), _seed1);
	uint64_t i = 1;
	uint64_t bin1 = hash & _mask;
	int empty_spot_in_bin1 = -1;
	int empty_spot_in_bin2 = -1;
    
	//const __m256i wl_item = _mm256_set1_epi32(_water_level_plus_1);
	const __m256i item = _mm256_set1_epi32((int)key);

	const __m256i match1 = _mm256_cmpeq_epi32(item, *((__m256i*) _keys[bin1]));
	const int     mask1 = _mm256_movemask_epi8(match1);
	if (mask1 != 0) {
		int tz1 = _tzcnt_u32(mask1);
		_values[bin1][tz1 >> 2] += value;
		return;
	}

	uint64_t bin2 = (hash >> 32) & _mask;
	const __m256i match2 = _mm256_cmpeq_epi32(item, *((__m256i*) _keys[bin2]));
	const int     mask2 = _mm256_movemask_epi8(match2);
	if (mask2 != 0) {
		int tz2 = _tzcnt_u32(mask2);
		_values[bin2][tz2 >> 2] += value;
		return;
	}

	if (!--_insertions_till_maintenance) {
		maintenance();
	}

	//const __m256i wl_match1 = _mm256_cmple_epi32(*((__m256i*) _values[bin1]), wl_item);
	const __m256i wl_match1 = _mm256_cmpgt_epi32(_wl_item, *((__m256i*) _values[bin1]));
	const int     wl_mask1 = _mm256_movemask_epi8(wl_match1);

	// Should we add deltas?
	if (wl_mask1 != 0) {
		int tz1 = _tzcnt_u32(wl_mask1);
		empty_spot_in_bin1 = tz1 >> 2;
		_keys[bin1][empty_spot_in_bin1] = key;
		_values[bin1][empty_spot_in_bin1] = _water_level+value;
		return;
	}


	//uint64_t bin2 = (hash >> 32) & _mask;


	const __m256i wl_match2 = _mm256_cmpgt_epi32(_wl_item, *((__m256i*) _values[bin2]));
	const int     wl_mask2 = _mm256_movemask_epi8(wl_match2);
	if (wl_mask2 != 0) {
		int tz2 = _tzcnt_u32(wl_mask2);
		empty_spot_in_bin2 = tz2 >> 2;
		_keys[bin2][empty_spot_in_bin2] = key;
		_values[bin2][empty_spot_in_bin2] = _water_level + value;
		return;
	}
	K kicked_key;
	V kicked_value;
	uint64_t kicked_from_bin;
	int coinFlip = rand();
	if (coinFlip & 0b1) {
		int kicked_index = (coinFlip >> 1) & 0b111;
		kicked_key = _keys[bin1][kicked_index];
		kicked_value = _values[bin1][kicked_index];

		_keys[bin1][kicked_index] = key;
		_values[bin1][kicked_index] = _water_level + value;
		kicked_from_bin = bin1;
		//insert(kicked_key, kicked_value);
	}
	else {
		int kicked_index = (coinFlip >> 1) & 0b111;
		kicked_key = _keys[bin2][kicked_index];
		kicked_value = _values[bin2][kicked_index];
		_keys[bin2][kicked_index] = key;
		_values[bin2][kicked_index] = _water_level + value;
		kicked_from_bin = bin2;
	}
	uint64_t kicked_hash = XXH64((void *) &kicked_key, sizeof(kicked_key), _seed1);
	uint64_t kicked_bin_1 = kicked_hash & _mask;

	if (kicked_bin_1 != kicked_from_bin) {
		move_to_bin(kicked_key, kicked_value, kicked_bin_1);
	}
	else {
		uint64_t kicked_bin_2 = (kicked_hash >> 32) & _mask;
		move_to_bin(kicked_key, kicked_value, kicked_bin_2);
	}
}

template<class K, class V>
void Cuckoo_waterLevel_HH_no_FP_SIMD_256<K, V>::move_to_bin(K key, V value, int bin )
{
    
//     counter++;
//     if(counter >=5){
//         maintenance();
//         counter=0;
//     }
    
    
    
	const __m256i wl_match = _mm256_cmpgt_epi32(_wl_item, *((__m256i*) _values[bin]));
	const int     wl_mask = _mm256_movemask_epi8(wl_match);
	if (wl_mask != 0) {
		int tz = _tzcnt_u32(wl_mask);
		int empty_spot_in_bin = tz >> 2;
		_keys[bin][empty_spot_in_bin] = key;
		_values[bin][empty_spot_in_bin] = value;
//         counter = 0;
		return;
	}

	int kicked_index = rand() & 0b111;
	K kicked_key = _keys[bin][kicked_index];
	V kicked_value = _values[bin][kicked_index];
	_keys[bin][kicked_index] = key;
	_values[bin][kicked_index] = value;

	uint64_t kicked_hash = XXH64((void *) &kicked_key, sizeof(kicked_key), _seed1);
	uint64_t kicked_bin_1 = kicked_hash & _mask;

	uint64_t kicked_from_bin = bin;

	if (kicked_bin_1 != kicked_from_bin) {
		move_to_bin(kicked_key, kicked_value, kicked_bin_1);
	}
	else {
		uint64_t kicked_bin_2 = (kicked_hash >> 32) & _mask;
		move_to_bin(kicked_key, kicked_value, kicked_bin_2);
	}
}

template<class K, class V>
V Cuckoo_waterLevel_HH_no_FP_SIMD_256<K, V>::query(K key)
{
	uint64_t hash = XXH64((void *) &key, sizeof(key), _seed1);
	uint64_t bin1 = hash & _mask;
	const __m256i item = _mm256_set1_epi32((int)key);
	const __m256i match1 = _mm256_cmpeq_epi32(item, *((__m256i*) _keys[bin1]));
	const int mask = _mm256_movemask_epi8(match1);
	if (mask != 0) {
		int tz1 = _tzcnt_u32(mask);
		return _values[bin1][tz1 >> 2];
	}
	uint64_t bin2 = (hash >> 32) & _mask;
	const __m256i match2 = _mm256_cmpeq_epi32(item, *((__m256i*) _keys[bin2]));
	const int mask2 = _mm256_movemask_epi8(match2);
	if (mask2 != 0) {
		//cout << tz1 << endl;
		int tz2 = _tzcnt_u32(mask2);
		return _values[bin2][tz2 >> 2];
	}
	return _water_level;
}


template<class K, class V>
inline void Cuckoo_waterLevel_HH_no_FP_SIMD_256<K, V>::set_water_level(V wl)
{
	_water_level = wl;
	_water_level_plus_1 = wl + 1;
	_wl_item = _mm256_set1_epi32(_water_level_plus_1);
}





template<class K, class V>
void Cuckoo_waterLevel_HH_no_FP_SIMD_256<K, V>::get8samples() {
    int indx = 0;
    eight_samples.ymm = _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
    __m256i t1 = _mm256_mul_epu32(eight_samples.ymm, m1);                           // 32x32->64 bit unsigned multiplication of even ele/*ment*/s of a
    __m256i t2 = _mm256_srli_epi64(t1, 32);                         // high dword of even numbered results
    __m256i t3 = _mm256_srli_epi64(eight_samples.ymm, 32);                          // get odd elements of a into position for multiplication
    __m256i t4 = _mm256_mul_epu32(t3, m1);                          // 32x32->64 bit unsigned multiplication of odd elements
    __m256i t5 = _mm256_set_epi32(-1, 0, -1, 0, -1, 0, -1, 0);      // mask for odd elements
    __m256i t7 = _mm256_blendv_epi8(t2, t4, t5);                    // blend two results
    __m256i t8 = _mm256_sub_epi32(eight_samples.ymm, t7);                           // subtract
    __m256i t9 = _mm256_srli_epi32(t8, shf1);                       // shift right logical
    __m256i t10 = _mm256_add_epi32(t7, t9);                         // add
    __m256i t11 = _mm256_srli_epi32(t10, shf2);                     // shift right logical 
    __m256i t12 = _mm256_mullo_epi32(t11, d1);                      // multiply quotient with divisor
    eight_samples.ymm = _mm256_sub_epi32(eight_samples.ymm, t12);                   // subtract
}

template<class K, class V>
void Cuckoo_waterLevel_HH_no_FP_SIMD_256<K, V>::get8samples1() {
    int indx = 0;
    eight_samples1.ymm = _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
    __m256i t1 = _mm256_mul_epu32(eight_samples1.ymm, m11);                           // 32x32->64 bit unsigned multiplication of even elements of a
    __m256i t2 = _mm256_srli_epi64(t1, 32);                         // high dword of even numbered results
    __m256i t3 = _mm256_srli_epi64(eight_samples1.ymm, 32);                          // get odd elements of a into position for multiplication
    __m256i t4 = _mm256_mul_epu32(t3, m11);                          // 32x32->64 bit unsigned multiplication of odd elements
    __m256i t5 = _mm256_set_epi32(-1, 0, -1, 0, -1, 0, -1, 0);      // mask for odd elements
    __m256i t7 = _mm256_blendv_epi8(t2, t4, t5);                    // blend two results
    __m256i t8 = _mm256_sub_epi32(eight_samples1.ymm, t7);                           // subtract
    __m256i t9 = _mm256_srli_epi32(t8, shf11);                       // shift right logical
    __m256i t10 = _mm256_add_epi32(t7, t9);                         // add
    __m256i t11 = _mm256_srli_epi32(t10, shf12);                     // shift right logical 
    __m256i t12 = _mm256_mullo_epi32(t11, d11);                      // multiply quotient with divisor
    eight_samples1.ymm = _mm256_sub_epi32(eight_samples1.ymm, t12);                   // subtract
}





template<class K, class V>
inline void Cuckoo_waterLevel_HH_no_FP_SIMD_256<K, V>::maintenance()
{
//     cout<< "maintenance" << endl;
// 	int size = hh_nr_bins_over_two * 8;
// 	V copy_of_values[hh_nr_bins_over_two][8];
// 	memcpy(copy_of_values, _values, size * sizeof(V));
// 	V* one_array = (V*)copy_of_values;
// 	sort(one_array, one_array + size);
// 	set_water_level(one_array[size >> 1]);
// // 	cout << "water_level = " << _water_level << endl;
// 	_insertions_till_maintenance = _insertions_till_maintenance_restart;
    
    
    
    
    std::priority_queue <int> p;
    int left_to_fill = _k;
    while (left_to_fill > 0) {
        get8samples();
        get8samples1();
        for (int i = 0; i < 8; ++i) {
            
            p.push(_values[eight_samples.i[i]][eight_samples1.i[i]]);
        }
        left_to_fill -= 8;
    }

    int left_to_sample = _Z-_k+ left_to_fill;
    while (left_to_fill++ < 0) {
        p.pop();
    }
    int top = p.top();

    while (left_to_sample > 0) {
        get8samples();
        get8samples1();
        for (int i = 0; i < 8; ++i) {
            if (top > _values[eight_samples.i[i]][eight_samples1.i[i]]) {
                p.pop();
                p.push(_values[eight_samples.i[i]][eight_samples1.i[i]]);
                top = p.top();
            }
        }
        left_to_sample -= 8;
    }
    
    updateZK();
    
    set_water_level(top);
	_insertions_till_maintenance = _insertions_till_maintenance_restart;
}





template<class K, class V>
void Cuckoo_waterLevel_HH_no_FP_SIMD_256<K, V>::updateZK()
{
    _delta=_delta/2.0;
    _k=ceil(((_alpha*_gamma*(2+_gamma - _alpha*_gamma)) / (pow(_gamma -_alpha*_gamma,2))) * log(1/_delta));
    _Z = (int)( (_k*(1+_gamma)) / (_alpha*_gamma) );
    if (_Z & 0b111) // integer multiple of 8
        _Z = _Z - (_Z & 0b111) + 8;
    _k = _Z * (_alpha * _gamma) / (1 + _gamma)+0.5;
//     cout<< "@@@@@@@@@@@@@@@"<<endl;
//     cout<< "delta "<< _delta <<endl;
//     cout<< "Z " <<_Z <<endl;
//     cout<< "k " <<_k <<endl;
//     cout<< "insersions " <<counter <<endl;
//     counter=0;
}



template<class K, class V>
void Cuckoo_waterLevel_HH_no_FP_SIMD_256<K, V>::test_correctness()
{
	unordered_map<K, V> true_map;
	uint64_t S = 0;
	//float eps = 0.01;

	for (int i = 0; i < (1 << 21); ++i) {
		uint32_t id = rand();
		//uint32_t val = (rand() << 15) | rand();
		uint32_t val =  rand();
		true_map[id] += val;
		S += val;
		insert(id, val);
		//cout << "id = " << id << " val = " << val << " query(id) = " << query(id) << endl;
	}
	int nr_counters = hh_nr_bins_over_two * 8;
	float eps = 2.0 / nr_counters;
	for (auto it = true_map.begin(); it != true_map.end(); ++it) {
		if ((it->second > query(it->first)) || (it->second < query(it->first) - S*eps)) {
			cout << it->first << " " << it->second << " " << query(it->first) << endl;
			exit(1);
		}
	}
}

template<class K, class V>
void Cuckoo_waterLevel_HH_no_FP_SIMD_256<K, V>::test_speed()
{
	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < 10; ++i)
	{
		for (int64_t j = 0; j < 3.6 * hh_nr_bins; ++j) {
			insert(j, i + j);
		}
		//for (int64_t j = 0; j < 3.6 * hh_nr_bins; ++j) {
		//	try_remove(j);
		//}
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;
}
#endif //CUCKOO_HH_H
