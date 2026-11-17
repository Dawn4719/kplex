#pragma once

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <new>
#include <utility>
#include <immintrin.h>

class CuckooHash {
	using int32 = int32_t;
	static const int32 unfilled = -1;
	static const int32 buff_size = 4;

	int32 capacity = 0;
	int32 mask = 0;
	int32 size = 0;
	int32 *hashtable = nullptr;

	static int32 *allocate(int32 cap) {
		void *ptr = nullptr;
		if(posix_memalign(&ptr, 16, sizeof(int32)*size_t(cap)) != 0) throw std::bad_alloc();
		memset(ptr, unfilled, sizeof(int32)*size_t(cap));
		return static_cast<int32 *>(ptr);
	}

	int32 hash1(const int32 x) const { return x & mask; }
	int32 hash2(const int32 x) const { return ~x & mask; }

	void rehash(int32 **table) {
		int32 old_capacity = capacity;
		mask = mask == 0? 1 : ((mask << 1) | 1);
		capacity = (mask+1)*buff_size;
		int32 *new_table = allocate(capacity);
		for(int32 i = 0;i < old_capacity;++ i) {
			if((*table)[i] != unfilled) insert((*table)[i], &new_table);
		}
		std::swap(*table, new_table);
		free(new_table);
	}

	void insert(const int32 &value, int32 **table) {
		int32 hs = hash1(value);
		for(int32 i = 0;i < buff_size;++ i) {
			if((*table)[hs*buff_size+i] == unfilled) { (*table)[hs*buff_size+i] = value; return; }
		}
		hs = hash2(value);
		for(int32 i = 0;i < buff_size;++ i) {
			if((*table)[hs*buff_size+i] == unfilled) { (*table)[hs*buff_size+i] = value; return; }
		}

		bool use_hash1 = true;
		int32 u = value;
		for(int32 i = 0;i < mask;++ i) {
			int32 replaced;
			hs = use_hash1? hash1(u) : hash2(u);
			int32 j = 0;
			for(;j < buff_size;++ j) if((*table)[hs*buff_size+j] == unfilled) break;
			if(j == buff_size) {
				replaced = (*table)[hs*buff_size];
				for(j = 1;j < buff_size;++ j) (*table)[hs*buff_size+j-1] = (*table)[hs*buff_size+j];
				(*table)[hs*buff_size+j-1] = u;
			}
			else {
				replaced = (*table)[hs*buff_size+j];
				(*table)[hs*buff_size+j] = u;
			}
			use_hash1 = hs == hash2(replaced);
			u = replaced;
			if(u == unfilled) return;
		}
		rehash(table);
		insert(u, table);
	}

public:
	CuckooHash() {}
	CuckooHash(const CuckooHash &) = delete;
	CuckooHash &operator=(const CuckooHash &) = delete;
	CuckooHash(CuckooHash &&other) noexcept: capacity(other.capacity), mask(other.mask), size(other.size), hashtable(other.hashtable) {
		other.hashtable = nullptr; other.capacity = other.mask = other.size = 0;
	}
	~CuckooHash() { if(hashtable != nullptr) free(hashtable); }

	void reserve(int32 expected) {
		if(hashtable != nullptr&&capacity >= expected) return;
		mask = mask == 0? 1 : ((mask << 1) | 1);
		while(expected >= mask*buff_size) mask = (mask << 1) | 1;
		capacity = (mask+1)*buff_size;
		if(hashtable != nullptr) free(hashtable);
		hashtable = allocate(capacity);
	}

	void insert(const int32 &value) {
		if(find(value)) return;
		insert(value, &hashtable);
		++ size;
	}

	bool find(const int32 &value) const {
		__m128i cmp = _mm_set1_epi32(value);
		__m128i b1 = _mm_load_si128((const __m128i *)&hashtable[buff_size*hash1(value)]);
		__m128i b2 = _mm_load_si128((const __m128i *)&hashtable[buff_size*hash2(value)]);
		__m128i flag = _mm_or_si128(_mm_cmpeq_epi32(cmp, b1), _mm_cmpeq_epi32(cmp, b2));
		return _mm_movemask_epi8(flag) != 0;
	}
};

