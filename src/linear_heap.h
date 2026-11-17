#pragma once


#include <vector>
#include "graph.h"

class ListLinearHeap {
	ui n, key_cap, min_key;
	std::vector<ui> key, head, pre, next;
public:
	ListLinearHeap(ui n_, ui key_cap_, ui cnt, const ui *ids, const ui *keys):
		n(n_), key_cap(key_cap_), min_key(key_cap_), key(n_), head(key_cap_+1, n_), pre(n_), next(n_) {
		for(ui i = 0;i < cnt;i ++) {
			const ui id = ids[i], x = keys[id];
			key[id] = x; pre[id] = n; next[id] = head[x];
			if(head[x] != n) pre[head[x]] = id;
			head[x] = id;
			if(x < min_key) min_key = x;
		}
	}
	ui get_key(ui id) const { return key[id]; }
	bool get_min(ui &id, ui &x) {
		while(min_key <= key_cap&&head[min_key] == n) ++ min_key;
		if(min_key > key_cap) return false;
		id = head[min_key]; x = min_key;
		return true;
	}
	bool pop_min(ui &id, ui &x) {
		if(!get_min(id, x)) return false;
		key[id] = key_cap+1;
		head[x] = next[id];
		if(head[x] != n) pre[head[x]] = n;
		return true;
	}
	void decrement(ui id, ui dec) {
		if(key[id] > key_cap) return ;
		if(pre[id] == n) { head[key[id]] = next[id]; if(next[id] != n) pre[next[id]] = n; }
		else { next[pre[id]] = next[id]; if(next[id] != n) pre[next[id]] = pre[id]; }
		const ui x = key[id] -= dec;
		pre[id] = n; next[id] = head[x];
		if(head[x] != n) pre[head[x]] = id;
		head[x] = id;
		if(x < min_key) min_key = x;
	}
};

