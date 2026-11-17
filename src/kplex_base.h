#pragma once


#include <algorithm>
#include <memory>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include "graph.h"
#include "linear_heap.h"
#include "cuckoo_hash.h"

class KPlexBase {
public:
	typedef unsigned long long Count;
	Count answer = 0;
	KPlexBase(const Graph &g, int k_, int q_): n(g.n), k(k_), q(q_), in_pstart(g.pstart.data()), in_edges(g.edges.data()) {}
	virtual ~KPlexBase() {}
	void run();

protected:
	typedef uint64_t u64;
	ui n;
	int k, q;
	const ept *in_pstart;
	const ui *in_edges;

	std::vector<ept> st;
	std::vector<ui> ed;
	std::vector<CuckooHash> hash;
	void reorder();
	void build_hash();
	void peel(ui root, ui &active_n, std::vector<char> &active, ListLinearHeap &heap);
	void compact(std::vector<char> &active, ui &active_n);
	void alloc_vertex_arrays();
public:
	Graph *release_input = nullptr;
protected:

	std::vector<char> tag, grp;
	std::vector<int> cn, cn2;
	std::vector<ui> N1, N2, T2, queue_buf, adjA, adjB, radj, S2, e22, a22, a22start, cursor, pos2, e11;
	std::vector<size_t> seg_begin, seg_end, bseg_begin, bseg_end;
	std::vector<ui> hc, touched;
	std::vector<ui> L, gid, lstart;
	std::vector<int> lid, ladj;
	void process_root(ui r, const std::vector<char> &active);
	bool is_onehop(size_t i) const { return tag[gid[i]] == 1; }
	virtual void run_root(size_t s) = 0;

	std::vector<int> M;

	int resid[64];
	std::vector<u64> capm;
	std::vector<int> capg;
	int pool_ncap = 0, pool_r = 0;
	Count pool_free = 0;
	static const int PT_MAX = 16;
	u64 pt_mask[PT_MAX];
	int pt_size[PT_MAX], pt_cnt[PT_MAX], pt_n = 0;
	Count pf_n[PT_MAX+2];
	int pf_r[PT_MAX+2], pf_cnt = 0;
	void comb_prod();
	void prod_enum(int i);
	void prod_comb(int i, Count start, int left);
	void type_prod(int rem);
	void comb_list(Count list_size, Count start, Count picked, Count r);
	void comb(Count list_size, Count r);
	void pool_types(int ncap, u64 one);
	void pt_counts(int g, int used);
	void pool_rec_g(int start, int t, u64 ex, u64 one);
};


inline void KPlexBase::reorder() {
	std::vector<ui> ids(n), degree(n), order(n), rank(n);
	std::vector<char> alive(n, 1);
	for(ui i = 0;i < n;i ++) { ids[i] = i; degree[i] = in_pstart[i+1]-in_pstart[i]; }
	ListLinearHeap heap(n, n-1, n, ids.data(), degree.data());
	for(ui i = 0;i < n;i ++) {
		ui u = 0, key = 0;
		heap.pop_min(u, key);
		alive[u] = 0; order[i] = u; rank[u] = i;
		for(ept j = in_pstart[u];j < in_pstart[u+1];j ++) if(alive[in_edges[j]]) heap.decrement(in_edges[j], 1);
	}
	st.assign(n+1, 0); ed.resize(in_pstart[n]);
	for(ui i = 0;i < n;i ++) { const ui u = order[i]; st[i+1] = st[i]+(in_pstart[u+1]-in_pstart[u]); }
	std::vector<ept> pos(st.begin(), st.end()-1);
	for(ui i = n;i > 0;i --) {
		const ui u = order[i-1];
		for(ept j = in_pstart[u];j < in_pstart[u+1];j ++) ed[pos[rank[in_edges[j]]]++] = i-1;
	}
}

inline void KPlexBase::build_hash() {
	hash.clear(); hash.resize(n);
	for(ui u = 0;u < n;u ++) {
		hash[u].reserve(int32_t(st[u+1]-st[u]));
		for(ept j = st[u];j < st[u+1];j ++) hash[u].insert(int32_t(ed[j]));
	}
}

inline void KPlexBase::peel(ui root, ui &active_n, std::vector<char> &active, ListLinearHeap &heap) {
	if(root < n) heap.decrement(root, heap.get_key(root));
	ui u, key;
	while(heap.get_min(u, key)&&(int(key) < q-1-k||u == root)) {
		heap.pop_min(u, key); active[u] = 0; -- active_n;
		for(ept j = st[u];j < st[u+1];j ++) if(active[ed[j]]) heap.decrement(ed[j], 1);
	}
}

inline void KPlexBase::compact(std::vector<char> &active, ui &active_n) {
	std::vector<ui> nid(n, ui(-1));
	ui c = 0;
	for(ui v = 0;v < n;v ++) if(active[v]) nid[v] = c++;
	ept tot = 0;
	for(ui v = 0;v < n;v ++) if(active[v]) for(ept j = st[v];j < st[v+1];j ++) tot += active[ed[j]] != 0;
	std::vector<ept> st2(c+1, 0);
	std::vector<ui> ed2(tot);
	ept p = 0;
	for(ui v = 0;v < n;v ++) if(active[v]) {
		for(ept j = st[v];j < st[v+1];j ++) if(active[ed[j]]) ed2[p++] = nid[ed[j]];
		st2[nid[v]+1] = p;
	}
	st.swap(st2); ed.swap(ed2);
	n = c; active.assign(n, 1); active_n = n;
}

inline void KPlexBase::alloc_vertex_arrays() {
	tag.assign(n, 0); grp.assign(n, 0); cn.assign(n, 0); cn2.assign(n, 0); lid.assign(n, -1); hc.assign(n, 0);
	seg_begin.assign(n, 0); seg_end.assign(n, 0); bseg_begin.assign(n, 0); bseg_end.assign(n, 0); pos2.assign(n, 0);
}

inline void KPlexBase::run() {
	if(q < 1||k < 0||q < 2*k+1) throw std::invalid_argument("requires q >= 1, k >= 0 and q >= 2k+1");
	if(q > 64) throw std::invalid_argument("q > 64 is not supported (prefix bitmask)");
	answer = 0;
	if(n < ui(q)) return ;
	if(q == 1) { answer = n; return ; }
	reorder();
	if(release_input) { std::vector<ept>().swap(release_input->pstart); std::vector<ui>().swap(release_input->edges); in_pstart = nullptr; in_edges = nullptr; }
	std::vector<char> active(n, 1);
	std::unique_ptr<ListLinearHeap> heap;
	auto make_heap = [&]() {
		std::vector<ui> ids(n), degree(n);
		for(ui i = 0;i < n;i ++) { ids[i] = i; degree[i] = st[i+1]-st[i]; }
		heap.reset(new ListLinearHeap(n, n == 0? 0 : n-1, n, ids.data(), degree.data()));
	};
	make_heap();
	ui active_n = n;
	peel(n, active_n, active, *heap);
	if(active_n < n) { heap.reset(); compact(active, active_n); make_heap(); }
	if(k >= 2) build_hash();
	alloc_vertex_arrays();
	M.assign(q, 0);
	for(ui r = 0;r < n&&active_n >= ui(q);r ++) if(active[r]) {
		process_root(r, active);
		peel(r, active_n, active, *heap);
	}
}


inline void KPlexBase::process_root(ui r, const std::vector<char> &active) {
	const int min_deg = q-1-k, cap1 = k, cap2 = k-1, min_cn2 = q-2*k;
	N1.clear(); N2.clear();
	for(ept j = st[r];j < st[r+1]&&ed[j] > r;j ++) if(active[ed[j]]) { const ui v = ed[j]; N1.push_back(v); tag[v] = 1; grp[v] = 1; cn[v] = 0; cn2[v] = 0; }
	if(int(N1.size()) < min_deg) { for(ui v : N1) tag[v] = 0; return ; }
	auto clear_tags = [&]() { for(ui v : N1) tag[v] = 0; for(ui v : N2) tag[v] = 0; };
	auto root_ok = [&]() { int d = 0; for(ui v : N1) if(tag[v] == 1) ++ d; return d >= min_deg; };
	const size_t n1 = N1.size();
	for(ui v : N1) { const ui *e = ed.data()+st[v]; __builtin_prefetch(e); __builtin_prefetch(e+16); }

	e11.clear();
	for(size_t i = 1;i < n1;i ++) {
		const ui v = N1[i], lo = N1[i-1];
		for(ept j = st[v];j < st[v+1]&&ed[j] >= lo;j ++) {
			const ui x = ed[j];
			if(tag[x] == 1) { e11.push_back(v); e11.push_back(x); ++ cn[v]; ++ cn[x]; }
		}
	}
	size_t off = 0;
	for(ui v : N1) { seg_begin[v] = seg_end[v] = off; off += cn[v]; }
	adjA.resize(off);
	for(size_t i = 0;i < e11.size();i += 2) { const ui a = e11[i], b = e11[i+1]; adjA[seg_end[a]++] = b; adjA[seg_end[b]++] = a; }
	if(1+cap1 < min_deg) {
		queue_buf.clear();
		for(ui v : N1) if(1+cn[v]+cap1 < min_deg) { tag[v] = 3; queue_buf.push_back(v); }
		for(size_t i = 0;i < queue_buf.size();i ++) {
			const ui x = queue_buf[i];
			for(size_t j = seg_begin[x];j < seg_end[x];j ++) { const ui w = adjA[j]; if(tag[w] == 1&&1+(-- cn[w])+cap1 < min_deg) { tag[w] = 3; queue_buf.push_back(w); } }
		}
		if(!root_ok()) { clear_tags(); return ; }
	}

	size_t tot = 0;
	if(k > 0) for(ui v : N1) if(tag[v] == 1) tot += st[v+1]-st[v];
	if(touched.size() < tot) touched.resize(tot);
	ui *HC = hc.data(), *TO = touched.data();
	size_t nt = 0;
	if(k > 0) for(ui v : N1) {
		if(tag[v] != 1) continue;
		for(ept j = st[v];j < st[v+1]&&ed[j] > r;j ++) {
			const ui w = ed[j];
			if(tag[w] == 1||tag[w] == 3) continue;
			const ui c = HC[w];
			HC[w] = c+1; TO[nt] = w; nt += c == 0;
		}
	}
	off = 0;
	for(size_t i = 0;i < nt;i ++) {
		const ui w = TO[i], c = HC[w];
		HC[w] = 0;
		if(int(c) >= min_cn2&&active[w]) { N2.push_back(w); tag[w] = 2; grp[w] = 2; cn[w] = int(c); cn2[w] = 0; seg_begin[w] = seg_end[w] = off; off += c; }
	}
	radj.resize(off);
	if(!N2.empty()) for(ui v : N1) {
		if(tag[v] != 1) continue;
		for(ept j = st[v];j < st[v+1]&&ed[j] > r;j ++) { const ui w = ed[j]; if(tag[w] == 2) { radj[seg_end[w]++] = v; ++ cn2[v]; } }
	}
	off = 0;
	for(ui v : N1) { bseg_begin[v] = bseg_end[v] = off; off += cn2[v]; }
	adjB.resize(off);
	for(ui w : N2) for(size_t j = seg_begin[w];j < seg_end[w];j ++) { const ui x = radj[j]; adjB[bseg_end[x]++] = w; }

	bool have22 = false;
	auto bound = [&](ui v) -> int {
		if(grp[v] == 1) return 1+cn[v]+std::min(cap1, cn2[v]);
		return cn[v]+std::min(cap2, have22? cn2[v] : cap2);
	};
	auto drop = [&](ui v) { if((tag[v] == 1||tag[v] == 2)&&bound(v) < min_deg) { tag[v] = 3; queue_buf.push_back(v); } };
	auto peel_root = [&]() {
		queue_buf.clear();
		for(ui v : N1) drop(v);
		for(ui w : N2) drop(w);
		for(size_t i = 0;i < queue_buf.size();i ++) {
			const ui x = queue_buf[i];
			if(grp[x] == 1) {
				for(size_t j = seg_begin[x];j < seg_end[x];j ++) { const ui w = adjA[j]; if(tag[w] == 1) { -- cn[w]; drop(w); } }
				for(size_t j = bseg_begin[x];j < bseg_end[x];j ++) { const ui w = adjB[j]; if(tag[w] == 2) { -- cn[w]; drop(w); } }
			}
			else {
				for(size_t j = seg_begin[x];j < seg_end[x];j ++) { const ui v = radj[j]; if(tag[v] == 1) { -- cn2[v]; drop(v); } }
				if(have22) for(ui t = a22start[pos2[x]];t < a22start[pos2[x]+1];t ++) { const ui y = S2[a22[t]]; if(tag[y] == 2) { -- cn2[y]; drop(y); } }
			}
		}
	};
	peel_root();
	if(!root_ok()) { clear_tags(); return ; }

	S2.clear();
	for(ui w : N2) if(tag[w] == 2) { pos2[w] = S2.size(); S2.push_back(w); }
	if(cap2 > 0&&S2.size() >= 2) {
		e22.clear();
		for(size_t i = 0;i+1 < S2.size();i ++) {
			const ui w = S2[i];
			if(S2.size()-i-1 < size_t(st[w+1]-st[w])) {
				for(size_t t = i+1;t < S2.size();t ++) if(hash[w].find(int32_t(S2[t]))) { e22.push_back(i); e22.push_back(t); }
			}
			else {
				for(ept j = st[w];j < st[w+1]&&ed[j] > r;j ++) { const ui x = ed[j]; if(tag[x] == 2&&pos2[x] > i) { e22.push_back(i); e22.push_back(pos2[x]); } }
			}
		}
		a22start.assign(S2.size()+1, 0);
		for(ui x : e22) ++ a22start[x+1];
		for(size_t i = 0;i < S2.size();i ++) a22start[i+1] += a22start[i];
		a22.resize(e22.size());
		cursor.assign(a22start.begin(), a22start.end()-1);
		for(size_t i = 0;i < e22.size();i += 2) { a22[cursor[e22[i]]++] = e22[i+1]; a22[cursor[e22[i+1]]++] = e22[i]; }
		for(size_t i = 0;i < S2.size();i ++) cn2[S2[i]] = a22start[i+1]-a22start[i];
		have22 = true;
		peel_root();
		if(!root_ok()) { clear_tags(); return ; }
	}

	L.clear();
	long long e12 = 0;
	for(ui w : S2) if(tag[w] == 2) { L.push_back(w); e12 += cn[w]; }
	std::sort(L.begin(), L.end());
	const size_t s2n = L.size();
	const bool onehop_first = 2*e12 <= (long long)(q+2)*(long long)s2n;
	if(onehop_first) {
		T2.assign(L.begin(), L.end());
		L.clear();
		for(size_t i = N1.size();i > 0;i --) if(tag[N1[i-1]] == 1) L.push_back(N1[i-1]);
		L.insert(L.end(), T2.begin(), T2.end());
	}
	else for(size_t i = N1.size();i > 0;i --) if(tag[N1[i-1]] == 1) L.push_back(N1[i-1]);
	const size_t s = L.size();
	if(s+1 < size_t(q)) { clear_tags(); return ; }

	if(gid.size() < s) { gid.resize(s); capm.resize(s); capg.resize(s); }
	for(size_t i = 0;i < s;i ++) { lid[L[i]] = i; gid[i] = L[i]; }
	lstart.resize(s+1); ladj.clear(); lstart[0] = 0;
	for(size_t i = 0;i < s;i ++) {
		const ui v = L[i];
		if(grp[v] == 1) {
			for(size_t j = seg_begin[v];j < seg_end[v];j ++) if(lid[adjA[j]] >= 0) ladj.push_back(lid[adjA[j]]);
			for(size_t j = bseg_begin[v];j < bseg_end[v];j ++) if(lid[adjB[j]] >= 0) ladj.push_back(lid[adjB[j]]);
		}
		else {
			for(size_t j = seg_begin[v];j < seg_end[v];j ++) if(lid[radj[j]] >= 0) ladj.push_back(lid[radj[j]]);
			if(have22) for(ui t = a22start[pos2[v]];t < a22start[pos2[v]+1];t ++) { const ui y = S2[a22[t]]; if(lid[y] >= 0) ladj.push_back(lid[y]); }
		}
		lstart[i+1] = ladj.size();
	}
	for(size_t i = 0;i < s;i ++) lid[L[i]] = -1;
	run_root(s);
	clear_tags();
}


__attribute__((noinline)) inline void KPlexBase::comb_list(Count list_size, Count start, Count picked, Count r) {
	const int m = int(r-picked);
	if(start+Count(m) > list_size) return ;
	const Count n = list_size;
	if(m <= 3) {
		if(m == 0) ++ answer;
		else if(m == 1) { for(Count i = start;i < n;i ++) ++ answer; }
		else if(m == 2) { for(Count i = start;i < n-1;i ++) for(Count j = i+1;j < n;j ++) ++ answer; }
		else for(Count i = start;i < n-2;i ++) for(Count j = i+1;j < n-1;j ++) for(Count l = j+1;l < n;l ++) ++ answer;
		return ;
	}
	const int L = m-3;
	Count idx[64];
	int t = 0;
	idx[0] = start;
	for(;;) {
		if(idx[t]+Count(m-t) > n) {
			if(t == 0) break;
			++ idx[-- t];
			continue;
		}
		if(t+1 < L) { idx[t+1] = idx[t]+1; ++ t; continue; }
		for(Count i = idx[t]+1;i < n-2;i ++) for(Count j = i+1;j < n-1;j ++) for(Count l = j+1;l < n;l ++) ++ answer;
		++ idx[t];
	}
}
__attribute__((always_inline)) inline void KPlexBase::comb(Count list_size, Count r) {
	if(r > list_size) return ;
	if(list_size-r < r) r = list_size-r;
	if(r == 0) { ++ answer; return ; }
	Count cnt = 0;
	if(r == 1) { for(Count i = 0;i < list_size;i ++) ++ cnt; answer += cnt; return ; }
	if(r == 2) { for(Count i = 0;i+1 < list_size;i ++) for(Count j = i+1;j < list_size;j ++) ++ cnt; answer += cnt; return ; }
	comb_list(list_size, 0, 0, r);
}

__attribute__((noinline)) inline void KPlexBase::pool_types(int ncap, u64 one) {
	int nt = 0;
	for(int i = 0;i < ncap;i ++) {
		const u64 m = capm[i];
		int t = 0;
		while(t < nt&&pt_mask[t] != m) ++ t;
		if(t == nt) {
			if(nt == PT_MAX) {
				for(int x = 0;x < ncap;x ++) capg[x] = x+1;
				pool_rec_g(0, 0, 0, one);
				return ;
			}
			pt_mask[nt] = m; pt_size[nt] = 0; ++ nt;
		}
		++ pt_size[t];
	}
	if(nt == 1||ncap >= 3*nt) { pt_n = nt; pt_counts(0, 0); return ; }
	int off = 0;
	for(int t = 0;t < nt;t ++) {
		const int e = off+pt_size[t];
		for(int i = off;i < e;i ++) { capm[i] = pt_mask[t]; capg[i] = e; }
		off = e;
	}
	pool_rec_g(0, 0, 0, one);
}
inline void KPlexBase::pt_counts(int g, int used) {
	if(g == pt_n) {
		const int rem = pool_r-used;
		if(Count(rem) <= pool_free) type_prod(rem);
		return ;
	}
	const u64 m = pt_mask[g];
	int lim = std::min(pt_size[g], pool_r-used);
	for(u64 b = m;b;b &= b-1) lim = std::min(lim, resid[__builtin_ctzll(b)]);
	for(int c = 0;c <= lim;c ++) {
		pt_cnt[g] = c;
		for(u64 b = m;b;b &= b-1) resid[__builtin_ctzll(b)] -= c;
		pt_counts(g+1, used+c);
		for(u64 b = m;b;b &= b-1) resid[__builtin_ctzll(b)] += c;
	}
	pt_cnt[g] = 0;
}
inline void KPlexBase::type_prod(int rem) {
	pf_cnt = 0;
	for(int t = 0;t < pt_n;t ++) if(pt_cnt[t]) { pf_n[pf_cnt] = Count(pt_size[t]); pf_r[pf_cnt ++] = pt_cnt[t]; }
	if(rem) { pf_n[pf_cnt] = pool_free; pf_r[pf_cnt ++] = rem; }
	comb_prod();
}
inline void KPlexBase::comb_prod() {
	if(pf_cnt == 0) { ++ answer; return ; }
	int bi = 0;
	double best = -1;
	for(int i = 0;i < pf_cnt;i ++) {
		if(Count(pf_r[i]) > pf_n[i]) return ;
		double c = 1;
		for(int x = 0;x < pf_r[i];x ++) c = c*double(pf_n[i]-Count(x))/double(x+1);
		if(c > best) { best = c; bi = i; }
	}
	std::swap(pf_n[bi], pf_n[pf_cnt-1]); std::swap(pf_r[bi], pf_r[pf_cnt-1]);
	prod_enum(0);
}
inline void KPlexBase::prod_enum(int i) {
	if(i == pf_cnt-1) { comb(pf_n[i], Count(pf_r[i])); return ; }
	const int r = int(std::min(Count(pf_r[i]), pf_n[i]-Count(pf_r[i])));
	if(r == 0) prod_enum(i+1);
	else if(r == 1&&i == pf_cnt-2) { const Count n = pf_n[i+1]; const int rr = pf_r[i+1]; for(Count x = 0;x < pf_n[i];x ++) comb(n, Count(rr)); }
	else prod_comb(i, 0, r);
}
inline void KPlexBase::prod_comb(int i, Count start, int left) {
	for(Count x = start;x+Count(left) <= pf_n[i];x ++) {
		if(left == 1) prod_enum(i+1);
		else prod_comb(i, x+1, left-1);
	}
}
inline void KPlexBase::pool_rec_g(int start, int t, u64 ex, u64 one) {
	const int left = pool_r-t;
	comb(pool_free, Count(left));
	if(left == 0) return ;
	if(left == 1) {
		Count cnt = 0;
		for(int i = start;i < pool_ncap;) { const int e = capg[i]; if(capm[i]&ex) { i = e; continue; } for(;i < e;i ++) ++ cnt; }
		answer += cnt;
		return ;
	}
	if(left == 2) {
		Count cnt = 0;
		for(int i = start;i < pool_ncap;) {
			const int e = capg[i]; const u64 m = capm[i];
			if(m&ex) { i = e; continue; }
			const u64 ex2 = ex|(m&one);
			for(;i < e;i ++) {
				comb(pool_free, 1);
				for(int j = i+1;j < pool_ncap;) { const int f = capg[j]; if(capm[j]&ex2) { j = f; continue; } for(;j < f;j ++) ++ cnt; }
			}
		}
		answer += cnt;
		return ;
	}
	for(int i = start;i < pool_ncap&&Count(pool_ncap-i)+pool_free >= Count(left);) {
		const int e = capg[i]; const u64 m = capm[i];
		if(m&ex) { i = e; continue; }
		u64 nex = ex, none = one;
		for(u64 b = m;b;b &= b-1) {
			const int j = __builtin_ctzll(b);
			if(-- resid[j] == 0) { nex |= u64(1) << j; none &= ~(u64(1) << j); }
			else if(resid[j] == 1) none |= u64(1) << j;
		}
		for(;i < e&&Count(pool_ncap-i)+pool_free >= Count(left);i ++) pool_rec_g(i+1, t+1, nex, none);
		for(u64 b = m;b;b &= b-1) ++ resid[__builtin_ctzll(b)];
		i = e;
	}
}

