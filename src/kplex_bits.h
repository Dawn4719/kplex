#pragma once


#include <immintrin.h>
#include "kplex_base.h"

class KPlexBits : public KPlexBase {
public:
	KPlexBits(const Graph &g, int k_, int q_): KPlexBase(g, k_, q_) {}

private:
	size_t W = 0;
	std::vector<u64> rows, root_row, CB, BC;
	std::vector<int> csize;
	std::vector<const u64 *> PR;
	Count lz = 0;
	std::vector<u64> PLB;

	const u64 *row_of(int v) const { return rows.data()+size_t(v)*W; }
	static bool bit(const u64 *b, int v) { return (b[v>>6] >> (v&63))&1; }
	u64 *cls(int need, int t) { return BC.data()+(size_t(need)*(k+1)+t)*W; }
	u64 *poolb_at(int need) { return PLB.data()+size_t(need)*W; }

	size_t h0 = 0; int n1v = 0; u64 h1mask = 0; bool rview = false;
	std::vector<u64> H1;
	u64 hop1(const u64 *b) const {
		const size_t w = h0 >> 6; const unsigned sh = unsigned(h0&63);
		u64 x = b[w] >> sh;
		if(sh&&w+1 < W) x |= b[w+1] << (64-sh);
		return x&h1mask;
	}

	void run_root(size_t s) override;
	int extract_safe(int need);
	void pool_list(int r, int dd, int slot, Count zc, u64 cand = ~u64(0));
	void pool_single(int slot, const u64 *row, u64 conflict);
	int pool_filter(int src, int n, const u64 *row, u64 rows_mask, int dst);
	int pool_adj(int slot, const u64 *row);
	void listing(int need, int d, u64 sat, u64 s1);
	struct HopView { bool ready; u64 C, F, K[32], P[64]; };
	void view_prep(int need, int d, HopView &hv);
	void view_child(int need, int d, int u, u64 mu, u64 sat, u64 s1, HopView &hv);
	int add_prefix(u64 mu, int d, u64 &ns1) {
		ns1 &= ~mu;
		for(u64 b = mu;b;b &= b-1) { const int j = __builtin_ctzll(b); if(++ M[j] == k-1) ns1 |= u64(1) << j; }
		return M[d] = __builtin_popcountll(mu);
	}
	void list_pairs(int d, u64 s1);

	u64 R1buf[64];
	const u64 *R1 = R1buf;
	u64 P1[64];
	u64 C1[64], F1[64];
	std::vector<u64> K1;
	u64 &k1(int need, int t) { return K1[size_t(need)*(k+1)+t]; }
	void enter_compact(int nl, int dd, u64 nsat, u64 ns1);
	void tail_compact(int need, int d, size_t wu, int u);
	void listing1(int need, int d, u64 sat, u64 s1);
	void branch1(int need, int d, u64 sat, u64 s1, u64 C, u64 F);
	void list_pairs1(int d, u64 s1);
	void pool_list1(int r, int dd, u64 F, u64 cand = ~u64(0));
	struct Rc { u64 p0, p1, p2, p3, p4; };
	u64 NB1[64], EX1[64], MM1[64];
	int r1 = 0;
	Count fz1 = 0;
	template<bool MULTI> void pool_rec1(u64 cand, int t, Rc rc);
};

inline void KPlexBits::run_root(size_t s) {
	const int top = q-1;
	if(K1.empty()) K1.assign(size_t(q)*(k+1), 0);
	if(s <= 64) {
		u64 one = 0;
		for(size_t i = 0;i < s;i ++) {
			u64 r = 0;
			for(ui l = lstart[i];l < lstart[i+1];l ++) r |= u64(1) << ladj[l];
			R1buf[i] = r;
			if(is_onehop(i)) one |= u64(1) << i;
		}
		const u64 all = s == 64? ~u64(0) : (u64(1) << s)-1;
		for(int t = 0;t <= k;t ++) k1(top, t) = 0;
		k1(top, 0) = one;
		if(k >= 1) k1(top, 1) = all&~one;
		C1[top] = all; F1[top] = 0; P1[0] = one; M[0] = 0; R1 = R1buf;
		listing1(top, 1, k == 0? 1 : 0, k == 1? 1 : 0);
		return ;
	}
	W = (s+63)/64;
	rows.assign(s*W, 0); root_row.assign(W, 0);
	if(CB.size() < size_t(q)*W) CB.resize(size_t(q)*W);
	if(BC.size() < size_t(q)*(k+1)*W) BC.resize(size_t(q)*(k+1)*W);
	for(int t = 0;t <= k;t ++) std::fill(cls(top, t), cls(top, t)+W, 0);
	csize.assign(q, 0); PR.assign(q, nullptr);
	if(PLB.size() < size_t(q)*W) PLB.resize(size_t(q)*W);
	std::fill(poolb_at(top), poolb_at(top)+W, 0);
	u64 *C = &CB[size_t(top)*W];
	std::fill(C, C+W, 0);
	n1v = 0; h0 = s;
	for(size_t i = 0;i < s;i ++) {
		u64 *r = &rows[i*W];
		for(ui l = lstart[i];l < lstart[i+1];l ++) r[ladj[l]>>6] |= u64(1) << (ladj[l]&63);
		C[i>>6] |= u64(1) << (i&63);
		if(is_onehop(i)) { root_row[i>>6] |= u64(1) << (i&63); if(n1v ++ == 0) h0 = i; }
		cls(top, is_onehop(i)? 0 : 1)[i>>6] |= u64(1) << (i&63);
	}
	csize[top] = s;
	rview = n1v > 0&&n1v <= 64;
	if(rview) {
		h1mask = n1v == 64? ~u64(0) : (u64(1) << n1v)-1;
		H1.resize(n1v);
		for(int i = 0;i < n1v;i ++) H1[i] = hop1(row_of(int(h0)+i));
	}
	M[0] = 0; lz = 0; PR[0] = root_row.data();
	listing(top, 1, k == 0? 1 : 0, k == 1? 1 : 0);
}


inline int KPlexBits::extract_safe(int need) {
	u64 *C = &CB[size_t(need)*W], *PB = poolb_at(need);
	const u64 *Bk = cls(need, k);
	queue_buf.clear();
	for(int t = 0;t <= k;t ++) {
		const u64 *Bt = cls(need, t);
		const int limit = k-t;
		for(size_t w = 0;w < W;w ++) for(u64 b = C[w]&Bt[w];b;b &= b-1) {
			const int v = int(w*64+__builtin_ctzll(b));
			if(need-1 <= limit) { queue_buf.push_back(v); continue; }
			const u64 *r = row_of(v);
			int nd = -1;
			for(size_t x = 0;x < W&&nd <= limit;x ++) nd += __builtin_popcountll(((t < k? C[x]&~Bk[x] : C[x])|PB[x])&~r[x]);
			if(nd <= limit) queue_buf.push_back(v);
		}
	}
	const int extracted = queue_buf.size();
	for(int i = 0;i < extracted;i ++) {
		const int v = queue_buf[i];
		C[v>>6] &= ~(u64(1) << (v&63));
		PB[v>>6] |= u64(1) << (v&63);
	}
	csize[need] -= extracted;
	return extracted;
}

inline void KPlexBits::pool_list(int r, int dd, int slot, Count zc, u64 cand) {
	if(Count(r) > zc) return ;
	const u64 *ZB = poolb_at(slot);
	u64 tight = 0, one = 0;
	if(dd < 64) cand &= (u64(1) << dd)-1;
	if(r >= 2) for(;cand;cand &= cand-1) {
		const int j = __builtin_ctzll(cand);
		if(M[j] >= k||k-M[j] >= r) continue;
		const u64 *pr = PR[j];
		int c = 0;
		for(size_t w = 0;w < W;w ++) c += __builtin_popcountll(ZB[w]&~pr[w]);
		if(c > k-M[j]) { tight |= u64(1) << j; resid[j] = k-M[j]; if(resid[j] == 1) one |= u64(1) << j; }
	}
	if(!tight) { comb(zc, r); return ; }
	int ncap = 0;
	auto add = [&](int v) {
		u64 m = 0;
		for(u64 t = tight;t;t &= t-1) { const int j = __builtin_ctzll(t); if(!bit(PR[j], v)) m |= u64(1) << j; }
		if(m) capm[ncap++] = m;
	};
	for(size_t w = 0;w < W;w ++) for(u64 b = ZB[w];b;b &= b-1) add(int(w*64+__builtin_ctzll(b)));
	pool_ncap = ncap; pool_r = r; pool_free = zc-Count(ncap);
	pool_types(ncap, one);
}

inline void KPlexBits::pool_single(int slot, const u64 *row, u64 conflict) {
	const u64 *Z = poolb_at(slot);
	for(size_t w = 0;w < W;w ++) {
		u64 x = Z[w];
		if(row) x &= row[w];
		for(u64 t = conflict;t&&x;t &= t-1) x &= PR[__builtin_ctzll(t)][w];
		for(;x;x &= x-1) ++ answer;
	}
}

__attribute__((noinline)) inline int KPlexBits::pool_adj(int slot, const u64 *row) {
	const u64 *Z = poolb_at(slot);
	int c = 0;
	for(size_t w = 0;w < W;w ++) c += __builtin_popcountll(Z[w]&row[w]);
	return c;
}

inline int KPlexBits::pool_filter(int src, int n, const u64 *row, u64 rows_mask, int dst) {
	const u64 *Z = poolb_at(src);
	u64 *O = poolb_at(dst);
	int c = 0;
	for(size_t w = 0;w < W;w ++) {
		u64 x = n? Z[w] : 0;
		if(row) x &= row[w];
		for(u64 t = rows_mask;t&&x;t &= t-1) x &= PR[__builtin_ctzll(t)][w];
		O[w] = x; c += __builtin_popcountll(x);
	}
	return c;
}


inline void KPlexBits::listing(int need, int d, u64 sat, u64 s1) {
	u64 *C = &CB[size_t(need)*W];
	if(need == 1) {
		comb(lz, 1);
		for(size_t w = 0;w < W;w ++) for(u64 b = C[w];b;b &= b-1) ++ answer;
		return ;
	}
	if(need == 2) { list_pairs(d, s1); return ; }
	int extracted = 0;
	if(csize[need] > 0) extracted = extract_safe(need);
	lz += extracted;
	pool_list(need, d, need, lz);
	const int cs = csize[need];
	const Count lz_state = lz;
	u64 *NX = &CB[size_t(need-1)*W];
	const u64 *B0 = cls(need, 0), *Bk = cls(need, k);
	const u64 ubit = u64(1) << d;
	HopView hv; hv.ready = false;
	bool pool2 = true;
	if(rview&&h0 > 0) { pool2 = false; const u64 *Z = poolb_at(need); for(size_t w = 0;w < W&&!pool2;w ++) pool2 = (Z[w]&~root_row[w]) != 0; }
	int i = 0, f1 = -1;
	for(size_t wu = 0;wu < W;wu ++) for(u64 bu = C[wu];bu&&Count(i)+need <= Count(cs)+lz;bu &= bu-1, ++ i) {
		const int u = int(wu*64+__builtin_ctzll(bu));
		if(h0 == 0&&size_t(u) >= size_t(n1v)) {
			if(f1 < 0) { f1 = 0; const u64 *Z = poolb_at(need); for(size_t w = 0;w < W;w ++) f1 += __builtin_popcountll(Z[w]&root_row[w]); }
			if(f1 < need-(k-M[0])) { wu = W; break; }
		}
		if(Count(cs-i)+lz <= 64) {
			tail_compact(need, d, wu, u);
			branch1(need, d, sat, s1, C1[need], F1[need]);
			wu = W; break;
		}
		u64 mu = 0;
		for(int j = 0;j < d;j ++) if(!bit(PR[j], u)) mu |= u64(1) << j;
		const bool usat = int(__builtin_popcountll(mu)) == k;
		const u64 newly = mu&s1;
		if(rview&&((newly&1)||(!pool2&&size_t(u) >= h0))) { view_child(need, d, u, mu, sat, s1, hv); continue; }
		const u64 *ru = row_of(u);
		lz = Count(pool_filter(need, int(lz_state), usat? ru : nullptr, newly, need-1));
		int cnt = 0, adj_all = 0;
		size_t hi = wu;
		for(size_t w = 0;w < wu;w ++) NX[w] = 0;
		for(size_t w = wu;w < W;w ++) {
			u64 x = C[w]&(usat? ru[w] : (ru[w]|~Bk[w]));
			if(w == wu) x &= (u&63) == 63? 0 : ~u64(0) << ((u&63)+1);
			for(u64 t = newly;t&&x;t &= t-1) x &= PR[__builtin_ctzll(t)][w];
			NX[w] = x;
			if(x) { hi = w; cnt += __builtin_popcountll(x); adj_all += __builtin_popcountll(x&B0[w]&ru[w]); }
		}
		csize[need-1] = cnt;
		if(Count(cnt)+lz >= Count(need-1)) {
			u64 nsat = sat|newly, ns1 = s1;
			const int mu_n = add_prefix(mu, d, ns1);
			if(usat) nsat |= ubit;
			if(mu_n == k-1) ns1 |= ubit;
			bool alive = true;
			const int rmax = need-1-adj_all;
			for(int j = 0;j <= d&&alive&&rmax > 1;j ++) {
				const int r = k-M[j];
				if(r == 0||r >= rmax) continue;
				const int req = need-1-r;
				const u64 *pr = j == d? ru : PR[j];
				int a = 0;
				for(size_t w = wu;w <= hi;w ++) a += __builtin_popcountll(NX[w]&pr[w]);
				if(a < req&&(Count(a)+lz < Count(req)||a+pool_adj(need-1, pr) < req)) { alive = false; }
			}
			if(alive&&Count(cnt)+lz <= 64) { PR[d] = ru; enter_compact(need-1, d+1, nsat, ns1); }
			else if(alive) {
				PR[d] = ru;
				for(int t = 0;t <= k;t ++) {
					const u64 *Bt = cls(need, t), *Bp = t? cls(need, t-1) : nullptr;
					u64 *Ct = cls(need-1, t);
					for(size_t w = wu;w <= hi;w ++) Ct[w] = NX[w]&((Bt[w]&ru[w])|(Bp? Bp[w]&~ru[w] : 0));
				}
				listing(need-1, d+1, nsat, ns1);
			}
			for(u64 b = mu;b;b &= b-1) -- M[__builtin_ctzll(b)];
		}
		lz = lz_state;
	}
	lz -= extracted;
}

inline void KPlexBits::view_prep(int need, int d, HopView &hv) {
	hv.C = hop1(&CB[size_t(need)*W]); hv.F = hop1(poolb_at(need));
	for(int t = 0;t <= k;t ++) hv.K[t] = hop1(cls(need, t));
	for(int j = 0;j < d;j ++) hv.P[j] = hop1(PR[j]);
	hv.ready = true;
}

__attribute__((noinline)) inline void KPlexBits::view_child(int need, int d, int u, u64 mu, u64 sat, u64 s1, HopView &hv) {
	if(!hv.ready) view_prep(need, d, hv);
	const bool usat = int(__builtin_popcountll(mu)) == k;
	const u64 newly = mu&s1, ru = hop1(row_of(u));
	u64 keep = usat? ru : ~u64(0);
	for(u64 t = newly&~u64(1);t;t &= t-1) keep &= hv.P[__builtin_ctzll(t)];
	const u64 later = size_t(u) < h0? h1mask : size_t(u) >= h0+size_t(n1v)? 0 : h1mask&~((u64(2) << (size_t(u)-h0))-1);
	const u64 X = hv.C&later&keep&(usat? ru : (ru|~hv.K[k]));
	const u64 Fn = hv.F&keep;
	const int cnt = __builtin_popcountll(X), lzn = __builtin_popcountll(Fn);
	if(cnt+lzn < need-1) return ;
	u64 nsat = sat|newly, ns1 = s1;
	const int mu_n = add_prefix(mu, d, ns1);
	if(usat) nsat |= u64(1) << d;
	if(mu_n == k-1) ns1 |= u64(1) << d;
	bool alive = true;
	const int rmax = need-1-__builtin_popcountll(X&hv.K[0]&ru);
	for(int j = 0;j <= d&&alive&&rmax > 1;j ++) {
		const int r = k-M[j];
		if(r == 0||r >= rmax) continue;
		const int req = need-1-r;
		const u64 pr = j == d? ru : hv.P[j];
		const int a = __builtin_popcountll(X&pr);
		if(a < req&&(a+lzn < req||a+__builtin_popcountll(Fn&pr) < req)) { alive = false; }
	}
	if(alive) {
		for(int j = 0;j < d;j ++) P1[j] = hv.P[j];
		P1[d] = ru;
		for(int t = 0;t <= k;t ++) k1(need-1, t) = X&((hv.K[t]&ru)|(t? hv.K[t-1]&~ru : 0));
		C1[need-1] = X; F1[need-1] = Fn; R1 = H1.data();
		listing1(need-1, d+1, nsat, ns1);
	}
	for(u64 b = mu;b;b &= b-1) -- M[__builtin_ctzll(b)];
}

inline void KPlexBits::list_pairs(int d, u64 s1) {
	const u64 *C = &CB[size_t(2)*W], *Bk = cls(2, k);
	pool_list(2, d, 2, lz, s1);
	const int np = int(lz);
	for(size_t wu = 0;wu < W;wu ++) for(u64 bu = C[wu];bu;bu &= bu-1) {
		const int u = int(wu*64+__builtin_ctzll(bu));
		const u64 *ru = row_of(u);
		const bool full_u = bit(Bk, u);
		u64 conflict = 0;
		for(u64 t = s1;t;t &= t-1) { const int j = __builtin_ctzll(t); if(!bit(PR[j], u)) conflict |= u64(1) << j; }
		if(np) pool_single(2, full_u? ru : nullptr, conflict);
		for(size_t w = wu;w < W;w ++) {
			u64 x = C[w]&(full_u? ru[w] : (ru[w]|~Bk[w]));
			if(w == wu) x &= (u&63) == 63? 0 : ~u64(0) << ((u&63)+1);
			for(u64 t = conflict;t&&x;t &= t-1) x &= PR[__builtin_ctzll(t)][w];
			for(;x;x &= x-1) ++ answer;
		}
	}
}


__attribute__((target("bmi2"))) static inline uint64_t compact_bits(const uint64_t *b, const size_t *wid, const uint64_t *mk, const int *off, size_t nw) {
	uint64_t x = 0;
	for(size_t i = 0;i < nw;i ++) x |= _pext_u64(b[wid[i]], mk[i]) << off[i];
	return x;
}

__attribute__((target("bmi2"), noinline)) inline void KPlexBits::enter_compact(int nl, int dd, u64 nsat, u64 ns1) {
	const u64 *NX = &CB[size_t(nl)*W], *ru = PR[dd-1];
	const u64 *PB = poolb_at(nl);
	u64 mk[64]; int off[64];
	int tot = 0;
	size_t nw = 0, wid[64];
	for(size_t w = 0;w < W;w ++) { const u64 m = NX[w]|PB[w]; if(m) { wid[nw] = w; mk[nw] = m; off[nw] = tot; tot += __builtin_popcountll(m); ++ nw; } }
	auto cp = [&](const u64 *b) -> u64 { return compact_bits(b, wid, mk, off, nw); };
	int c = 0;
	for(size_t i = 0;i < nw;i ++) for(u64 b = mk[i];b;b &= b-1, ++ c) {
		const int v = int(wid[i]*64+__builtin_ctzll(b));
		if((NX[v>>6] >> (v&63))&1) R1buf[c] = cp(row_of(v));
	}
	for(int j = 0;j < dd;j ++) P1[j] = cp(PR[j]);
	C1[nl] = cp(NX); F1[nl] = cp(PB); R1 = R1buf;
	for(int t = 0;t <= k;t ++) {
		const u64 *Bt = cls(nl+1, t), *Bp = t? cls(nl+1, t-1) : nullptr;
		u64 x = 0;
		for(size_t i = 0;i < nw;i ++) { const size_t w = wid[i]; x |= _pext_u64(NX[w]&((Bt[w]&ru[w])|(Bp? Bp[w]&~ru[w] : 0)), mk[i]) << off[i]; }
		k1(nl, t) = x;
	}
	listing1(nl, dd, nsat, ns1);
}

__attribute__((target("bmi2"), noinline)) inline void KPlexBits::tail_compact(int need, int d, size_t wu, int u) {
	const u64 *C = &CB[size_t(need)*W], *PB = poolb_at(need);
	u64 mk[64]; int off[64];
	int tot = 0;
	size_t nw = 0, wid[64];
	for(size_t w = 0;w < W;w ++) {
		const u64 m = (w < wu? 0 : w == wu? C[w]&(~u64(0) << (u&63)) : C[w])|PB[w];
		if(m) { wid[nw] = w; mk[nw] = m; off[nw] = tot; tot += __builtin_popcountll(m); ++ nw; }
	}
	auto cp = [&](const u64 *b) -> u64 { return compact_bits(b, wid, mk, off, nw); };
	int c = 0;
	for(size_t i = 0;i < nw;i ++) for(u64 b = mk[i];b;b &= b-1, ++ c) {
		const int v = int(wid[i]*64+__builtin_ctzll(b));
		if((C[v>>6] >> (v&63))&1) R1buf[c] = cp(row_of(v));
	}
	for(int j = 0;j < d;j ++) P1[j] = cp(PR[j]);
	C1[need] = cp(C); F1[need] = cp(PB); R1 = R1buf;
	for(int t = 0;t <= k;t ++) k1(need, t) = cp(cls(need, t))&C1[need];
}

inline void KPlexBits::pool_list1(int r, int dd, u64 F, u64 cand) {
	const Count zc = __builtin_popcountll(F);
	if(Count(r) > zc) return ;
	u64 tight = 0, rs = 0;
	if(dd < 64) cand &= (u64(1) << dd)-1;
	if(r >= 2) for(;cand;cand &= cand-1) {
		const int j = __builtin_ctzll(cand);
		if(M[j] >= k||k-M[j] >= r) continue;
		const u64 nb = F&~P1[j];
		if(__builtin_popcountll(nb) > k-M[j]) { tight |= u64(1) << j; rs |= nb; NB1[j] = nb; resid[j] = k-M[j]; }
	}
	if(!tight) { comb(zc, r); return ; }
	const int ncap = __builtin_popcountll(rs), nj = __builtin_popcountll(tight);
	{
		int sres = 0;
		for(u64 t = tight;t;t &= t-1) sres += resid[__builtin_ctzll(t)];
		if(zc-Count(ncap)+Count(std::min(ncap, sres)) < Count(r)) return ;
	}
	if(nj <= 2) {
		const int j0 = __builtin_ctzll(tight), j1 = nj == 2? __builtin_ctzll(tight&(tight-1)) : j0;
		if(resid[j0] == 1&&resid[j1] == 1) {
			const Count fz = zc-Count(ncap);
			comb(fz, r);
			pf_n[0] = Count(ncap); pf_r[0] = 1; pf_n[1] = fz; pf_r[1] = r-1; pf_cnt = 2;
			comb_prod();
			if(nj == 2) {
				pf_n[0] = Count(__builtin_popcountll(NB1[j0]&~NB1[j1])); pf_r[0] = 1;
				pf_n[1] = Count(__builtin_popcountll(NB1[j1]&~NB1[j0])); pf_r[1] = 1;
				pf_n[2] = fz; pf_r[2] = r-2; pf_cnt = r > 2? 3 : 2;
				comb_prod();
			}
			return ;
		}
	}
	bool all1 = true;
	for(u64 t = tight;t;t &= t-1) if(resid[__builtin_ctzll(t)] != 1) { all1 = false; break; }
	if(!all1&&nj <= 4) {
		int nt = 0;
		if(ncap < 6) {
			bool one_type = true;
			for(u64 t = tight;t;t &= t-1) if(NB1[__builtin_ctzll(t)] != rs) { one_type = false; break; }
			if(one_type) { pt_mask[0] = tight; pt_size[0] = ncap; nt = 1; }
		}
		else {
			int jl[4], c = 0;
			for(u64 t = tight;t;t &= t-1) jl[c++] = __builtin_ctzll(t);
			for(int m = 1;m < (1 << nj)&&3*nt <= ncap;m ++) {
				u64 cell = rs, pm = 0;
				for(int b = 0;b < nj;b ++) { if((m >> b)&1) { cell &= NB1[jl[b]]; pm |= u64(1) << jl[b]; } else cell &= ~NB1[jl[b]]; }
				if(cell) { pt_mask[nt] = pm; pt_size[nt] = __builtin_popcountll(cell); ++ nt; }
			}
		}
		if(nt == 1||(nt > 1&&ncap >= 3*nt)) { pool_r = r; pool_free = zc-Count(ncap); pt_n = nt; pt_counts(0, 0); return ; }
	}
	for(u64 b = rs;b;b &= b-1) { const int v = __builtin_ctzll(b); EX1[v] = 0; MM1[v] = 0; }
	Rc rc = {0, 0, 0, 0, 0};
	bool multi = false;
	for(u64 t = tight;t;t &= t-1) {
		const int j = __builtin_ctzll(t);
		const u64 nb = NB1[j], jb = u64(1) << j;
		if(resid[j] == 1) { for(u64 b = nb;b;b &= b-1) EX1[__builtin_ctzll(b)] |= nb; }
		else {
			multi = true;
			for(u64 b = nb;b;b &= b-1) MM1[__builtin_ctzll(b)] |= jb;
			const unsigned x = unsigned(resid[j]);
			if(x&1) rc.p0 |= jb;
			if(x&2) rc.p1 |= jb;
			if(x&4) rc.p2 |= jb;
			if(x&8) rc.p3 |= jb;
			if(x&16) rc.p4 |= jb;
		}
	}
	r1 = r; fz1 = zc-Count(ncap);
	if(multi) pool_rec1<true>(rs, 0, rc); else pool_rec1<false>(rs, 0, rc);
}

template<bool MULTI> inline void KPlexBits::pool_rec1(u64 cand, int t, Rc rc) {
	const int left = r1-t;
	comb(fz1, Count(left));
	if(left == 1) { answer += __builtin_popcountll(cand); return ; }
	if(Count(__builtin_popcountll(cand))+fz1 < Count(left)) return ;
	const u64 one = MULTI? rc.p0&~(rc.p1|rc.p2|rc.p3|rc.p4) : 0;
	if(left == 2) {
		for(u64 b = cand;b;b &= b-1) {
			const int v = __builtin_ctzll(b);
			u64 nc = (b&(b-1))&~EX1[v];
			if(MULTI) for(u64 m = MM1[v]&one;m;m &= m-1) nc &= ~NB1[__builtin_ctzll(m)];
			answer += fz1;
			answer += __builtin_popcountll(nc);
		}
		return ;
	}
	for(u64 b = cand;b;b &= b-1) {
		const int v = __builtin_ctzll(b);
		u64 nc = (b&(b-1))&~EX1[v];
		const u64 S = MULTI? MM1[v] : 0;
		Rc nr = rc;
		if(MULTI&&S) {
			for(u64 m = S&one;m;m &= m-1) nc &= ~NB1[__builtin_ctzll(m)];
			u64 br = S;
			nr.p0 = rc.p0^br; br &= ~rc.p0;
			nr.p1 = rc.p1^br; br &= ~rc.p1;
			nr.p2 = rc.p2^br; br &= ~rc.p2;
			nr.p3 = rc.p3^br; br &= ~rc.p3;
			nr.p4 = rc.p4^br;
		}
		if(!nc) comb(fz1, Count(left-1));
		else if(Count(__builtin_popcountll(nc))+fz1 >= Count(left-1)) pool_rec1<MULTI>(nc, t+1, nr);
	}
}

inline void KPlexBits::listing1(int need, int d, u64 sat, u64 s1) {
	u64 C = C1[need], F = F1[need];
	if(need == 1) { answer += __builtin_popcountll(F); answer += __builtin_popcountll(C); return ; }
	if(need == 2) { list_pairs1(d, s1); return ; }
	if(C) {
		u64 ex = 0;
		const u64 U = C|F, full = C&k1(need, k);
		for(int t = 0;t <= k;t ++) {
			const int limit = k-t;
			const u64 Ut = t < k? U&~full : U;
			for(u64 b = C&k1(need, t);b;b &= b-1) {
				const int v = __builtin_ctzll(b);
				u64 N = Ut&~R1[v]&~(u64(1) << v);
				int nd = __builtin_popcountll(N);
				if(std::min(need-1, nd) <= limit) ex |= u64(1) << v;
			}
		}
		C &= ~ex; F |= ex;
	}
	pool_list1(need, d, F);
	branch1(need, d, sat, s1, C, F);
}

inline void KPlexBits::branch1(int need, int d, u64 sat, u64 s1, u64 C, u64 F) {
	const int cs = __builtin_popcountll(C), lzs = __builtin_popcountll(F);
	const u64 B0 = k1(need, 0), Bk = k1(need, k);
	const u64 ubit = u64(1) << d;
	int i = 0;
	for(u64 bu = C;bu&&i+need <= cs+lzs;bu &= bu-1, ++ i) {
		const int u = __builtin_ctzll(bu);
		u64 mu = 0;
		for(int j = 0;j < d;j ++) if(!((P1[j] >> u)&1)) mu |= u64(1) << j;
		const int mu_n = __builtin_popcountll(mu);
		const bool usat = mu_n == k;
		const u64 newly = mu&s1, ru = R1[u];
		u64 keep = usat? ru : ~u64(0);
		for(u64 t = newly;t;t &= t-1) keep &= P1[__builtin_ctzll(t)];
		const u64 Fn = F&keep;
		const u64 X = C&(u == 63? 0 : ~u64(0) << (u+1))&keep&(ru|~Bk);
		const int cnt = __builtin_popcountll(X), lzn = __builtin_popcountll(Fn);
		if(cnt+lzn < need-1) continue;
		const int adj_all = __builtin_popcountll(X&B0&ru);
		u64 nsat = sat|newly, ns1 = s1;
		add_prefix(mu, d, ns1);
		if(usat) nsat |= ubit;
		if(mu_n == k-1) ns1 |= ubit;
		bool alive = true;
		const int rmax = need-1-adj_all;
		for(int j = 0;j <= d&&alive&&rmax > 1;j ++) {
			const int r = k-M[j];
			if(r == 0||r >= rmax) continue;
			const int req = need-1-r;
			const u64 pr = j == d? ru : P1[j];
			const int a = __builtin_popcountll(X&pr);
			if(a < req&&(a+lzn < req||a+__builtin_popcountll(Fn&pr) < req)) { alive = false; }
		}
		if(alive) {
			P1[d] = ru;
			for(int t = 0;t <= k;t ++) k1(need-1, t) = X&((k1(need, t)&ru)|(t? k1(need, t-1)&~ru : 0));
			C1[need-1] = X; F1[need-1] = Fn;
			listing1(need-1, d+1, nsat, ns1);
		}
		for(u64 b = mu;b;b &= b-1) -- M[__builtin_ctzll(b)];
	}
}

inline void KPlexBits::list_pairs1(int d, u64 s1) {
	const u64 C = C1[2], F = F1[2], Bk = k1(2, k);
	pool_list1(2, d, F, s1);
	for(u64 bu = C;bu;bu &= bu-1) {
		const int u = __builtin_ctzll(bu);
		const u64 ru = R1[u];
		const bool full_u = (Bk >> u)&1;
		u64 keep = ~u64(0);
		for(u64 t = s1;t;t &= t-1) { const int j = __builtin_ctzll(t); if(!((P1[j] >> u)&1)) keep &= P1[j]; }
		if(F) answer += __builtin_popcountll(F&keep&(full_u? ru : ~u64(0)));
		answer += __builtin_popcountll(C&(u == 63? 0 : ~u64(0) << (u+1))&keep&(full_u? ru : (ru|~Bk)));
	}
}

