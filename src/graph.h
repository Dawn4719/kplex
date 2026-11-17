#pragma once


#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

typedef unsigned int ui;
typedef unsigned int ept;

struct Graph {
	ui n = 0;
	std::vector<ept> pstart;
	std::vector<ui> edges;

	static FILE *open_file(const std::string &name, const char *mode) {
		FILE *f = fopen(name.c_str(), mode);
		if(f == nullptr) { fprintf(stderr, "Can not open file: %s\n", name.c_str()); exit(1); }
		return f;
	}

	void read_text(const std::string &dir) {
		FILE *f = open_file(dir+"/edges.txt", "r");
		ui m = 0;
		if(fscanf(f, "%u%u", &n, &m) != 2) { fprintf(stderr, "bad edges.txt\n"); exit(1); }
		std::vector<std::pair<ui, ui> > e;
		e.reserve(size_t(m)*2);
		for(ui i = 0;i < m;i ++) {
			ui a, b;
			if(fscanf(f, "%u%u", &a, &b) != 2||a >= n||b >= n) { fprintf(stderr, "bad edge in edges.txt\n"); exit(1); }
			if(a == b) continue;
			e.push_back(std::make_pair(a, b)); e.push_back(std::make_pair(b, a));
		}
		fclose(f);
		std::sort(e.begin(), e.end());
		e.erase(std::unique(e.begin(), e.end()), e.end());
		pstart.assign(n+1, 0); edges.resize(e.size());
		for(size_t i = 0;i < e.size();i ++) { ++ pstart[e[i].first+1]; edges[i] = e[i].second; }
		for(ui v = 0;v < n;v ++) pstart[v+1] += pstart[v];
	}

	void read_binary(const std::string &dir) {
		FILE *f = open_file(dir+"/b_degree.bin", "rb");
		ui int_size = 0, m = 0;
		if(fread(&int_size, sizeof(ui), 1, f) != 1||int_size != sizeof(ui)) { fprintf(stderr, "bad b_degree.bin\n"); exit(1); }
		if(fread(&n, sizeof(ui), 1, f) != 1||fread(&m, sizeof(ui), 1, f) != 1) { fprintf(stderr, "bad b_degree.bin\n"); exit(1); }
		std::vector<ui> degree(n);
		if(fread(degree.data(), sizeof(ui), n, f) != n) { fprintf(stderr, "bad b_degree.bin\n"); exit(1); }
		fclose(f);
		f = open_file(dir+"/b_adj.bin", "rb");
		pstart.assign(n+1, 0); edges.resize(m);
		for(ui v = 0;v < n;v ++) {
			ui *buf = edges.data()+pstart[v], d = 0;
			if(degree[v] > 0&&fread(buf, sizeof(ui), degree[v], f) != degree[v]) { fprintf(stderr, "bad b_adj.bin\n"); exit(1); }
			std::sort(buf, buf+degree[v]);
			for(ui j = 0;j < degree[v];j ++) if(buf[j] != v&&(j == 0||buf[j] != buf[j-1])) buf[d++] = buf[j];
			pstart[v+1] = pstart[v]+d;
		}
		fclose(f);
		edges.resize(pstart[n]);
	}
};

