#ifndef _CLASSES_H
#define _CLASSES_H

#include "utils.h"

class snp {
public:
	int chr;
	string id;
	int pos;
	string a0, a1;
	int idx;
	bool done;

	snp(int chr, string & id, int pos, string & a0, string & a1, int idx) {
		this->chr = chr;
		this->id = id;
		this->pos = pos;
		this->a0 = a0;
		this->a1 = a1;
		this->idx = idx;
		done = false;
	}
};

class snp_set {
public:
	vector < snp * > vec_snp;
	multimap < int, snp *> map_snp;

	//C/D
	snp_set() {};
	~snp_set() {};

	//M
	int size() {
		return vec_snp.size();
	}

	void push(snp * s) {
		vec_snp.push_back(s);
		map_snp.insert ( std::pair< int, snp * >( s->pos , s));
	}

	snp * get(int p, string & r0, string & r1) {
		pair < multimap < int, snp *>::iterator, multimap < int, snp *>::iterator > seqM = map_snp.equal_range(p);
		for (multimap < int, snp *>::iterator itM = seqM.first; itM != seqM.second; ++itM)
			if (r0 == itM->second->a0 && r1 == itM->second->a1) return itM->second;
		return NULL;
	}
};


class genotype {
public:
	string id;
	int idx;

	vector < float > L;
	vector < float > P;
	vector < char > G;
	vector < bool > h1;
	vector < bool > h2;

	//C/D
	genotype(string & id, int idx){
		this->id = id;
		this->idx = idx;
	}

	~genotype () {
		L.clear(); P.clear(); G.clear(); h1.clear(); h2.clear();
	}

	int call(float threshold) {
		int n = 0;
		for (int l = 0 ; l < h1.size() ; l ++) {
			G[l] = -1;
			if (P[3 * l + 0] >= threshold) G[l] = 0;
			if (P[3 * l + 1] >= threshold) G[l] = 1;
			if (P[3 * l + 2] >= threshold) G[l] = 2;
			if (G[l] >= 0) n ++;
		}
		return n;
	}
};

class chunk {
public:
	string id;
	vector < int > mappingS;
	vector < vector < bool > > H;
	vector < vector < float > > P;
	vector < bool > switches;

	chunk(char *) ;
	~chunk();

	bool operator < (const chunk & c) const;

	void assemble(chunk & pc);
};

class data {
public:
	vector < genotype * > G;
	snp_set M;
	vector < chunk * > C;

	data();
	~data();

	void readLikelihoods(char *);
	void readPosteriors(char *);
	void assemble();
	void call(double);
	void writeHaplotypes(char *);
	void writeGenotypes(char *);
};

#endif
