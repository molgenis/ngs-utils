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
	int chunk_id;
	int chunk_idx;

	snp(int chr, string & id, int pos, string & a0, string & a1, int idx) {
		this->chr = chr;
		this->id = id;
		this->pos = pos;
		this->a0 = a0;
		this->a1 = a1;
		this->idx = idx;
		chunk_id = -1;
		chunk_idx = -1;
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

class chunk {
public:
	string id;
	vector < int > mappingS;
	vector < vector < bool > > H;
	vector < bool > switches;

	chunk(char *) ;
	~chunk();

	void assemble(chunk *);
};

class data {
public:
	string id;
	double threshold;
	snp_set M;
	vector < chunk * > C;
	vector < string > I;
	vector < bool > males;

	data(string, double threshold);
	~data();

	void readLikelihoodsMap();
	void readHaploidGuys(char *);
	void readPosteriorsMapAndHap(char *);
	void assembleMap();
	void writeHaplotypes(char *);
	void writeGenotypes(char *);
};

#endif
