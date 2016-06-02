#include "classes.h"


int hamming (vector < bool > & h1, vector < bool > & h2) {
	assert(h1.size() == h2.size());
	int c = 0;
	for(int l = 0 ; l < h1.size(); l++) if (h1[l] != h2[l]) c++;
	return c;
}

chunk::chunk(char * id) {
	this->id = string(id);
}

chunk::~chunk() {
	mappingS.clear();
	H.clear();
	P.clear();
	switches.clear();
}

bool chunk::operator < (const chunk & c) const {
	return (mappingS[0] < c.mappingS[0]);
}

void chunk::assemble(chunk & pc) {
	//mapping starting point of the overlap
	int b1 = -1;
	int e1 = pc.mappingS.size() - 1;
	for (int l = 0 ; l < pc.mappingS.size() && b1 < 0 ; l ++ ) if (pc.mappingS[l] == mappingS[0]) b1 = l;
	if (b1 < 0) {cout << "Chunk [" << id << "] does not overlap chunk [" << pc.id << "], code 1" << endl; exit(1);}

	//mapping starting point of the overlap
	int b2 = 0;
	int e2 = -1;
	for (int l = 0 ; l < mappingS.size() && e2 < 0 ; l ++ ) if (mappingS[l] == pc.mappingS.back()) e2 = l;
	if (e2 < 0) {cout << "Chunk [" << id << "] does not overlap chunk [" << pc.id << "], code 2" << endl; exit(1);}

	//overlap size
	int osize = e1 - b1 + 1;
	int osize2 = e2 - b2 + 1;
	if (osize != osize2) {cout << "Overlaps are differents between [" << id << "] does not overlap chunk [" << pc.id << "] o1=" << osize << " o2=" << osize2 <<endl; exit(1);}
	else cout << "  * overlap=" << osize << endl;

	//blanking posteriors
	int cpt = 0;
	int i1, i2;
	for (i1 = b1, i2 = b2; i1 <= e1 && i2 <= e2; i1++, i2 ++) {
		if (cpt < (osize / 2)) mappingS[i2] = -1;
		else pc.mappingS[i1] = -1;
	}

	//switching haplotypes
	int n_switches = 0;
	for (int i = 0 ; i < switches.size() ; i ++) {
		vector < bool > h11 = vector < bool > (pc.H[2 * i + 0].begin() + b1, pc.H[2 * i + 0].end());
		vector < bool > h12 = vector < bool > (pc.H[2 * i + 1].begin() + b1, pc.H[2 * i + 1].end());
		vector < bool > h21 = vector < bool > (H[2 * i + 0].begin(), pc.H[2 * i + 0].begin() + e2 + 1);
		vector < bool > h22 = vector < bool > (H[2 * i + 1].begin(), pc.H[2 * i + 1].begin() + e2 + 1);
		int score0 = 0, score1 = 0;
		if (!pc.switches[i]) {
			score0 = hamming(h11, h21) + hamming(h12, h22);
			score1 = hamming(h11, h22) + hamming(h12, h21);
		} else {
			score0 = hamming(h12, h21) + hamming(h11, h22);
			score1 = hamming(h12, h22) + hamming(h11, h21);
		}

		if (score0 <= score1) switches[i] = false;
		else {
			switches[i] = true;
			n_switches ++;
		}
	}
	cout << "  * #switches=" <<n_switches << " out of " << switches.size() << " possible" << endl;
}
