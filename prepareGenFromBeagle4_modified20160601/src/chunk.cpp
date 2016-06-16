#include "classes.h"


int hamming (vector < bool > & h1, vector < bool > & h2) {
	assert(h1.size() == h2.size());
	int c = 0;
	for(int l = 0 ; l < h1.size(); l++) if (h1[l] != h2[l]) c++;
	return c;
}

int discordance (vector < bool > & h11, vector < bool > & h12, vector < bool > & h21, vector < bool > & h22) {
	int c = 0;
	for(int l = 0 ; l < h11.size(); l++) {
		int g1 = h11[l] + h12[l];
		int g2 = h21[l] + h22[l];
		if (g1 != g2) c++;
	}
	return c;
}


chunk::chunk(char * id) {
	this->id = string(id);
}

chunk::~chunk() {
	mappingS.clear();
}

void chunk::assemble(chunk * pc) {
	//mapping starting point of the overlap
	int b1 = -1;
	int e1 = pc->mappingS.size() - 1;
	for (int l = 0 ; l < pc->mappingS.size() && b1 < 0 ; l ++ ) if (pc->mappingS[l] == mappingS[0]) b1 = l;
	if (b1 < 0) {cout << "Chunk [" << id << "] does not overlap chunk [" << pc->id << "], code 1" << endl; exit(1);}

	//mapping starting point of the overlap
	int b2 = 0;
	int e2 = -1;
	for (int l = 0 ; l < mappingS.size() && e2 < 0 ; l ++ ) if (mappingS[l] == pc->mappingS.back()) e2 = l;
	if (e2 < 0) {cout << "Chunk [" << id << "] does not overlap chunk [" << pc->id << "], code 2" << endl; exit(1);}

	//overlap size
	int osize = e1 - b1 + 1;
	int osize2 = e2 - b2 + 1;
	if (osize != osize2) {cout << "Overlaps are differents between [" << id << "] does not overlap chunk [" << pc->id << "] o1=" << osize << " o2=" << osize2 <<endl; exit(1);}
	else cout << "  * overlap=" << osize << endl;

	//blanking overlap
	cout << "  * prev chunk end = " << pc->mappingS[b1 + ((int)ceil(osize * 1.0 / 2)) - 1] << endl;
	cout << "  * curr chunk beg = " << mappingS[(e2 - osize / 2 + 1)] << endl;
	for (int l1 = b1 + ((int)ceil(osize * 1.0 / 2)); l1 <= e1 ; l1 ++) pc->mappingS[l1] = -1;
	for (int l2 = b2; l2 <= (e2 - osize / 2) ; l2 ++) mappingS[l2] = -1;


	//switching haplotypes
	int n_switches = 0;
	for (int i = 0 ; i < switches.size() ; i ++) {
		vector < bool > h11, h12, h21, h22;
		for (int l1 = b1 ; l1 <= e1; l1 ++) {
			h11.push_back(pc->H[l1][2*i+0]);
			h12.push_back(pc->H[l1][2*i+1]);
		}
		for (int l2= b2 ; l2 <= e2; l2 ++) {
			h21.push_back(H[l2][2*i+0]);
			h22.push_back(H[l2][2*i+1]);
		}

		int score0 = 0, score1 = 0;

		if (!pc->switches[i]) {
			score0 = hamming(h11, h21) + hamming(h12, h22);
			score1 = hamming(h11, h22) + hamming(h12, h21);
		} else {
			score0 = hamming(h12, h21) + hamming(h11, h22);
			score1 = hamming(h12, h22) + hamming(h11, h21);
		}

		//cout << i << " d=" << discordance(h11, h12, h21, h22) << " s0=" << score0 << " s1=" << score1 << endl;
		if (score0 <= score1) switches[i] = false;
		else {
			switches[i] = true;
			n_switches ++;
		}
	}
	cout << "  * #switches=" <<n_switches << " out of " << switches.size() << " possible" << endl;

}
