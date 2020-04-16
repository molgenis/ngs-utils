#include "utils.h"

/******************************************************/
/*                  UTILS STATISTICS                  */
/******************************************************/
long seed = -123456789;

namespace putils {
	void initRandom(long s) {
		if (s == -1) seed = - time(NULL);
		else seed = - s;
	}

	double getRandom() {
		int j;
		long k;
		static long iy=0;
		static long iv[NTAB];
		double temp;
		if (seed <= 0 || !iy) {
			if (-(seed) < 1) seed=1;
			else seed = -(seed);
			for (j=NTAB+7;j>=0;j--) {
				k=(seed)/IQ;
				seed=IA*(seed-k*IQ)-IR*k;
				if (seed < 0) seed += IM;
				if (j < NTAB) iv[j] = seed;
			}
			iy=iv[0];
		}
		k=(seed)/IQ;
		seed=IA*(seed-k*IQ)-IR*k;
		if (seed < 0) seed += IM;
		j=iy/NDIV;
		iy=iv[j];
		iv[j] = seed;
		if ((temp=AM*iy) > RNMX) return RNMX;
		else return temp;
	}

	long getSeed() {
		return seed;
	}

	int getRandom(int n) {
		return (int)floor(getRandom() * n);
	}

	void normalise(vector < double > & v) {
		double sum = 0.0;
		for (int i = 0 ; i < v.size() ; i++) sum += v[i];
		if (sum != 0.0) for (int i = 0 ; i < v.size() ; i++) v[i] /= sum;
	}

	int sample(vector< double > & v, double sum) {
		double csum = v[0];
		double u = getRandom() * sum;
		for (int i = 0; i < v.size() - 1; ++i) {
			if ( u < csum ) return i;
			csum += v[i+1];
		}
		return v.size() - 1;
	}

	double median (vector < double> & V) {
		sort(V.begin(), V.end());
		return V[V.size() / 2];
	}

	double mean (vector < double> & V) {
		double m  = 0.0;
		for (int i = 0 ; i < V.size() ; i ++) m += V[i];
		return m / V.size();
	}
};

/******************************************************/
/*                  UTILS ALGORITHM                   */
/******************************************************/
namespace autils {
	int max(vector < double > & v) {
		double max = -1e300;
		int index_max = 0;
		for (int i = 0 ; i < v.size() ; i ++)
			if (v[i] > max) {
				max = v[i];
				index_max = i;
				}
		return index_max;
	}

	int max(vector < int > & v) {
		int max = -1000000000;
		int index_max = 0;
		for (int i = 0 ; i < v.size() ; i ++)
		if (v[i] > max) {
			max = v[i];
			index_max = i;
		}
		return index_max;
	}

	int count(vector < bool > &v) {
		int c=0;
		for (int i = 0 ; i < v.size() ; i ++) if (v[i]) c++;
		return c;
	}

	void decompose(int min, vector < vector < int > > & B, vector < vector < vector < int > > > & BB) {
		if (B.size() < 2 * min || B.size() == 2) {
			BB.push_back(B);
			return;
		}

		int l = putils::getRandom(B.size() - 2 * min) + min;
		vector < vector < int > > LB = vector < vector < int > > (B.begin(), B.begin() + l + 1);
		vector < vector < int > > RB = vector < vector < int > > (B.begin() + l - 1, B.end());
		vector < vector < vector < int > > > LBB;
		vector < vector < vector < int > > > RBB;
		decompose(min, LB, LBB);
		decompose(min, RB, RBB);
		BB = LBB;
		BB.insert(BB.end(), RBB.begin(), RBB.end());
	}

	int checkDuo (int pa1, int pa2, int ca1, int ca2) {
		if (pa1 == 0 && pa2 == 0 && ca1 == 0 && ca2 == 0) { return 0; }
		if (pa1 == 0 && pa2 == 0 && ca1 == 0 && ca2 == 1) { return 0; }
		if (pa1 == 0 && pa2 == 0 && ca1 == 1 && ca2 == 0) { return 0; }
		if (pa1 == 0 && pa2 == 0 && ca1 == 1 && ca2 == 1) { return -1; }
		if (pa1 == 0 && pa2 == 1 && ca1 == 0 && ca2 == 0) { return 1; }
		if (pa1 == 0 && pa2 == 1 && ca1 == 0 && ca2 == 1) { return 0; }
		if (pa1 == 0 && pa2 == 1 && ca1 == 1 && ca2 == 0) { return 0; }
		if (pa1 == 0 && pa2 == 1 && ca1 == 1 && ca2 == 1) { return 1; }
		if (pa1 == 1 && pa2 == 0 && ca1 == 0 && ca2 == 0) { return 1; }
		if (pa1 == 1 && pa2 == 0 && ca1 == 0 && ca2 == 1) { return 0; }
		if (pa1 == 1 && pa2 == 0 && ca1 == 1 && ca2 == 0) { return 0; }
		if (pa1 == 1 && pa2 == 0 && ca1 == 1 && ca2 == 1) { return 1; }
		if (pa1 == 1 && pa2 == 1 && ca1 == 0 && ca2 == 0) { return -1; }
		if (pa1 == 1 && pa2 == 1 && ca1 == 0 && ca2 == 1) { return 0; }
		if (pa1 == 1 && pa2 == 1 && ca1 == 1 && ca2 == 0) { return 0; }
		if (pa1 == 1 && pa2 == 1 && ca1 == 1 && ca2 == 1) { return 0; }
		return 0;
	}

	int checkTrio (int fa1, int fa2, int ma1, int ma2, int ca1, int ca2) {
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 0; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return -1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return -1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return 1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return -1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 0; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 0; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return -1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return 2; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 0; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 0; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 2; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 2; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 0; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 0; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return 2; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return 2; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 0; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 0; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 2; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 2; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 0; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 0; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return 2; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return -1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 0; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 0; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return -1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return -1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return -1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 0; }
		return 0;
	}


	int solveDuoF(int & fa1, int & fa2, int & ca1, int & ca2) {
		int gf = fa1 + fa2;
		int gc = ca1 + ca2;

		if (gf == 0 && gc == 2) { return -1;};
		if (gf == 2 && gc == 0) { return -1;};

		if (gf == 0 && gc == 1) {fa1 = 0; fa2 = 0; ca1 = 0; ca2 = 1; return 1;};
		if (gf == 1 && gc == 0) {fa1 = 0; fa2 = 1; ca1 = 0; ca2 = 0; return 1;};
		if (gf == 1 && gc == 2) {fa1 = 1; fa2 = 0; ca1 = 1; ca2 = 1; return 1;};
		if (gf == 2 && gc == 1) {fa1 = 1; fa2 = 1; ca1 = 1; ca2 = 0; return 1;};
		return 0;
	}

	int solveDuoM(int & ma1, int & ma2, int & ca1, int & ca2) {
		int gm = ma1 + ma2;
		int gc = ca1 + ca2;

		if (gm == 0 && gc == 2) { return -1; };
		if (gm == 2 && gc == 0) { return -1; };

		if (gm == 0 && gc == 1) {ma1 = 0; ma2 = 0; ca1 = 1; ca2 = 0; return 1;};
		if (gm == 1 && gc == 0) {ma1 = 1; ma2 = 0; ca1 = 0; ca2 = 0; return 1;};
		if (gm == 1 && gc == 2) {ma1 = 0; ma2 = 1; ca1 = 1; ca2 = 1; return 1;};
		if (gm == 2 && gc == 1) {ma1 = 1; ma2 = 1; ca1 = 0; ca2 = 1; return 1;};
		return 0;
	}


	int solveTrio (int & fa1, int & fa2, int & ma1, int & ma2, int & ca1, int & ca2) {
		int gf = fa1 + fa2;
		int gm = ma1 + ma2;
		int gc = ca1 + ca2;

		if (gf == 0 && gm == 0 && gc == 1) { return -1;};
		if (gf == 0 && gm == 0 && gc == 2) { return -1;};
		if (gf == 0 && gm == 1 && gc == 2) { return -1;};
		if (gf == 0 && gm == 2 && gc == 0) { return -1;};
		if (gf == 0 && gm == 2 && gc == 2) { return -1;};
		if (gf == 1 && gm == 0 && gc == 2) { return -1;};
		if (gf == 1 && gm == 2 && gc == 0) { return -1;};
		if (gf == 2 && gm == 0 && gc == 0) { return -1;};
		if (gf == 2 && gm == 0 && gc == 2) { return -1;};
		if (gf == 2 && gm == 1 && gc == 0) { return -1;};
		if (gf == 2 && gm == 2 && gc == 0) { return -1;};
		if (gf == 2 && gm == 2 && gc == 1) { return -1;};

		if (gf == 0 && gm == 1 && gc == 0) {fa1 = 0; fa2 = 0; ma1 = 1; ma2 = 0; ca1 = 0; ca2 = 0; return 1;};
		if (gf == 0 && gm == 1 && gc == 1) {fa1 = 0; fa2 = 0; ma1 = 0; ma2 = 1; ca1 = 0; ca2 = 1; return 1;};
		if (gf == 0 && gm == 2 && gc == 1) {fa1 = 0; fa2 = 0; ma1 = 1; ma2 = 1; ca1 = 0; ca2 = 1; return 1;};
		if (gf == 1 && gm == 0 && gc == 0) {fa1 = 0; fa2 = 1; ma1 = 0; ma2 = 0; ca1 = 0; ca2 = 0; return 1;};
		if (gf == 1 && gm == 0 && gc == 1) {fa1 = 1; fa2 = 0; ma1 = 0; ma2 = 0; ca1 = 1; ca2 = 0; return 1;};
		if (gf == 1 && gm == 1 && gc == 0) {fa1 = 0; fa2 = 1; ma1 = 1; ma2 = 0; ca1 = 0; ca2 = 0; return 1;};
		if (gf == 1 && gm == 1 && gc == 2) {fa1 = 1; fa2 = 0; ma1 = 0; ma2 = 1; ca1 = 1; ca2 = 1; return 1;};
		if (gf == 1 && gm == 2 && gc == 1) {fa1 = 0; fa2 = 1; ma1 = 1; ma2 = 1; ca1 = 0; ca2 = 1; return 1;};
		if (gf == 1 && gm == 2 && gc == 2) {fa1 = 1; fa2 = 0; ma1 = 1; ma2 = 1; ca1 = 1; ca2 = 1; return 1;};
		if (gf == 2 && gm == 0 && gc == 1) {fa1 = 1; fa2 = 1; ma1 = 0; ma2 = 0; ca1 = 1; ca2 = 0; return 1;};
		if (gf == 2 && gm == 1 && gc == 1) {fa1 = 1; fa2 = 1; ma1 = 1; ma2 = 0; ca1 = 1; ca2 = 0; return 1;};
		if (gf == 2 && gm == 1 && gc == 2) {fa1 = 1; fa2 = 1; ma1 = 0; ma2 = 1; ca1 = 1; ca2 = 1; return 1;};
		return 0;
	}

	int switch_error(vector < bool > & t1, vector < bool > & t2, vector < bool > & m12, vector < bool > & h1, vector < bool > & h2) {
		int error = 0;

		//cerr << utils::bool2str(t1) << endl;
		//cerr << utils::bool2str(t2) << endl;
		//cerr << utils::bool2str(h1) << endl;
		//cerr << utils::bool2str(h2) << endl;

		assert(t1.size() == h1.size());
		assert(t1.size() == h2.size());
		assert(t1.size() == t2.size());
		assert(t1.size() == m12.size());

		int locus = -1;

		for (int l = 0 ; l < t1.size() && locus < 0 ; l ++) if (t1[l] != t2[l] && m12[l]) locus = l;

		if (locus < 0) return -1;

		bool swap = (t1[locus] == h2[locus]);
		for (int l = locus + 1 ; l < t1.size() ; l ++)
			if (m12[l]) {
				if (!swap && t1[l] != h1[l]) {
					swap = true;
					error++;
				} else if (swap && t1[l] != h2[l]) {
					swap = false;
					error++;
				}
			}

		return error;
	}

};

/******************************************************/
/*                  UTILS STRING                      */
/******************************************************/
namespace sutils {
	vector<string> tokenize(string & str,string d) {
		vector<string> tokens;
		string::size_type lastPos = str.find_first_not_of(d, 0);
		string::size_type pos = str.find_first_of(d, lastPos);
		while (string::npos != pos || string::npos != lastPos) {
			tokens.push_back(str.substr(lastPos, pos - lastPos));
			lastPos = str.find_first_not_of(d, pos);
			pos = str.find_first_of(d, lastPos);
		}
		return tokens;
	}

	vector<string> tokenize(string & str,string d, int n) {
		vector<string> tokens;
		int cpt = 0;
		string::size_type lastPos = str.find_first_not_of(d, 0);
		string::size_type pos = str.find_first_of(d, lastPos);
		while ((string::npos != pos || string::npos != lastPos) && cpt < n) {
			tokens.push_back(str.substr(lastPos, pos - lastPos));
			lastPos = str.find_first_not_of(d, pos);
			pos = str.find_first_of(d, lastPos);
			cpt++;
		}
		return tokens;
	}

	string int2str(int n) {
		ostringstream s2( stringstream::out );
		s2 << n;
		return s2.str();
	}

	string int2str(vector < int > & v) {
		ostringstream s2( stringstream::out );
		for (int l = 0 ; l < v.size() ; l++) {
			if (v[l] < 0) s2 << "-";
			else s2 << v[l] ;
		}
		return s2.str();
	}

	string long2str(long int n) {
		ostringstream s2( stringstream::out );
		s2 << n;
		return s2.str();
	}

	string double2str(double n, int prc) {
		ostringstream s2;
		s2 << setiosflags( ios::fixed );
		if ( prc > 0 ) s2.precision(prc);
		s2 << n;
		return s2.str();
	}

	string double2str(vector < double > &v, int prc) {
		ostringstream s2;
		s2 << setiosflags( ios::fixed );
		if ( prc >= 0 ) s2.precision(prc);
		for (int l = 0 ; l < v.size() ; l++) {
			s2 << v[l] << " ";
			}
		return s2.str();
	}

	string bool2str(vector<bool> & v) {
		ostringstream s2( stringstream::out );
		for (int l = 0 ; l < v.size() ; l++) {
			if (v[l]) s2 << "1";
			else s2 << "0";
		}
		return s2.str();
	}

	string date2str(time_t * t, string format) {
		struct tm * timeinfo = localtime(t);
		char buffer[128];
		strftime(buffer, 128, format.c_str(), timeinfo);
		ostringstream s2( stringstream::out );
		s2 << buffer;
		return s2.str();
	}
};

/******************************************************/
/*                  UTILS FILE                        */
/******************************************************/
namespace futils {
	bool isFile(string f) {
		ifstream inp;
		inp.open(f.c_str(), ifstream::in);
		if(inp.fail()) {
			inp.clear(ios::failbit);
			inp.close();
			return false;
		}
		inp.close();
		return true;
	}

	bool createFile(string f) {
		ofstream out;
		out.open(f.c_str(), ofstream::out);
		if(out.fail()) {
			out.clear(ios::failbit);
			out.close();
			return false;
		}
		out.close();
		return true;
	}

	string extensionFile(string & filename) {
		if (filename.find_last_of(".") != string::npos)
			return filename.substr(filename.find_last_of(".") + 1);
		return "";
	}


	void bool2binary(vector < bool > & V, ostream &fd) {
		int nb = V.size();
		fd.write((char*)&nb, sizeof(int));
		int cpt_byte = 0;
		int cpt_bin = 0;
		int nb_byte = (int)ceil( (V.size() * 1.0) / 8);
		while (cpt_byte < nb_byte) {
			bitset<8> byte_bitset;
			for (int l = 0; l < 8 && cpt_bin < V.size() ; l++) {
				byte_bitset[l] = V[cpt_bin];
				cpt_bin ++;
			}
			char byte_char = (char)byte_bitset.to_ulong();
			fd.write(&byte_char, 1);
			cpt_byte++;
		}
	}

	bool binary2bool(vector < bool > & V, istream & fd) {
		int nb;
		fd.read((char*)&nb, sizeof(int));
		if (!fd) return false;
		int cpt_byte = 0;
		int cpt_bin = 0;
		int nb_byte = (int)ceil( (nb * 1.0) / 8);
		V = vector < bool >(nb);
		while (cpt_byte < nb_byte) {
			char byte_char;
			fd.read(&byte_char, 1);
			if (!fd) return false;
			bitset<8> byte_bitset = byte_char;
			for (int l = 0; l < 8 && cpt_bin < nb ; l++) {
				V[cpt_bin] = byte_bitset[l];
				cpt_bin++;
			}
			cpt_byte++;
		}
		return true;
	}
};

/******************************************************/
/*                  INPUT FILE                        */
/******************************************************/
ifile::ifile() {
}

ifile::ifile(string filename , bool binary) {
	open(filename, binary);
}

ifile::~ifile() {
	close();
}

string ifile::name() {
	return file;
}

bool ifile::open(string filename, bool binary) {
	file = filename;
	string ext = futils::extensionFile(filename);
	if (ext == "gz") {
		fd.open(file.c_str(), ios::in | ios::binary);
		push(gzip_decompressor());
	} else if (ext == "bz2") {
		fd.open(file.c_str(), ios::in | ios::binary);
		push(bzip2_decompressor());
	} else if (binary) {
		fd.open(file.c_str(), ios::in | ios::binary);
	} else  {
		fd.open(file.c_str());
	}
	if (fd.fail()) return false;
	push(fd);
	return true;
}

void ifile::close() {
	 if (!empty()) reset();
	 fd.close();
}

/******************************************************/
/*                  OUTPUT FILE                       */
/******************************************************/
ofile::ofile() {
}

ofile::ofile(string filename , bool binary) {
	open(filename, binary);
}

ofile::~ofile() {
	close();
}

string ofile::name() {
	return file;
}

bool ofile::open(string filename, bool binary) {
	file = filename;
	string ext = futils::extensionFile(filename);
	if (ext == "gz") {
		fd.open(file.c_str(), ios::out | ios::binary);
		push(gzip_compressor());
	} else if (ext == "bz2") {
		fd.open(file.c_str(), ios::out | ios::binary);
		push(bzip2_compressor());
	} else if (binary) {
		fd.open(file.c_str(), ios::out | ios::binary);
	} else  {
		fd.open(file.c_str());
	}
	if (fd.fail()) return false;
	push(fd);
	return true;
}

void ofile::close() {
	 if (!empty()) reset();
	 fd.close();
}

/******************************************************/
/*                  LOG FILE                          */
/******************************************************/
lfile::lfile() {
	verboseC = true;
	verboseL = true;
}

lfile::~lfile() {
	close();
}

string lfile::name() {
	return file;
}

bool lfile::open(string filename) {
	file = filename;
	if (futils::extensionFile(file) != "log") file += ".log";
	fd.open(file.c_str());
	if (fd.fail()) return false;
	return true;
}

void lfile::close() {
	 fd.close();
}

string lfile::getPrefix() {
	return file.substr(0, file.find_last_of("."));
}

void lfile::muteL() {
	verboseL = false;
}

void lfile::unmuteL() {
	verboseL = true;
}

void lfile::muteC() {
	verboseC = false;
}

void lfile::unmuteC() {
	verboseC = true;
}


void lfile::print(string s) {
	if (verboseL) { fd << s; fd.flush(); }
	if (verboseC) { cerr << s; cerr.flush(); }
}

void lfile::printC(string s) {
	if (verboseC) { cerr << s; cerr.flush(); }
}

void lfile::printL(string s) {
	if (verboseL) { fd << s; fd.flush(); }
}

void lfile::println(string s) {
	if (verboseL) { fd << s << endl; }
	if (verboseC) { cerr << s << endl; }
}

void lfile::printlnC(string s) {
	if (verboseC) { cerr << s << endl; }
}

void lfile::printlnL(string s) {
	if (verboseL) { fd << s << endl; }
}

void lfile::warning(string s) {
	cerr << "\33[33mWARNING:\33[0m " << s << endl;
	if (verboseL) fd << "WARNING: " << s << endl;
}

void lfile::error(string s) {
	cerr << "\33[33mERROR:\33[0m " << s << endl;
	if (verboseL) fd << "ERROR: " << s << endl;
	close();
	exit(1);
}
