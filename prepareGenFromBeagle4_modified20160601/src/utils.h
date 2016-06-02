#ifndef _UTILS_H
#define _UTILS_H

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI 3.14159265358979323846

#include <string>
#include <vector>
#include <queue>
#include <map>
#include <bitset>
#include <list>
#include <tr1/unordered_map>
#include <bitset>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <pthread.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

using namespace std;
using namespace boost::iostreams;

/******************************************************/
/*                  UTILS STATISTICS                  */
/******************************************************/
namespace putils {
	void initRandom(long s);
	double getRandom();
	int getRandom(int);
	long getSeed();
	void normalise(vector < double > & v);
	int sample(vector< double > & v, double sum);
	double median (vector < double> &);
	double mean (vector < double> &);
};

/******************************************************/
/*                  UTILS ALGORITHM                   */
/******************************************************/
namespace autils {
	int max(vector < double > & v);
	int max(vector < int > & v);
	int count(vector < bool > &);
	void decompose(int min, vector < vector < int > > & B, vector < vector < vector < int > > > & BB);
	int checkDuo (int pa1, int pa2, int ca1, int ca2);
	int checkTrio (int fa1, int fa2, int ma1, int ma2, int ca1, int ca2) ;

	int solveDuoF(int & fa1, int & fa2, int & ca1, int & ca2);
	int solveDuoM(int & ma1, int & ma2, int & ca1, int & ca2);
	int solveTrio (int & fa1, int & fa2, int & ma1, int & ma2, int & ca1, int & ca2);
	int switch_error(vector < bool > & t1, vector < bool > & t2, vector < bool > & m12, vector < bool > & h1, vector < bool > & h2);


};

/******************************************************/
/*                  UTILS STRING                      */
/******************************************************/
namespace sutils {
	vector<string> tokenize(string & str,string d);
	vector<string> tokenize(string & str,string d, int n);
	string int2str(int n);
	string int2str(vector < int > & v);
	string long2str(long int n);
	string double2str(double n, int prc = 4);
	string double2str(vector < double > &v, int prc = 4);
	string bool2str(vector<bool> & v);
	string date2str(time_t * t, string format);
};

/******************************************************/
/*                  UTILS FILE                        */
/******************************************************/
namespace futils {
	bool isFile(string f);
	bool createFile(string f);
	string extensionFile(string & filename);
	void bool2binary(vector < bool > & V, ostream &fd);
	bool binary2bool(vector < bool > & V, istream & fd);
};

/******************************************************/
/*                  INPUT FILE                        */
/******************************************************/
class ifile : public filtering_istream {
private:
	string file;
	ifstream fd;

public:
	ifile();
	ifile(string filename , bool binary = false);
	~ifile();
	string name();
	bool open(string filename, bool binary = false);
	void close();
};

/******************************************************/
/*                  OUTPUT FILE                       */
/******************************************************/
class ofile : public filtering_ostream {
private:
	string file;
	ofstream fd;

public:
	ofile();
	ofile(string filename , bool binary = false);
	~ofile();
	string name();
	bool open(string filename, bool binary = false);
	void close();
};

/******************************************************/
/*                  LOG FILE                          */
/******************************************************/
class lfile {
private:
	string file;
	ofstream fd;
	bool verboseC;
	bool verboseL;

public:
	lfile();
	~lfile();
	string name();
	bool open(string filename = "file.log");
	void close();
	string getPrefix();
	void muteL();
	void unmuteL();
	void muteC();
	void unmuteC();
	void print(string s);
	void printC(string s);
	void printL(string s);
	void println(string s);
	void printlnC(string s);
	void printlnL(string s);
	void warning(string s);
	void error(string s);
};

#endif
