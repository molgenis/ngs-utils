#include "classes.h"

int main(int argc, char ** argv) {
	string buffer;
	vector < string > tok;

	char * likelihoods = NULL;
	vector < char * > posteriors;
	char * haplotypes = NULL;
	char * genotypes = NULL;
	double threshold = 0.995;

	//Reading parameters
	for (int a =1 ; a < argc ; a++) {
		if (strcmp(argv[a], "--likelihoods") == 0) {likelihoods = argv[a+1];}
		if (strcmp(argv[a], "--haplotypes") == 0) {haplotypes = argv[a+1];}
		if (strcmp(argv[a], "--genotypes") == 0) {genotypes = argv[a+1];}
		if (strcmp(argv[a], "--threshold") == 0) {threshold= atoi(argv[a+1]);}
		if (strcmp(argv[a], "--posteriors") == 0) {
			int cpt = a + 1;
			while (cpt >= 0) {
				if (argv[cpt][0] == '-') cpt = -1;
				else {
					posteriors.push_back(argv[cpt]);
					cpt ++;
				}
			}
		}
	}

	//Checking parameters
	if (likelihoods == NULL) { cout << "Missing --likelihoods arguments!" << endl; exit(1);}
	else cout << "Likelihoods:\t[" << likelihoods << "]" << endl;

	if (posteriors.size() == 0) { cout << "Missing --posteriors arguments!" << endl; exit(1);}
	else {
		cout << "Posteriors:\n" << endl;
		for (int p = 0 ; p < posteriors.size() ; p++) cout << "\t[" << posteriors[p] << "]" << endl;
	}

	if (haplotypes == NULL) { cout << "Missing --haplotypes arguments!" << endl; exit(1);}
	else cout << "Haplotypes:\t[" << haplotypes << "]" << endl;
	if (genotypes == NULL) { cout << "Missing --haplotypes arguments!" << endl; exit(1);}
	else cout << "Genotypes:\t[" << genotypes << "]" << endl;
	cout << "Threshold:\t" << threshold << endl;

	cout << endl;
	data D;
	D.readLikelihoods(likelihoods);
	for (int p = 0 ; p < posteriors.size() ; p++) D.readPosteriors(posteriors[p]);
	D.assemble();
	D.call(threshold);
	D.writeHaplotypes(haplotypes);
	D.writeGenotypes(genotypes);
}
