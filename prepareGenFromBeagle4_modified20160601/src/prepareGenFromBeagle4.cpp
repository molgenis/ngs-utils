#include "classes.h"

int main(int argc, char ** argv) {
	string buffer;
	vector < string > tok;

	char * likelihoods = NULL;
	vector < char * > posteriors;
	char * output = NULL;
	double threshold = 0.995;
	char * haploid = NULL;

	//Reading parameters
	for (int a =1 ; a < argc ; a++) {
		if (strcmp(argv[a], "--likelihoods") == 0) {likelihoods = argv[a+1];}
		if (strcmp(argv[a], "--output") == 0) {output = argv[a+1];}
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
		if (strcmp(argv[a], "--haploid") == 0) {haploid = argv[a+1];}
	}


	//Checking parameters
	if (likelihoods == NULL) { cout << "Missing --likelihoods argument!" << endl; exit(1);}
	else cout << "Likelihoods:\t[" << likelihoods << "]" << endl;

	if (posteriors.size() == 0) { cout << "Missing --posteriors argument!" << endl; exit(1);}
	else cout << "Posteriors:\t" << posteriors.size() << " VCF files" << endl;

	if (output == NULL) { cout << "Missing --output argument!" << endl; exit(1);}
	else cout << "Output Prefix:\t[" << output << "]" << endl;
	cout << "Threshold:\t" << threshold << endl;

	cout << endl;
	data D = data(string(likelihoods), threshold);
	D.readLikelihoodsMap();
	if (haploid != NULL) D.readHaploidGuys(haploid);
	for (int p = 0 ; p < posteriors.size() ; p++) D.readPosteriorsMapAndHap(posteriors[p]);
	D.assembleMap();
	D.writeHaplotypes(output);
	D.writeGenotypes(output);
}
