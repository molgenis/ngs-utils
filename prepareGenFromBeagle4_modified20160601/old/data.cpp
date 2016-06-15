#include "classes.h"

data::data() {
}

data::~data() {
	for (int i = 0 ; i < G.size() ; i ++) delete G[i];
	G.clear();
	for (int i = 0 ; i < C.size() ; i ++) delete C[i];
	C.clear();
}

void data::readLikelihoods (char * vcf) {
	string buffer;
	vector < string > tok, field, prob;

	cout << "Reading likelihoods from [" << vcf << "]" << endl;
	int line = 0;
	ifile fd(vcf);
	while(getline(fd, buffer, '\n')) {
		tok = sutils::tokenize(buffer, "\t");
		if (tok[0] == "#CHROM") {
			for (int i = 9 ; i < tok.size() ; i ++) G.push_back(new genotype(tok[i], G.size()));
			cout << "  * #inds = " << G.size() << endl;
		} else if (tok[0][0] != '#') {
			M.push(new snp (atoi(tok[0].c_str()), tok[2], atoi(tok[1].c_str()), tok[3], tok[4], M.size()));
			int idx_gl = -1;
			field = sutils::tokenize(tok[8], ":");
			for(int f = 0 ; f < field.size() && idx_gl < 0 ; f++) if (field[f] == "GL") idx_gl = f;
			if (idx_gl < 0) { cout << "No GL field in VCF file at line = " << line << endl; exit(1); }
			for (int i = 9 ; i < tok.size() ; i ++) {
				field = sutils::tokenize(tok[i], ":");
				prob = sutils::tokenize(field[idx_gl], ",");
				double l0 = atof(prob[0].c_str());
				double l1 = atof(prob[1].c_str());
				double l2 = atof(prob[2].c_str());
				G[i - 9]->L.push_back(l0);
				G[i - 9]->L.push_back(l1);
				G[i - 9]->L.push_back(l2);
				G[i - 9]->P.push_back(0.0);
				G[i - 9]->P.push_back(0.0);
				G[i - 9]->P.push_back(0.0);
				G[i - 9]->G.push_back(-1);
				G[i - 9]->h1.push_back(false);
				G[i - 9]->h2.push_back(false);
			}
		}

		if (line % 1000 == 0) { cout << "\r" << line ;cout.flush();}
		line ++;
	}
	cout << "\r  * #sites = " << M.size() << endl;
}

void data::readPosteriors (char * vcf) {
	string buffer;
	vector < string > tok, field, prob;
	chunk * c = new chunk(vcf);

	cout << "Reading posteriors from [" << vcf << "]" << endl;
	int line = 0;
	ifile fd(vcf);
	while(getline(fd, buffer, '\n')) {
		sutils::tokenize(buffer, "\t");
		if (tok[0] == "#CHROM") {
			for (int i = 9 ; i < tok.size() ; i ++)
				if (G[i-9]->id != tok[i]) { cout << "Sample inconsistency like=" << G[i-9]->id << " post=" << tok[i] << endl;exit(1); }
			c->switches = vector < bool > (G.size(), false);
			cout << "  * All " << G.size() << " samples found" << endl;
		} else if (tok[0][0] != '#') {
			snp * s = M.get(atoi(tok[1].c_str()), tok[3], tok[4]);
			if (s == NULL) { cout << "Site unfound in likelihoods file! pos = " << atoi(tok[1].c_str()) << endl; exit(1); }
			c->mappingS.push_back(s->idx);
			c->H.push_back(vector < bool > (2 * G.size(), false));
			c->P.push_back(vector < float > (3 * G.size(), 0.0));
			for (int i = 9 ; i < tok.size() ; i ++) {
				field = sutils::tokenize(tok[i], ":");
				prob = sutils::tokenize(field[2], ",");
				c->H.back()[(i-9) * 2 + 0] = (tok[i][0] == '1');
				c->H.back()[(i-9) * 2 + 1] = (tok[i][2] == '1');
				c->P.back()[(i-9) * 3 + 0] = atof(prob[0].c_str());
				c->P.back()[(i-9) * 3 + 1] = atof(prob[1].c_str());
				c->P.back()[(i-9) * 3 + 2] = atof(prob[2].c_str());
			}
		}

		if (line % 1000 == 0) { cout << "\r" << line ;cout.flush();}
		line ++;
	}
	cout << "\r  * #sites = " << c->mappingS.size() << endl;
	C.push_back(c);
}

void data::assemble() {
	cout << "Sorting the " << C.size() << " chunks of data" << endl;
	sort(C.begin(), C.end());
	for (int c = 1 ; c < C.size() ; c++) {
		cout << "Assembling chunk " << c -1 << " witch chunk " << c << endl;
		C[c]->assemble(*C[c-1]);
	}

	cout << "Copying chunk genotype & haplotype data" << endl;
	for (int c = 0 ; c < C.size() ; c++) {
		for (int l = 0 ; l < C[c]->mappingS.size() ; l ++) {
			if (C[c]->mappingS[l] >= 0) {
				//copying posteriors
				for (int i = 0 ; i < G.size() ; i ++) {
					G[i]->P[3 * C[c]->mappingS[l] + 0] = C[c]->P[l][3 * i + 0];
					G[i]->P[3 * C[c]->mappingS[l] + 1] = C[c]->P[l][3 * i + 1];
					G[i]->P[3 * C[c]->mappingS[l] + 2] = C[c]->P[l][3 * i + 2];
				}
				//copying haplotypes
				for (int i = 0 ; i < G.size() ; i ++) {
					if (!C[c]->switches[i]) {
						G[i]->h1[C[c]->mappingS[l]] = C[c]->H[l][2 * i + 0];
						G[i]->h2[C[c]->mappingS[l]] = C[c]->H[l][2 * i + 1];
					} else {
						G[i]->h1[C[c]->mappingS[l]] = C[c]->H[l][2 * i + 1];
						G[i]->h2[C[c]->mappingS[l]] = C[c]->H[l][2 * i + 0];
					}
				}
			}
		}
	}
}

void data::call(double threshold) {
	int n_called = 0;
	int n_total = 0;
	cout << "Calling genotypes from posteriors" << endl;
	for (int i = 0 ; i < G.size() ; i ++) {
		n_called += G[i]->call(threshold);
		n_total += M.size();
	}
	cout << "  * #called = " << sutils::double2str(n_called * 100 / n_total, 1) << "%" << endl;
}

void data::writeHaplotypes(char * hap) {
	string fhap = string(hap) + ".hap.gz";
	string fsam = string(hap) + ".hap.sample";

	cout << "Writing samples in [" << fsam << "]" << endl;
	ofile fds (fhap);
	fds << "ID_1 ID_2 missing" << endl << "0 0 0" << endl;
	for (int i = 0 ; i < G.size() ; i ++) fds << G[i]->id << " " << G[i]->id << " 0" << endl;
	fds.close();

	cout << "Writing haplotypes in [" << fhap << "]" << endl;
	ofile fdh (fsam);
	for (int s = 0 ; s < M.size() ; s ++) {
		fdh << M.vec_snp[s]->chr << " " << M.vec_snp[s]->id << " " << M.vec_snp[s]->pos << " " << M.vec_snp[s]->a0 << " " << M.vec_snp[s]->a1;
		for (int i = 0 ; i < G.size() ; i ++) fdh << " " << G[i]->h1[s] << " " << G[i]->h2[s];
		fdh << endl;
	}
	fdh.close();
}

void data::writeGenotypes(char * gen) {
	string fgen = string(gen) + ".gen.gz";
	string fsam = string(gen) + ".gen.sample";

	cout << "Writing samples in [" << fsam << "]" << endl;
	ofile fds (fsam);
	fds << "ID_1 ID_2 missing" << endl << "0 0 0" << endl;
	for (int i = 0 ; i < G.size() ; i ++) fds << G[i]->id << " " << G[i]->id << " 0" << endl;
	fds.close();

	cout << "Writing genotypes in [" << fgen << "]" << endl;
	ofile fdg (fgen);
	for (int s = 0 ; s < M.size() ; s ++) {
		fdg << M.vec_snp[s]->chr << " " << M.vec_snp[s]->id << " " << M.vec_snp[s]->pos << " " << M.vec_snp[s]->a0 << " " << M.vec_snp[s]->a1;
		for (int i = 0 ; i < G.size() ; i ++) {
			switch (G[i]->G[s]) {
			case 0:	fdg << " 1 0 0"; break;
			case 1:	fdg << " 0 1 0"; break;
			case 2:	fdg << " 0 0 1"; break;
			default: fdg << " " << G[i]->P[3 * i + 0] << " " << G[i]->P[3 * i + 1] << " " << G[i]->P[3 * i + 2]; break;
			}
		}
		fdg << endl;
	}
	fdg.close();
}
