#include "classes.h"
#include <iostream>

data::data(string id, double threshold) {
	this->id = id;
	this->threshold = threshold;
}

data::~data() {
	for (int i = 0 ; i < C.size() ; i ++) delete C[i];
	C.clear();
}

void data::readLikelihoodsMap () {
	string buffer;
	vector < string > tok, field, prob;

	cout << "Reading Likelihoods Map in [" << id << "]" << endl;
	int line = 0;
	ifile fd(id);
	while(getline(fd, buffer, '\n')) {
		if (buffer[0] == '#' && buffer[1] == 'C' && buffer[2] == 'H' && buffer[3] == 'R' && buffer[4] == 'O') {
			tok = sutils::tokenize(buffer, "\t");
			for (int i = 9 ; i < tok.size() ; i ++) {
				I.push_back(tok[i]);
				males.push_back(false);
			}
			cout << "  * #inds = " << I.size() << endl;
		} else if (buffer[0] != '#') {
			tok = sutils::tokenize(buffer, "\t", 5);
			vector < string > alleles_alt = sutils::tokenize(tok[4], ",");
			if (alleles_alt[0] != "." && alleles_alt[0] != "0" && alleles_alt.size() == 1)
				M.push(new snp (atoi(tok[0].c_str()), tok[2], atoi(tok[1].c_str()), tok[3], tok[4], M.size()));
		}

		if (line > 0 && line % 1000 == 0) { cout << "\r" << line ;cout.flush();}
		line ++;
	}
	cout << "\r  * #sites = " << M.size() << endl;
}

void data::readHaploidGuys(char * haploid) {
	string buffer;
	vector < string > tok;
	ifile fd(haploid);
	int n_set = 0;
	while(getline(fd, buffer, '\n')) {
		bool done = false;
		for (int i = 0 ; i  < I.size() && !done ; i ++) if (I[i] == buffer) {
			males[i] = true;
			done = true;
			n_set++;
		}
	}
	cerr << "  * #males=" << n_set << endl;
	fd.close();
}

void data::readPosteriorsMapAndHap (char * vcf) {
	string buffer;
	vector < string > tok, field, prob;
	chunk * c = new chunk(vcf);

	cout << "Reading Posteriors Map in [" << vcf << "]";
	ifile fd(vcf);
	while(getline(fd, buffer)) {
		if (buffer[0] == '#' && buffer[1] == 'C' && buffer[2] == 'H' && buffer[3] == 'R' && buffer[4] == 'O') {
			tok = sutils::tokenize(buffer, "\t");
			for (int i = 9 ; i < tok.size() ; i ++) if (I[i-9] != tok[i]) { cout << "Sample inconsistency like=" << I[i-9] << " post=" << tok[i] << endl;exit(1); }
			c->switches = vector < bool > (I.size(), false);
		} else if (buffer[0] != '#') {
			tok = sutils::tokenize(buffer, "\t", 5);
			vector < string > alleles_alt = sutils::tokenize(tok[4], ",");
			if (alleles_alt[0] != "." && alleles_alt[0] != "0" && alleles_alt.size() == 1) {
				snp * s = M.get(atoi(tok[1].c_str()), tok[3], tok[4]);
				if (s == NULL) { cout << "\nSite unfound in likelihoods file! pos = " << atoi(tok[1].c_str()) << endl; exit(1); }
				c->mappingS.push_back(s->idx);

				int n_tok = 9 + I.size();
				int n_cha = buffer.size();
				int i_cha = 0;
				c->H.push_back(vector < bool > (2 * I.size(), false));
				for (int t = 0 ; t < n_tok ; t ++) {
					if (t >= 9) {
						if (buffer[i_cha+0] == '0') c->H.back()[(t-9) * 2 + 0] = false;
						else c->H.back()[(t-9) * 2 + 0] = true;
						if (buffer[i_cha+2] == '0') c->H.back()[(t-9) * 2 + 1] = false;
						else c->H.back()[(t-9) * 2 + 1] = true;
					}
					for (; buffer[i_cha] != '\t' && i_cha < n_cha ; i_cha ++);
					i_cha ++;
				}
			}
		}
	}
	cout << "\t[#sites = " << c->mappingS.size() << ", first=" << M.vec_snp[c->mappingS[0]]->pos << ", last =" << M.vec_snp[c->mappingS.back()]->pos << "]" << endl;
	C.push_back(c);
}

struct less_than_chunk {
	inline bool operator() (const chunk * c1, const chunk * c2) {
		return (c1->mappingS[0] < c2->mappingS[0]);
    }
};

void data::assembleMap() {
	cout << "Sorting the " << C.size() << " chunks of data" << endl;
	sort(C.begin(), C.end(), less_than_chunk());
	for (int c = 1 ; c < C.size() ; c++) {
		cout << "Assembling [" << C[c -1]->id << "] with [" << C[c]->id << "]" << endl;
		C[c]->assemble(C[c-1]);
	}
	cout << "Mapping sites to chunks" << endl;
	for (int c = 0 ; c < C.size() ; c++) {
		for (int l = 0 ; l < C[c]->mappingS.size() ; l ++) {
			if (C[c]->mappingS[l] >= 0) {
				M.vec_snp[C[c]->mappingS[l]]->chunk_id = c;
				M.vec_snp[C[c]->mappingS[l]]->chunk_idx = l;
			}
		}
	}
	cout << "Checking mapping sites to chunks" << endl;
	for (int l = 0 ; l < M.vec_snp.size() ; l ++)
		if (M.vec_snp[l]->chunk_id < 0)
			cout << "Site " << l << " with pos=" << M.vec_snp[l]->pos << " is unmapped" << endl;
}


void data::writeHaplotypes(char * output) {
	string fhap = string(output) + ".hap.gz";
	string fsam = string(output) + ".hap.sample";

	cout << "Writing samples in [" << fsam << "]" << endl;
	ofile fds (fsam);
	fds << "ID_1 ID_2 missing" << endl << "0 0 0" << endl;
	for (int i = 0 ; i < I.size() ; i ++) fds << I[i] << " " << I[i] << " 0" << endl;
	fds.close();

	cout << "Writing haplotypes in [" << fhap << "]" << endl;
	ofile fdh (fhap);
	for (int s = 0 ; s < M.size() ; s ++) {
		fdh << M.vec_snp[s]->chr << " " << M.vec_snp[s]->id << " " << M.vec_snp[s]->pos << " " << M.vec_snp[s]->a0 << " " << M.vec_snp[s]->a1;
		for (int i = 0 ; i < I.size() ; i ++) {
			int chunk_id = M.vec_snp[s]->chunk_id;
			int chunk_idx = M.vec_snp[s]->chunk_idx;
			if (C[chunk_id]->switches[i]) fdh << " " << C[chunk_id]->H[chunk_idx][2 * i + 1] << " " << C[chunk_id]->H[chunk_idx][2 * i + 0];
			else fdh << " " << C[chunk_id]->H[chunk_idx][2 * i + 0] << " " << C[chunk_id]->H[chunk_idx][2 * i + 1];
		}
		fdh << endl;
	}
	fdh.close();
}

void data::writeGenotypes(char * output) {
	string fgen = string(output) + ".gen.gz";
	string fsam = string(output) + ".gen.sample";

	cout << "Writing samples in [" << fsam << "]" << endl;
	ofile fds (fsam);
	fds << "ID_1 ID_2 missing" << endl << "0 0 0" << endl;
	for (int i = 0 ; i < I.size() ; i ++) fds << I[i] << " " << I[i] << " 0" << endl;
	fds.close();

	cout << "Writing genotypes in [" << fgen << "]" << endl;
	string like_buffer = "#";
	int prev_chunk = -1;
	long geno_post = 0;
	long geno_like = 0;

	ofile fdg (fgen);
	ifile fdl (id);
	ifile fdc;

	for (int s = 0 ; s < M.size();) {
		getline(fdl, like_buffer, '\n');
		if (like_buffer[0] != '#') {
			vector < string > like_tok = sutils::tokenize(like_buffer, "\t");
			vector < string > alleles_alt = sutils::tokenize(like_tok[4], ",");

//cout << like_tok[1] << "\n";
			if (alleles_alt[0] != "0" && alleles_alt[0] != "." && alleles_alt.size() == 1) {
				//get current site
				snp * site = M.vec_snp[s];

				//get corresponding chunk
				int curr_chunk = site->chunk_id;
				int curr_chunk_idx = site->chunk_idx;

				//if chunk not opened, open it
				string chunk_buffer = "#";
				if (curr_chunk != prev_chunk) {
					if (prev_chunk >= 0) fdc.close();
					fdc.open(C[curr_chunk]->id, false);
					prev_chunk = curr_chunk;
				}

				//parse chunk to relevant site
				bool done = false;
				while (! done) {
					getline(fdc, chunk_buffer, '\n');
					vector < string > chunk_tok2 = sutils::tokenize(chunk_buffer, "\t", 5);
					if (chunk_buffer[0] != '#' && atoi(chunk_tok2[1].c_str()) == site->pos && chunk_tok2[3] == site->a0 && chunk_tok2[4] == site->a1) done = true;
				}

				//tokenize line of chunk
				vector < string > chunk_tok = sutils::tokenize(chunk_buffer, "\t");

				assert(chunk_tok[1] == like_tok[1]);
				assert(chunk_tok[3] == like_tok[3]);
				assert(chunk_tok[4] == like_tok[4]);
				assert(chunk_tok.size() == like_tok.size());
				assert(atoi(like_tok[1].c_str()) == site->pos);
				assert(like_tok[3] == site->a0);
				assert(like_tok[4] == site->a1);

				//find GL field in like
				int like_idx_gl = -1;
				vector < string > like_field = sutils::tokenize(like_tok[8], ":");
				for(int f = 0 ; f < like_field.size() && like_idx_gl < 0 ; f++) if (like_field[f] == "GL") like_idx_gl = f;
				assert(like_idx_gl >= 0);

				//find GP field in chunk
				int chunk_idx_gp = -1;
				vector < string > chunk_field = sutils::tokenize(chunk_tok[8], ":");
				for(int f = 0 ; f < chunk_field.size() && chunk_idx_gp < 0 ; f++) if (chunk_field[f] == "GP") chunk_idx_gp = f;
				assert(chunk_idx_gp >= 0);

				//write site to output
				fdg << M.vec_snp[s]->chr << " " << M.vec_snp[s]->id << " " << M.vec_snp[s]->pos << " " << M.vec_snp[s]->a0 << " " << M.vec_snp[s]->a1;

				//parse data
				for (int i = 9 ; i < like_tok.size() ; i ++) {
					chunk_field = sutils::tokenize(chunk_tok[i], ":");

					vector < string > chunk_prob = sutils::tokenize(chunk_field[chunk_idx_gp], ",");

					double p0 = atof(chunk_prob[0].c_str());
					double p1 = atof(chunk_prob[1].c_str());
					double p2 = atof(chunk_prob[2].c_str());
					
					char str[80];					
					int g = -1;
					threshold = 0.995;
					
//					printf("%f    %f",p2,threshold);
//puts(str);
					if (p0 >= threshold) g = 0;
					if (p1 >= threshold && !males[i-9]) g = 1;
					if (p2 >= threshold) g = 2;


					if (g >= 0) {
						switch (g) {
						case 0 : 	fdg << " 1 0 0"; break;
						case 1 : 	fdg << " 0 1 0"; break;
						case 2 : 	fdg << " 0 0 1"; break;
						}
						geno_post ++;
					} else if (like_tok[i] == ".") {
						fdg << " 0.333 0.333 0.333";
						geno_like++;

						
					} else {
						like_field = sutils::tokenize(like_tok[i], ":");

						vector < string > like_prob = sutils::tokenize(like_field[like_idx_gl], ",");          
//sprintf(str,"test2");
//cout << like_field[like_idx_gl] << "\n";
				
//puts(str);
						double l0 = pow(10, atof(like_prob[0].c_str()));

						double l1 = pow(10, atof(like_prob[1].c_str()));

						double l2 = pow(10, atof(like_prob[2].c_str()));
//printf("%f   $f",l0,l2);
//puts(str);
						double sum = l0 + l1 + l2;
//	                                        sprintf(str,"test6");
//puts(str);
					fdg << " " << sutils::double2str(l0/sum, 5) << " " << sutils::double2str(l1/sum, 5) << " " << sutils::double2str(l2/sum, 5);

//                                        sprintf(str,"test7");
//puts(str);
						geno_like++;
					}
				}

				fdg << endl;
				s++;

			}
			if (s % 100 == 0) { cout << "\r" << s << " / " << M.size(); cout.flush(); }
		}
	}
	cout << "\r  * #fixd = " << geno_post << " (" << sutils::double2str(geno_post * 100.0 / (geno_post + geno_like), 2) << "%)" << endl;
	cout << "  * #miss = " << geno_like << " (" << sutils::double2str(geno_like * 100.0 / (geno_post + geno_like), 2) << "%)" << endl;


	fdl.close();
	fdg.close();
	fdc.close();
}

