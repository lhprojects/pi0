#include <vector>
#include <stdio.h>
#include "blossom5-v2.05/GEOM/GeomPerfectMatching.h"
#include "blossom5-v2.05/PerfectMatching.h"

#include <TTree.h>
#include <TFile.h>


// control printf
int verbose = 0;
// if something wrong, you got seen something
bool check_perfect_matching = true;
// mass of pi0
double const mpi0 = 0.135;

//a good parameter for A0.15 C0.01
double chi2_threshould = 100;
double startalpha = 10;

// a good parameter for A0.03 C0.01
//double chi2_threshould = 100000;
//double startalpha = 20;

//double chi2_threshould = 50;
//double startalpha = 10;

void OnEvent(int nGamma, double const *en, double const *px, double const *py, double const *pz,
 double a, double c, std::vector<std::pair<int,int> > &pairs) {

		// select photon
		std::vector<bool> gammas_maybe(nGamma);

		for (int g1 = 0; g1 < nGamma; g1++) {
			for (int g2 = g1+ 1; g2 < nGamma; g2++) {

				double const tpx = px[g1] + px[g2];
				double const tpy = py[g1] + py[g2];
				double const tpz = pz[g1] + pz[g2];
				double const ten = en[g1] + en[g2];

				if(en[g1] < 1E-3) continue;
				if(en[g2] < 1E-3) continue;
				double const e1 = std::max(en[g1], 1E-6);
				double const e2 = std::max(en[g2], 1E-6);

				double mass = sqrt(fabs(ten * ten - (tpx * tpx + tpy * tpy + tpz * tpz)));
				// dmass2 = 0.25 (dE1/E oplus dE2/E2)
				// dmass2 = 0.25 (c**2 E1**2 + c**2 E2**2 + a**2*E1 + a**2*E2)
				// dmass2 = 0.25 (c**2 (E1**2 + E2**2) + a**2(E1 + E2))
				double dmass = 0.5 * sqrt(c*c + a*a/e1 + c*c + a*a/e2)*mass;
				double chi2 =  (mass - mpi0)*(mass - mpi0)/(dmass*dmass);
				if(verbose > 0) {
					printf("%2d[%f] %2d[%f] %f %f %f\n", g1, en[g1], g2, en[g2],  mass, dmass, chi2);
				}
				if(chi2 < chi2_threshould) {				
					gammas_maybe[g1] = true;
					gammas_maybe[g2] = true;
				}
			}
		}

		std::vector<int> gammas_pi0_candidates;
		std::vector<int> gammas_single;
		for(int i = 0; i < (int)gammas_maybe.size(); ++i) {
			if(gammas_maybe[i]) {
				gammas_pi0_candidates.push_back(i);
			} else{
				gammas_single.push_back(i);
			}
		}

		if(verbose > 0) {
			printf("%2d never be candidate\n", (int)gammas_single.size());
			printf("%2d candidates\n", (int)gammas_pi0_candidates.size());
		}




		// use only candidates
		int n_candidates = gammas_pi0_candidates.size();
		PerfectMatching::REAL weights[100*100];
		int edges[2*100*100];
		int edge_num = 0;
		int node_num = 2*n_candidates;

		for (int g1 = 0; g1 < n_candidates; g1++) {
			for (int g2 = g1 + 1; g2 < n_candidates; g2++) {

				int g1_0 = gammas_pi0_candidates[g1];
				int g2_0 = gammas_pi0_candidates[g2];

				double tpx = px[g1_0] + px[g2_0];
				double tpy = py[g1_0] + py[g2_0];
				double tpz = pz[g1_0] + pz[g2_0];
				double ten = en[g1_0] + en[g2_0];

				double const e1 = std::max(en[g1_0], 0.);
				double const e2 = std::max(en[g2_0], 0.);

				double mass = sqrt(fabs(ten * ten - (tpx * tpx + tpy * tpy + tpz * tpz)));
				// dmass2 = 0.25 (dE1/E oplus dE2/E2)
				// dmass2 = 0.25 (c**2 E1**2 + c**2 E2**2 + a**2*E1 + a**2*E2)
				// dmass2 = 0.25 (c**2 (E1**2 + E2**2) + a**2(E1 + E2))
				double dmass = 0.5 * sqrt(c*c + a*a/e1 + c*c + a*a/e2)*mass;
				double chi2 =  (mass - mpi0)*(mass - mpi0)/(dmass*dmass);

				if(chi2 < chi2_threshould) {

					double chi2_wieght = chi2;

					edges[2*edge_num] = 2 * g1;
					edges[2*edge_num+1] = 2 * g2;
					weights[edge_num] = chi2_wieght;
					edge_num += 1;

					edges[2*edge_num] = 2 * g1 + 1;
					edges[2*edge_num+1] = 2 * g2 + 1;
					weights[edge_num] = chi2_wieght;
					edge_num += 1;
				}
			}
		}

		int const link_edge_start = edge_num;		
		double dalpha = 0.1;
		double endalpha = startalpha +1E-6;

		for(double alpha = startalpha; alpha < endalpha;  alpha += dalpha) {

			edge_num = link_edge_start;
			for(int g = 0; g < n_candidates; ++g) {
				edges[2*edge_num] = 2 * g;
				edges[2*edge_num+1] = 2 * g + 1;
				weights[edge_num] = 2*alpha;
				edge_num += 1;
			}

			struct PerfectMatching::Options options;
			options.verbose = 0;
			PerfectMatching* pm = new PerfectMatching(node_num, edge_num);

			for (int e = 0; e < edge_num; e++) {
				pm->AddEdge(edges[2*e], edges[2*e+1], weights[e]);
			}

			pm->options = options;


			pm->Solve();

			if (check_perfect_matching) {
				int res = CheckPerfectMatchingOptimality(node_num, edge_num, edges, weights, pm);
				if (res != 0) {
					if(verbose > 0) {
						printf("check optimality: res=%d (%s)\n", res, (res == 0) ? "ok" : ((res == 1) ? "error" : "fatal error"));
					}
				}
			}

			double cost = ComputePerfectMatchingCost(node_num, edge_num, edges, weights, pm);
			if(cost < 0) {
				printf("cost = %f\n", cost);
			}



			////////////////////////////////////////////
			// read out the result
			////////////////////////////////////////////

			if(verbose > 0) {
				for(int g = 0; g < (int)gammas_single.size(); ++g) {
					printf("single gamma %d (not a candidate)\n", gammas_single.at(g));
				}

				for(int g = 0; g < n_candidates; ++g) {
					bool v = pm->GetSolution(link_edge_start + g);
					if(v) {
						printf("single gamma %d\n", gammas_pi0_candidates.at(g));
					}
				}

				for(int e = 0; e < link_edge_start; e += 2) {
					bool v = pm->GetSolution(e);
					if(v) {
						int e1 = edges[2*e]/2;
						int e2 = edges[2*e+1]/2;
						int e1_0 = gammas_pi0_candidates.at(e1);
						int e2_0 = gammas_pi0_candidates.at(e2);
						printf("paired gamma %d %d\n", e1_0, e2_0);
					}
				}
			}



			for(int e = 0; e < link_edge_start; e += 2) {
				bool v = pm->GetSolution(e);
				if(v) {
					int e1 = edges[2*e]/2;
					int e2 = edges[2*e+1]/2;
					int e1_0 = gammas_pi0_candidates.at(e1);
					int e2_0 = gammas_pi0_candidates.at(e2);
					pairs.push_back(std::make_pair(e1_0, e2_0));
				}
			}
		
			std::sort(pairs.begin(), pairs.end());
			delete pm;
		}


}

int main(int argc, char* argv[])
{

	double a, c;

	char file[1000];
	if(argc < 4) {
		printf("pi0 a c rootfile\n");
		exit(1);
	} else {
		sscanf(argv[1], "%lf", &a);
		sscanf(argv[2], "%lf", &c);
		sscanf(argv[3], "%s", &file[0]);
	}

	TFile *f  = new TFile(file, "r");
	TTree * tree = (TTree*)f->Get("EvtTree");

	int Phhoton_n;
	std::vector<double> *Photon_Smear_E = NULL;
	std::vector<double> *Photon_Smear_Px = NULL;
	std::vector<double> *Photon_Smear_Py = NULL;
	std::vector<double> *Photon_Smear_Pz = NULL;
	
	int Pi0_n;
	std::vector<int> *Photon_Index = NULL;
	tree->SetBranchAddress("Photon_n",  &Phhoton_n);
	tree->SetBranchAddress("Photon_Smear_E",  &Photon_Smear_E);
	tree->SetBranchAddress("Photon_Smear_Px",  &Photon_Smear_Px);
	tree->SetBranchAddress("Photon_Smear_Py",  &Photon_Smear_Py);
	tree->SetBranchAddress("Photon_Smear_Pz",  &Photon_Smear_Pz);
	tree->SetBranchAddress("Pi0_n",  &Pi0_n);
	tree->SetBranchAddress("Photon_Index",  &Photon_Index);


	double correct_pair = 0;
	double total_pair = 0;
	double total_reconstructed_pair = 0;
	printf("%d\n", (int)tree->GetEntries());
    for (int evt = 0; evt < tree->GetEntries(); ++evt) {

		tree->GetEntry(evt);
		int const nGamma = Phhoton_n;

		if(Photon_Index->size() != 2 * Pi0_n) {
			printf("haah\n");
			exit(1);
		}

		double const *en = Photon_Smear_E->data();
		double const *px = Photon_Smear_Px->data();
		double const *py = Photon_Smear_Py->data();
		double const *pz = Photon_Smear_Pz->data();


		std::vector<std::pair<int,int> > pairs;
		OnEvent(nGamma, en, px, py, pz, a, c, pairs);



		///////////////////////////////
		/// print truth
		///////////////////////////////

		for(int g = 0; g < (int)pairs.size(); ++g) {
			for(int g2 = 0; g2 < Pi0_n; ++g2) {
				int a = pairs[g].first;
				int b = pairs[g].second;
				int c  = Photon_Index->at(2*g2);
				int d  = Photon_Index->at(2*g2 + 1);
				if((a == c && b == d) || (a == d && b == c)) {
					correct_pair += 1;
					break;
				}
			}
		}
		total_pair += Pi0_n;
		total_reconstructed_pair += pairs.size();

		if(verbose > 0) {
			for(int g = 0; g < Pi0_n; ++g) {
				printf(" truth: paired gamma %d %d\n", Photon_Index->at(2*g), Photon_Index->at(2*g+1));
			}
		}

		if(verbose) break;

	}

	printf("eff. %f\n", correct_pair/total_pair);
	printf("pur. %f\n", correct_pair/total_reconstructed_pair);
	printf("eff. x pur. %f\n", correct_pair/total_pair * correct_pair/total_reconstructed_pair);

    return 0;
}
