#ifndef SIMUL_LIB_H_TEST
#define SIMUL_LIB_H_TEST

#include <vector>
#include <iostream>
#include <fstream>
#include <boost/dynamic_bitset.hpp>
#include <cmath>
#include "ranbin.h"
#include "MersenneTwister.h"
#include <omp.h>
using namespace std;
extern MTRand rnd;

struct chr_diplo
{
	boost::dynamic_bitset<> chr1; // selected loci
	boost::dynamic_bitset<> chr2; // selected loci
};


inline void IntializePop(chr_diplo ** &pop,double ** &Fik, int &s_num,int &N,int &n,int &loc)
{
	int i,j,k;

	for (i=0 ; i < loc ; i++)
	{
		for(j=0 ; j< N ; j++)
		{
			for (k=0; k < s_num+1; k++)
			{
				pop[i][j].chr1[k]=1;
				pop[i][j].chr2[k]=1;
			}
		}
	}
	for (i=loc ; i < n ; i++)
	{
		for(j=0 ; j< N ; j++)
		{
			for (k=0; k < s_num+1; k++)
			{
				pop[i][j].chr1[k]=0;
				pop[i][j].chr2[k]=0;
			}
		}
	}
	for(i=0; i < n ; i++)
	{
		for(k=0; k < (s_num+1) ; k++)
			Fik[i][k]=0.0;
	}
};

inline void ComputeFreq(chr_diplo ** &pop, double ** &Fik, double ** &TmpFik, int &n,int &N, int &nbS)
{
	int i,j,k;
	for (i=0; i< n ; i++)
	{
		for (k=0 ; k < nbS ; k++)
		{
			TmpFik[i][k]=Fik[i][k];
			Fik[i][k]=0.0;
			for (j=0; j < N; j++)
			{
				Fik[i][k] += pop[i][j].chr1[k];
				Fik[i][k] += pop[i][j].chr2[k];
			}
			Fik[i][k]/=2*N;
		}
	}
};

inline void StoreEq1(ofstream &global_results,double ** &Fik, int &n, int &nbS,double &s, double &m)
{


	for (int k=0;k<nbS-1;k++)
	{
		{
			cout <<  Fik[0][k] << "\t";
			global_results <<  m << "\t" << s << "\t" << 2*m/s << "\t";
			global_results << 1- Fik[0][k] << endl;
		}
	}

}
inline void StoreFST(ofstream &global_results,double ** &Fik, int &n, int &nbS,double &s, double &m, double n_r, int N)
{


			//	for (int k=0;k<nbS-1;k++)
			//	{
			//		{
			//			cout <<  Fik[0][k] << "\t";
			//			global_results <<  m << "\t" << s << "\t" << 2*m/s << "\t";
			//			global_results << 1- Fik[0][k] << endl;
			//		}
			//	}
				double me = m * n_r / (s + n_r);
				double FSTtheo =  1 / (1 + 8 * N * me);


				global_results << N << "\t" << s << "\t" << m << "\t" << n_r << "\t" << me << "\t"<< FSTtheo  << "\t" << Fik[0][1] << "\t" << Fik[1][1] << "\t";

}
inline void AfficheStoreFreq(ofstream &fichierResults,double ** &Fik, int &n, int &nbS,int &gen)
{
//	cout << "ca passe ou pas ?" << endl;
	cout << "gen: " << gen << "\t";

	for (int i=0;i<n;i++)
	{
		fichierResults << gen << "\t";
		fichierResults << i << "\t";
		for(int k=0;k<nbS;k++)
		{
			cout <<  Fik[i][k] << "\t";
			fichierResults << Fik[i][k] << "\t";
		}
		fichierResults << endl;
		cout << "|\t";
	}
	cout << endl;
};

inline void AfficheStoreFST(ofstream &fichierFST, double * &FSTwright, int &nbS, int &gen, double &theovalue)
{

//	cout << "FST at gen: " << gen << "\t|\t";
//	cout << "Theoretical value:" << theovalue << "\t|\t";
	fichierFST << gen << "\t";
	fichierFST << theovalue << "\t";

	for(int k = 0 ; k < nbS ; k++)
	{
		fichierFST << FSTwright[k] << "\t";
//		cout << FSTwright[k] << "\t";
	}
	fichierFST << endl;
//	cout << endl;
}



inline void AfficheFitness(chr_diplo ** &pop, double ** Wij,double * wbar, int &n, int &N)
{
	int i,j;
	for (i=0;i < n;i++)
	{
		cout << "deme: " << i << endl;
		for (j=0;j < N;j++)
		{
			cout << "ind:" << j  << " deme: " << i << endl;
			cout << "chr1:" << pop[i][j].chr1 << endl;
			cout << "chr2:" << pop[i][j].chr2 << endl;
			cout << "ind:" << j << " w:"  << Wij[i][j] << endl;
		}
		cout << "wbar: " << wbar[i] << endl;
	}
};

inline void CompareFreqAdap(double ** &Fik,double ** &TmpFik,int &n,int &N,int &s_num, bool &EqAdap)
{
	int twoN= 2*N;
	EqAdap= true;
	for (int i=0; i<n ; i++)
	{
		for (int k=0;k<s_num;k++)
		{
//			cout << fabs( Fik[i][k] - TmpFik[i][k] ) << endl;
//			cout << sqrt(( TmpFik[i][k] * (1 - TmpFik[i][k]) ) / twoN)  << endl;
			if ( fabs( Fik[i][k] - TmpFik[i][k] ) >  sqrt(( TmpFik[i][k] * (1 - TmpFik[i][k]) ) / twoN) )
				EqAdap=false;
		}
	}
};
inline void CompareFreqNeutre(double ** &Fik,double ** &TmpFik,int &n,int &N,int &s_num, bool &eqNeutre)
{
	int twoN= 2*N;
	eqNeutre= true;
	for (int i=0; i<n ; i++)
	{
		if ( fabs( Fik[i][s_num] - TmpFik[i][s_num] ) >  sqrt(( TmpFik[i][s_num] * (1 - TmpFik[i][s_num]) ) / twoN) )
			eqNeutre=false;
	}
};

inline void ComputeWeirFst(double ** Fik,double &FSTWeir,int &N,int &n,int &nbS)
{
    int k,i;
//	int nloci;
//	nloci = loci.size();


	double * pbar = new double [nbS];	// fréquence allèle moyenne (métapop)
	double ss,ss2,sumss,nc;
	double num,den;
	double SSI,SSP,MSI,MSP;

	for (i=0;i<nbS;i++)
		pbar[i] = 0.0;

	// Attention peut-être à caster les variables suivantes ?

	ss = 2 * N; 	// sample size per deme
	sumss = ss * n; // total sample size (across demes)
	ss2 = 4 * pow((double)N,2) * n;
    nc = (sumss - ss2 / sumss) / (n - 1.0);

    num = 0.0;
    den = 0.0;

	for (k = 0; k < nbS; k++) {

		for (i = 0; i < n; i++) {
	        pbar[k] += Fik[i][k];  // Fik -> fréquece dans dème i du locus s_num

		}
	    pbar[k] /= n;
		for (i = 0; i < n; i++) {

			SSI = ss * (Fik[i][k] - pow(Fik[i][k],2));

			SSP = ss * pow((Fik[i][k] - pbar[k]),2);
			MSI = SSI / (sumss - n);
			MSP = SSP / (n - 1.0);
	//		FST = (MSP - MSI)  / (MSP + (nc - 1) * MSI); // ca c'est un FST monolocus...
			num += (MSP - MSI);
			den += (MSP + (nc - 1) * MSI);
		}
    }
//	cout << "coucou2 " << num << endl;
	if (den > 0.0) {
		FSTWeir = num / den;
	} else {
		FSTWeir = 999.999; // ou bien indiquer FST = 999.999 ???
	}

//    if (sp!=0) FstWright = spi / sp;        // Division by 0 if fixed in both deme
//    else FstWright=1;
//    cout << "FstWright: " << FstWright << endl;
};



inline void rec(double &a_r,double &n_r, int &N_1, int &s_num, int &nbS, chr_diplo &res, chr_diplo &c1, chr_diplo &c2)
{
	int nbCo;
	vector<int> pos;
	int j;
	boost::dynamic_bitset<> rec;
	boost::dynamic_bitset<> off1;
	boost::dynamic_bitset<> off2;
	res.chr1.clear();
	res.chr2.clear();

	// First parent (c1)

	// number and positions of cross-overs
	//int nbCo = int(poisdev(param.Get_L()));
	//for (j = 0; j < nbCo; j++)
	//	pos.push_back(rnd.randInt(nS_2));
	//sort(pos.begin(), pos.end());
	nbCo=0;
	for (j=0; j< s_num; j++)
	{
		if ( rnd.rand() < a_r )
		{
			pos.push_back(j);
			nbCo++;
		}
	}
	if (rnd.rand() < n_r)
	{
		pos.push_back(s_num);
		nbCo++;
	}
	// recombination mask:

	for (j = 0; j < nbCo; j++)
		rec.resize(pos[j], (j % 2) == 0 ? 0 : 1);
	rec.resize(nbS, (nbCo % 2) == 0 ? 0 : 1);

	if (rnd.rand()>0.5) rec.flip();

	off1 = (c1.chr1 & rec);
	rec.flip();
	off2 = (c1.chr2 & rec);
	res.chr1 = (off1 | off2);

	// Second Parent (c2)

	// number and positions of cross-overs
	pos.clear();
	rec.clear();
	off1.clear();
	off2.clear();
	nbCo=0;
	for (j=0; j< s_num; j++)
	{
		if ( rnd.rand() < a_r )
		{
			pos.push_back(j);
			nbCo++;
		}
	}
	if (rnd.rand() < n_r )
	{
		pos.push_back(s_num);
		nbCo++;
	}
	// recombination mask:
	for (j = 0; j < nbCo; j++)
		rec.resize(pos[j], (j % 2) == 0 ? 0 : 1);
	rec.resize(nbS, (nbCo % 2) == 0 ? 0 : 1);
	if (rnd.rand()>0.5) rec.flip();
	off1 = (c2.chr1 & rec);
	rec.flip();
	off2 = (c2.chr2 & rec);

	res.chr2 = (off1 | off2);

};


inline void mutation_all_site(chr_diplo ** &pop, int nbS, int n, int N,double U)
{
	int i,j,k,mut;
	int nbS_1=nbS-1;

	for (i = 0; i < n; i++)
	{
		for (j=0; j< N ; j++)
		{
			mut = int(poisdev(U));
			for (k = 0; k < mut; k++)
			{
				int l = rnd.randInt(nbS_1);
				if (rnd.rand()>0.5)		// mutation fall on chr1
				{
					if (pop[i][j].chr1[l]==0)
						pop[i][j].chr1[l]=1;
					else
						pop[i][j].chr1[l]=0;
				}
				else					// mutation fall on chr2
				{
					if (pop[i][j].chr2[l]==0)
						pop[i][j].chr2[l]=1;
					else
						pop[i][j].chr2[l]=0;
				}
			}
		}
	}
}

inline void fitness(double &hom_00_0, double &hom_11_0,double &het_10_0,double &hom_00_1, double &hom_11_1,double &het_10_1, int &s_num, int &loc, int &n, int &N, boost::dynamic_bitset<> &mask, double ** &Wij,double * &wbar,double * &wmax,chr_diplo ** &pop)
{
	int i,j;
	double w;
	int cpt_hom_11,cpt_hom_00,cpt_het_10;


//	cout << "hello?" << endl;
	// loop for first habitat / deme : Allele 1 are favored
	for (i = 0; i < loc; i++)
	{
//		cout << "deme: " << i << " loc: " << loc << endl;

		wbar[i] = 0;
		wmax[i] = 0;
		for (j=0; j < N ; j++)
		{
//			cout << "indiv: " << j << endl;
//			cout << "chr1:" << pop[i][j].chr1 << endl;
//			cout << "chr2:" << pop[i][j].chr1 << endl;

			w = 1.0;
			mask = pop[i][j].chr1 | pop[i][j].chr2;			// count number of 11,10 and 01 genotypes in chromosomes

//			cout << "mask:" << mask << endl;

			mask.resize(s_num,0);							// exclude neutral marker

//			cout << "resized mask: " << mask << endl;

			mask.flip();									// count number of 00 genotypes
//			cout << "flipped mask: " << mask << endl;

			cpt_hom_00=mask.count();
//			cout << "count of 00:" << cpt << endl;

			if (cpt_hom_00!=0)
				w *= pow(hom_00_0,cpt_hom_00);						// multiplicative effect on fitness of 00 genotypes

			mask = pop[i][j].chr1 ^ pop[i][j].chr2;			// count number of 10 genotypes in chromosomes
			mask.resize(s_num,0);
			cpt_het_10=mask.count();
			if (cpt_het_10!=0)
				w *= pow(het_10_0,cpt_het_10);						// multiplicative effect on fitness of 11 genotypes

			cpt_hom_11 = s_num - cpt_hom_00 - cpt_het_10; // number of hom11 genotypes

			if (cpt_hom_11!=0)
				w *= pow(hom_11_0,cpt_hom_11);

			mask.clear();
			Wij[i][j] = w;									// store fitness of indiv
			wbar[i] += w;
			if (wmax[i] < w)
				wmax[i] = w;								// keep trace of maximal fitness in deme
		}
		wbar[i]= wbar[i] / N;								// mean fitness in deme
	}
	// Loop for habitat/deme 1, where alleles 00 are favored
	for (i = loc; i < n; i++)
	{
//		cout << "deme: " << i << " loc: " << param.Get_loc() << endl;

		wbar[i] = 0;
		wmax[i] = 0;
		for (j=0; j < N ; j++)
		{
			w = 1.0;
			mask = pop[i][j].chr1 & pop[i][j].chr2; // number of 11 genotypes
			mask.resize(s_num,0);					// exclude neutral marker of count
			cpt_hom_11 = mask.count();					// Number of 11 genotypes
			if (cpt_hom_11 != 0)
				w *= pow(hom_11_1,cpt_hom_11);				// multiplicative effect on fitness


			mask.clear();

			mask = pop[i][j].chr1 ^ pop[i][j].chr2;		// number of 1/0 and 0/1 genotypes
			mask.resize(s_num,0);
			cpt_het_10=mask.count();
			if (cpt_het_10 != 0)
				w *= pow(het_10_1,cpt_het_10);

			mask.clear();

			cpt_hom_00 = s_num - cpt_hom_11 - cpt_het_10;
			if (cpt_hom_00 != 0)
				w*= pow(hom_00_1,cpt_hom_00);

			Wij[i][j] = w;						// same as above, mean and max fitness calculation
			wbar[i] += w;
			if (wmax[i] < w)
				wmax[i] = w;
		}
		wbar[i]= wbar[i] / N;
	}

}


inline void twoDemeMigration(int &n, int &N,int &N_1, double &m, double &a_r,double &n_r, int &s_num, int &nbS, double ** &Wij,  double * &wmax, chr_diplo ** &temp, chr_diplo ** &pop)
{
	int i,j,ind;
	int par1,par2;
	// Deme 0:
	for (ind = 0;  ind < N ; ind++)
	{
		// Phillopatric

		if (rnd.rand()>m)
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[0][par1] / wmax[0]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[0][par2] / wmax[0]);
			// recombination
			rec(a_r, n_r, N_1, s_num, nbS, temp[0][ind], pop[0][par1], pop[0][par2]);
		}
		else  // migrants:
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[1][par1] / wmax[1]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[1][par2] / wmax[1]);
			// recombination:
			rec(a_r,n_r,N_1,s_num, nbS,temp[0][ind], pop[1][par1], pop[1][par2]);
		}
	}
	// Deme 1:
	for (ind = 0;  ind < N ; ind++)
	{

		// Phillopatric:
		if (rnd.rand()>m)
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[1][par1] / wmax[1]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[1][par2] / wmax[1]);
			// recombination:
			rec(a_r,n_r,N_1,s_num, nbS,temp[1][ind], pop[1][par1], pop[1][par2]);
		}
		else   // migrants:
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[0][par1] / wmax[0]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[0][par2] / wmax[0]);
			// recombination:
			rec(a_r,n_r,N_1,s_num, nbS,temp[1][ind], pop[0][par1], pop[0][par2]);
		}
	}
	for (i=0; i < n ; i++)
	{
		for (j=0; j< N ; j++)
		{
			pop[i][j]=temp[i][j];
		}
	}

};


inline void twoDemeMigrationWithNeutralMig(int &n, int &N,int &N_1, double &m, double &a_r,double &n_r, int &s_num, int &nbS, double ** &Wij,  double * &wmax, chr_diplo ** &temp, chr_diplo ** &pop)
{
	int i,j,ind;
	int par1,par2;
	// Deme 0:
	for (ind = 0;  ind < N ; ind++)
	{
		// Phillopatric

		if (rnd.rand()>m)
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[0][par1] / wmax[0]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[0][par2] / wmax[0]);
			// recombination
			rec(a_r, n_r, N_1, s_num, nbS, temp[0][ind], pop[0][par1], pop[0][par2]);
		}
		else  // migrants:
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[1][par1] / wmax[1]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[1][par2] / wmax[1]);
			// recombination:
			rec(a_r,n_r,N_1,s_num, nbS,temp[0][ind], pop[1][par1], pop[1][par2]);
		}
	}
	// Deme 1:

	for (ind = 0;  ind < N ; ind++)
	{

		// Phillopatric:
		if (rnd.rand()>m)
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[1][par1] / wmax[1]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[1][par2] / wmax[1]);
			// recombination:
			rec(a_r,n_r,N_1,s_num, nbS,temp[1][ind], pop[1][par1], pop[1][par2]);
		}
		else   // migrants:
		{
			do
			{
				par1 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[0][par1] / wmax[0]);

			do
			{
				par2 = rnd.randInt(N_1);
			}
			while (rnd.rand() > Wij[0][par2] / wmax[0]);
			// recombination:
			rec(a_r,n_r,N_1,s_num, nbS,temp[1][ind], pop[0][par1], pop[0][par2]);
			temp[1][ind].chr1[s_num]=1;
			temp[1][ind].chr2[s_num]=1;
		}
	}
	#pragma omp parallel for
	for (i=0; i < n ; i++)
	{
		for (j=0; j< N ; j++)
		{
			pop[i][j]=temp[i][j];
		}
	}

};


inline void ComputeResultValues(double ** &Fik, int N,int &nbS, double n_r, double &me,double &Fst)
{
	cout << "Frequency of neutral allele in deme 2: " << Fik[1][nbS-1] << endl;
	me = 2* Fik[1][nbS-1];
	double piTS = 1/(8*N*me);
	double piT = (1 - Fik[0][nbS-1]) + 0.5*N*n_r + piTS;
	Fst = piTS / piT ;

};
inline void ComputeWrightFst(double ** &Fik,double * &FstWright,int &N,int &n,int &s_num)
{
	int i,k;
	double p = 0;
	double spi = 0;
	double sp = 0;
	for (k=0 ; k < s_num+1; k++)
	{

		for (i=0;i<n;i++)
		{
			p+=Fik[i][k];
			spi+=pow(Fik[i][k],2);
		}
		p/=n;
		spi/=n;
		spi-=pow(p,2);
		sp = p*(1-p);
		if (sp!=0) FstWright[k] = spi / sp;		// Division by 0 if fixed in both deme
		else
			FstWright[k]=0;
//		cout << "FstWright: " << FstWright << endl;
		p=sp=spi=0;
	}



};

#endif
