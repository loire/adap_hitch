/*
 * main.cpp
 *
 *  Created on: 25 janv. 2013
 *      Author: etienne
 */
#include <iostream>
#include <fstream>
#include <string>
#include <boost/dynamic_bitset.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <omp.h>
#include <iomanip>
#include <libconfig.h++>
#include <iterator>
#include "simul_lib.h"

namespace po = boost::program_options;
using namespace std;
using namespace libconfig;
MTRand rnd;
ofstream fichierResults;
ofstream fichierFST;
ifstream fichierConfig;

ofstream global_results;

int main(int ac, char* av[])
{
	//////////////////////////////////////////////////////
	// Parsing command line options:
	//////////////////////////////////////////////////////
	string Config_filename;
	string Output_filename;
	cout << "hello" << endl;
	try
	{
		po::options_description desc("Allowed options");

		desc.add_options()
        				("help,h", "produce help message")
        				("configFile,c", po::value<string>()->default_value("parametres.txt"), "configFile")
						("goutFile,g",po::value<string>()->default_value("global_results.txt"), "goutFile")
						("fstFile,s", po::value<string>()->default_value("FST.txt"), "fstFile");
		po::variables_map vm;
		po::store(po::parse_command_line(ac , av, desc), vm);

		po::notify(vm);

		if (vm.count("help"))
		{
			cout << desc << "\n";
			return 1;
		}
//		if (vm.count("fstFile"))
//		{
//			cout << "fstFile was set to "
//					<< vm["fstFile"].as <string> () << "\n";
//			Config_filename = vm["fstFile"].as <string> ();
//		}
//		else
//		{
//			cout << "fstFile was not set, using default\n";
//			Config_filename = vm["fstFile"].as <string> ();
//		}
		if (vm.count("configFile"))
		{
			cout << "configFile was set to "
					<< vm["configFile"].as <string> () << "\n";
			Config_filename = vm["configFile"].as <string> ();
		}
		else
		{
			cout << "configFile was not set, using default\n";
			Config_filename = vm["configFile"].as <string> ();
		}

		///////////////// Set global outfile //////////////////
		if (vm.count("goutFile"))
		{
			cout << "goutFile was set to "
					<< vm["goutFile"].as <string> () << "\n";
			Output_filename = vm["goutFile"].as <string> ();
		}
		else
		{
			cout << "goutFile was not set, using default\n";
			Output_filename = vm["goutFile"].as <string> ();
		}


	}
	catch(exception& e)
	{
		cerr << "error: " << e.what() << "\n";
		return 1;
	}
	catch(...)
	{
		cerr << "Exception of unknown type!\n";
	}

	////////////////////////////////////////////////////////////
	//////////// Parsing Config file //////////////////////////
	///////////////////////////////////////////////////////////
	Config cfg;
	cfg.setAutoConvert(true);
	try
	{
		cfg.readFile(Config_filename.c_str());
	}
	catch(const FileIOException &fioex)
	{
		cout << "I/O error while reading file." << endl;
		return(1);
	}
	catch(const ParseException &pex)
	{
		std::cout << "Parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << std::endl;
		return(1);
	}
	////////////////////////////////////////////////////////////////////////////////////////
	///////////// Initalizing simulation parameters ////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////

//	 Trying to parse:
	try
	{
		cfg.lookup("metapop.n");			// Number of deme
		cfg.lookup("metapop.N");	// Size of deme
		cfg.lookup("metapop.loc");		// Postion of frontier between habitats
		cfg.lookup("metapop.m");		// Migration rate
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No metapop setting in configuration file." << endl;
		return(1);
	}

	try
	{
		cfg.lookup("chrom.structure_genetic");
		cfg.lookup("chrom.recombination_map");
		cfg.lookup("chrom.adap_s");
		cfg.lookup("chrom.adap_h");
		cfg.lookup("chrom.assort_strength");
		cfg.lookup("chrom.assort_h");
		cfg.lookup("chrom.traits_cost");
		cfg.lookup("chrom.traits_h");
		cfg.lookup("chrom.mutation_rate");
	}
	catch(const SettingNotFoundException &nfex)
	{
		cout << "No chrom setting in configuration file." << endl;
		return(1);
	}

	try
	{
		cfg.lookup("simul.NbGen_max");
		cfg.lookup("simul.pas");
		cfg.lookup("simul.mut");

	}
	catch(const SettingNotFoundException &nfex)
	{
		cout << "No simul setting in configuration file." << endl;
		return(1);
	}
	/////////////// Single Parameter values //////////////////
	int  n = cfg.lookup("metapop.n");			// Number of deme
	int loc = cfg.lookup("metapop.loc");		// Postion of frontier between habitats
	int NbGen = cfg.lookup("simul.NbGen_max");
	int pas = cfg.lookup("simul.pas");
	bool mut = cfg.lookup("simul.mut");

	/////////////// Array of parameter values ////////////////


	Setting &Genstruct = cfg.lookup("chrom.structure_genetic");
	Setting &Recmap = cfg.lookup("chrom.recombination_map");
	Setting &N_vals = cfg.lookup("metapop.N");
	Setting &m_vals = cfg.lookup("metapop.m");
	Setting &s_vals = cfg.lookup("chrom.adap_s");
	Setting &as_vals = cfg.lookup("chrom.assort_strength");
	Setting &as_h_vals = cfg.lookup("chrom.assort_h");
	Setting &ts_vals = cfg.lookup("chrom.traits_cost");
	Setting &ts_h_vals = cfg.lookup("chrom.traits_h");
	Setting &U_vals = cfg.lookup("chrom.mutation_rate");
	Setting &h_vals = cfg.lookup("chrom.adap_h");



	/////////////////////////

	global_results.open(Output_filename.c_str(),ios::out | ios::app);
//	global_results << "m\t s \t 2m/s \t q1" << endl;

	///////////////////////////////
	// don't loop on Parameters values yet
	//////////////////////////////
	int N = N_vals[0];			// Size of deme
	double m = m_vals[0];		// Migration rate
	double U = U_vals[0];		// Mutation rate

	double s = s_vals[0];		// adaptation locus strength
	double h = h_vals[0];		// dominance of adaptation locus
	double as = as_vals[0];		// Assortment strength
	double as_h = as_h_vals[0];
	double ts = ts_vals[0];
	double ts_h = ts_h_vals[0];


	////// First value for each parameter is choosen //////////


	/////////////////////////////////////////////////////////////
	/////////// Print paramaters values /////////////////////////
	/////////////////////////////////////////////////////////////
	cout << "Number of demes: " << n << endl;
	cout << "size of demes: " << N << endl;
	cout << "migration rate: " << m << endl;
	cout << "Max number of generation: " << NbGen << endl;
	cout << "frontier of habitat at deme: " << loc << endl;
	cout << "mutation ? " << mut << endl;
	if (mut) cout << "mutation rate: " << U << endl;
	cout << "pas :" << pas << endl;


	cout << "__________________________________" << endl;
	cout << "selective values of adaptative locus: " << s << endl;
	cout << "dominance of adapative locus: " << h << endl;
	cout << "Strength of assortment locus: " << as << endl;
	cout << "Dominance of assortment locus: " << as_h << endl;
	cout << "Traits cost: " << ts << endl;
	cout << "Dominance of traits cost:" << ts_h << endl;
	cout << "__________________________________" << endl;

	/////////////////////////////////////////////////////
	//////////// Declaring working variables ////////////
	////////////////////////////////////////////////////
//	int i,j;
//	int N_1=N-1;
//	int nbS = s_num+1;
//	int gen=0;
//	bool EqAdap = false;
//	bool eqNeutre = false;
//	boost::dynamic_bitset<> mask;
//	mask.resize(nbS,0);
//
//	cout << "number of site: " << nbS << endl;
//
//
//	///////////////////////////////// Fitness function ///////////////////////////////////////////
//	//////////////////////////////////////////////////////////////////////////////////////////////
//	////////////////////////// May differ according to theoretical values to be tested ///////////
//	//////////////////////////////////////////////////////////////////////////////////////////////
//	//////////////////////////////////////////////////////////////////////////////////////////////
//
//
////	double hom_00_0 = 1 - s;								// selective values of homozygotes  00 in deme 0
////	double hom_11_0 = 1 + s;								// selective values of homozygotes 11 in deme 0
////	double het_10_0 = hom_00_0 + (hom_11_0-hom_00_0)*h;		// selective values of heterozygotes in deme 0
////
////	double hom_00_1 = 1 + s;								// selective values of homozygotes  00 in deme 1
////	double hom_11_1 = 1 - s;								// selective values of homozygotes 11 in deme 1
////	double het_10_1 = hom_11_1 + (hom_00_1-hom_11_1)*h;		// selective values of heterozygotes in deme 1
//
//
//	double hom_00_0 = 1 - s;								// selective values of homozygotes  00 in deme 0
//	double hom_11_0 = 1;								// selective values of homozygotes 11 in deme 0
//	double het_10_0 = 1-h*s;		// selective values of heterozygotes in deme 0
//
//	double hom_00_1 = 1;								// selective values of homozygotes  00 in deme 1
//	double hom_11_1 = 1 - s;								// selective values of homozygotes 11 in deme 1
//	double het_10_1 = 1-h*s;		// selective values of heterozygotes in deme 1
//
//
//
//	cout << "selective values of homozygotes 00 in deme 0: " << hom_00_0 << endl;
//	cout << "selective values of homozygotes 11 in deme 0: " << hom_11_0 << endl;
//	cout << "selective values of heterozygotes in deme 0: " <<  het_10_0 << endl;
//	cout << endl;
//	cout << "selective values of homozygotes 00 in deme 1: " << hom_00_1 << endl;
//	cout << "selective values of homozygotes 11 in deme 1: " << hom_11_1 << endl;
//	cout << "selective values of heterozygotes in deme 1: " <<  het_10_1 << endl;
//	cout << endl;
//
//
//
//
//	cout.precision(4);
//	cout.setf(ios::fixed,ios::floatfield);
//	/////////////////////////////// Outputfilename //////////////////////
//	string ofilename("Results_n"+boost::lexical_cast<string>(n));
//	ofilename+="_N"+boost::lexical_cast<string>(N);
//	ofilename+="_nbS"+boost::lexical_cast<string>(nbS);
//	ofilename+="_m"+boost::lexical_cast<string>(m);
//	ofilename+="_s"+boost::lexical_cast<string>(s);
//	ofilename+="_h"+boost::lexical_cast<string>(h);
//	ofilename+="_ar"+boost::lexical_cast<string>(a_r);
//	ofilename+="_nr"+boost::lexical_cast<string>(n_r);
//	if (mut) ofilename+="_U"+boost::lexical_cast<string>(U);
//	ofilename+=".txt";
//
//	string FSTfilename = "FST_of_" + ofilename;
//
//	cout << "nom du fichier de sortie pour les frequence alleliques:" << ofilename << endl;
//	cout << "nom du fichier de sortie pour les FST : " << FSTfilename << endl;
//
//	fichierResults.open(ofilename.c_str(),ios::out);
//	fichierFST.open(FSTfilename.c_str(),ios::out);
//
//
//	////////////////////////////////////////////////////////////////
//	//////////////// Trying to allocate memory for population //////
//	////////////////////////////////////////////////////////////////
//
//	try
//	{
//		chr_diplo ** pop = new chr_diplo *[n];
//		chr_diplo ** temp = new chr_diplo *[n];
//		double ** Wij = new double *[n];
//		double ** Fik = new double *[n];
//		double ** TmpFik = new double *[n];
//		for ( i=0; i < n;i++)
//		{
//			pop[i] = new chr_diplo[N]; // ALlocate memory to store chr_diplo in each deme
//			temp[i] = new chr_diplo[N]; // copy of the other one
//			Wij[i] = new double[N]; // Allocate memory to store fitnesses in each deme
//			Fik[i] = new double[nbS];
//			TmpFik[i] = new double[nbS];
//		}
//
//		double * wbar = new double [n];
//		double * wmax = new double [n];
//		///////////////////////////////////////////////////////////////////////////
//		// Allocate memory for fitness of locus
//		///////////////////////////////////////////////////////////////////////////
//		for (i = 0; i < n; i++) // for each deme
//		{
//			for (j=0; j < N ; j++)	// for each chrom
//			{
//				pop[i][j].chr1.resize(nbS);			// create chrom (dynamic biteset) with nbSv sites set to 0
//				pop[i][j].chr2.resize(nbS);			// create chrom (dynamic biteset) with nbSv sites set to 0
//			}
//		}
//	}
//	catch(exception& e)
//	{
//		cout << "Memory problem, try smaller parameter values" << endl;
//	}
//
//	///////////////////////////////////////////////////////////////
//	//////// Now allocating memory for real //////////////////////
//	//////////////////////////////////////////////////////////////
//
//	chr_diplo ** pop = new chr_diplo *[n];
//	chr_diplo ** temp = new chr_diplo *[n];
//	double ** Wij = new double *[n];
//	double ** Fik = new double *[n];
//	double ** TmpFik = new double *[n];
//
//	for ( i=0; i < n;i++)
//	{
//		pop[i] = new chr_diplo[N]; // ALlocate memory to store chr_diplo in each deme
//		temp[i] = new chr_diplo[N]; // copy of the other one
//		Wij[i] = new double[N]; // Allocate memory to store fitnesses in each deme
//		Fik[i] =  new double[nbS];
//		TmpFik[i] = new double[nbS];
//	}
//
//	double * wbar = new double [n];
//	double * wmax = new double [n];
//
//	///////////////////////////////////////////////////////////////////////////
//	// Allocate memory chromosomes
//	///////////////////////////////////////////////////////////////////////////
//	for (i = 0; i < n; i++) // for each deme
//	{
//		for (j=0; j < N ; j++)	// for each chrom
//		{
//			pop[i][j].chr1.resize(nbS,0);			// create chrom (dynamic biteset) with nbSv sites set to 0
//			pop[i][j].chr2.resize(nbS,0);			// create chrom (dynamic biteset) with nbSv sites set to 0
//		}
//	}
//
//
//	double me,FstFN, FSTWeir;
//	double * FSTwright= new double [nbS];
//
//
//	//////////////////// Equilibrium values of FST for differents parameters values /////
//	//////////// First: 2 locus, 2 demes, no mutation, symetric migration, selection
//	double m_e = m*(n_r/(s+n_r));
//	double expectedFST_migration_selection_drift = 1 / (1+ 8 * (N * m_e) );
//	cout << "expectedFST_migration_selection_drift: "  << expectedFST_migration_selection_drift << endl;
//
//	cout.precision(10);
//	cout.setf(ios::fixed,ios::floatfield);
//	//////// Neutral locus with two allele in two deme.
//	double mu_carre = pow( (1 - 2 * m ),2) * pow( (1 - 2 * U ),2) ;
//	cout << mu_carre << endl;
//	cout << (2*N*mu_carre + mu_carre) << endl ;
//	cout << mu_carre / (2*N*mu_carre + mu_carre) << endl;
//	double expectedFST_mutation_migration_drift = mu_carre / (2*N*mu_carre + mu_carre);
//	cout << "expectedFST_migration_mutation_drift: " << expectedFST_mutation_migration_drift << endl;
//
//	double theovalue;
//	cout.precision(4);
//	cout.setf(ios::fixed,ios::floatfield);
//	///////////////////////////////////////////////
//	//// initialize population ////////////////////
//	/////////////////////////////////////////////
//
//	if (not mut)
//
//	{
//		cout << "intialize populations" << endl;
//		IntializePop(pop,Fik,s_num,N,n,loc);
//		theovalue = expectedFST_migration_selection_drift;
//	}
//	else
//		theovalue = expectedFST_mutation_migration_drift;
//
//	ComputeFreq(pop,Fik,TmpFik,n,N,nbS);
//	AfficheStoreFreq(fichierResults,Fik,n,nbS,gen);
//	AfficheStoreFST(fichierFST,FSTwright,nbS,gen,theovalue);
//	////////////////////////////////////////////////////////////////////////////////////////
//	///// Burning loop (Waiting for adaptation locus to reach equilibrium /////////////////
//	///////////////////////////////////////////////////////////////////////////////////////
//
//	gen=1;
////				cout << "Burning loop (Waiting for adaptation locus to reach equilibrium)..." << endl;
////				do
////				{
////
////
////					if (mut) mutation_all_site(pop, nbS,  n, N, U);
////			//		cout << "before fitness" << endl;
////					fitness(hom_00_0,hom_11_0,het_10_0,hom_00_1,hom_11_1,het_10_1, s_num, loc,n, N, mask, Wij, wbar,wmax,pop);
////					twoDemeMigration(n, N, N_1, m, a_r,n_r,s_num, nbS, Wij, wmax, temp, pop);
////					if (gen % pas == 0)
////					{
////						ComputeFreq(pop,Fik,TmpFik,n,N,nbS);
////				//		ComputeWrightFst(Fik,FSTwright,N,n,s_num);
////						AfficheStoreFreq(fichierResults,Fik,n,nbS,gen);
////			//			AfficheStoreFST(fichierFST,FSTwright,nbS,gen,theovalue);
////			//			AfficheFitness(pop,Wij,wbar,n,N);
////					//	CompareFreqNeutre(Fik,TmpFik,n,N,s_num,eqNeutre);
////					}
////					gen++;
////				} while((EqAdap == false) and (gen < NbGen));
////			//	} while(gen < NbGen);
////				cout << "Burning took " << gen << " generations" << endl;
//	//StoreEq1(global_results,Fik, n, nbS,s, m);
//
//
//	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	////////////////  Make a generation with all migrants from deme 1 to deme 2 fixed for the neutral mutation ////////
//	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////	cout << "Pulse chase (fixed neutral mutation in migrants" << endl;
////
////
////	fitness(hom_00_0,hom_11_0,het_10_0,hom_00_1,hom_11_1,het_10_1, s_num, loc,n, N, mask, Wij, wbar,wmax,pop);
////	twoDemeMigrationWithNeutralMig(n, N, N_1, m, a_r,n_r,s_num, nbS, Wij, wmax, temp, pop);
////	ComputeFreq(pop,Fik,TmpFik,n,N,nbS);
////	gen++;
////	AfficheStoreFreq(fichierResults,Fik,n,nbS,gen);
////	ComputeWrightFst(Fik,FSTwright,N,n,s_num);
//
//	///////////////////////// Fix neutral mutation in deme //////////////////////////////////////////////////////
////	if (not mut)
////	{
////		cout << "Fixing mutation in deme 0" << endl;
////		for (j=0;j<N;j++)
////		{
////			pop[0][j].chr1[s_num]=1;
////			pop[0][j].chr2[s_num]=1;
////		}
////	}
////
////
////	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
////	////////////// Second burning loop (Waiting for equilibrium in neutral mutation ////////////////////////////////
////	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////	cout << "second burning loop (waiting for neutral mutation to reach equilibrium)" << endl;
//	do{
//		if (mut) mutation_all_site(pop, nbS,  n, N, U);
//		fitness(hom_00_0,hom_11_0,het_10_0,hom_00_1,hom_11_1,het_10_1, s_num, loc,n, N, mask, Wij, wbar,wmax,pop);
////					AfficheStoreFreq(fichierResults,Fik,n,nbS,gen);
//		twoDemeMigration(n, N, N_1, m, a_r,n_r,s_num, nbS, Wij, wmax, temp, pop);
//		if (gen % pas == 0)
//		{
//			ComputeFreq(pop,Fik,TmpFik,n,N,nbS);
//			AfficheStoreFreq(fichierResults,Fik,n,nbS,gen);
//			CompareFreqNeutre(Fik,TmpFik,n,N,s_num,eqNeutre);
//			ComputeWeirFst(Fik,FSTWeir, N, n, s_num);
//			cout << "FSTWeir multilocus: " << setprecision(10) << FSTWeir << " Expected value: " << expectedFST_mutation_migration_drift <<  endl;
//		}
//		gen++;
////	} while((eqNeutre == false) and (gen < NbGen+NbGen));
//	}	while(gen < NbGen+NbGen);
//		cout << "equilibrium reached in " << gen << " generations" << endl;
//
//
//	///////////////////////////////////////////////////////////////////
//	////////////// Variables of results	///////////////////////////////
//	///////////////////////////////////////////////////////////////////
//
//
//	ComputeFreq(pop,Fik,TmpFik,n,N,nbS);
//	ComputeResultValues(Fik,N,nbS,n_r,me,FstFN);
//	ComputeWrightFst(Fik,FSTwright,N,n,s_num);
//	cout << "expectedFST_migration_selection_drift: "  << expectedFST_migration_selection_drift << endl;
//
//
////	cout << "Effective Migration Rate: " << me << endl;
////	cout << "Fst (feder & Nosil): " << FstFN << endl;
//
//	fichierResults.close();
//
//	///////////////////////////////////////////////////////////
//	///// End of program, de-allocate memory like a clean boy ///
//	//////////////////////////////////////////////////////////
//
//
//	for (i=0;i< n ; i++)
//	{
//		delete [] pop[i];
//		delete [] temp[i];
//		delete [] Wij[i];
//	}
//	delete [] pop;
//	delete [] temp;
//	delete [] Wij;
//	delete [] wbar;
//	delete [] wmax;


	return 0;
}
