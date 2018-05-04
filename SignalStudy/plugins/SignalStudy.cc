// -*- C++ -*-
//
// Package:    Plots/SignalStudy
// Class:      SignalStudy
// 
/**\class SignalStudy SignalStudy.cc Plots/SignalStudy/plugins/SignalStudy.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  local user
//         Created:  Sun, 29 Apr 2018 10:05:57 GMT
//
//

#include <memory>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <map>
#include <sys/stat.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "IIHEAnalysis.h"

#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#endif
using namespace std;
const float m_el = 0.000511 ;
const float mu_mass = 0.10565837;

template <typename T>
struct SortByPt
{
	bool operator () (const T& a, const T& b) const {
		return a.first.Pt() > b.first.Pt();
	}

};

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SignalStudy : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit SignalStudy(const edm::ParameterSet&);
		~SignalStudy();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------
		string InputFile;
		IIHEAnalysis* tree;
		ifstream file_db1;
		char datafile[2000];
		edm::Service<TFileService> fs;

		TH1D *h_Fill_SE_Trigger, *h_Fill_SP_Trigger, *h_Fill_Ele115_Trigger, *h_Fill_Combo_Trigger, *h_Fill_Filters, *h_Fill_ElePtEta, *h_Fill_EleVID;
		TH1D *h_Fill_MuonIso; TH1D *h_Fill_MuonID, *h_Fill_MuPtEta;
		TH1D *h_Fill_TauPtEta, *h_Fill_DMF;

		TH1D *h_Fill_DMF_MuLoose;
		TH1D *h_Fill_DMF_MuTight;
		TH1D *h_Fill_DMF_MuTight_EleLoose;
		TH1D *h_Fill_DMF_MuTight_EleMedium;
		TH1D *h_Fill_DMF_MuTight_EleTight;
		TH1D *h_Fill_DMF_MuTight_EleLoose_LooseMVA;
		TH1D *h_Fill_DMF_MuTight_EleLoose_MediumMVA;
		TH1D *h_Fill_DMF_MuTight_EleLoose_TightMVA;

		TH1D *h_Fill_MuLoose;
		TH1D *h_Fill_MuTight;
		TH1D *h_Fill_MuTight_EleLoose;
		TH1D *h_Fill_MuTight_EleMedium;
		TH1D *h_Fill_MuTight_EleTight;

		TH1D *h_Fill_MuTight_EleLoose_LooseMVA;
		TH1D *h_Fill_MuTight_EleLoose_MediumMVA;
		TH1D *h_Fill_MuTight_EleLoose_TightMVA;
		bool MatchingToGenTaus(TLorentzVector tau, int &genindex);
		bool GenFilter();
		TH1D *h_Fill_GenTauFilter, *h_Fill_Events;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SignalStudy::SignalStudy(const edm::ParameterSet& iConfig)

{
	//now do what ever initialization is needed
	InputFile = iConfig.getParameter<string>("InputFile");


}


SignalStudy::~SignalStudy()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
SignalStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;


	while(!file_db1.eof()) {
		file_db1>>datafile;
		if(file_db1.eof()) break;
		TFile *file_in=TFile::Open(datafile);

		TTree* treePtr = (TTree*) file_in->Get("IIHEAnalysis");
		tree = new IIHEAnalysis (treePtr);


		bool isTrigger;
		for (int iEntry = 0; iEntry < tree->GetEntries(); iEntry++)
		{

			isTrigger=false;
			tree->GetEntry(iEntry);
			// trigger
			// MET filters + PV
			// e 
			// tau
			h_Fill_Events->Fill(1.);
			std::cout<<"== start event==="<<std::endl;
			std::cout<<"gen filter:"<< GenFilter() << std::endl;
			if(!(GenFilter())) continue;
			//	std::cout<<"== event end==="<<std::endl;


			//			h_Fill_GenTauFilter->Fill(1.);

			/*
			   bool isMuon = false;
			   for (unsigned int iGen = 0; iGen < tree->mc_px->size(); iGen++){

			   TLorentzVector gen_part, gen_part2;
			   gen_part.SetPxPyPzE(tree->mc_px->at(iGen),tree->mc_py->at(iGen),tree->mc_pz->at(iGen),tree->mc_energy->at(iGen));
			   isMuon = abs(tree->mc_pdgId->at(iGen))==13  ? true : false ;
			   std::cout<<"status:"<< tree->mc_status->at(iGen) << std::endl;
			   if(isMuon) {
			   std::cout<<"muon hai====="<< std::endl;
			// find mother
			if(tree->mc_mother_index->at(iGen).size() > 0) { for (unsigned int j = 0; j < tree->mc_mother_index->at(iGen).size() ; j++) { std::cout<<"mother pdgID"<< fabs(tree->mc_mother_pdgId->at(iGen).at(j)) << std::endl;  }
			}
			}

			}
			if(!isMuon) continue;
			std::cout<<"== event end==="<<std::endl;
			*/
			h_Fill_GenTauFilter->Fill(1.);

                        if(tree->trig_HLT_Ele27_WPTight_Gsf_accept)  { h_Fill_SE_Trigger->Fill(1.);}
                        if(tree->trig_HLT_Photon175_accept) { h_Fill_SP_Trigger->Fill(1.);}
                        if(tree->trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept) { h_Fill_Ele115_Trigger->Fill(1.);}

                        if( tree->trig_HLT_Ele27_WPTight_Gsf_accept || tree->trig_HLT_Photon175_accept || tree->trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept ) {isTrigger=true; h_Fill_Combo_Trigger->Fill(1);}

			if(!isTrigger) continue;
			if(!tree->trig_Flag_goodVertices_accept) continue;
			if(!tree->trig_Flag_globalTightHalo2016Filter_accept) continue;
			if(!tree->trig_Flag_HBHENoiseFilter_accept) continue;
			if(!tree->trig_Flag_HBHENoiseIsoFilter_accept) continue;
			if(!tree->trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept) continue;
			if(!tree->trig_Flag_BadPFMuonFilter_accept) continue;
			if(!tree->trig_Flag_BadChargedCandidateFilter_accept) continue;

			h_Fill_Filters->Fill(1);

			// electron 

			vector<pair<TLorentzVector, unsigned int>> ElePairs;
			ElePairs.clear();

			for (unsigned int dau1index = 0; dau1index < tree->gsf_et->size(); dau1index++){
				float ET1 = tree->gsf_caloEnergy->at(dau1index)*sin(2.*atan(exp(-1.*tree->gsf_eta->at(dau1index)))) ;
				if( ET1 < 0 ) continue;
				TLorentzVector DauEle1 ;
				DauEle1.SetPtEtaPhiM(ET1, tree->gsf_eta->at(dau1index), tree->gsf_phi->at(dau1index),m_el);
				pair<TLorentzVector, unsigned int> tmpPair;
				tmpPair.first = DauEle1; tmpPair.second = dau1index;
				ElePairs.push_back(tmpPair);
			}
			std::sort(ElePairs.begin(), ElePairs.end(),SortByPt<pair<TLorentzVector,unsigned int>>());


			////// cuts
			if( ElePairs.size() == 0) continue;
			int mu_count1, mu_count2, mu_count3;
			mu_count1 = mu_count2 = mu_count3 = 0;

			for(unsigned int k = 0; k < ElePairs.size() ; k++) {
				bool isMuPtEta, isHighPtMuon, isMuonIso;
				isMuPtEta = isHighPtMuon = isMuonIso= false;
				if(ElePairs.at(k).first.Pt() > 50.  && fabs(tree->gsf_sc_eta->at(ElePairs.at(k).second)) < 2.5 ) { isMuPtEta=true;}				

				if((tree->gsf_VIDHEEP7->at(ElePairs.at(k).second))) { isHighPtMuon=true; isMuonIso=true;}

				if(isMuPtEta) mu_count1++;
				if(isMuPtEta && isHighPtMuon) mu_count2++;
				if(isMuPtEta && isHighPtMuon && isMuonIso) mu_count3++;
			}

			if(mu_count1 > 0 ) h_Fill_MuPtEta->Fill(1);
			if(mu_count2 > 0 ) h_Fill_MuonID->Fill(1);
			if(mu_count3 > 0 ) h_Fill_MuonIso->Fill(1);

			if( ! (mu_count1 > 0 && mu_count2 > 0  && mu_count3 > 0)) continue;

			// DR separation with taus and picking highest pt tau
			vector<pair<TLorentzVector, unsigned int>> TauPairs;
			TauPairs.clear();
			//std::cout<<"matched gen index:"<< tree->mu_mc_bestDR->at(ElePairs.at(0).second) << std::endl;
			for (unsigned int dau2index = 0; dau2index < tree->tau_px->size(); dau2index++){

				if(tree->tau_pt->at(dau2index) <= 0. ) continue;

				TLorentzVector DauTau1;
				DauTau1.SetPxPyPzE(tree->tau_px->at(dau2index), tree->tau_py->at(dau2index), tree->tau_pz->at(dau2index),tree->tau_energy->at(dau2index));
				if(!(ElePairs.at(0).first.DeltaR(DauTau1) > 0.5)) continue;
				pair<TLorentzVector, unsigned int> tmpPair1;
				tmpPair1.first = DauTau1; tmpPair1.second = dau2index;
				TauPairs.push_back(tmpPair1);
			}
			std::sort(TauPairs.begin(), TauPairs.end(),SortByPt<pair<TLorentzVector,unsigned int>>());
			if( TauPairs.size() == 0) continue;


			int counter1, counter2, counter2a, counter3, counter3a, counter4, counter5, counter6, counter7, counter8, counter8a,  counter9, counter10, counter11, counter12, counter13, counter14, counter15;
			counter1 = counter2 = counter3 = counter4 = counter5 = counter6 = counter7 = counter8 = counter9 = counter10 = counter11 = counter12 = counter13 = counter14 = counter15 = 0;
			counter2a = counter3a = counter8a = 0;
			for(unsigned int j = 0 ; j < TauPairs.size() ; j++ ) { 

				bool ispteta;
				ispteta = false;
				bool isDMF, isEleLoose, isEleTight, isMuonLoose, isMuonTight, isEleMedium;
				isDMF = isEleLoose = isEleTight = isMuonLoose = isMuonTight = isEleMedium = false;

				bool isLooseMVA(false), isMediumMVA(false), isTightMVA(false);
				bool isLooseMVANew(false), isMediumMVANew(false), isTightMVANew(false);



				if((TauPairs.at(j).first.Pt() > 30. && fabs(TauPairs.at(j).first.Eta()) < 2.3)) { ispteta = true;}


				if( tree->tau_decayModeFinding->at( TauPairs.at(j).second) > 0.5 ) {isDMF = true;}
				// discrimination against e/mu
				if(tree->tau_againstMuonLoose3->at(TauPairs.at(j).second) > 0.5) { isMuonLoose = true; }
				if(tree->tau_againstMuonTight3->at(TauPairs.at(j).second) > 0.5) { isMuonTight = true; }  

				if((tree->tau_againstElectronTightMVA6->at(TauPairs.at(j).second) > 0.5)) {isEleTight = true;}
				if((tree->tau_againstElectronMediumMVA6->at(TauPairs.at(j).second) > 0.5)) {isEleMedium= true;}
				if((tree->tau_againstElectronVLooseMVA6->at(TauPairs.at(j).second) > 0.5)) {isEleLoose = true;}

				if(tree->tau_byLooseIsolationMVArun2v1DBoldDMwLT->at(TauPairs.at(j).second) > 0.5) { isLooseMVA = true;}
				if(tree->tau_byMediumIsolationMVArun2v1DBoldDMwLT->at(TauPairs.at(j).second) > 0.5) {isMediumMVA = true;}
				if(tree->tau_byTightIsolationMVArun2v1DBoldDMwLT->at(TauPairs.at(j).second) > 0.5) {isTightMVA = true;}

				if(tree->tau_byLooseIsolationMVArun2v1DBnewDMwLT->at(TauPairs.at(j).second) > 0.5) { isLooseMVANew = true;}
				if(tree->tau_byMediumIsolationMVArun2v1DBnewDMwLT->at(TauPairs.at(j).second) > 0.5) {isMediumMVANew = true;}
				if(tree->tau_byTightIsolationMVArun2v1DBnewDMwLT->at(TauPairs.at(j).second) > 0.5) {isTightMVANew = true;}


				if(ispteta) counter1++;

				if(ispteta && isDMF )  counter2++;


				if(ispteta && isDMF && isMuonLoose ) counter2a++;

				if(ispteta && isDMF && isMuonTight ) counter3++;

				if(ispteta && isDMF && isMuonTight && isEleLoose ) counter3a++;

				if(ispteta && isDMF && isMuonTight && isEleMedium ) counter4++;

				if(ispteta && isDMF && isMuonTight && isEleTight ) counter5++;

				if(ispteta && isDMF && isMuonTight && isEleLoose && isLooseMVA ) counter6++; 

				if(ispteta && isDMF && isMuonTight && isEleLoose && isMediumMVA ) counter7++;
				if(ispteta && isDMF && isMuonTight && isEleLoose && isTightMVA ) counter8++;


				if(ispteta  && isMuonLoose ) counter8a++;

				if(ispteta  && isMuonTight ) counter9++;


				if(ispteta && isMuonTight && isEleLoose ) counter10++;

				if(ispteta && isMuonTight && isEleMedium ) counter11++;

				if(ispteta && isMuonTight && isEleTight ) counter12++;

				if(ispteta && isMuonTight && isEleLoose && isLooseMVANew ) counter13++;

				if(ispteta && isMuonTight && isEleLoose && isMediumMVANew ) counter14++;
				if(ispteta && isMuonTight && isEleLoose && isTightMVANew ) counter15++;

			}

			if(counter1 > 0 ) {h_Fill_TauPtEta->Fill(1.); }
			if(counter2 > 0 ) { h_Fill_DMF->Fill(1.);  }
			if(counter2a > 0 ) { h_Fill_DMF_MuLoose->Fill(1.);  }
			if(counter3 > 0 ) { h_Fill_DMF_MuTight->Fill(1.); }
			if(counter3a > 0 ) {h_Fill_DMF_MuTight_EleLoose->Fill(1.); }
			if(counter4 > 0 ) {h_Fill_DMF_MuTight_EleMedium->Fill(1.); }
			if(counter5 > 0 ) {h_Fill_DMF_MuTight_EleTight->Fill(1.); }
			if(counter6 > 0 ) {h_Fill_DMF_MuTight_EleLoose_LooseMVA->Fill(1.); }
			if(counter7 > 0 ) {h_Fill_DMF_MuTight_EleLoose_MediumMVA->Fill(1.); }
			if(counter8 > 0 ) {h_Fill_DMF_MuTight_EleLoose_TightMVA->Fill(1.); }

			if(counter8a > 0 ) {h_Fill_MuLoose->Fill(1.); }

			if(counter9 > 0 ) {h_Fill_MuTight->Fill(1.); }
			if(counter10 > 0 ) {h_Fill_MuTight_EleLoose->Fill(1.); }
			if(counter11 > 0 ) {h_Fill_MuTight_EleMedium->Fill(1.); }
			if(counter12 > 0 ) {h_Fill_MuTight_EleTight->Fill(1.); }
			if(counter13 > 0 ) {h_Fill_MuTight_EleLoose_LooseMVA->Fill(1.); }
			if(counter14 > 0 ) {h_Fill_MuTight_EleLoose_MediumMVA->Fill(1.);  }
			if(counter15 > 0 ) { h_Fill_MuTight_EleLoose_TightMVA->Fill(1.);  }


		}

		file_in->Close();

	}
	file_db1.close();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
	void 
SignalStudy::beginJob()
{
	file_db1.open(InputFile);   
	h_Fill_GenTauFilter = fs->make<TH1D>("h_Fill_GenTauFilter","h_Fill_GenTauFilter",2,0,2);
	h_Fill_Events = fs->make<TH1D>("h_Fill_Events","h_Fill_Events", 2,0,2);                                     
	h_Fill_SE_Trigger = fs->make<TH1D>("h_Fill_SE_Trigger","h_Fill_SE_Trigger",2,0,2);
	h_Fill_SP_Trigger = fs->make<TH1D>("h_Fill_SP_Trigger","h_Fill_SP_Trigger",2,0,2);
	h_Fill_Ele115_Trigger = fs->make<TH1D>("h_Fill_Ele115_Trigger","h_Fill_Ele115_Trigger",2,0,2);
	h_Fill_Combo_Trigger  = fs->make<TH1D>("h_Fill_Combo_Trigger","h_Fill_Combo_Trigger",2,0,2);
	h_Fill_Filters = fs->make<TH1D>("h_Fill_Filters","h_Fill_Filters",2,0,2);
	h_Fill_MuPtEta = fs->make<TH1D>("h_Fill_MuPtEta","h_Fill_MuPtEta",2,0,2);
	h_Fill_MuonID = fs->make<TH1D>("h_Fill_MuonID","h_Fill_MuonID",2,0,2);
	h_Fill_MuonIso = fs->make<TH1D>("h_Fill_MuonIso","h_Fill_MuonIso",2,0,2);

	h_Fill_TauPtEta = fs->make<TH1D>("h_Fill_TauPtEta","h_Fill_TauPtEta",2,0,2);
	h_Fill_DMF  = fs->make<TH1D>("h_Fill_DMF","h_Fill_DMF",2,0,2);
	h_Fill_DMF_MuLoose = fs->make<TH1D>("h_Fill_DMF_MuLoose","h_Fill_DMF_MuLoose",2,0,2);
	h_Fill_DMF_MuTight = fs->make<TH1D>("h_Fill_DMF_MuTight","h_Fill_DMF_MuTight",2,0,2);
	h_Fill_DMF_MuTight_EleLoose = fs->make<TH1D>("h_Fill_DMF_MuTight_EleLoose","h_Fill_DMF_MuTight_EleLoose",2,0,2);
	h_Fill_DMF_MuTight_EleMedium = fs->make<TH1D>("h_Fill_DMF_MuTight_EleMedium","h_Fill_DMF_MuTight_EleMedium",2,0,2);
	h_Fill_DMF_MuTight_EleTight = fs->make<TH1D>("h_Fill_DMF_MuTight_EleTight","h_Fill_DMF_MuTight_EleTight",2,0,2);
	h_Fill_DMF_MuTight_EleLoose_LooseMVA = fs->make<TH1D>("h_Fill_DMF_MuTight_EleLoose_LooseMVA","h_Fill_DMF_MuTight_EleLoose_LooseMVA",2,0,2);
	h_Fill_DMF_MuTight_EleLoose_MediumMVA  = fs->make<TH1D>("h_Fill_DMF_MuTight_EleLoose_MediumMVA","h_Fill_DMF_MuTight_EleLoose_MediumMVA",2,0,2);
	h_Fill_DMF_MuTight_EleLoose_TightMVA = fs->make<TH1D>("h_Fill_DMF_MuTight_EleLoose_TightMVA","h_Fill_DMF_MuTight_EleLoose_TightMVA", 2,0,2);

	h_Fill_MuLoose = fs->make<TH1D>("h_Fill_MuLoose","h_Fill_MuLoose",2,0,2);
	h_Fill_MuTight = fs->make<TH1D>("h_Fill_MuTight","h_Fill_MuTight",2,0,2);
	h_Fill_MuTight_EleLoose = fs->make<TH1D>("h_Fill_MuTight_EleLoose","h_Fill_MuTight_EleLoose",2,0,2);
	h_Fill_MuTight_EleMedium = fs->make<TH1D>("h_Fill_MuTight_EleMedium","h_Fill_MuTight_EleMedium",2,0,2);
	h_Fill_MuTight_EleTight  = fs->make<TH1D>("h_Fill_MuTight_EleTight","h_Fill_MuTight_EleTight",2,0,2);

	h_Fill_MuTight_EleLoose_LooseMVA = fs->make<TH1D>("h_Fill_MuTight_EleLoose_LooseMVA","h_Fill_MuTight_EleLoose_LooseMVA",2,0,2);
	h_Fill_MuTight_EleLoose_MediumMVA = fs->make<TH1D>("h_Fill_MuTight_EleLoose_MediumMVA","h_Fill_MuTight_EleLoose_MediumMVA",2,0,2);
	h_Fill_MuTight_EleLoose_TightMVA = fs->make<TH1D>("h_Fill_MuTight_EleLoose_TightMVA","h_Fill_MuTight_EleLoose_TightMVA",2,0,2);
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
SignalStudy::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SignalStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

bool SignalStudy::MatchingToGenTaus(TLorentzVector tau, int &genindex){
	double DRact = 0.5;
	genindex=-1;
	bool ismatched=false;
	for (unsigned int iGen = 0; iGen < tree->mc_px->size(); iGen++){

		TLorentzVector gen_part, gen_part2;
		gen_part.SetPxPyPzE(tree->mc_px->at(iGen),tree->mc_py->at(iGen),tree->mc_pz->at(iGen),tree->mc_energy->at(iGen));
		bool isMuon = abs(tree->mc_pdgId->at(iGen))==15  ? true : false ;
		unsigned int  moth_ind = tree->mc_mother_index->at(iGen).at(0);
		if( isMuon) { // &&   moth_ind > 0 ) { 
			//      if(isSignalZ) { if ( fabs(tree->mc_pdgId->at(moth_ind)) < 32 ) continue;}
			// to check if it is hadronic decay

			bool ishadronicdecay(false);
			// find a neutrino

			int neutrino = 0;
			for (unsigned int iGen2 = 0; iGen2 < tree->mc_px->size(); iGen2++){

				gen_part2.SetPxPyPzE(tree->mc_px->at(iGen2),tree->mc_py->at(iGen2),tree->mc_pz->at(iGen2),tree->mc_energy->at(iGen2));
				if(fabs(tree->mc_pdgId->at(iGen2))==16 || fabs(tree->mc_pdgId->at(iGen2))==14 || fabs(tree->mc_pdgId->at(iGen2))== 12) {
					if((tree->mc_mother_index->at(iGen2).at(0)) > 0) {
						if(fabs(tree->mc_pdgId->at(tree->mc_mother_index->at(iGen2).at(0))) == 15 && int(tree->mc_mother_index->at(iGen2).at(0)) == int(iGen) ) {
							neutrino++;
						}
					}
				}
			}
			// its a hadronic decay as only one neutrino is involved in it
			if(neutrino == 1) {
				gen_part = gen_part - gen_part2; // subtracting neutrino 4 momentum

				if(tau.DeltaR(gen_part) < DRact) {

					DRact=tau.DeltaR(gen_part);
					genindex=iGen;
					ismatched=true;

				}

			}
		}
		}
		return ismatched;
	}

	bool SignalStudy::GenFilter(){
		bool hadtau_present=false;
		for (unsigned int iGen = 0; iGen < tree->mc_px->size(); iGen++){

			TLorentzVector gen_part, gen_part2;
			gen_part.SetPxPyPzE(tree->mc_px->at(iGen),tree->mc_py->at(iGen),tree->mc_pz->at(iGen),tree->mc_energy->at(iGen));
			bool isMuon = abs(tree->mc_pdgId->at(iGen))==15  ? true : false ;
			unsigned int  moth_ind = tree->mc_mother_index->at(iGen).at(0);
			if( isMuon) { // &&   moth_ind > 0 ) { 

				int neutrino = 0;
				for (unsigned int iGen2 = 0; iGen2 < tree->mc_px->size(); iGen2++){

					gen_part2.SetPxPyPzE(tree->mc_px->at(iGen2),tree->mc_py->at(iGen2),tree->mc_pz->at(iGen2),tree->mc_energy->at(iGen2));
					if(fabs(tree->mc_pdgId->at(iGen2))==16 || fabs(tree->mc_pdgId->at(iGen2))==14 || fabs(tree->mc_pdgId->at(iGen2))== 12) {
						if((tree->mc_mother_index->at(iGen2).at(0)) > 0) {
							if(fabs(tree->mc_pdgId->at(tree->mc_mother_index->at(iGen2).at(0))) == 15 && int(tree->mc_mother_index->at(iGen2).at(0)) == int(iGen) ) {
								neutrino++;
							}
						}
					}
				}
				std::cout<<"neutrino:"<<neutrino << std::endl;
				// its a hadronic decay as only one neutrino is involved in it
				if(neutrino == 1) { hadtau_present = true;}
			}
			}

			return hadtau_present;

		}


		//define this as a plug-in
		DEFINE_FWK_MODULE(SignalStudy);
