// -*- C++ -*-
//
// Package:    trackPAT
// Class:      trackPAT
// 
//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  Fryday Anuary 20 2017        |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
//=================================================

// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/FWLite/interface/EventBase.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include "TLorentzVector.h"
#include <utility>
#include <string>

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

//
// class decleration
//

class trackPAT : public edm::EDAnalyzer {
public:
  explicit trackPAT(const edm::ParameterSet&);
  ~trackPAT();
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);
  int const getMuCat(reco::Muon const& muon) const;
  bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);

  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

 
    // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label;
  

  TTree*      tree_;
 

  // *************************************
  unsigned int             nB;
  std::vector<float>       *B_phi_px1_track, *B_phi_py1_track, *B_phi_pz1_track;
  std::vector<float>       *B_phi_px2_track, *B_phi_py2_track, *B_phi_pz2_track;

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
trackPAT::trackPAT(const edm::ParameterSet& iConfig)
  :
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
 
  tree_(0),
  
  nB(0),
  B_phi_px1_track(0), B_phi_py1_track(0), B_phi_pz1_track(0), 
  B_phi_px2_track(0), B_phi_py2_track(0), B_phi_pz2_track(0)
 
{
   //now do what ever initialization is needed
}

trackPAT::~trackPAT()
{

}

//
// member functions
//

// ------------ method called to for each event  ------------
void trackPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  //*********************************
  // Get event content information
  //*********************************  

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);
 
  for(View<pat::PackedCandidate>::const_iterator iTrack1 = thePATTrackHandle->begin();
      iTrack1 != thePATTrackHandle->end(); ++iTrack1 )
    {
      /*for (unsigned int i = 0, n = thePATTrackHandle.size(); i<n ; ++i )
	{
	const pat::PackedCandidate *iTrack1 =  (thePATTrackHandle)[i];*/
      
      
      if(iTrack1->charge()==0) continue;
      if(fabs(iTrack1->pdgId())!=211) continue;
      if(iTrack1->pt()<0.5) continue;// probaly the best is 0.95. see MiniAOD link

      
      for(View<pat::PackedCandidate>::const_iterator iTrack2 = iTrack1+1;
	  iTrack2 != thePATTrackHandle->end(); ++iTrack2 ) 
	{
	  
	  if(iTrack1==iTrack2) continue;
	  if(iTrack2->charge()==0) continue;
	  if(fabs(iTrack2->pdgId())!=211) continue;
	  //if(iTrack1->charge() == iTrack2->charge()) continue;
	  if(iTrack2->pt()<0.5) continue;

	  B_phi_px1_track->push_back(iTrack1->px());
	  B_phi_py1_track->push_back(iTrack1->py());
	  B_phi_pz1_track->push_back(iTrack1->pz());
	  
	  
	  B_phi_px2_track->push_back(iTrack2->px());
	  B_phi_py2_track->push_back(iTrack2->py());
	  B_phi_pz2_track->push_back(iTrack2->pz());

	   nB++;
	  	  
	}
    }
  
  
   //fill the tree and clear the vectors
   if (nB > 0 ) 
     {

       //std::cout << "filling tree" << endl;
       tree_->Fill();
     }
   // *********

   nB = 0; 

   B_phi_px1_track->clear(); B_phi_py1_track->clear(); B_phi_pz1_track->clear();
   B_phi_px2_track->clear(); B_phi_py2_track->clear(); B_phi_pz2_track->clear();
 
}


// ------------ method called once each job just before starting event loop  ------------

void 
trackPAT::beginJob()
{

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","track ntuple");

  tree_->Branch("nB",&nB,"nB/i");

  tree_->Branch("B_phi_px1_track", &B_phi_px1_track);
  tree_->Branch("B_phi_py1_track", &B_phi_py1_track);
  tree_->Branch("B_phi_pz1_track", &B_phi_pz1_track);
  
  tree_->Branch("B_phi_px2_track", &B_phi_px2_track);
  tree_->Branch("B_phi_py2_track", &B_phi_py2_track);
  tree_->Branch("B_phi_pz2_track", &B_phi_pz2_track);
  
}


// ------------ method called once each job just after ending the event loop  ------------
void trackPAT::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(trackPAT);

