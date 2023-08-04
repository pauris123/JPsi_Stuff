
#include <memory>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <vector>
#include <algorithm>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "DataFormats/Common/interface/ValidHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/Common/interface/ValueMap.h"


#include "/afs/cern.ch/user/n/nstrautn/CMSSW_12_5_0/src/MiniAOD_ElectronNtuplizer/ElectronAnalyzer/interface/ElectronEffArea.h"
//#include "/interface/ElectronEffArea.h"

float Ele_Eff_Area(float SCeta_stuff);

float zero = 0.0;

class ElectronAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ElectronAnalyzer(const edm::ParameterSet&);
      ~ElectronAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::ElectronCollection> elecCollToken;

      edm::EDGetTokenT<reco::VertexCollection> PV_Token_;
      edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> SV_Token_;

      edm::InputTag elecSrc_;
      
      edm::EDGetTokenT<double> theRhoToken;
      //edm::EDGetTokenT<float> theRhoToken;
      edm::InputTag rhoSrc_;
      typedef std::vector<PileupSummaryInfo> PileupSummaryInfoCollection;
      edm::EDGetTokenT<PileupSummaryInfoCollection> pileupToken;
      edm::InputTag pileupSrc_;
      edm::InputTag TracksSrc_;

      edm::InputTag PV_Src_;
      edm::InputTag SV_Src_;

      edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
      edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
      edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
      edm::Handle<edm::ValueMap<bool> > wp80_iso_id_decisions;
      edm::Handle<edm::ValueMap<bool> > wp90_iso_id_decisions;
      edm::Handle<edm::ValueMap<bool> > wp80_noiso_id_decisions;
      edm::Handle<edm::ValueMap<bool> > wp90_noiso_id_decisions;

      edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > elewp80IsoIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > elewp90IsoIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > elewp80noIsoIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > elewp90noIsoIdMapToken_;

      TTree *electron_tree;
      std::vector<float> ele_pt,ele_eta,ele_phi,scl_eta,ele_oldsigmaietaieta,ele_oldsigmaiphiiphi,ele_oldcircularity,ele_oldr9,ele_scletawidth,ele_sclphiwidth,ele_he,ele_oldhe,
      ele_kfchi2,ele_gsfchi2,ele_fbrem,ele_ep,ele_eelepout,ele_IoEmIop,ele_deltaetain,ele_deltaphiin,ele_deltaetaseed, ele_vx, ele_vy, ele_vz,
      ele_psEoverEraw,ele_pfPhotonIso,ele_pfChargedHadIso,ele_pfNeutralHadIso,ele_PFPUIso,ElectronMVAEstimatorRun3Winter22IsoV1Values,ElectronMVAEstimatorRun3Winter22NoIsoV1Values,ElectronMVAEstimatorRun2Fall17IsoV2Values,ElectronMVAEstimatorRun2Fall17IsoV1Values,
      ElectronMVAEstimatorRun2Fall17NoIsoV1Values,ElectronMVAEstimatorRun2Fall17NoIsoV2Values,
      gen_pt, gen_eta, gen_phi, mother_pt, mother_eta, mother_phi, relISO_a, relISO_a_corr, ele_ip3D, ele_ip3D_dB, ele_ip3D_dB_error, ele_dz, ele_dzError, ele_dxy, ele_dxyError, track_pt, track_eta, track_phi, track_charge, track_dz, track_dxy;
      std::vector<int> ele_kfhits, ele_chi2_hits,ele_gsfhits, ele_expected_inner_hits,ele_charge,ele_mother,ele_gmother,ele_ggmother;
      float Diele_mass, rho, dR, relISO_0, relISO_1,relISO_0_corr, relISO_1_corr, Diele_pt;
      int numele, PFnumele, numTrack, ele_sameVertex, ele_ip3D_match;
      TLorentzVector P,P0,P1;
      std::vector<bool> ele_isPF, cutBasedElectronID_Winter22_122X_V1_veto, cutBasedElectronID_Winter22_122X_V1_loose, cutBasedElectronID_Winter22_122X_V1_medium, cutBasedElectronID_Winter22_122X_V1_tight,
      mvaEleID_Winter22_iso_V1_wp90, mvaEleID_Winter22_iso_V1_wp80, mvaEleID_Winter22_noIso_V1_wp90, mvaEleID_Winter22_noIso_V1_wp80,track_isEle,cutBasedElectronID_Fall17_94X_V2_veto, cutBasedElectronID_Fall17_94X_V2_loose, cutBasedElectronID_Fall17_94X_V2_medium, cutBasedElectronID_Fall17_94X_V2_tight,
      mvaEleID_Fall17_iso_V2_wp90, mvaEleID_Fall17_iso_V2_wp80, mvaEleID_Fall17_noIso_V2_wp90, mvaEleID_Fall17_noIso_V2_wp80,mvaEleID_Fall17_noIso_V2_wpLoose_unsopported,
      mvaEleID_Fall17_iso_V2_wpHZZ_unsopported;

};

ElectronAnalyzer::ElectronAnalyzer(const edm::ParameterSet& iConfig):
elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("elecSrc")),
rhoSrc_(iConfig.getUntrackedParameter<edm::InputTag>("rhoSrc")),
pileupSrc_(iConfig.getUntrackedParameter<edm::InputTag>("pileupSrc")),

PV_Src_(iConfig.getUntrackedParameter<edm::InputTag>("PV_Src")),
SV_Src_(iConfig.getUntrackedParameter<edm::InputTag>("SV_Src")),

eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLoose"))),
eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMedium"))),
eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTight"))),
elewp80IsoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIsoWp80"))),
elewp90IsoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIsoWp90"))),
elewp80noIsoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleNoIsoWp80"))),
elewp90noIsoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleNoIsoWp90")))

// eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("egmGsfElectronIDs:cutBasedElectronID-Winter22-122X-V1-loose"))),
// eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("egmGsfElectronIDs:cutBasedElectronID-Winter22-122X-V1-medium"))),
// eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("egmGsfElectronIDs:cutBasedElectronID-Winter22-122X-V2-tight"))),
// elewp80IsoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("egmGsfElectronIDs:mvaEleID-RunIIIWinter22-iso-V1-wp80"))),
// elewp90IsoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("egmGsfElectronIDs:mvaEleID-RunIIIWinter22-iso-V1-wp90"))),
// elewp80noIsoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("mvaEleID-RunIIIWinter22-noIso-V1-wp80"))),
// elewp90noIsoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("mvaEleID-RunIIIWinter22-noIso-V1-wp90")))

{
   elecCollToken = consumes<pat::ElectronCollection>(elecSrc_);

   PV_Token_ = consumes<reco::VertexCollection>(PV_Src_);
   SV_Token_ = consumes<reco::VertexCompositePtrCandidateCollection>(SV_Src_);

   theRhoToken = consumes<double>(rhoSrc_);
   pileupToken = consumes<PileupSummaryInfoCollection>(pileupSrc_);

   edm::Service<TFileService> fs;
   electron_tree = fs->make<TTree>("Events", "Events");

   electron_tree->Branch("numele",&numele);
   electron_tree->Branch("PFnumele",&PFnumele);
   electron_tree->Branch("ele_sameVertex",&ele_sameVertex);
   electron_tree->Branch("ele_ip3D_match",&ele_ip3D_match);
   electron_tree->Branch("ele_eta",&ele_eta);
   electron_tree->Branch("ele_phi",&ele_phi);
   electron_tree->Branch("ele_isPF",&ele_isPF);
   electron_tree->Branch("ele_vx",&ele_vx);
   electron_tree->Branch("ele_vy",&ele_vy);
   electron_tree->Branch("ele_vz",&ele_vz);
   electron_tree->Branch("Diele_mass",&Diele_mass);
   electron_tree->Branch("Diele_pt",&Diele_pt);

   electron_tree->Branch("ele_ip3D",&ele_ip3D);
   electron_tree->Branch("ele_ip3D_dB",&ele_ip3D_dB);
   electron_tree->Branch("ele_ip3D_dB_error",&ele_ip3D_dB_error);
   electron_tree->Branch("ele_dz",&ele_dz);
   electron_tree->Branch("ele_dzError",&ele_dzError);
   electron_tree->Branch("ele_dxy",&ele_dxy);
   electron_tree->Branch("ele_dxyError",&ele_dxyError);

   //BDT ID branches
   electron_tree->Branch("scl_eta",&scl_eta);
   electron_tree->Branch("ele_pt",&ele_pt);
   electron_tree->Branch("ele_oldsigmaietaieta",&ele_oldsigmaietaieta);
   electron_tree->Branch("ele_oldsigmaiphiiphi",&ele_oldsigmaiphiiphi);
   electron_tree->Branch("ele_oldcircularity",&ele_oldcircularity);
   electron_tree->Branch("ele_oldr9",&ele_oldr9 );
   electron_tree->Branch("ele_scletawidth",&ele_scletawidth);
   electron_tree->Branch("ele_sclphiwidth",&ele_sclphiwidth);
   electron_tree->Branch("ele_he",&ele_he);
   electron_tree->Branch("ele_oldhe",&ele_oldhe);
   electron_tree->Branch("ele_kfhits",&ele_kfhits);
   electron_tree->Branch("ele_kfchi2",&ele_kfchi2 );
   electron_tree->Branch("ele_gsfchi2",&ele_gsfchi2);
   electron_tree->Branch("ele_chi2_hits",&ele_chi2_hits);
   electron_tree->Branch("ele_fbrem",&ele_fbrem);
   electron_tree->Branch("ele_gsfhits",&ele_gsfhits);
   electron_tree->Branch("ele_expected_inner_hits",&ele_expected_inner_hits);
   electron_tree->Branch("ele_ep",&ele_ep);
   electron_tree->Branch("ele_eelepout",&ele_eelepout);
   electron_tree->Branch("ele_IoEmIop",&ele_IoEmIop);
   electron_tree->Branch("ele_deltaetain",&ele_deltaetain);
   electron_tree->Branch("ele_deltaphiin",&ele_deltaphiin);
   electron_tree->Branch("ele_deltaetaseed",&ele_deltaetaseed);
   electron_tree->Branch("ele_psEoverEraw",&ele_psEoverEraw);
   electron_tree->Branch("rho",&rho);
   electron_tree->Branch("ele_charge",&ele_charge);

   //isolation variables
   electron_tree->Branch("ele_pfPhotonIso",&ele_pfPhotonIso);
   electron_tree->Branch("ele_pfChargedHadIso",&ele_pfChargedHadIso);
   electron_tree->Branch("ele_pfNeutralHadIso",&ele_pfNeutralHadIso);
   electron_tree->Branch("ele_PFPUIso",&ele_PFPUIso);
   electron_tree->Branch("relISO_0",&relISO_0);
   electron_tree->Branch("relISO_1",&relISO_1);
   electron_tree->Branch("relISO_a",&relISO_a);
   electron_tree->Branch("relISO_0_corr",&relISO_0_corr);
   electron_tree->Branch("relISO_1_corr",&relISO_1_corr);
   electron_tree->Branch("relISO_a_corr",&relISO_a_corr);

   //Already implemented ID
   electron_tree->Branch("cutBasedElectronID_Winter22_122X_V1_veto",&cutBasedElectronID_Winter22_122X_V1_veto);
   electron_tree->Branch("cutBasedElectronID_Winter22_122X_V1_loose",&cutBasedElectronID_Winter22_122X_V1_loose);
   electron_tree->Branch("cutBasedElectronID_Winter22_122X_V1_medium",&cutBasedElectronID_Winter22_122X_V1_medium);
   electron_tree->Branch("cutBasedElectronID_Winter22_122X_V1_tight",&cutBasedElectronID_Winter22_122X_V1_tight);
   electron_tree->Branch("mvaEleID_Winter22_iso_V1_wp90",&mvaEleID_Winter22_iso_V1_wp90);
   electron_tree->Branch("mvaEleID_Winter22_iso_V1_wp80",&mvaEleID_Winter22_iso_V1_wp80);
   electron_tree->Branch("mvaEleID_Winter22_noIso_V1_wp90",&mvaEleID_Winter22_noIso_V1_wp90);
   electron_tree->Branch("mvaEleID_Winter22_noIso_V1_wp80",&mvaEleID_Winter22_noIso_V1_wp80);

   electron_tree->Branch("cutBasedElectronID_Fall17_94X_V2_veto",&cutBasedElectronID_Fall17_94X_V2_veto);
   electron_tree->Branch("cutBasedElectronID_Fall17_94X_V2_loose",&cutBasedElectronID_Fall17_94X_V2_loose);
   electron_tree->Branch("cutBasedElectronID_Fall17_94X_V2_medium",&cutBasedElectronID_Fall17_94X_V2_medium);
   electron_tree->Branch("cutBasedElectronID_Fall17_94X_V2_tight",&cutBasedElectronID_Fall17_94X_V2_tight);
   electron_tree->Branch("mvaEleID_Fall17_iso_V2_wp90",&mvaEleID_Fall17_iso_V2_wp90);
   electron_tree->Branch("mvaEleID_Fall17_iso_V2_wp80",&mvaEleID_Fall17_iso_V2_wp80);
   electron_tree->Branch("mvaEleID_Fall17_noIso_V2_wp90",&mvaEleID_Fall17_noIso_V2_wp90);
   electron_tree->Branch("mvaEleID_Fall17_noIso_V2_wp80",&mvaEleID_Fall17_noIso_V2_wp80);

   //V2 and V1 MVA scores
   electron_tree->Branch("ElectronMVAEstimatorRun3Winter22IsoV1Values",&ElectronMVAEstimatorRun3Winter22IsoV1Values);
   electron_tree->Branch("ElectronMVAEstimatorRun3Winter22NoIsoV1Values",&ElectronMVAEstimatorRun3Winter22NoIsoV1Values);

   electron_tree->Branch("ElectronMVAEstimatorRun2Fall17IsoV2Values",&ElectronMVAEstimatorRun2Fall17IsoV2Values);
   electron_tree->Branch("ElectronMVAEstimatorRun2Fall17IsoV1Values",&ElectronMVAEstimatorRun2Fall17IsoV1Values);
   electron_tree->Branch("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",&ElectronMVAEstimatorRun2Fall17NoIsoV2Values);
   electron_tree->Branch("ElectronMVAEstimatorRun2Fall17NoIsoV1Values",&ElectronMVAEstimatorRun2Fall17NoIsoV1Values);

   //unsopported IDs

   electron_tree->Branch("mvaEleID_Fall17_noIso_V2_wpLoose_unsopported",&mvaEleID_Fall17_noIso_V2_wpLoose_unsopported);
   electron_tree->Branch("mvaEleID_Fall17_iso_V2_wpHZZ_unsopported",&mvaEleID_Fall17_iso_V2_wpHZZ_unsopported);

   //Gen Info
   electron_tree->Branch("gen_pt",&gen_pt);
   electron_tree->Branch("gen_eta",&gen_eta);
   electron_tree->Branch("gen_phi",&gen_phi);
   electron_tree->Branch("mother_pt",&mother_pt);
   electron_tree->Branch("mother_eta",&mother_eta);
   electron_tree->Branch("mother_phi",&mother_phi);
   electron_tree->Branch("ele_mother",&ele_mother);

   electron_tree->Branch("ele_gmother",&ele_gmother);
   electron_tree->Branch("ele_ggmother",&ele_ggmother);

   //deltaR
   electron_tree->Branch("dR",&dR);
   
   
}


ElectronAnalyzer::~ElectronAnalyzer()
{


}

void
ElectronAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   edm::Handle< std::vector<pat::Electron>> electrons;
   iEvent.getByToken(elecCollToken, electrons);

   edm::Handle<double> rhoHandle;
   iEvent.getByToken(theRhoToken, rhoHandle);

   edm::Handle<PileupSummaryInfoCollection> pileupSummaryInfos;
   iEvent.getByToken(pileupToken, pileupSummaryInfos);

   edm::Handle< std::vector<reco::Vertex>> PV_verticies;
   iEvent.getByToken(PV_Token_, PV_verticies);

   edm::Handle< std::vector<reco::VertexCompositePtrCandidate>> SV_verticies;
   iEvent.getByToken(SV_Token_, SV_verticies);

   iEvent.getByToken(eleLooseIdMapToken_, loose_id_decisions);
   iEvent.getByToken(eleMediumIdMapToken_, medium_id_decisions);
   iEvent.getByToken(eleTightIdMapToken_ , tight_id_decisions);

   iEvent.getByToken(elewp80IsoIdMapToken_ , wp80_iso_id_decisions);
   iEvent.getByToken(elewp90IsoIdMapToken_ , wp90_iso_id_decisions);
   iEvent.getByToken(elewp80noIsoIdMapToken_ , wp80_noiso_id_decisions);
   iEvent.getByToken(elewp90noIsoIdMapToken_ , wp90_noiso_id_decisions);


   //edm::Handle<std::vector<reco::GsfTrack>> track_handle; // track stuff
   //iEvent.getByToken(tracksToken,track_handle);

   //auto track_handle = iEvent.getHandle(tracksToken); // Does not work in cmssw_10_x_x
   //reco::TrackCollection eleColl = *(track_handle.product());

   numTrack=0;
   numele=0;
   PFnumele=0;
   Diele_mass=0;
   Diele_pt=0;
   ele_sameVertex=0;
   ele_ip3D_match=0;
   rho=0;
   dR=0;
   relISO_0=0;
   relISO_1=0;
   relISO_0_corr=0;
   relISO_1_corr=0;
   relISO_a.clear();
   relISO_a_corr.clear();
   ele_eta.clear();
   ele_vx.clear();
   ele_vy.clear();
   ele_vz.clear();
   ele_phi.clear();
   ele_isPF.clear();
   scl_eta.clear();
   ele_pt.clear();
   ele_ip3D.clear();
   ele_ip3D_dB.clear();
   ele_ip3D_dB_error.clear();
   ele_dz.clear();
   ele_dzError.clear();
   ele_dxy.clear();
   ele_dxyError.clear();
   ele_mother.clear();
   ele_oldsigmaietaieta.clear();
   ele_oldsigmaiphiiphi.clear();
   ele_oldcircularity.clear();
   ele_oldr9.clear();
   ele_scletawidth.clear();
   ele_sclphiwidth.clear();
   ele_he.clear();
   ele_oldhe.clear();
   ele_kfhits.clear();
   ele_kfchi2.clear();
   ele_gsfchi2.clear();
   ele_chi2_hits.clear();
   ele_fbrem.clear();
   ele_gsfhits.clear();
   ele_expected_inner_hits.clear();
   ele_ep.clear();
   ele_eelepout.clear();
   ele_IoEmIop.clear();
   ele_deltaetain.clear();
   ele_deltaphiin.clear();
   ele_deltaetaseed.clear();
   ele_psEoverEraw.clear();
   ele_pfPhotonIso.clear();
   ele_pfChargedHadIso.clear();
   ele_pfNeutralHadIso.clear();
   ele_PFPUIso.clear();
   ele_charge.clear();
   cutBasedElectronID_Winter22_122X_V1_veto.clear();
   cutBasedElectronID_Winter22_122X_V1_loose.clear();
   cutBasedElectronID_Winter22_122X_V1_medium.clear();
   cutBasedElectronID_Winter22_122X_V1_tight.clear();
   mvaEleID_Winter22_iso_V1_wp80.clear();
   mvaEleID_Winter22_iso_V1_wp90.clear();
   mvaEleID_Winter22_noIso_V1_wp80.clear();
   mvaEleID_Winter22_noIso_V1_wp90.clear();
   ElectronMVAEstimatorRun3Winter22IsoV1Values.clear();
   ElectronMVAEstimatorRun3Winter22NoIsoV1Values.clear();

   cutBasedElectronID_Fall17_94X_V2_veto.clear();
   cutBasedElectronID_Fall17_94X_V2_loose.clear();
   cutBasedElectronID_Fall17_94X_V2_medium.clear();
   cutBasedElectronID_Fall17_94X_V2_tight.clear();
   mvaEleID_Fall17_iso_V2_wp80.clear();
   mvaEleID_Fall17_iso_V2_wp90.clear();
   mvaEleID_Fall17_noIso_V2_wp80.clear();
   mvaEleID_Fall17_noIso_V2_wp90.clear();
   ElectronMVAEstimatorRun2Fall17IsoV1Values.clear();
   ElectronMVAEstimatorRun2Fall17IsoV2Values.clear();
   ElectronMVAEstimatorRun2Fall17NoIsoV1Values.clear();
   ElectronMVAEstimatorRun2Fall17NoIsoV2Values.clear();
   mvaEleID_Fall17_noIso_V2_wpLoose_unsopported.clear();
   mvaEleID_Fall17_iso_V2_wpHZZ_unsopported.clear();
   gen_eta.clear();
   gen_phi.clear();
   gen_pt.clear();
   mother_eta.clear();
   mother_phi.clear();
   mother_pt.clear();

   ele_gmother.clear();
   ele_ggmother.clear();

   
   // Track stuff
   

   //for (PileupSummaryInfoCollection::const_iterator pileupSummaryInfo = pileupSummaryInfos->begin(); pileupSummaryInfo != pileupSummaryInfos->end(); ++pileupSummaryInfo)
    //{
     // if ( pileupSummaryInfo->getBunchCrossing() == 0 )
      //{
	      //std::cout<<pileupSummaryInfo->getTrueNumInteractions()<<std::endl;  //getTrueNumInteractions, getPU_NumInteractions
      //}
    //}

   for (auto it = electrons->cbegin(); it != electrons->cend(); ++it)
   {
     numele++;
     ;
     if(it->isPF())
     {
     pat::ElectronRef electronRef(electrons,numele);
     PFnumele++;
     ele_eta.push_back(it->eta());
     ele_phi.push_back(it->phi());
     ele_isPF.push_back(it->isPF());
     scl_eta.push_back(it->superCluster()->eta());
     ele_pt.push_back(it->pt());
     ele_vx.push_back((it->closestCtfTrackRef().isAvailable() && it->closestCtfTrackRef().isNonnull()) ? it->closestCtfTrackRef()->vx() : -99);
     ele_vy.push_back((it->closestCtfTrackRef().isAvailable() && it->closestCtfTrackRef().isNonnull()) ? it->closestCtfTrackRef()->vy() : -99);
     ele_vz.push_back((it->closestCtfTrackRef().isAvailable() && it->closestCtfTrackRef().isNonnull()) ? it->closestCtfTrackRef()->vz() : -99);
     ele_ip3D.push_back(it->ip3d());
     ele_ip3D_dB.push_back(it->dB(pat::Electron::PV3D)); // fabs(it->dB(pat::Electron::PV3D))
     ele_ip3D_dB_error.push_back(it->edB(pat::Electron::PV3D)); // fabs(it->dB(pat::Electron::PV3D))
     ele_dz.push_back(it->dB(pat::Electron::PVDZ)); // fabs(it->dB(pat::Electron::PVDZ))
     ele_dzError.push_back(it->edB(pat::Electron::PVDZ)); // fabs(it->edB(pat::Electron::PVDZ))
     ele_dxy.push_back(it->dB(pat::Electron::PV2D)); // fabs(it->dB(pat::Electron::PV2D))
     ele_dxyError.push_back(it->edB(pat::Electron::PV2D)); // fabs(it->edB(pat::Electron::PV2D))
     ele_oldsigmaietaieta.push_back(it->full5x5_sigmaIetaIeta());
     ele_oldsigmaiphiiphi.push_back(it->full5x5_sigmaIphiIphi());
     ele_oldcircularity.push_back(1.0-(it->full5x5_e1x5())/(it->full5x5_e5x5()));
     ele_oldr9.push_back(it->full5x5_r9());
     ele_scletawidth.push_back(it->superCluster()->etaWidth());
     ele_sclphiwidth.push_back(it->superCluster()->phiWidth());
     ele_he.push_back(it->hadronicOverEm());
     ele_oldhe.push_back(it->full5x5_hcalOverEcal());
     ele_kfhits.push_back((it->closestCtfTrackRef().isAvailable() && it->closestCtfTrackRef().isNonnull()) ? it->closestCtfTrackRef()->hitPattern().trackerLayersWithMeasurement() : -1);
     ele_kfchi2.push_back((it->closestCtfTrackRef().isAvailable() && it->closestCtfTrackRef().isNonnull()) ? it->closestCtfTrackRef()->normalizedChi2() : 0);
     ele_gsfchi2.push_back(it->gsfTrack()->normalizedChi2());
     ele_chi2_hits.push_back(it->gsfTrack()->normalizedChi2());
     ele_fbrem.push_back(it->fbrem());
     ele_gsfhits.push_back(it->gsfTrack()->hitPattern().trackerLayersWithMeasurement());
     ele_expected_inner_hits.push_back(it->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
     //ele_conversionVertexFitProbability.push_back(it->convVtxFitProb()) ;
     ele_ep.push_back(it->eSuperClusterOverP());
     ele_eelepout.push_back(it->eEleClusterOverPout());
     ele_IoEmIop.push_back(1.0/(it->ecalEnergy())-1.0/(it->trackMomentumAtVtx().R()));
     ele_deltaetain.push_back(it->deltaEtaSuperClusterTrackAtVtx());
     ele_deltaphiin.push_back(it->deltaPhiSuperClusterTrackAtVtx());
     ele_deltaetaseed.push_back(it->deltaEtaSeedClusterTrackAtCalo());
     ele_psEoverEraw.push_back((it->superCluster()->preshowerEnergy())/(it->superCluster()->rawEnergy()));
     ele_pfPhotonIso.push_back(it->pfIsolationVariables().sumPhotonEt);
     ele_pfChargedHadIso.push_back(it->pfIsolationVariables().sumChargedHadronPt);
     ele_pfNeutralHadIso.push_back(it->pfIsolationVariables().sumNeutralHadronEt);
     ele_PFPUIso.push_back(it->pfIsolationVariables().sumPUPt);
     ele_charge.push_back(it->charge());
    //  cutBasedElectronID_Winter22_122X_V1_veto.push_back(it->electronID("egmGsfElectronIDs:cutBasedElectronID-Winter22-122X-V1-veto"));
     cutBasedElectronID_Winter22_122X_V1_loose.push_back((*loose_id_decisions)[electronRef]);
     cutBasedElectronID_Winter22_122X_V1_medium.push_back((*medium_id_decisions)[electronRef]);
     cutBasedElectronID_Winter22_122X_V1_tight.push_back((*tight_id_decisions)[electronRef]);
     mvaEleID_Winter22_iso_V1_wp90.push_back((*wp80_iso_id_decisions)[electronRef]);
     mvaEleID_Winter22_iso_V1_wp80.push_back((*wp90_iso_id_decisions)[electronRef]);
     mvaEleID_Winter22_noIso_V1_wp90.push_back((*wp80_noiso_id_decisions)[electronRef]);
     mvaEleID_Winter22_noIso_V1_wp80.push_back((*wp90_noiso_id_decisions)[electronRef]);
     ElectronMVAEstimatorRun3Winter22IsoV1Values.push_back(0);
     ElectronMVAEstimatorRun3Winter22NoIsoV1Values.push_back(0);

     cutBasedElectronID_Fall17_94X_V2_veto.push_back(it->electronID("cutBasedElectronID-Fall17-94X-V2-veto"));
     cutBasedElectronID_Fall17_94X_V2_loose.push_back(it->electronID("cutBasedElectronID-Fall17-94X-V2-loose"));
     cutBasedElectronID_Fall17_94X_V2_medium.push_back(it->electronID("cutBasedElectronID-Fall17-94X-V2-medium"));
     cutBasedElectronID_Fall17_94X_V2_tight.push_back(it->electronID("cutBasedElectronID-Fall17-94X-V2-tight"));
     mvaEleID_Fall17_iso_V2_wp90.push_back(it->electronID("mvaEleID-Fall17-iso-V2-wp90"));
     mvaEleID_Fall17_iso_V2_wp80.push_back(it->electronID("mvaEleID-Fall17-iso-V2-wp80"));
     mvaEleID_Fall17_noIso_V2_wp90.push_back(it->electronID("mvaEleID-Fall17-noIso-V2-wp90"));
     mvaEleID_Fall17_noIso_V2_wp80.push_back(it->electronID("mvaEleID-Fall17-noIso-V2-wp80"));
     ElectronMVAEstimatorRun2Fall17IsoV1Values.push_back(it->userFloat("ElectronMVAEstimatorRun2Fall17IsoV1Values"));
     ElectronMVAEstimatorRun2Fall17IsoV2Values.push_back(it->userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values"));
     ElectronMVAEstimatorRun2Fall17NoIsoV1Values.push_back(it->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV1Values"));
     ElectronMVAEstimatorRun2Fall17NoIsoV2Values.push_back(it->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV2Values"));
     mvaEleID_Fall17_noIso_V2_wpLoose_unsopported.push_back(it->electronID("mvaEleID-Fall17-noIso-V2-wpLoose"));
     mvaEleID_Fall17_iso_V2_wpHZZ_unsopported.push_back(it->electronID("mvaEleID-Fall17-iso-V2-wpHZZ"));

     if ((it->genParticleRef ()).isNonnull ())
     {
       ele_mother.push_back(it->genParticle()->mother(0)->pdgId());
       gen_pt.push_back(it->genParticle()->pt());
       gen_eta.push_back(it->genParticle()->eta());
       gen_phi.push_back(it->genParticle()->phi());
       mother_pt.push_back(it->genParticle()->mother(0)->pt());
       mother_eta.push_back(it->genParticle()->mother(0)->eta());
       mother_phi.push_back(it->genParticle()->mother(0)->phi());

       auto *mom = it->genParticle()->mother();
	     int gmom=0;
	     int ggmom=0;
       if(mom->mother()!=0){
		    gmom=mom->mother()->pdgId();
		    if(mom->mother()->mother()!=0){
			   ggmom=mom->mother()->mother()->pdgId();
		    }
	     }
       
	     ele_gmother.push_back(gmom);
        ele_ggmother.push_back(ggmom);
     }
     }
   }

  // track stuff
  
   //for (auto it = TrackHandle->cbegin(); it != TrackHandle->cend(); ++it)

   if(PFnumele==2 && ele_charge[0]==(-1.0)*ele_charge[1])
   {
     P0.SetPtEtaPhiM(ele_pt[0],ele_eta[0],ele_phi[0],0.00051);
     P1.SetPtEtaPhiM(ele_pt[1],ele_eta[1],ele_phi[1],0.00051);
     P=P0+P1;
     Diele_mass=P.M();
     Diele_pt=P.Pt();
     rho=*(rhoHandle.product());
     dR=deltaR(ele_eta[0],ele_phi[0],ele_eta[1],ele_phi[1]);
     relISO_0 = (ele_pfPhotonIso[0]+ele_pfChargedHadIso[0]+ele_pfNeutralHadIso[0])/ele_pt[0];
     relISO_1 = (ele_pfPhotonIso[1]+ele_pfChargedHadIso[1]+ele_pfNeutralHadIso[1])/ele_pt[1];
     relISO_a = {relISO_0,relISO_1};
     
     relISO_0_corr = (ele_pfChargedHadIso[0]+std::max(zero,ele_pfPhotonIso[0]+ele_pfNeutralHadIso[0]-rho*Ele_Eff_Area(scl_eta[0])))/ele_pt[0];
     relISO_1_corr = (ele_pfChargedHadIso[1]+std::max(zero,ele_pfPhotonIso[1]+ele_pfNeutralHadIso[1]-rho*Ele_Eff_Area(scl_eta[1])))/ele_pt[1];
     relISO_a_corr = {relISO_0_corr,relISO_1_corr};

   
     if(ele_vx[0] == ele_vx[1] && ele_vy[0] == ele_vy[1] && ele_vz[0] == ele_vz[1]){
      ele_sameVertex = 1;
     }else{
      ele_sameVertex = 0;
     }

     if((ele_ip3D_dB[0] >= (ele_ip3D_dB[1]-ele_ip3D_dB_error[1]) && ele_ip3D_dB[0] <= (ele_ip3D_dB[1]+ele_ip3D_dB_error[1])) || (ele_ip3D_dB[1] >= (ele_ip3D_dB[0]-ele_ip3D_dB_error[0]) && ele_ip3D_dB[1] <= (ele_ip3D_dB[0]+ele_ip3D_dB_error[0]))){
      ele_ip3D_match = 1;
     }else{
      ele_ip3D_match = 0;
     }

     
     
     
     //if(dR < 0.3){

      //relISO_0_woEle = (ele_pfPhotonIso[0]+ele_pfChargedHadIso[0]+ele_pfNeutralHadIso[0]-ele_pt[1])/ele_pt[0];
      //relISO_1_woEle = (ele_pfPhotonIso[1]+ele_pfChargedHadIso[1]+ele_pfNeutralHadIso[1]-ele_pt[0])/ele_pt[1];
      //relISO_a_woEle = {relISO_0_woEle,relISO_1_woEle};

     //}
     //else{

      //relISO_0_woEle = (ele_pfPhotonIso[0]+ele_pfChargedHadIso[0]+ele_pfNeutralHadIso[0])/ele_pt[0];
      //relISO_1_woEle = (ele_pfPhotonIso[1]+ele_pfChargedHadIso[1]+ele_pfNeutralHadIso[1])/ele_pt[1];
      //relISO_a_woEle = {relISO_0_woEle,relISO_1_woEle};

     //}


     if(Diele_mass>2 && Diele_mass<4)
     electron_tree->Fill();
   }

}

void
ElectronAnalyzer::beginJob()
{
}

void
ElectronAnalyzer::endJob()
{
}

void
ElectronAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

DEFINE_FWK_MODULE(ElectronAnalyzer);

