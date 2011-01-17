// -*- C++ -*-
/**\class electronCandidateFilter electronCandidateFilter.cc EWKSoftware/EDMTupleSkimmerFilter/src/electronCandidateFilter.cc

 Description: <one line class summary>

 Implementation:

*/

#ifndef electronCandidateFilter_H
#define electronCandidateFilter_H
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
#include <vector>
#include <iostream>
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "DataFormats/Math/interface/deltaR.h"

//
#include "TString.h"
#include "TMath.h"

//
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"

//
// class declaration
//

class ElectronCandidateFilter : public edm::EDFilter {
public:
  explicit ElectronCandidateFilter(const edm::ParameterSet&);
  ~ElectronCandidateFilter();
  
  static float scEt(const pat::Electron& ele) { return ele.caloEnergy()/TMath::CosH(ele.gsfTrack()->eta()); } 
  
  struct GreaterBySCEt {
    bool operator()( const pat::Electron& t1, const pat::Electron& t2 ) const 
    {
      float et1=ElectronCandidateFilter::scEt(t1);
      float et2=ElectronCandidateFilter::scEt(t2);
      return et1 > et2;
    }
  };

private:
  virtual Bool_t filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  Bool_t isInFiducial(Double_t eta);
      // ----------member data ---------------------------
  GreaterBySCEt scEtComparator_;

  Double_t ETCut_;

  edm::InputTag triggerEvent_;
  std::string hltpath_;

  edm::InputTag electronCollection_;

  edm::InputTag ebRecHits_;
  edm::InputTag eeRecHits_;
  //  edm::InputTag PrimaryVerticesCollection_;
  Double_t BarrelMaxEta_;
  Double_t EndCapMaxEta_;
  Double_t EndCapMinEta_;

  Bool_t useEcalDrivenElectrons_;
  Bool_t useSpikeRejection_;
  Double_t spikeCleaningSwissCrossCut_;


};
#endif

//
// constructors and destructor
//
ElectronCandidateFilter::ElectronCandidateFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // I N P U T      P A R A M E T E R S  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // Cuts
  ETCut_ = iConfig.getUntrackedParameter<double>("ETCut");

  useEcalDrivenElectrons_ = iConfig.getUntrackedParameter<Bool_t>("useEcalDrivenElectrons", false);

  Double_t BarrelMaxEta_D = 1.4442;
  Double_t EndCapMinEta_D = 1.566;
  Double_t EndCapMaxEta_D = 2.5;
  BarrelMaxEta_ = iConfig.getUntrackedParameter<double>("BarrelMaxEta", BarrelMaxEta_D);
  EndCapMaxEta_ = iConfig.getUntrackedParameter<double>("EndCapMaxEta", EndCapMaxEta_D);
  EndCapMinEta_ = iConfig.getUntrackedParameter<double>("EndCapMinEta", EndCapMinEta_D);
  // trigger related
  triggerEvent_=iConfig.getUntrackedParameter<edm::InputTag>("triggerEvent");
  hltpath_=iConfig.getUntrackedParameter<std::string>("hltpath");

  electronCollection_=iConfig.getUntrackedParameter<edm::InputTag>("electronCollection");
  //
  ebRecHits_ =  iConfig.getUntrackedParameter<edm::InputTag>("ebRecHits");
  eeRecHits_ =  iConfig.getUntrackedParameter<edm::InputTag>("eeRecHits");

  // spike cleaning:
  useSpikeRejection_ = iConfig.getUntrackedParameter<Bool_t>("useSpikeRejection");
  if (useSpikeRejection_) {
    spikeCleaningSwissCrossCut_= iConfig.getUntrackedParameter<Double_t>("spikeCleaningSwissCrossCut");
  }

  // now print a summary with what exactly the filter does:
  std::cout << "ElectronCandidateFilter..." << std::endl;
  std::cout << "ElectronCandidateFilter: HLT Path   " << hltpath_ << std::endl;
  std::cout << "ElectronCandidateFilter: ET  > " << ETCut_ << std::endl;

  if (useEcalDrivenElectrons_) {
    std::cout << "ElectronCandidateFilter: Electron Candidate is required to"
	      << " be ecal driven" << std::endl;
  }

  if (useSpikeRejection_) {
    std::cout << "ElectronCandidateFilter: Spike Cleaning will be done with the Swiss Cross Criterion"
	      << " cutting at " << spikeCleaningSwissCrossCut_ << std::endl;
  }

  std::cout << "ElectronCandidateFilter: Fiducial Cut: " << std::endl;
  std::cout << "ElectronCandidateFilter:    BarrelMax: "<<BarrelMaxEta_<<std::endl;
  std::cout << "ElectronCandidateFilter:    EndcapMin: " << EndCapMinEta_
	    << "  EndcapMax: " << EndCapMaxEta_
	    <<std::endl;

  produces< pat::ElectronCollection > 
    ("filteredPATElectronCandidates").setBranchAlias("filteredPATElectronCandidates");

}


ElectronCandidateFilter::~ElectronCandidateFilter()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ElectronCandidateFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace pat;


   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   // TRIGGER REQUIREMENT 
   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   // the event should pass an electron trigger

   edm::Handle< pat::TriggerEvent > triggerEvent;
   iEvent.getByLabel( triggerEvent_, triggerEvent );

   // ask for trigger accept; otherwise we don't even start
   if(!(triggerEvent->path(hltpath_)->wasRun() && triggerEvent->path(hltpath_)->wasAccept())){
     return false;
   }

   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   // Electron collection
   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   edm::Handle<pat::ElectronCollection> patElectron;
   iEvent.getByLabel(electronCollection_, patElectron);
   const pat::ElectronCollection *pElecs = patElectron.product();

   pat::ElectronCollection::const_iterator elec;
   // check how many electrons there are in the event
   if (   pElecs->size() == 0) {
     //std::cout << "No electrons found in this event" << std::endl;
     return false; // RETURN if no elecs in the event
   }

   //Using an auxiliary collection to manipulate electrons
   pat::ElectronCollection myElectrons;
   myElectrons.resize(pElecs->size());
   std::copy(pElecs->begin(),pElecs->end(),myElectrons.begin());
   std::sort(myElectrons.begin(), myElectrons.end(), scEtComparator_);

   
   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   // spike rejection only in EB. Already done in 39X but can be useful for other checks
   // for example on timing of central rechit or adding removal of double spikes (spikes contamination
   // is normally low for electrons due to track matching requirements)
   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

   if (useSpikeRejection_)  {
     for (pat::ElectronCollection::iterator iEle=myElectrons.begin();iEle!=myElectrons.end();++iEle)
       {
	 if (!iEle->isEB())
	   continue;
	 edm::Handle<EcalRecHitCollection> recHits;
	 iEvent.getByLabel(ebRecHits_, recHits);
	 const EcalRecHitCollection *myRecHits = recHits.product();     
	 const DetId seedId = iEle->superCluster()->seed()->seed();

	 EcalSeverityLevelAlgo severity;
	 Double_t swissCross = severity.swissCross(seedId, *myRecHits);

	 if (swissCross > spikeCleaningSwissCrossCut_) {
	   //removing spikes from electron collection
	   myElectrons.erase(iEle);
	 }
       }
   }

   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   // ET CUT: at least one electron in the event with ET>ETCut_-*-*-*-*-*
   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   if (scEt(*(myElectrons.begin())) < ETCut_) {
     return false; // RETURN if the highest ET elec has ET< ETcut
   }

   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   // demand that highest Et electron it is in fiducial:
   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   if (!isInFiducial(myElectrons.begin()->caloPosition().eta())) {
       return false; // RETURN highest ET electron is not in fiducial
   }
   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   // demand also that it is ecal driven
   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   if (useEcalDrivenElectrons_) {
     if (!myElectrons.begin()->ecalDrivenSeed()) {
       return false; // RETURN highest ET electron is not ecal driven
     }
   }

   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   // Demand also that leading electron is trigger matched
   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   if (myElectrons.begin()->triggerObjectMatchesByPath(hltpath_).size()==0)
     return false; // RETURN highest ET electron is not trigger matched

   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   // Add some user informations to electrons. Mostly to show functionality... 
   // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   for (pat::ElectronCollection::iterator iEle=myElectrons.begin();iEle!=myElectrons.end();++iEle)
       {
	 int firstPXBHit = (int) iEle->gsfTrack()->hitPattern().hasValidHitInFirstPixelBarrel();
	 iEle->addUserInt("PassValidFirstPXBHit",firstPXBHit);

	 int numberOfInnerHits = (int) iEle->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
	 iEle->addUserInt("NumberOfExpectedMissingHits",numberOfInnerHits);

	 float matched_dr_distance=-1;
	 if (myElectrons.begin()->triggerObjectMatchByPath(hltpath_))
	   matched_dr_distance=deltaR((*iEle),(*(myElectrons.begin()->triggerObjectMatchByPath(hltpath_))));
	 iEle->addUserFloat("HLTMatchingDR", matched_dr_distance);
       }

  

   //Put new collection in the event
   auto_ptr<pat::ElectronCollection> filteredElectronCandidates(new pat::ElectronCollection);
   filteredElectronCandidates->resize(myElectrons.size());
   std::copy(myElectrons.begin(),myElectrons.end(),filteredElectronCandidates->begin());

   iEvent.put( filteredElectronCandidates, "filteredPATElectronCandidates");

   return true;
}
   
// ------------ method called once each job just after ending the event loop  -
void 
ElectronCandidateFilter::endJob() {
}

Bool_t ElectronCandidateFilter::isInFiducial(Double_t eta)
{
  if (TMath::Abs(eta) < BarrelMaxEta_) return true;
  else if (TMath::Abs(eta) < EndCapMaxEta_ && TMath::Abs(eta) > EndCapMinEta_)
    return true;
  return false;

}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronCandidateFilter);
