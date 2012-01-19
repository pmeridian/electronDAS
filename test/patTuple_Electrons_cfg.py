# PAT example for electrons where we add custom electronIds to electrons, trigger matching,
# We also filter the pat electron collection and add user variables to pat::Electrons
# This runs from RECO or AOD and creates patTuples for further processing

from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *

## global tag for data
process.GlobalTag.globaltag = cms.string('GR_R_42_V23::All') 

# add pf met
from PhysicsTools.PatAlgos.tools.metTools import *
removeMCMatching(process, ['All'])
removeAllPATObjectsBut(process, ['Electrons','METs','Jets'])
addPfMET(process, 'PF')

# Add PF jets
from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = False,
                 doBTagging   = False,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])),
                 doType1MET   = False,
                 doL1Cleaning = False,
                 doL1Counters = False,
                 doJetID      = False
                 )

#for isolation correction
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJetsForRhoCorrection = process.kt6PFJets.clone(doRhoFastjet = True)
process.kt6PFJetsForRhoCorrection.Rho_EtaMax = cms.double(2.5)

##
process.load("Analysis.electronDAS.simpleEleIdSequence_cff")
##Customize electron producers
process.patElectrons.isoDeposits = cms.PSet()
process.patElectrons.userIsolation = cms.PSet()
process.patElectrons.addElectronID = cms.bool(True)
process.patElectrons.electronIDSources = cms.PSet(
    simpleEleId95relIso= cms.InputTag("simpleEleId95relIso"),
    simpleEleId90relIso= cms.InputTag("simpleEleId90relIso"),
    simpleEleId85relIso= cms.InputTag("simpleEleId85relIso"),
    simpleEleId80relIso= cms.InputTag("simpleEleId80relIso"),
    simpleEleId70relIso= cms.InputTag("simpleEleId70relIso"),
    simpleEleId60relIso= cms.InputTag("simpleEleId60relIso"),
    simpleEleId95cIso= cms.InputTag("simpleEleId95cIso"),
    simpleEleId90cIso= cms.InputTag("simpleEleId90cIso"),
    simpleEleId85cIso= cms.InputTag("simpleEleId85cIso"),
    simpleEleId80cIso= cms.InputTag("simpleEleId80cIso"),
    simpleEleId70cIso= cms.InputTag("simpleEleId70cIso"),
    simpleEleId60cIso= cms.InputTag("simpleEleId60cIso"),    
    )
process.patElectrons.addGenMatch = cms.bool(False)
process.patElectrons.embedGenMatch = cms.bool(False)
process.patElectrons.usePV = cms.bool(False)

#
process.patElectronIDs = cms.Sequence(process.simpleEleIdSequence)

# switch on PAT trigger for trigger matching and embedding
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )
process.patTrigger.addL1Algos = cms.bool( True )

# adding trigger matching
# firing trigger objects used in succeeding HLT path 'HLT_Ele17_SW_TighterEleIdIsol_L1R_v2'
process.electronTriggerMatchHLT = cms.EDProducer(
  "PATTriggerMatcherDRDPtLessByR"                 # match by DeltaR only, best match by DeltaR
  , src     = cms.InputTag( "selectedPatElectrons" )
  , matched = cms.InputTag( "patTrigger" )          # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
  , matchedCuts = cms.string( 'path( "HLT_Ele*" )' )          # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
  #, andOr                      = cms.bool( False )  # AND
  #, filterIdsEnum              = cms.vstring( 'TriggerElectron' ) # wildcard, overlaps with 'filterIds'
  #, filterIds                  = cms.vint32( 0 )    # wildcard, overlaps with 'filterIdsEnum'
  #, filterLabels               = cms.vstring( '*' ) # wildcard
  #, pathNames                  = cms.vstring(
  #    'HLT_Ele17_SW_TighterEleIdIsol_L1R_v3'
  #  )
  , pathL3FilterAccepted = cms.bool( True ) # select only trigger objects used in trigger decision
  , pathLastFilterAcceptedOnly = cms.bool( False )   
  , collectionTags             = cms.vstring( '*' ) # wildcard
  , maxDPtRel = cms.double( 0.5 )
  , maxDeltaR = cms.double( 0.5 )
  , resolveAmbiguities    = cms.bool( True )        # only one match per trigger object
  , resolveByMatchQuality = cms.bool( True )        # take best match found per reco object: by DeltaR here (s. above)
  )

#switchOnTriggerMatching(
#    process,
#    triggerMatchers = [ 'electronTriggerMatchHLT' ]
#    )
# embedding trigger matching in electron objects

switchOnTriggerMatchEmbedding(
    process,
    triggerMatchers = [ 'electronTriggerMatchHLT' ]
 )

# Select jets
process.selectedPatJets.cut = cms.string('pt > 30')
process.selectedPatJetsAK5PF.cut = cms.string('pt > 30')

# Add the files
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
    'file:///cmsrm/pc24_2/meridian/data/ZSkim_RAWRECO.root'
    ] );
process.source.fileNames = readFiles

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
        throw = cms.bool(False),
            HLTPaths = ["HLT_Ele*"]
        )

# the electron filter. We start from selectedPatElectrons with triggerMatch (no selection by default)
# and add additional selection and user informations
process.electronPATFilter = cms.EDFilter(
    'ElectronCandidateFilter',
    ### the input collections needed:
    electronCollection = cms.untracked.InputTag("selectedPatElectronsTriggerMatch"),
    ebRecHits = cms.untracked.InputTag("reducedEcalRecHitsEB"),
    eeRecHits = cms.untracked.InputTag("reducedEcalRecHitsEE"),
    ### here the preselection is applied
    # fiducial cuts:
    BarrelMaxEta = cms.untracked.double(1.4442),
    EndCapMinEta = cms.untracked.double(1.566),
    EndCapMaxEta = cms.untracked.double(2.5),
    # demand ecal driven electron:
    useEcalDrivenElectrons = cms.untracked.bool(True),
    ETCut = cms.untracked.double(20.)                                  
    )

###########################
# Z (tag & probe) filter
##########################

ELECTRON_ET_CUT_MIN = 20.0
TAG_ELECTRON_ET_CUT_MIN = 20.0
ELECTRON_CUTS = "(abs(superCluster.eta)<2.5) && (ecalEnergy*sin(superClusterPosition.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"
MASS_CUT_MIN = 40.
PRESCALE = 1

process.goodElectrons = cms.EDFilter("GsfElectronRefSelector",
                                     src = cms.InputTag( 'gsfElectrons' ),
                                     cut = cms.string( ELECTRON_CUTS )
                                     )

process.PassingWP90 = cms.EDFilter("GsfElectronRefSelector",
                                    src = cms.InputTag( 'gsfElectrons' ),
                                    cut = cms.string(
    str(ELECTRON_CUTS)
    + " && (gsfTrack.trackerExpectedHitsInner.numberOfHits<=1 && !(-0.02<convDist<0.02 && -0.02<convDcot<0.02))" #wrt std WP90 allowing 1 numberOfMissingExpectedHits
    + " && (ecalEnergy*sin(superClusterPosition.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"
    + " && ((isEB"
    + " && ( dr03TkSumPt/p4.Pt <0.12 && dr03EcalRecHitSumEt/p4.Pt < 0.09 && dr03HcalTowerSumEt/p4.Pt  < 0.1 )"
    + " && (sigmaIetaIeta<0.01)"
    + " && ( -0.8<deltaPhiSuperClusterTrackAtVtx<0.8 )"
    + " && ( -0.007<deltaEtaSuperClusterTrackAtVtx<0.007 )"
    + " && (hadronicOverEm<0.12)"
    +                                                   ")"
    + " || (isEE"
    + " && ( dr03TkSumPt/p4.Pt <0.07 && dr03EcalRecHitSumEt/p4.Pt < 0.07 && dr03HcalTowerSumEt/p4.Pt  < 0.07 )"
    + " && (sigmaIetaIeta<0.03)"
    + " && ( -0.7<deltaPhiSuperClusterTrackAtVtx<0.7 )"
    + " && ( -0.009<deltaEtaSuperClusterTrackAtVtx<0.009 )"
    + " && (hadronicOverEm<0.1) "
    + "))"
    )
      )

process.Zele_sequence = cms.Sequence(
    process.goodElectrons+
    process.PassingWP90
    )
process.tagGsf = cms.EDProducer("CandViewShallowCloneCombiner",
                                                                decay = cms.string("PassingWP90 goodElectrons"),
                                                                checkCharge = cms.bool(False),
                                                                cut   = cms.string("mass > " + str(MASS_CUT_MIN))
                                                                )
process.tagGsfCounter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("tagGsf"),
                                     minNumber = cms.uint32(1)
                                     )

process.tagGsfFilter = cms.Sequence(process.tagGsf * process.tagGsfCounter)
process.tagGsfSeq = cms.Sequence( process.Zele_sequence * process.tagGsfFilter )

process.prescaler = cms.EDFilter("Prescaler",
                                 prescaleFactor = cms.int32(PRESCALE),
                                 prescaleOffset = cms.int32(0)
                                 )

process.filter = cms.Sequence( process.prescaler * process.tagGsfSeq )

        

# let it run. We write only events that pass the final electron filter
process.p = cms.Path(
    process.hltFilter
    *process.filter
    *process.kt6PFJetsForRhoCorrection
    *process.patElectronIDs
    *process.patDefaultSequence
    *process.electronPATFilter
    )

# rename output file
process.out.fileName = 'electronsPATTuple.root'            ##  (e.g. 'myTuple.root')
process.out.outputCommands += [
        'keep *_electronPATFilter_*_*',
        'keep *_gsfElectronCores_*_*', #important needed to dereference up to the core which keep seeds information
        'keep *_offlineBeamSpot_*_*',
        'keep *_offlinePrimaryVertices*_*_*',
        'keep edmTriggerResults_TriggerResults*_*_*',
        'keep *_hltTriggerSummaryAOD_*_*',
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
        'keep edmConditionsIn*Block_conditionsInEdm_*_*',
        'keep *_kt6PFJetsForRhoCorrection_rho_*',
        'keep *_kt6PFJetsForRhoCorrection_sigma_*'
        ]

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# process all the events
process.maxEvents.input = 500
process.options.wantSummary = True


