# PAT example for electrons where we add custom electronIds to electrons, trigger matching,
# We also filter the pat electron collection and add user variables to pat::Electrons
# This runs from RECO or AOD and creates patTuples for further processing

from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *

## global tag for data
process.GlobalTag.globaltag = cms.string('START42_V17::All') 


removeAllPATObjectsBut(process, ['Electrons','METs','Jets'])

# add pf met
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

# Add PF jets
from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = False,
                 doBTagging   = False,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute'])),
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
process.patElectrons.addGenMatch = cms.bool(True)
process.patElectrons.embedGenMatch = cms.bool(True)
process.patElectrons.usePV = cms.bool(False)

#
process.patElectronIDs = cms.Sequence(process.simpleEleIdSequence)

# switch on PAT trigger for trigger matching and embedding
from PhysicsTools.PatAlgos.tools.trigTools import *
#REDIGI39X for 39X MC
#switchOnTrigger( process , 'patTrigger', 'patTriggerEvent', 'patDefaultSequence' )
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
    'file:///cmsrm/pc24_2/meridian/data/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_AODSIM_PU_S4_START42_V11-v1.root'
    ] );
process.source.fileNames = readFiles

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
        throw = cms.bool(False),
            HLTPaths = ["HLT_Ele*"]
        )

# The electron filter. We start from selectedPatElectrons with triggerMatch (no selection by default)
# and add additional selection and some user informations
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



# let it run. We write only events that pass the final electron filter
process.p = cms.Path(
    process.hltFilter
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
process.maxEvents.input = 10000
process.options.wantSummary = True


