# PAT example for electrons where we add custom electronIds to electrons, trigger matching,
# We also filter the pat electron collection and add user variables to pat::Electrons
# This runs from RECO or AOD and creates patTuples for further processing

from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

## global tag for data
process.GlobalTag.globaltag = cms.string('GR_R_39X_V5::All') 

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
, andOr                      = cms.bool( False )  # AND
, filterIdsEnum              = cms.vstring( 'TriggerElectron' ) # wildcard, overlaps with 'filterIds'
, filterIds                  = cms.vint32( 0 )    # wildcard, overlaps with 'filterIdsEnum'
, filterLabels               = cms.vstring( '*' ) # wildcard
, pathNames                  = cms.vstring(
    'HLT_Ele17_SW_TighterEleIdIsol_L1R_v2'
  )
, pathLastFilterAcceptedOnly = cms.bool( True )   # select only trigger objects used in last filters of succeeding paths
, collectionTags             = cms.vstring( '*' ) # wildcard
, maxDPtRel = cms.double( 0.5 )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True )        # take best match found per reco object: by DeltaR here (s. above)
)

switchOnTriggerMatching(
    process,
    triggerMatchers = [ 'electronTriggerMatchHLT' ]
    )
# embedding trigger matching in electron objects
switchOnTriggerMatchEmbedding(
    process,
    triggerMatchers = [ 'electronTriggerMatchHLT' ]
    )


# Select jets
process.selectedPatJets.cut = cms.string('pt > 10')
process.selectedPatJetsAK5PF.cut = cms.string('pt > 10')

# Add the files
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
    'file:///cmsrm/pc24_2/meridian/D6B89C71-4B12-E011-8F7D-001A92971B5E.root'
    ] );
process.source.fileNames = readFiles

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

# the electron filter. We start from selectedPatElectrons with triggerMatch (no selection by default)
# and add additional selection and user informations
process.electronPATFilter = cms.EDFilter(
    'ElectronCandidateFilter',
    ### the input collections needed:
    electronCollection = cms.untracked.InputTag("selectedPatElectronsTriggerMatch","","PAT"),
    triggerEvent = cms.untracked.InputTag("patTriggerEvent","","PAT"),
    hltpath = cms.untracked.string("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2"), 
    ebRecHits = cms.untracked.InputTag("reducedEcalRecHitsEB"),
    eeRecHits = cms.untracked.InputTag("reducedEcalRecHitsEE"),
    ### here the preselection is applied
    # fiducial cuts:
    BarrelMaxEta = cms.untracked.double(1.4442),
    EndCapMinEta = cms.untracked.double(1.566),
    EndCapMaxEta = cms.untracked.double(2.5),
    # demand ecal driven electron:
    useEcalDrivenElectrons = cms.untracked.bool(True),
    # demand offline spike cleaning with the Swiss Cross criterion:
    useSpikeRejection = cms.untracked.bool(False),
    spikeCleaningSwissCrossCut = cms.untracked.double(0.95),
    # demand geometrically matched to an HLT object 
    # ET Cut in the SC
    ETCut = cms.untracked.double(20.),                                  
    )

# let it run. We write only events that pass the final electron filter
process.p = cms.Path(
    process.patElectronIDs*
    process.patDefaultSequence*
    process.electronPATFilter
    )

# rename output file
process.out.fileName = 'electronsPATTuple.root'            ##  (e.g. 'myTuple.root')
process.out.outputCommands += [
        'keep *_electronPATFilter_*_*'
        ]
process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# process all the events
process.maxEvents.input = -1
process.options.wantSummary = True


