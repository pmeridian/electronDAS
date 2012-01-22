#! /usr/bin/env python
import sys
import math
import time
from array import array
from ROOT import *
print "loading FWLite libraries... ",
from DataFormats.FWLite import Events, Handle
print "done."

def printInfo(i,event):
    print "processing entry # " + str(i) + " from Run "+ str(event.eventAuxiliary().id().run()) + " lumi "+str(event.eventAuxiliary().id().luminosityBlock()) + " @ " + time.asctime(time.localtime(time.time()))

def checkFullElectronId(ele,cut):
    # Requiring the full ElectronId (value of float is 7)
    if not ele.ecalDrivenSeed():
        return false
    if cut == "":
        return true
    else:
        if ele.isElectronIDAvailable(cut):
            #level 7 is used since we want the full id 
            if fabs(ele.electronID(cut)-7)<0.1:
                return true
            else:
                return false
        else:
            #returning false if electronId not calculated
            return false

def checkElectronIdBit(ele,cut,bit):
    # Check each single bit entering in the final electronId value
    # Bit 0: ElectronId
    # Bit 1: Isolation
    # Bit 2: Conversion Rejection
    if not ele.ecalDrivenSeed():
        return false
    if cut == "":
        return true
    else:
        if ele.isElectronIDAvailable(cut):
            #bit 0 is used for electronId, bit1 for isolation, bit2 for conversionRejection
            if int(ele.electronID(cut)) & (1<<(bit)):
                return true
            else:
                return false
        else:
            return true
    
def main():
    # if you want to run in batch mode
    ROOT.gROOT.SetBatch()
    # maximum number of events. -1 run over all
    maxNevents = 30000
    rand= TRandom3(0)

    #collections 
    electronCollection = "electronPATFilter"
    electronCollectionInstance = "filteredPATElectronCandidates"
    metCollection ="patMETsPF"
    vertexCollection = "offlinePrimaryVerticesWithBS"
    rhoCollection = "kt6PFJetsForRhoCorrection"
    rhoCollectionInstance = "rho"
    
    #change this according to your needs
    outfilename = "electronDistributionsMC.root"
    outputroot = TFile( outfilename, "RECREATE")
        
#    prefixFnal = 'dcache:/pnfs/cms/WAX/11/store/user/meridian/electronDAS/'
#    prefixCern = 'rfio:/castor/cern.ch/user/m/meridian/electronDAS/'
    prefixPisa = '/gpfs/gpfsddn/srm/cms/store/user/cmsdas/2012/ElectronShortExercise/'
    prefixLocal = './'

    prefix = prefixPisa

    # Kinamatic cuts used in the analysis
    ptCut = 20.
#    metCut = 20.
#    mtCut = 40.

    # PAT ntuples with electronCollection
    files = [
        # DATA
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_153_1_af3.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_98_1_wVK.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_62_1_ZUG.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_97_1_TTG.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_87_1_f4r.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_57_1_xAL.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_49_1_OPG.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_95_1_G3u.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_91_1_zVW.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_61_1_Pyp.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_86_1_LrO.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_90_1_d6O.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_94_1_LXl.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_53_1_3jO.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_88_1_hLL.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_84_1_EQI.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_48_1_A1Y.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_5_1_0RS.root",
#        "DATA/DoubleElectron-Run2011B-ZElectronSkimOnTheFly-30Nov2011-v2/electronsPATTuple_68_1_j8Z.root",
#    
#       MC DY files
       "MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/electronsPATTuple_10_1_HJi.root",
       "MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/electronsPATTuple_1_1_ZaP.root",
       "MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/electronsPATTuple_2_1_UJp.root",
       "MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/electronsPATTuple_3_1_0e6.root",
       "MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/electronsPATTuple_4_1_6M3.root",
       "MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/electronsPATTuple_5_1_Zzn.root",
       "MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/electronsPATTuple_6_1_VwQ.root",
       "MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/electronsPATTuple_7_1_54n.root",
       "MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/electronsPATTuple_8_2_pzv.root",
       "MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/electronsPATTuple_9_1_jHf.root",
    
        ]
    
    fullpath_files = []
    
    for afile in files:
        fullpath_files.append( prefix+afile )

    print 'Processing files '
    for ifile in fullpath_files:
            print ifile
            
    events = Events ( fullpath_files , maxEvents=maxNevents )

    handleElectron = Handle ("vector<pat::Electron>")
    handleVertices = Handle ("vector<reco::Vertex>")
    handleRho = Handle ("double")
#    handleMET = Handle ("vector<pat::MET>")

    electronIdLevels= [ "","simpleEleId90cIso","simpleEleId80cIso","simpleEleId70cIso" ]

    # ****************************************************************
    # Booking histograms
    # ****************************************************************
    histogram = {}
    histogram["nEle"] = TH1F("nEle","nEle", 20, -0.5, 19.5)


    for id in electronIdLevels:
        # general  kinematics carachteristics at various level of ID requirements
        histogram["nVertices_" + id] = TH1F("nVertices_" + id, "nVertices", 50, -0.5, 49.5)
        histogram["rho_" + id ] = TH1F("rho_" + id, "rho", 200, 0., 50.)
        histogram["ele_pt_" + id] = TH1F("ele_pt_"+id,"Ele p_{T} [GeV/c]", 50, 0, 150)
        histogram["ele_scpt_" + id] = TH1F("ele_scpt_"+id,"Ele p_{T} (from SC) [GeV/c]", 50, 0, 150)
        histogram["ele_eta_" + id] = TH1F("ele_eta_"+id,"ele #eta", 50, -2.5, 2.5)
        histogram["ele_phi_"+ id] = TH1F("ele_phi_"+id,"ele #phi", 50, -math.pi, math.pi)

        # Z candidates mass plots general at various level of ID requirements
#        histogram["met_"+ id] = TH1F("met_"+id,"MET [GeV/c]", 50, 0. , 100.)
#        histogram["mt_"+ id] = TH1F("mt_"+id,"Transverse Mass [GeV/c^{2}]", 50, 0. , 200.)
        histogram["mee_"+ id] = TH1F("mee_"+id,"Invariant Mass [GeV/c^{2}]", 240, 30. , 150.)
        histogram["mee60_"+ id] = TH1F("mee60_"+id,"Invariant Mass [GeV/c^{2}]", 240, 30. , 150.)


        #Id variables just after preselection for EB and EE
        if id == "": 
            histogram["ele_EB_sigmaIetaIeta"+ id] = TH1F("ele_EB_sigmaIetaIeta"+id,"#sigma_{i#etai#eta} (EB)", 100, 0., 0.03)
            histogram["ele_EB_HOverE"+ id] = TH1F("ele_EB_HOverE"+id,"H/E (EB)", 100, 0., 0.15)
            histogram["ele_EB_CombinedIso"+ id] = TH1F("ele_EB_CombinedIso"+id,"CombinedIso (EB)", 100, 0., 0.15)
            histogram["ele_EB_ExpMissHits"+ id] = TH1F("ele_EB_ExpMissHits"+id,"Exp Miss. Hits (EB)", 10, -0.5, 9.5)
            histogram["ele_EE_sigmaIetaIeta"+ id] = TH1F("ele_EE_sigmaIetaIeta"+id,"#sigma_{i#etai#eta} (EE)", 100, 0.015, 0.05)
            histogram["ele_EE_HOverE"+ id] = TH1F("ele_EE_HOverE"+id,"H/E (EE)", 100, 0., 0.15)
            histogram["ele_EE_CombinedIso"+ id] = TH1F("ele_EE_CombinedIso"+id,"CombinedIso (EE)", 100, 0., 0.15)
            histogram["ele_EE_ExpMissHits"+ id] = TH1F("ele_EE_ExpMissHits"+id,"Exp Miss. Hits (EE)", 10, -0.5, 9.5)

            #Bkg enriched distributions
            histogram["ele_EB_sigmaIetaIetaBg_"+ id] = TH1F("ele_EB_sigmaIetaIetaBg_"+id,"#sigma_{i#etai#eta} (EB)", 100, 0., 0.03)
            histogram["ele_EB_HOverEBg_"+ id] = TH1F("ele_EB_HOverEBg_"+id,"H/E (EB)", 100, 0., 0.15)
            histogram["ele_EB_CombinedIsoBg_"+ id] = TH1F("ele_EB_CombinedIsoBg_"+id,"CombinedIso (EB)", 100, 0., 0.15)
            histogram["ele_EB_ExpMissHitsBg_"+ id] = TH1F("ele_EB_ExpMissHitsBg_"+id,"Exp Miss. Hits (EB)", 10, -0.5, 9.5)
            histogram["ele_EE_sigmaIetaIetaBg_"+ id] = TH1F("ele_EE_sigmaIetaIetaBg_"+id,"#sigma_{i#etai#eta} (EE)", 100, 0.015, 0.05)
            histogram["ele_EE_HOverEBg_"+ id] = TH1F("ele_EE_HOverEBg_"+id,"H/E (EE)", 100, 0., 0.15)
            histogram["ele_EE_CombinedIsoBg_"+ id] = TH1F("ele_EE_CombinedIsoBg_"+id,"CombinedIso (EE)", 100, 0., 0.15)
            histogram["ele_EE_ExpMissHitsBg_"+ id] = TH1F("ele_EE_ExpMissHitsBg_"+id,"Exp Miss. Hits (EE)", 10, -0.5, 9.5)

        # Some N-1 plots for EB and EE
        if id != "": # not booking them when no selection is applied

            # Mass plots for Tag&Probe 
            histogram["meetp_"+ id] = TH1F("meetp_"+id,"Invariant Mass [GeV/c^{2}]", 240, 30. , 150.)
            histogram["mee_EB_tp_"+ id] = TH1F("mee_EB_tp_"+id,"Invariant Mass [GeV/c^{2}]", 240, 30. , 150.)
            histogram["mee_EE_tp_"+ id] = TH1F("mee_EE_tp_"+id,"Invariant Mass [GeV/c^{2}]", 240, 30. , 150.)
            for id1 in electronIdLevels:
                if id1 != "": # not booking them when no selection is applied
                    histogram["meetp_pass_t"+ id +"_p"+ id1] = TH1F("meetp_pass_t"+ id+"_p"+ id1,"Invariant Mass Pass [GeV/c^{2}]", 240, 30. , 150.)
                    histogram["meetp_fail_t"+ id +"_p"+ id1] = TH1F("meetp_fail_t"+ id+"_p"+ id1,"Invariant Mass Pass [GeV/c^{2}]", 240, 30. , 150.)
                    histogram["mee_EB_tp_pass_t"+ id+"_p"+ id1] = TH1F("mee_EB_tp_pass_t"+ id+"_p"+ id1,"Invariant Mass Pass [GeV/c^{2}]", 240, 30. , 150.)
                    histogram["mee_EB_tp_fail_t"+ id+"_p"+ id1] = TH1F("mee_EB_tp_fail_t"+ id+"_p"+ id1,"Invariant Mass Pass [GeV/c^{2}]", 240, 30. , 150.)
                    histogram["mee_EE_tp_pass_t"+ id+"_p"+ id1] = TH1F("mee_EE_tp_pass_t"+ id+"_p"+ id1,"Invariant Mass Pass [GeV/c^{2}]", 240, 30. , 150.)
                    histogram["mee_EE_tp_fail_t"+ id+"_p"+ id1] = TH1F("mee_EE_tp_fail_t"+ id+"_p"+ id1,"Invariant Mass Pass [GeV/c^{2}]", 240, 30. , 150.)

            # Probe ID variables plots
            histogram["ele_EB_sigmaIetaIetaProbe_"+ id] = TH1F("ele_EB_sigmaIetaIetaProbe_"+id,"#sigma_{i#etai#eta} (EB)", 100, 0., 0.03)
            histogram["ele_EB_HOverEProbe_"+ id] = TH1F("ele_EB_HOverEProbe_"+id,"H/E (EB)", 100, 0., 0.15)
            histogram["ele_EB_CombinedIsoProbe_"+ id] = TH1F("ele_EB_CombinedIsoProbe_"+id,"CombinedIso (EB)", 100, 0., 0.15)
            histogram["ele_EB_ExpMissHitsProbe_"+ id] = TH1F("ele_EB_ExpMissHitsProbe_"+id,"Exp Miss. Hits (EB)", 10, -0.5, 9.5)
            histogram["ele_EE_sigmaIetaIetaProbe_"+ id] = TH1F("ele_EE_sigmaIetaIetaProbe_"+id,"#sigma_{i#etai#eta} (EE)", 100, 0.015, 0.05)
            histogram["ele_EE_HOverEProbe_"+ id] = TH1F("ele_EE_HOverEProbe_"+id,"H/E (EE)", 100, 0., 0.15)
            histogram["ele_EE_CombinedIsoProbe_"+ id] = TH1F("ele_EE_CombinedIsoProbe_"+id,"CombinedIso (EE)", 100, 0., 0.15)
            histogram["ele_EE_ExpMissHitsProbe_"+ id] = TH1F("ele_EE_ExpMissHitsProbe_"+id,"Exp Miss. Hits (EE)", 10, -0.5, 9.5)

            # N-1 ID variables plots
            histogram["ele_EB_sigmaIetaIetaNMinusOne_"+ id] = TH1F("ele_EB_sigmaIetaIetaNMinusOne_"+id,"#sigma_{i#etai#eta} (EB)", 100, 0., 0.03)
            histogram["ele_EB_HOverENMinusOne_"+ id] = TH1F("ele_EB_HOverENMinusOne_"+id,"H/E (EB)", 100, 0., 0.15)
            histogram["ele_EB_CombinedIsoNMinusOne_"+ id] = TH1F("ele_EB_CombinedIsoNMinusOne_"+id,"CombinedIso (EB)", 100, 0., 0.15)
            histogram["ele_EB_ExpMissHitsNMinusOne_"+ id] = TH1F("ele_EB_ExpMissHitsNMinusOne_"+id,"Exp Miss. Hits (EB)", 10, -0.5, 9.5)
            histogram["ele_EE_sigmaIetaIetaNMinusOne_"+ id] = TH1F("ele_EE_sigmaIetaIetaNMinusOne_"+id,"#sigma_{i#etai#eta} (EE)", 100, 0.015, 0.05)
            histogram["ele_EE_HOverENMinusOne_"+ id] = TH1F("ele_EE_HOverENMinusOne_"+id,"H/E (EE)", 100, 0., 0.15)
            histogram["ele_EE_CombinedIsoNMinusOne_"+ id] = TH1F("ele_EE_CombinedIsoNMinusOne_"+id,"CombinedIso (EE)", 100, 0., 0.15)
            histogram["ele_EE_ExpMissHitsNMinusOne_"+ id] = TH1F("ele_EE_ExpMissHitsNMinusOne_"+id,"Exp Miss. Hits (EE)", 10, -0.5, 9.5)
        
    for ih in histogram.keys():
        histogram[ih].Sumw2()
        histogram[ih].SetXTitle( histogram[ih].GetTitle() )
    
    # loop over events
    i = 0 # event counter

    print events.size()
    for event in events:
        i= i + 1
        if i%1000 == 0:
            printInfo(i,event)
            
        # check if maximum number of events was asked
#        if maxNevents > 0 and maxNevents == i:
#           print "Maximum number of events read "+str(maxNevents) 
#           break
        
        event.getByLabel (electronCollection, electronCollectionInstance, handleElectron)
        event.getByLabel (vertexCollection, handleVertices)
        event.getByLabel (rhoCollection, rhoCollectionInstance, handleRho)
        electrons = handleElectron.product()
        vertices = handleVertices.product()
        rho = handleRho.product()
        
#        event.getByLabel (metCollection, handleMET)
#        met = handleMET.product()

        Nelectrons = electrons.size()
        histogram["nEle"].Fill( Nelectrons )

        Nvertices = vertices.size()

        
        # ***********************************************************************************************
        # first loop on electrons filling kinematic variables and id variables for the leading electron
        # ***********************************************************************************************
        for electron in electrons:
            for id in electronIdLevels:
                if checkFullElectronId(electron,id):
                    histogram["ele_pt_" + id].Fill( electron.pt() )
                    histogram["ele_scpt_" + id].Fill( electron.caloEnergy()/math.cosh(electron.gsfTrack().eta()) )
                    histogram["ele_eta_" + id].Fill( electron.eta() )
                    histogram["ele_phi_" + id].Fill( electron.phi() )
                    #fill W selection plots (applying only 2nd lepton veto without id and pt Cut)
                    elePt=electron.caloEnergy()/math.cosh(electron.gsfTrack().eta())
                    if elePt<ptCut:
                        continue
                    #                     histogram["met_" + id].Fill( met[0].pt() )
                    #                     a = TVector3()
                    #                     a.SetPtEtaPhi(elePt, 0. ,electron.phi())
                    #                     b = TVector3()
                    #                     b.SetPtEtaPhi( met[0].pt(), 0. , met[0].phi())
                    #                     mt = sqrt(2 * a.Mag() * b.Mag() * (1-math.cos(a.Angle(b))) )
                    #                     histogram["mt_" + id].Fill( mt )
                    if id == "": #filling some id variables just after preselection for the leading electron
                        tkIso = electron.dr03TkSumPt();
                        ecalIsoPed = (max([0.,electron.dr04EcalRecHitSumEt()-1.]) if electron.isEB() else electron.dr04EcalRecHitSumEt());
                        hcalIso = electron.dr04HcalTowerSumEt();
                        if (electron.isEB()):
                            histogram["ele_EB_CombinedIso"].Fill( (tkIso+ecalIsoPed+hcalIso)/electron.pt() )
                            histogram["ele_EB_sigmaIetaIeta"].Fill( electron.sigmaIetaIeta() )
                            histogram["ele_EB_HOverE"].Fill( electron.hadronicOverEm() )
                            histogram["ele_EB_ExpMissHits"].Fill( electron.gsfTrack().trackerExpectedHitsInner().numberOfHits() )
                        else:
                            histogram["ele_EE_CombinedIso"].Fill( (tkIso+ecalIsoPed+hcalIso)/electron.pt() )
                            histogram["ele_EE_sigmaIetaIeta"].Fill( electron.sigmaIetaIeta() )
                            histogram["ele_EE_HOverE"].Fill( electron.hadronicOverEm() )
                            histogram["ele_EE_ExpMissHits"].Fill( electron.gsfTrack().trackerExpectedHitsInner().numberOfHits() )

                        
        # ****************************************************************
        # fill plots for symmetric Z selection
        # ****************************************************************
        if Nelectrons>1:
            for id in electronIdLevels:
                #Avoid full combinatorics just use the first 2 highest pt electrons
                elePt1=electrons[0].caloEnergy()/math.cosh(electrons[0].gsfTrack().eta())
                elePt2=electrons[1].caloEnergy()/math.cosh(electrons[1].gsfTrack().eta())
                if (elePt1<ptCut or elePt2<ptCut): # requiring both electrons above ptCut
                    continue
                v1=TLorentzVector()
                v2=TLorentzVector()
                v1.SetPtEtaPhiM(elePt1,electrons[0].eta(),electrons[0].phi(),0.)
                v2.SetPtEtaPhiM(elePt2,electrons[1].eta(),electrons[1].phi(),0.)
                vZ=v1+v2
                
                # ************************************************************************************************************************************************************************************************
                # Rudimentary tag & probe selection for different tag selection. Probe here is a gsfElectron satistying the ptCut. Arbitration of the tag is done using a random number.
                # It will be useful only to study the ID cut efficiencies
                # ************************************************************************************************************************************************************************************************
                tagIsE1 = rand.Uniform()
                if ( (id != "") and ( ( (tagIsE1 <0.5) and checkFullElectronId(electrons[0],id) ) or ( (tagIsE1 >0.5) and checkFullElectronId(electrons[1],id) ) ) ): 
                    if (vZ.Mag()<60. or vZ.Mag()>120.):
                        continue

                    histogram["meetp_" + id].Fill( vZ.Mag() )

                    iele=1
                    if ( tagIsE1> 0.5):
                        iele=0

                    if (electrons[iele].isEB()):
                        histogram["mee_EB_tp_" + id].Fill( vZ.Mag() )
                    else:
                        histogram["mee_EE_tp_" + id].Fill( vZ.Mag() )

                    tkIso = electrons[iele].dr03TkSumPt()
                    ecalIsoPed = (max([0.,electrons[iele].dr04EcalRecHitSumEt()-1.]) if electrons[iele].isEB() else electrons[iele].dr04EcalRecHitSumEt())
                    hcalIso = elecAtrons[iele].dr04HcalTowerSumEt()

                    if (electrons[iele].isEB()):
                        histogram["ele_EB_CombinedIsoProbe_"+ id].Fill( (tkIso+ecalIsoPed+hcalIso)/electrons[iele].pt() )
                        histogram["ele_EB_sigmaIetaIetaProbe_"+ id].Fill( electrons[iele].sigmaIetaIeta() )
                        histogram["ele_EB_HOverEProbe_"+ id].Fill( electrons[iele].hadronicOverEm() )
                        histogram["ele_EB_ExpMissHitsProbe_"+ id].Fill( electrons[iele].gsfTrack().trackerExpectedHitsInner().numberOfHits() ) 
                    else:
                        histogram["ele_EE_CombinedIsoProbe_"+ id].Fill( (tkIso+ecalIsoPed+hcalIso)/electrons[iele].pt() )
                        histogram["ele_EE_sigmaIetaIetaProbe_"+ id].Fill( electrons[iele].sigmaIetaIeta() )
                        histogram["ele_EE_HOverEProbe_"+ id].Fill( electrons[iele].hadronicOverEm() )
                        histogram["ele_EE_ExpMissHitsProbe_"+ id].Fill( electrons[iele].gsfTrack().trackerExpectedHitsInner().numberOfHits() )

                    for id1 in electronIdLevels:
                        if id1=="":
                            continue
                        if (checkFullElectronId(electrons[iele],id1)):
                            histogram["meetp_pass_t" + id+"_p" + id1].Fill( vZ.Mag() )
                            if (electrons[iele].isEB()):
                                histogram["mee_EB_tp_pass_t" + id+"_p" + id1].Fill( vZ.Mag() )
                            else:
                                histogram["mee_EE_tp_pass_t" + id+"_p" + id1].Fill( vZ.Mag() )
                        else:
                            histogram["meetp_fail_t" + id+"_p" + id1].Fill( vZ.Mag() )
                            if (electrons[iele].isEB()):
                                histogram["mee_EB_tp_fail_t" + id+"_p" + id1].Fill( vZ.Mag() )
                            else:
                                histogram["mee_EE_tp_fail_t" + id+"_p" + id1].Fill( vZ.Mag() )
                                
                if (checkFullElectronId(electrons[0],id) and checkFullElectronId(electrons[1],id) ): # checking id for both electrons
                    histogram["mee_" + id].Fill( vZ.Mag() )
                    histogram["nVertices_" + id ].Fill( Nvertices )
                    histogram["rho_" + id ].Fill( rho[0] )

                    if (vZ.Mag()<60. or vZ.Mag()>120.):
                        continue

                    histogram["mee60_" + id].Fill( vZ.Mag() )

                    for iele in range(0,1):
                        tkIso = electrons[iele].dr03TkSumPt()
                        ecalIsoPed = (max([0.,electrons[iele].dr04EcalRecHitSumEt()-1.]) if electrons[iele].isEB() else electrons[iele].dr04EcalRecHitSumEt())
                        hcalIso = electrons[iele].dr04HcalTowerSumEt()
                        for id in range(1,len(electronIdLevels)): #avoid the first level without any selection
                            if checkElectronIdBit(electrons[iele],electronIdLevels[id],0) and checkElectronIdBit(electrons[iele],electronIdLevels[id],2): #requiring only ElectronID and Conv. Rejection
                                if (electrons[iele].isEB()):
                                    histogram["ele_EB_CombinedIsoNMinusOne_"+ electronIdLevels[id]].Fill( (tkIso+ecalIsoPed+hcalIso)/electrons[iele].pt() )
                                else:
                                    histogram["ele_EE_CombinedIsoNMinusOne_"+ electronIdLevels[id]].Fill( (tkIso+ecalIsoPed+hcalIso)/electrons[iele].pt() )
                                
                            if checkElectronIdBit(electrons[iele],electronIdLevels[id],1) and checkElectronIdBit(electrons[iele],electronIdLevels[id],2): #requiring only Isolation and Conv. Rejection
                                if (electrons[iele].isEB()):
                                    histogram["ele_EB_sigmaIetaIetaNMinusOne_"+ electronIdLevels[id]].Fill( electrons[iele].sigmaIetaIeta() )
                                    histogram["ele_EB_HOverENMinusOne_"+ electronIdLevels[id]].Fill( electrons[iele].hadronicOverEm() )
                                else:
                                    histogram["ele_EE_sigmaIetaIetaNMinusOne_"+ electronIdLevels[id]].Fill( electrons[iele].sigmaIetaIeta() )
                                    histogram["ele_EE_HOverENMinusOne_"+ electronIdLevels[id]].Fill( electrons[iele].hadronicOverEm() )
                                

                            if checkElectronIdBit(electrons[iele],electronIdLevels[id],0) and checkElectronIdBit(electrons[iele],electronIdLevels[id],1): #requiring only ElectronId and Isolation
                                if (electrons[iele].isEB()):
                                    histogram["ele_EB_ExpMissHitsNMinusOne_"+ electronIdLevels[id]].Fill( electrons[iele].gsfTrack().trackerExpectedHitsInner().numberOfHits() )
                                else:
                                    histogram["ele_EE_ExpMissHitsNMinusOne_"+ electronIdLevels[id]].Fill( electrons[iele].gsfTrack().trackerExpectedHitsInner().numberOfHits() )



                    
    # close loop over entries    
    outputroot.cd()

    # write histograms to file
    for key in histogram.keys():
        histogram[key].Write()

    outputroot.Close()

if __name__ == '__main__':
    main()
