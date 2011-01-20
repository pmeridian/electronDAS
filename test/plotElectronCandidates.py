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
#    if not ele.ecalDrivenSeed():
#        return false
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
#    if not ele.ecalDrivenSeed():
#        return false
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
    maxNevents = -1

    # Operating point
    electronCollection = "electronPATFilter"
    electronCollectionInstance = "filteredPATElectronCandidates"
    metCollection ="patMETsPF"
    
    #change this according to your needs
    outfilename = "electronDistributions.root"
    outputroot = TFile( outfilename, "RECREATE")
        
    prefixFnal = 'dcache:/pnfs/cms/WAX/11/store/user/meridian/electronDAS'
    prefixCern = 'rfio:/castor/cern.ch/user/m/meridian/electronDAS'
    prefixLocal = '/cmsrm/pc24_2/meridian'

    prefix = prefixLocal

    # Kinamatic cuts used in the analysis
    ptCut = 20.
    metCut = 20.
    mtCut = 40.

    # PAT ntuples with electronCollection
    files = [
        # DATA Run 148862-148864  
        "/Electron22DecPAT/electronsPATTuple_1_1_CY7.root",
        "/Electron22DecPAT/electronsPATTuple_2_1_6aR.root",
        "/Electron22DecPAT/electronsPATTuple_3_1_Kuc.root",
        "/Electron22DecPAT/electronsPATTuple_4_1_HCD.root",
        "/Electron22DecPAT/electronsPATTuple_5_1_Ga9.root",
        "/Electron22DecPAT/electronsPATTuple_6_1_u6B.root",
        "/Electron22DecPAT/electronsPATTuple_7_1_3Ax.root",
        "/Electron22DecPAT/electronsPATTuple_8_1_o8S.root",
        "/Electron22DecPAT/electronsPATTuple_9_1_XSq.root",
        "/Electron22DecPAT/electronsPATTuple_10_1_Jbe.root",
        "/Electron22DecPAT/electronsPATTuple_11_0_cKj.root",
        "/Electron22DecPAT/electronsPATTuple_12_0_jRu.root"

#        MC WEnu files 
#        "/WToENu_TuneZ2_7TeV-pythia6_v1/electronsPATTuple_10_1_Ldq.root",
#        "/WToENu_TuneZ2_7TeV-pythia6_v1/electronsPATTuple_11_1_uBy.root",
#        "/WToENu_TuneZ2_7TeV-pythia6_v1/electronsPATTuple_12_1_EO8.root",
#        "/WToENu_TuneZ2_7TeV-pythia6_v1/electronsPATTuple_13_1_q4l.root",
#        "/WToENu_TuneZ2_7TeV-pythia6_v1/electronsPATTuple_14_1_1dL.root",

#       MC DYtoEE files
#        "/DYToEE_M-20_TuneZ2_7TeV-pythia6_v1/electronsPATTuple_11_1_jTy.root",
#        "/DYToEE_M-20_TuneZ2_7TeV-pythia6_v1/electronsPATTuple_12_1_jeM.root",
#        "/DYToEE_M-20_TuneZ2_7TeV-pythia6_v1/electronsPATTuple_13_1_zsz.root",
#        "/DYToEE_M-20_TuneZ2_7TeV-pythia6_v1/electronsPATTuple_17_1_7Yk.root",
#        "/DYToEE_M-20_TuneZ2_7TeV-pythia6_v1/electronsPATTuple_18_1_m9E.root",

        ]
    
    fullpath_files = []
    
    for afile in files:
        fullpath_files.append( prefix+afile )

    print 'Processing files '
    for ifile in fullpath_files:
            print ifile
            
    events = Events ( fullpath_files )

    handleElectron = Handle ("vector<pat::Electron>")
    handleMET = Handle ("vector<pat::MET>")

    electronIdLevels= [ "","simpleEleId95cIso","simpleEleId80cIso","simpleEleId70cIso" ]

    # ****************************************************************
    # Booking histograms
    # ****************************************************************
    histogram = {}
    histogram["nEle"] = TH1F("nEle","nEle", 20, -0.5, 19.5)

    for id in electronIdLevels:
        # general electron kinematics carachteristics
        histogram["ele_pt_" + id] = TH1F("ele_pt_"+id,"Ele p_{T} [GeV/c]", 50, 0, 150)
        histogram["ele_scpt_" + id] = TH1F("ele_scpt_"+id,"Ele p_{T} (from SC) [GeV/c]", 50, 0, 150)
        histogram["ele_eta_" + id] = TH1F("ele_eta_"+id,"ele #eta", 50, -2.5, 2.5)
        histogram["ele_phi_"+ id] = TH1F("ele_phi_"+id,"ele #phi", 50, -math.pi, math.pi)

        # W and Z candidates plots
        histogram["met_"+ id] = TH1F("met_"+id,"MET [GeV/c]", 50, 0. , 100.)
        histogram["mt_"+ id] = TH1F("mt_"+id,"Transverse Mass [GeV/c^{2}]", 50, 0. , 200.)
        histogram["mee_"+ id] = TH1F("mee_"+id,"Invariant Mass [GeV/c^{2}]", 100, 0. , 200.)

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

        # Some N-1 plots for EB and EE
        if id != "": # not booking them when no selection is applied
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

    for event in events:
        i = i + 1
        if i%1000 == 0:
            printInfo(i,event)

        # check if maximum number of events was asked
        if maxNevents > 0 and maxNevents == i:
           print "Maximum number of events read "+str(maxNevents) 
           break
        
        event.getByLabel (electronCollection, electronCollectionInstance, handleElectron)
        electrons = handleElectron.product()
        event.getByLabel (metCollection, handleMET)
        met = handleMET.product()

        Nelectrons = electrons.size()
        histogram["nEle"].Fill( Nelectrons )

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

                    if Nelectrons == 1:
                        #fill W selection plots (applying only 2nd lepton veto without id and pt Cut)
                        elePt=electron.caloEnergy()/math.cosh(electron.gsfTrack().eta())
                        if elePt<ptCut:
                            continue
                        histogram["met_" + id].Fill( met[0].pt() )
                        a = TVector3()
                        a.SetPtEtaPhi(elePt, 0. ,electron.phi())
                        b = TVector3()
                        b.SetPtEtaPhi( met[0].pt(), 0. , met[0].phi())
                        mt = sqrt(2 * a.Mag() * b.Mag() * (1-math.cos(a.Angle(b))) )
                        histogram["mt_" + id].Fill( mt )
                        if id == "": #filling some id variables just after preselection for the leading electron
                            tkIso = electrons[0].dr03TkSumPt();
                            ecalIsoPed = (max([0.,electrons[0].dr04EcalRecHitSumEt()-1.]) if electrons[0].isEB() else electrons[0].dr04EcalRecHitSumEt());
                            hcalIso = electrons[0].dr04HcalTowerSumEt();
                            if (electrons[0].isEB()):
                                histogram["ele_EB_CombinedIso"].Fill( (tkIso+ecalIsoPed+hcalIso)/electrons[0].pt() )
                                histogram["ele_EB_sigmaIetaIeta"].Fill( electrons[0].sigmaIetaIeta() )
                                histogram["ele_EB_HOverE"].Fill( electrons[0].hadronicOverEm() )
                                histogram["ele_EB_ExpMissHits"].Fill( electrons[0].gsfTrack().trackerExpectedHitsInner().numberOfHits() )
                            else:
                                histogram["ele_EE_CombinedIso"].Fill( (tkIso+ecalIsoPed+hcalIso)/electrons[0].pt() )
                                histogram["ele_EE_sigmaIetaIeta"].Fill( electrons[0].sigmaIetaIeta() )
                                histogram["ele_EE_HOverE"].Fill( electrons[0].hadronicOverEm() )
                                histogram["ele_EE_ExpMissHits"].Fill( electrons[0].gsfTrack().trackerExpectedHitsInner().numberOfHits() )

                        
        # ****************************************************************
        # fill plots for symmetric Z selection
        # ****************************************************************
        if Nelectrons>1:
            for id in electronIdLevels:
                for x in range(len(electrons)-1):
                    for y in range(x+1,len(electrons)):
                        elePt1=electrons[x].caloEnergy()/math.cosh(electrons[x].gsfTrack().eta())
                        elePt2=electrons[y].caloEnergy()/math.cosh(electrons[y].gsfTrack().eta())
                        if (elePt1<ptCut or elePt2<ptCut): # requiring both electrons above ptCut
                            continue
                        if (checkFullElectronId(electrons[x],id) and checkFullElectronId(electrons[y],id) ): # checking id for both electrons
                            v1=TLorentzVector()
                            v2=TLorentzVector()
                            v1.SetPtEtaPhiM(elePt1,electrons[x].eta(),electrons[x].phi(),0.)
                            v2.SetPtEtaPhiM(elePt2,electrons[y].eta(),electrons[y].phi(),0.)
                            vZ=v1+v2
                            histogram["mee_" + id].Fill( vZ.Mag() )
                            
        # ****************************************************************                            
        # doing now n-1 plots for some Id variables for W like candidates
        # ****************************************************************
        else:
            elePt=electrons[0].caloEnergy()/math.cosh(electrons[0].gsfTrack().eta())
            a = TVector3()
            a.SetPtEtaPhi(electrons[0].caloEnergy()/math.cosh(electrons[0].gsfTrack().eta()), 0. ,electrons[0].phi())
            b = TVector3()
            b.SetPtEtaPhi( met[0].pt(), 0. , met[0].phi())
            mt = sqrt(2 * a.Mag() * b.Mag() * (1-math.cos(a.Angle(b))))
            if (elePt> ptCut and met[0].pt()> metCut and mt> mtCut): #MET and MT cuts to increase purity of W selection for N-1 plots
                tkIso = electrons[0].dr03TkSumPt();
                ecalIsoPed = (max([0.,electrons[0].dr04EcalRecHitSumEt()-1.]) if electrons[0].isEB() else electrons[0].dr04EcalRecHitSumEt());
                hcalIso = electrons[0].dr04HcalTowerSumEt();
                for id in range(1,len(electronIdLevels)): #avoid the first level without any selection
                    #Filling combined isolation N-1 plot
                    if checkElectronIdBit(electrons[0],electronIdLevels[id],0) and checkElectronIdBit(electrons[0],electronIdLevels[id],2): #requiring only ElectronID and Conv. Rejection
                        if (electrons[0].isEB()):
                            histogram["ele_EB_CombinedIsoNMinusOne_"+ electronIdLevels[id]].Fill( (tkIso+ecalIsoPed+hcalIso)/electrons[0].pt() )
                        else:
                            histogram["ele_EE_CombinedIsoNMinusOne_"+ electronIdLevels[id]].Fill( (tkIso+ecalIsoPed+hcalIso)/electrons[0].pt() )

                    #Filling SigmaIEtaIeta and HOverE N-1 plots
                    if checkElectronIdBit(electrons[0],electronIdLevels[id],1) and checkElectronIdBit(electrons[0],electronIdLevels[id],2): #requiring only Isolation and Conv. Rejection
                        if (electrons[0].isEB()):
                            histogram["ele_EB_sigmaIetaIetaNMinusOne_"+ electronIdLevels[id]].Fill( electrons[0].sigmaIetaIeta() )
                            histogram["ele_EB_HOverENMinusOne_"+ electronIdLevels[id]].Fill( electrons[0].hadronicOverEm() )
                        else:
                            histogram["ele_EE_sigmaIetaIetaNMinusOne_"+ electronIdLevels[id]].Fill( electrons[0].sigmaIetaIeta() )
                            histogram["ele_EE_HOverENMinusOne_"+ electronIdLevels[id]].Fill( electrons[0].hadronicOverEm() )

                    #Filling missHits variables N-1 Plots
                    if checkElectronIdBit(electrons[0],electronIdLevels[id],0) and checkElectronIdBit(electrons[0],electronIdLevels[id],1): #requiring only ElectronId and Isolation
                        if (electrons[0].isEB()):
                            histogram["ele_EB_ExpMissHitsNMinusOne_"+ electronIdLevels[id]].Fill( electrons[0].gsfTrack().trackerExpectedHitsInner().numberOfHits() )
                        else:
                            histogram["ele_EE_ExpMissHitsNMinusOne_"+ electronIdLevels[id]].Fill( electrons[0].gsfTrack().trackerExpectedHitsInner().numberOfHits() )
                    
    # close loop over entries    
    outputroot.cd()

    # write histograms to file
    for key in histogram.keys():
        histogram[key].Write()

    outputroot.Close()

if __name__ == '__main__':
    main()
