#! /usr/bin/env python


import sys
import math
from array import array
from ROOT import *
print "loading FWLite libraries... ",
from DataFormats.FWLite import Events, Handle

print "done."

def checkElectronId(ele,cut):
    if cut == "":
        return true
    else:
        if ele.isElectronIDAvailable(cut):
            if fabs(ele.electronID(cut)-7)<0.1:
                return true
            else:
                return false
        else:
            return false
    
def main():

    # if you want to run in batch mode
    ROOT.gROOT.SetBatch()
    # maximum number of events. -1 run over all
    maxNevents = -1

    # Operating point
    electronCollection = "electronPATFilter"
    electronCollectionInstance = "filteredPATElectronCandidates"
    metCollection ="patMETsPF"
    
    outfilename = "electronDistributions.root"
    outputroot = TFile( outfilename, "RECREATE")
        
    prefixFnal = 'dcache:/pnfs/cms/WAX/11'
    prefixCern = 'rfio:/castor/cern.ch/cms'
    prefixLocal = ''
    prefix = prefixLocal

    # PAT ntuples with electronCollection
    files = [
 #       '/cmsrm/pc24_2/meridian/D6B89C71-4B12-E011-8F7D-001A92971B5E.root'
        'electronsPATTuple.root'
        ]
    
    fullpath_files = []
    
    for afile in files:
        fullpath_files.append( prefix+afile )
    
    events = Events ( fullpath_files )

    handleElectron = Handle ("vector<pat::Electron>")
    handleMET = Handle ("vector<pat::MET>")

    electronIdLevels= [ "","simpleEleId95cIso","simpleEleId80cIso","simpleEleId70cIso" ]

    histogram = {}
    histogram["nEle"] = TH1F("nEle","nEle", 20, -0.5, 19.5)

    for id in electronIdLevels:
        histogram["ele_pt" + id] = TH1F("ele_pt"+id,"Ele p_{T} [GeV/c]", 50, 0, 300)
        histogram["ele_scpt" + id] = TH1F("ele_scpt"+id,"Ele p_{T} (from SC) [GeV/c]", 50, 0, 300)
        histogram["ele_eta" + id] = TH1F("ele_eta"+id,"ele #eta", 50, -2.5, 2.5)
        histogram["ele_phi"+ id] = TH1F("ele_phi"+id,"ele #phi", 50, -math.pi, math.pi)
        histogram["met"+ id] = TH1F("met"+id,"MET [GeV/c]", 50, 0. , 100.)
        histogram["mt"+ id] = TH1F("mt"+id,"Transverse Mass [GeV/c^{2}]", 50, 0. , 200.)
        histogram["mee"+ id] = TH1F("mee"+id,"Invariant Mass [GeV/c^{2}]", 100, 0. , 200.)

    for ih in histogram.keys():
        histogram[ih].Sumw2()
        histogram[ih].SetXTitle( histogram[ih].GetTitle() )
    
    # loop over events
    i = 0 # event counter

    for event in events:
        i = i + 1
        if i%100 == 0:
            print  "processing entry # " + str(i) + " from Run "+ str(event.eventAuxiliary().id().run()) + " lumi "+str(event.eventAuxiliary().id().luminosityBlock())

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
        
        for electron in electrons:
            for id in electronIdLevels:
                if checkElectronId(electron,id):
                    histogram["ele_pt" + id].Fill( electron.pt() )
                    histogram["ele_scpt" + id].Fill( electron.caloEnergy()/math.cosh(electron.gsfTrack().eta()) )
                    histogram["ele_eta" + id].Fill( electron.eta() )
                    histogram["ele_phi" + id].Fill( electron.phi() )
                    if Nelectrons == 1:
                        #fill also W selection plots (applying 2nd lepton veto without id)
                        histogram["met" + id].Fill( met[0].pt() )
                        a = TVector3()
                        a.SetPtEtaPhi(electron.caloEnergy()/math.cosh(electron.gsfTrack().eta()), 0. ,electron.phi())
                        b = TVector3()
                        b.SetPtEtaPhi( met[0].pt(), 0. , met[0].phi())
                        mt = sqrt(2 * a.Mag() * b.Mag() * (1-math.cos(a.Angle(b))) )
                        histogram["mt" + id].Fill( mt )

        if Nelectrons>1:
            for id in electronIdLevels:
                #fill plots for symmetric Z selection
                for x in range(len(electrons)-1):
                    for y in range(x+1,len(electrons)):
                        if (checkElectronId(electrons[x],id) and checkElectronId(electrons[y],id)):
                            v1=TLorentzVector()
                            v2=TLorentzVector()
                            v1.SetPtEtaPhiM(electrons[x].caloEnergy()/math.cosh(electrons[x].gsfTrack().eta()),electrons[x].eta(),electrons[x].phi(),0.)
                            v2.SetPtEtaPhiM(electrons[y].caloEnergy()/math.cosh(electrons[y].gsfTrack().eta()),electrons[y].eta(),electrons[y].phi(),0.)
                            vZ=v1+v2
                            histogram["mee" + id].Fill( vZ.Mag() )

    # close loop over entries    
    outputroot.cd()

    # write histograms to file
    for key in histogram.keys():
        histogram[key].Write()

    outputroot.Close()

if __name__ == '__main__':
    main()
