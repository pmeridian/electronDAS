#! /usr/bin/env python

from ROOT import *

import sys
import os
import math

gROOT.SetStyle("Plain")
gROOT.SetBatch()

inputfiles = {}
inputfiles["data"] = TFile("electronDistributionsDATA.root")
inputfiles["mc"] = TFile("electronDistributionsMC.root")

electronIdLevels= [ "","simpleEleId95cIso","simpleEleId80cIso","simpleEleId70cIso" ]
electronIdLevelsForNMinusOne= [ "simpleEleId70cIso","simpleEleId80cIso","simpleEleId95cIso","" ]

weight = {}
weight["data"] = 1.
weight["mc"] = 1. #possibity to reweight to match data luminosity

histogram = {}
histogram_MC = {}

canvas = TCanvas("c","c",1)

#Get all histograms
for input in inputfiles.keys():
    inputfiles[input].cd()
    for id in electronIdLevels:
        histogram[input,"ele_pt_" + id] = gDirectory.Get("ele_pt_" + id)
        histogram[input,"ele_scpt_" + id] = gDirectory.Get("ele_scpt_" + id)
        histogram[input,"ele_eta_" + id] = gDirectory.Get("ele_eta_" + id)
        histogram[input,"ele_phi_"+ id] = gDirectory.Get("ele_phi_" + id)
#        histogram[input,"met_"+ id] = gDirectory.Get("met_" + id)
#        histogram[input,"mt_"+ id] = gDirectory.Get("mt_" + id)
        histogram[input,"mee_"+ id] = gDirectory.Get("mee_" + id)
        histogram[input,"mee60_"+ id] = gDirectory.Get("mee60_" + id)
        if id != "":
            histogram[input,"meetp_"+ id] = gDirectory.Get("meetp_" + id)
            for id1 in electronIdLevels:
                histogram[input,"meetp_pass_t"+ id+"_p"+ id1] = gDirectory.Get("meetp_pass_t"+ id+"_p"+ id1)
                histogram[input,"meetp_fail_t"+ id+"_p"+ id1] = gDirectory.Get("meetp_fail_t"+ id+"_p"+ id1)
                
        for detector in ["EB","EE"]:
            if id != "":
                histogram[input,"ele_"+ detector + "_sigmaIetaIetaNMinusOne_"+ id] = gDirectory.Get("ele_"+ detector + "_sigmaIetaIetaNMinusOne_"+id)
                histogram[input,"ele_"+ detector + "_HOverENMinusOne_"+ id] = gDirectory.Get("ele_"+ detector + "_HOverENMinusOne_"+id)
                histogram[input,"ele_"+ detector + "_CombinedIsoNMinusOne_"+ id] = gDirectory.Get("ele_"+ detector + "_CombinedIsoNMinusOne_"+id)
                histogram[input,"ele_"+ detector + "_ExpMissHitsNMinusOne_"+ id] = gDirectory.Get("ele_"+ detector + "_ExpMissHitsNMinusOne_"+id)
                histogram[input,"ele_"+ detector + "_sigmaIetaIetaProbe_"+ id] = gDirectory.Get("ele_"+ detector + "_sigmaIetaIetaProbe_"+id)
                histogram[input,"ele_"+ detector + "_HOverEProbe_"+ id] = gDirectory.Get("ele_"+ detector + "_HOverEProbe_"+id)
                histogram[input,"ele_"+ detector + "_CombinedIsoProbe_"+ id] = gDirectory.Get("ele_"+ detector + "_CombinedIsoProbe_"+id)
                histogram[input,"ele_"+ detector + "_ExpMissHitsProbe_"+ id] = gDirectory.Get("ele_"+ detector + "_ExpMissHitsProbe_"+id)
            else:
                histogram[input,"ele_"+ detector + "_sigmaIetaIeta"] = gDirectory.Get("ele_"+ detector + "_sigmaIetaIeta")
                histogram[input,"ele_"+ detector + "_HOverE"] = gDirectory.Get("ele_"+ detector + "_HOverE")
                histogram[input,"ele_"+ detector + "_CombinedIso"] = gDirectory.Get("ele_"+ detector + "_CombinedIso")
                histogram[input,"ele_"+ detector + "_ExpMissHits"] = gDirectory.Get("ele_"+ detector + "_ExpMissHits")

os.system("mkdir -p plots")

# *----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----* 
# Plotting kinematic distributions for increasing tightness of electronId (linear and log) for data and MC
# *----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*
for input in inputfiles.keys():
    print "++++++++++++++++++ " + input + " ++++++++++++++++++++" 
    for key in ["ele_pt","ele_scpt","ele_eta","ele_phi","mee","mee60"]:
        print "Drawing " + key
        icolor=1
        for id in electronIdLevels:
            canvas.SetLogy(0)
            if key == "ele_pt":
                histogram[input,key + "_" + id].GetXaxis().SetRangeUser(0.,80.)
                
            histogram[input,key + "_" + id].Scale(weight[input])
            histogram[input,key + "_" + id].SetMarkerStyle(20)
            histogram[input,key + "_" + id].SetMarkerSize(1.)
            histogram[input,key + "_" + id].SetMarkerColor(icolor)
            #        histogram[input,key + "_" + id].SetLineWidth()
            histogram[input,key + "_" + id].SetLineColor(icolor)
            if id == "":
                histogram[input,key + "_" + id].SetMinimum(0.)
                histogram[input,key + "_" + id].Draw("E")
            else:
                histogram[input,key + "_" + id].Draw("ESAME")
            icolor=icolor+1
            canvas.Update()            
        canvas.SaveAs("plots/"+ key + "_"+ input + ".png")
            
        # Plotting also in Logarithmic scale
        icolor=1
        for id in electronIdLevels:
            if key == "ele_pt":
                histogram[input,key + "_" + id].GetXaxis().SetRangeUser(0.,80.)
            histogram[input,key + "_" + id].SetMinimum(0.1)
            histogram[input,key + "_" + id].SetMaximum(histogram[input,key + "_" + id].GetMaximum()*10.)
            canvas.SetLogy(1)
            if id == "":
                histogram[input,key + "_" + id].Draw("E")
            else:
                histogram[input,key + "_" + id].Draw("ESAME")
            icolor = icolor + 1
        canvas.SaveAs("plots/"+ key + "_"+ input + "_log.png")

# *----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*
# Now plotting also N-1 plot (linear and log) for data 
# *----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*

for key in ["sigmaIetaIeta","HOverE","CombinedIso","ExpMissHits"]:
    for type in ["NMinusOne_","Probe_"]:
        print "Drawing " + key
        for detector in ["EB","EE"]:
            for input in inputfiles.keys():
                icolor=4
                for id in electronIdLevelsForNMinusOne:
                    canvas.SetLogy(0)

                    if id != "":
                        histoKey="ele_"+ detector + "_"+ key + type + id
                    else:
                        histoKey="ele_"+ detector + "_"+ key
                
                        histogram[input,histoKey].SetMinimum(0.)
                        histogram[input,histoKey].SetMarkerStyle(20)
                        histogram[input,histoKey].SetMarkerSize(1.)
                        histogram[input,histoKey].SetMarkerColor(icolor)
                        #        histogram[input,histoKey].SetLineWidth()
                        histogram[input,histoKey].SetLineColor(icolor)
                        if id == "simpleEleId70cIso":
                            histogram[input,histoKey].DrawNormalized("E")
                        else:
                            histogram[input,histoKey].DrawNormalized("ESAME")
                        icolor=icolor-1
            canvas.SaveAs("plots/"+ key + "_" + detector + "_" + input+ "_"+ type+".png")

            for input in inputfiles.keys():
                icolor=4
                for id in electronIdLevelsForNMinusOne:
                    canvas.SetLogy(1)
                    
                    if id != "":
                        histoKey="ele_"+ detector + "_"+ key + type + id
                    else:
                        histoKey="ele_"+ detector + "_"+ key
                

                    histogram[input,histoKey].SetMarkerStyle(20)
                    histogram[input,histoKey].SetMarkerSize(1.)
                    histogram[input,histoKey].SetMarkerColor(icolor)
                    #        histogram[input,histoKey].SetLineWidth()
                    histogram[input,histoKey].SetLineColor(icolor)
                    if id == "simpleEleId70cIso" :
                        histogram[input,histoKey].DrawNormalized("E")
                    else:
                        histogram[input,histoKey].DrawNormalized("ESAME")
                    histogram[input,histoKey].SetMinimum(0.00001)
                    histogram[input,histoKey].SetMaximum(histogram[input,histoKey].GetMaximum()*10.)
                    canvas.Update()
                    icolor=icolor-1
                    canvas.SaveAs("plots/"+ key + "_" + detector + "_" + input+ "_"+ type+"log.png")

            
