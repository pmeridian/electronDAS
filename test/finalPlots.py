#! /usr/bin/env python

from ROOT import *

import sys
import os
import math

gROOT.SetStyle("Plain")
gROOT.SetBatch()

inputfiles = {}
inputfiles["mc"] = TFile("electronDistributionsMC.root")
inputfiles["data"] = TFile("electronDistributionsDATA.root")

electronIdLevels= [ "","simpleEleId90cIso","simpleEleId80cIso","simpleEleId70cIso" ]
electronIdLevelsForNMinusOne= [ "simpleEleId70cIso","" ]

histogram = {}

canvas = TCanvas("c","c",1)

#Get all histograms
for input in inputfiles.keys():
    inputfiles[input].cd()
    for id in electronIdLevels:
        for key in ["nVertices","rho","ele_pt","ele_scpt","ele_eta","ele_phi","mee","mee60"]:
            histogram[input,key+"_"+id] = gDirectory.Get(key+"_"+id)
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

#rescaling factor for MC (instead of lumi reweighting)
weight = {}
leg=TLegend(0.5,0.85,0.7,0.99)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.033)
for key in ["nVertices","rho","ele_pt","ele_scpt","ele_eta","ele_phi","mee","mee60"]:
    leg.Clear()
    print "++++++++++++++++++ " + key + " ++++++++++++++++++++"
    scaleFactor=histogram["data",key+"_simpleEleId90cIso"].Integral()/histogram["mc",key+"_simpleEleId90cIso"].Integral()
    maximum=max([histogram["data",key+"_"].GetMaximum(),histogram["data",key+"_"].GetMaximum()])*1.2
    weight["data"] = 1.
    weight["mc"] = scaleFactor 
    for input in ["mc","data"]:
        icolor=1
        for id in electronIdLevels:
            print icolor
            canvas.SetLogy(0)
            if key == "ele_pt":
                histogram[input,key + "_" + id].GetXaxis().SetRangeUser(0.,80.)
            if key == "mee":
                histogram[input,key + "_" + id].GetXaxis().SetRangeUser(40.,140.)
            if key == "mee60":
                histogram[input,key + "_" + id].GetXaxis().SetRangeUser(60.,120.)
                
            histogram[input,key + "_" + id].Scale(weight[input])
            histogram[input,key + "_" + id].SetMarkerStyle(20)
            histogram[input,key + "_" + id].SetMarkerSize(1.2)
            histogram[input,key + "_" + id].SetMarkerColor(icolor)
            histogram[input,key + "_" + id].SetLineWidth(2)
            histogram[input,key + "_" + id].SetLineColor(icolor)
            if id == "" and input == "mc":
                histogram[input,key + "_" + id].SetMinimum(0.)
                histogram[input,key + "_" + id].SetMaximum(maximum)
                histogram[input,key + "_" + id].Draw("HIST")
            elif input == "mc":
                histogram[input,key + "_" + id].Draw("HISTSAME")
            elif input == "data":
                print                 histogram[input,key + "_" + id].Integral()
                leg.AddEntry(histogram[input,key + "_" + id],id,"PL")
                histogram[input,key + "_" + id].Draw("ESAME")
            icolor=icolor+1
    canvas.Update()
    leg.Draw()
    canvas.SaveAs("plots/"+ key + ".png")
            
        # Plotting also in Logarithmic scale
    for input in ["mc","data"]:
        icolor=1
        for id in electronIdLevels:
#            if key == "ele_pt":
#                histogram[input,key + "_" + id].GetXaxis().SetRangeUser(0.,80.)
            histogram[input,key + "_" + id].SetMinimum(0.1)
            histogram[input,key + "_" + id].SetMaximum(histogram[input,key + "_" + id].GetMaximum()*30.)
            canvas.SetLogy(1)
            if id == "" and input == "mc":
                histogram[input,key + "_" + id].Draw("HIST")
            elif input=="mc":
                histogram[input,key + "_" + id].Draw("HISTSAME")
            elif input=="data":
                histogram[input,key + "_" + id].Draw("ESAME")
            icolor = icolor + 1
    leg.Draw()
    canvas.SaveAs("plots/"+ key + "_log.png")

# *----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*
# Now pleotting also N-1 plot (linear and log) for data and MC
# *----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*----*

canvas = TCanvas("c","c",1)
for key in ["sigmaIetaIeta","HOverE","CombinedIso","ExpMissHits"]:
    print "Drawing " + key
    for type in ["NMinusOne_","Probe_"]:
        for detector in ["EB","EE"]:
            for input in ["mc","data"]:
                icolor=len(electronIdLevelsForNMinusOne)
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
                    histogram[input,histoKey].SetLineWidth(2)
                    histogram[input,histoKey].SetLineColor(icolor)
                    if id == "simpleEleId70cIso" and input == "mc":
                        histogram[input,histoKey].DrawNormalized("HIST")
                    elif input == "mc":
                        histogram[input,histoKey].DrawNormalized("HISTSAME")
                    elif input == "data":
                        histogram[input,histoKey].DrawNormalized("ESAME")
                    canvas.Update()
                    icolor=icolor-1
            canvas.SaveAs("plots/"+ key + "_" + detector + "_"+ type+".png")

    for type in ["NMinusOne_","Probe_"]:
        for detector in ["EB","EE"]:
            for input in ["mc","data"]:
                icolor=len(electronIdLevelsForNMinusOne)
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
                    histogram[input,histoKey].SetMinimum(0.00001)
                    histogram[input,histoKey].SetMaximum(histogram[input,histoKey].GetMaximum()*10.)
                    histogram[input,histoKey].SetLineColor(icolor)
                    if id == "simpleEleId70cIso" and input == "mc":
                        histogram[input,histoKey].DrawNormalized("HIST")
                    elif input == "mc":
                        histogram[input,histoKey].DrawNormalized("HISTSAME")
                    elif input == "data":
                        histogram[input,histoKey].DrawNormalized("ESAME")

                    canvas.Update()
                    icolor=icolor-1
            canvas.SaveAs("plots/"+ key + "_" + detector + "_"+ type+ "_log.png")

            
