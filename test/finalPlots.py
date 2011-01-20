#! /usr/bin/env python

from ROOT import *

import sys
import os
import math

gROOT.SetStyle("Plain")
gROOT.SetBatch()

inputfile = TFile("electronDistributions.root")
inputfile.cd()

electronIdLevels= [ "","simpleEleId95cIso","simpleEleId80cIso","simpleEleId70cIso" ]
electronIdLevelsForNMinusOne= [ "simpleEleId70cIso","simpleEleId80cIso","simpleEleId95cIso","" ]

histogram = {}

canvas = TCanvas("c","c",1)

#Get all histograms
for id in electronIdLevels:
    histogram["ele_pt_" + id] = gDirectory.Get("ele_pt_" + id)
    histogram["ele_scpt_" + id] = gDirectory.Get("ele_scpt_" + id)
    histogram["ele_eta_" + id] = gDirectory.Get("ele_eta_" + id)
    histogram["ele_phi_"+ id] = gDirectory.Get("ele_phi_" + id)
    histogram["met_"+ id] = gDirectory.Get("met_" + id)
    histogram["mt_"+ id] = gDirectory.Get("mt_" + id)
    histogram["mee_"+ id] = gDirectory.Get("mee_" + id)
    for detector in ["EB","EE"]:
        if id != "":
            histogram["ele_"+ detector + "_sigmaIetaIetaNMinusOne_"+ id] = gDirectory.Get("ele_"+ detector + "_sigmaIetaIetaNMinusOne_"+id)
            histogram["ele_"+ detector + "_HOverENMinusOne_"+ id] = gDirectory.Get("ele_"+ detector + "_HOverENMinusOne_"+id)
            histogram["ele_"+ detector + "_CombinedIsoNMinusOne_"+ id] = gDirectory.Get("ele_"+ detector + "_CombinedIsoNMinusOne_"+id)
            histogram["ele_"+ detector + "_ExpMissHitsNMinusOne_"+ id] = gDirectory.Get("ele_"+ detector + "_ExpMissHitsNMinusOne_"+id)
        else:
            histogram["ele_"+ detector + "_sigmaIetaIeta"] = gDirectory.Get("ele_"+ detector + "_sigmaIetaIeta")
            histogram["ele_"+ detector + "_HOverE"] = gDirectory.Get("ele_"+ detector + "_HOverE")
            histogram["ele_"+ detector + "_CombinedIso"] = gDirectory.Get("ele_"+ detector + "_CombinedIso")
            histogram["ele_"+ detector + "_ExpMissHits"] = gDirectory.Get("ele_"+ detector + "_ExpMissHits")

os.system("mkdir -p plots")

for key in ["ele_pt","ele_scpt","ele_eta","ele_phi","met","mt","mee"]:
    print "Drawing " + key
    icolor=1
    for id in electronIdLevels:
        canvas.SetLogy(0)
        if key == "ele_pt":
            histogram[key + "_" + id].GetXaxis().SetRangeUser(0.,80.)
        histogram[key + "_" + id].SetMinimum(0.)
        histogram[key + "_" + id].SetMarkerStyle(20)
        histogram[key + "_" + id].SetMarkerSize(1.)
        histogram[key + "_" + id].SetMarkerColor(icolor)
#        histogram[key + "_" + id].SetLineWidth()
        histogram[key + "_" + id].SetLineColor(icolor)
        if id == "":
            histogram[key + "_" + id].Draw("E")
        else:
            histogram[key + "_" + id].Draw("ESAME")
        icolor=icolor+1
    canvas.SaveAs("plots/"+ key + ".png")

    # Plotting also in Logarithmic scale
    icolor=1
    for id in electronIdLevels:
        if key == "ele_pt":
            histogram[key + "_" + id].GetXaxis().SetRangeUser(0.,80.)
        histogram[key + "_" + id].SetMinimum(0.1)
        histogram[key + "_" + id].SetMaximum(histogram[key + "_" + id].GetMaximum()*10.)
        canvas.SetLogy(1)
        if id == "":
            histogram[key + "_" + id].Draw("E")
        else:
            histogram[key + "_" + id].Draw("ESAME")
        icolor = icolor + 1
    canvas.SaveAs("plots/"+ key + "_log.png")

for key in ["sigmaIetaIeta","HOverE","CombinedIso","ExpMissHits"]:
    print "Drawing " + key
    for detector in ["EB","EE"]:
        icolor=4
        for id in electronIdLevelsForNMinusOne:
            canvas.SetLogy(0)

            if id != "":
                histoKey="ele_"+ detector + "_"+ key + "NMinusOne_"+ id
            else:
                histoKey="ele_"+ detector + "_"+ key
                
            histogram[histoKey].SetMinimum(0.)
            histogram[histoKey].SetMarkerStyle(20)
            histogram[histoKey].SetMarkerSize(1.)
            histogram[histoKey].SetMarkerColor(icolor)
            #        histogram[histoKey].SetLineWidth()
            histogram[histoKey].SetLineColor(icolor)
            if id == "simpleEleId70cIso":
                histogram[histoKey].DrawNormalized("E")
            else:
                histogram[histoKey].DrawNormalized("ESAME")
            icolor=icolor-1
        canvas.SaveAs("plots/"+ key + "_" + detector + ".png")

        icolor=4
        for id in electronIdLevelsForNMinusOne:
            canvas.SetLogy(1)

            if id != "":
                histoKey="ele_"+ detector + "_"+ key + "NMinusOne_"+ id
            else:
                histoKey="ele_"+ detector + "_"+ key
                

            histogram[histoKey].SetMarkerStyle(20)
            histogram[histoKey].SetMarkerSize(1.)
            histogram[histoKey].SetMarkerColor(icolor)
            #        histogram[histoKey].SetLineWidth()
            histogram[histoKey].SetLineColor(icolor)
            if id == "simpleEleId70cIso":
                histogram[histoKey].DrawNormalized("E")
            else:
                histogram[histoKey].DrawNormalized("ESAME")
            histogram[histoKey].SetMinimum(0.00001)
            histogram[histoKey].SetMaximum(histogram[histoKey].GetMaximum()*10.)
            canvas.Update()
            icolor=icolor-1
        canvas.SaveAs("plots/"+ key + "_" + detector + "_log.png")

            
