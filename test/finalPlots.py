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

histogram = {}

canvas = TCanvas("c","c",1)

for id in electronIdLevels:
    histogram["ele_pt_" + id] = gDirectory.Get("ele_pt_" + id)
    histogram["ele_scpt_" + id] = gDirectory.Get("ele_scpt_" + id)
    histogram["ele_eta_" + id] = gDirectory.Get("ele_eta_" + id)
    histogram["ele_phi_"+ id] = gDirectory.Get("ele_phi_" + id)
    histogram["met_"+ id] = gDirectory.Get("met_" + id)
    histogram["mt_"+ id] = gDirectory.Get("mt_" + id)
    histogram["mee_"+ id] = gDirectory.Get("mee_" + id)

os.system("mkdir -p plots")


for key in histogram.keys():
    print "Drawing " + key
    canvas.SetLogy(0)
    histogram[key].SetMarkerStyle(20)
    histogram[key].SetMarkerSize(1.)
    histogram[key].Draw("E")
    canvas.SaveAs("plots/"+ key + ".png")
    histogram[key].SetMinimum(0.1)
    canvas.SetLogy(1)
    histogram[key].Draw("E")
    canvas.SaveAs("plots/"+ key + "_log.png")


