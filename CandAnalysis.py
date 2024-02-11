print("\n^^^^^^^^\tY(2S)->µµ + Phi->KK SPECTRUM ANALYSER\t^^^^^^^^\n\nImporting modules...")

import ROOT
ROOT.gROOT.SetBatch(True)

from time import time
import os
from sys import argv
from definitions import CandTreeDefinitions
import declarations
import pandas as pd
import matplotlib.pyplot as plt
import numpy


with open('Y2SPhiRun2List.txt') as f:
    allFiles = f.readlines()
    f.close()

for i in range(len(allFiles)):			
    allFiles[i] = allFiles[i].replace("\n", "")

#	TREES READING - using command line arguments
#	- no command line arguments: analysis of all data and store plots in a "test" directory
# 	- argv is an integer: analysis of first n root files, in a "test" directory
#	- argv is a string: analysis of all data and store plots in directory named <argv_value>


d = "test"
try: 
	sample = allFiles[:int(argv[1])] 
except IndexError:
	sample = allFiles[:]
except ValueError:
	sample = allFiles[:]
	d = argv[1]
											
if not os.path.isdir(f"./{d}/"):			
	os.system(f"mkdir {d}")							

# create a root file to store and edit plots, if necessary											
p = ROOT.TFile.Open(d+"/cand.root","RECREATE")		
												

#	PLOT TEMPLATE
def cprint (hist, name, opt="same", stats=False, f = p):
	title= "Candidate "+name
	c = ROOT.TCanvas(title, title)
	legend = ROOT.TLegend()
	multiple_hists = len(hist) > 1
		
	for i,h in enumerate(hist): 
		h.SetStats(stats)
		h.SetLineColor(i+1)
		h.Draw(opt)
		
		if multiple_hists:
			legend.AddEntry(h.GetPtr(), h.GetName(), "l")
	
	if multiple_hists: legend.Draw("same")
	
	c.SaveAs(d+"/"+name+".pdf")
	f.WriteObject(c, name)
	os.system(f"xdg-open {d}/"+name+".pdf")

# creating datasets for the analysis		
print("Creating Dataset...")

dataY2SKK = CandTreeDefinitions\
(ROOT.RDataFrame("rootuple/CandidateTree",sample))



#KK filters
dataKKrs = dataY2SKK.Filter("track1_charge * track2_charge < 0")
dataKKws = dataY2SKK.Filter("track1_charge * track2_charge > 0")

#############################################################
#					CANDIDATES KK DECAY						#
#############################################################


#	DEFINING DIFFERENT PLOTS
"""
#	MASS PLOT OF CANDIDATE	##############################
	
	#candidate mass
hist = dataY2SKK.Histo1D(("Candidate Mass", "Candidate Mass; m(K^{L})+m(K^{S})+ m(Y(2S))[GeV/c^{2}];Counts", 500, 11., 13.6), "candidate_vMass")
cprint(hist, d+"/candidateMass")

	#dimuon mass
hist = dataY2SKK.Histo1D(("Dimuon Mass", "Candidate Mass; m[GeV/c^{2}];Counts", 500, 9.4, 10.5), "dimuon_vMass")
cprint(hist, d+"/dimuonMass", stats=True)
"""

#	LEADING KAON PT PLOT

def kaonL_pt ():
	hist = dataKKrs.Histo1D(("Leading Kaon", "Leading Kaon p_{T};p_{T}(K^{L}) [GeV];Counts", 200, 0., 3.), "trackL_pT")
	cprint([hist], "LeadingK", stats=False)

# 	SOFT KAON PT PLOT
def kaonS_pt():
	hist = dataKKrs.Histo1D(("Soft Kaon", "Soft Kaon p_{T};p_{T}(K^{soft}) [GeV];Counts", 200, 0., 3.), "trackS_pT")
	cprint([hist], "SoftK", stats=False)
	
	
def prob():
	#dimuon = dataKKrs.Histo1D(("Dimuon vProb", "Dimuon vProb vs Candidate vProb;p;Counts", 200, 0., 1.), "dimuon_vProb")
	candidate = dataKKrs.Histo1D(("Candidate vProb", "Candidate vertex probability;P_{vtx}(Y(2S)KK);Counts", 500, 0., 1.), "candidate_vProb")
	
	cprint([candidate], "vProb")
	
	
	
	
	
	
	
#	KK INVARIANT MASS PLOT, WITH CUTS	#################
def m_kk():

	dataset = dataKKrs

		#borders definitions
	LValues = [0., 0.7, 0.8, 0.8, 0.9, 0.9, 1.0, 1.0, 1.1, 1.1, 1.2,
				1.2, 1.3, 1.3, 1.4, 1.4]
	SValues = [0., 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0, 
				1.0, 1.1, 1.1, 1.2, 1.2]
		
		#create array of histograms with all cuts  
	hists = []
	to_fit = []
	
	c0 = ROOT.TCanvas()

	legend = ROOT.TLegend(0.46, 0.15, 0.89, 0.50)
	
	for i,j in zip(LValues, SValues):
		dataset = dataset.Filter(f"trackL_pT > {i} & trackS_pT > {j}")
		
		hists.append(dataset\
		.Histo1D(("KK invariant mass", "#phi #rightarrow K^{+}K^{-};m(KK) [GeV];Counts", 100, 0.99, 1.06), "ditrack_mass"))
	
	for i, hist in enumerate(hists):
		if i%2 == 0: 
			hist.SetLineColor(i+1)	
			hist.GetYaxis().SetRangeUser(1000, 1000000)				
			legend.AddEntry(hist.GetPtr(), "p_{T} lead > "+str(LValues[i])+"; p_{T} soft > "+str(SValues[i]), "l")	
			hist.Draw("same")
	
	legend.Draw("")


	c0.Draw("")
	c0.SaveAs(d+"/abs_MassKK.pdf")
	p.WriteObject(c0, "MassKK_abs")
	os.system(f"xdg-open {d}/abs_MassKK.pdf")
	
	
	legend.Clear()
	c1 = ROOT.TCanvas()
	hists[0].SetStats(0)	#no stats

		#overlap plots with cuts
	for i, hist in enumerate(hists):
		
		hist.Scale(1./hist.Integral())	#normalization of plots
		hist.GetYaxis().SetRangeUser(0, 0.015)	#range of y axis
		hist.SetLineColor(i+1)					
		legend.AddEntry(hist.GetPtr(), "p_{T} lead > "+str(LValues[i])+"; p_{T} soft > "+str(SValues[i]), "l")	
		hist.Draw("same")	#non si distinguono molto bene i dati
		
	legend.Draw("")


	c1.Draw("")
	c1.SaveAs(d+"/MassKK.pdf")
	p.WriteObject(c1, "MassKK_norm")
	os.system(f"xdg-open {d}/MassKK.pdf")
		
		
		 #only last cut
	#hists[-1].SetTitle("#phi #rightarrow K^{+}K^{-}: p_{T}(K^{L}) > 1.0, p_{T}(K^{S}) > 0.8")
#	cprint(hists[-1], d+"/MassKK_lastcut")



def m_kk_quality ():
	
	dataset = dataKKrs
	
	hists = [ dataset.Histo1D(("KK invariant mass", "#phi #rightarrow K^{+}K^{-};m(KK) [GeV];Counts", 100, 0.99, 1.06), "ditrack_mass") ]

	
#	dfkk = dataKKrs.Filter(f"trackL_pT > {ptL} & trackS_pT > {ptS}")\
#	.Filter("ditrack_mass > 1.0 & ditrack_mass < 1.04")\
#	.Filter("track1_pvAssocQ + track2_pvAssocQ > 11")\
#	.AsNumpy(columns=["ditrack_mass"])

#	kkroodata = ROOT.RooDataSet.from_numpy({"ditrack_mass": dfkk["ditrack_mass"]}, [kkmass])
#	kkroohist = kkroodata.binnedClone()
	
	
	c0 = ROOT.TCanvas()

	legend = ROOT.TLegend(0.46, 0.15, 0.89, 0.50)
	
	hists[0].SetLineColor(1)
	legend.AddEntry(hists[0].GetPtr(), "no pv_Q cut", "l")	
	
	hists[0].Draw("")
	
	for i in range(9,14):
		dataset = dataset.Filter(f"track1_pvAssocQ + track2_pvAssocQ > {i}")
	
		hists.append(dataset\
		.Histo1D(("KK invariant mass", "#phi #rightarrow K^{+}K^{-};m(KK) [GeV];Counts", 100, 0.99, 1.06), "ditrack_mass"))
		
	for i, hist in enumerate(hists[1:], start=1):
		hist.SetLineColor(i+1)	
		#hist.GetYaxis().SetRangeUser(1000, 1000000)
						
		legend.AddEntry(hist.GetPtr(), "pv_Q_{track1} + pv_Q_{track2} > "+f"{i+8}", "l")	
		hist.Draw("same")

	
	legend.Draw("")
	
	c0.Draw("")
	c0.SaveAs(d+"/abs_MassKK_qual.pdf")
	p.WriteObject(c0, "qual_MassKK_abs")
	os.system(f"xdg-open {d}/abs_MassKK_qual.pdf")
	
	
	c1 = ROOT.TCanvas()
		#overlap plots with cuts
	for i, hist in enumerate(hists):
		
		hist.Scale(1./hist.Integral())	#normalization of plots
		#hist.GetYaxis().SetRangeUser(0, 0.015)	#range of y axis
		hist.SetLineColor(i+1)					
		#legend.AddEntry(hist.GetPtr(), "p_{T} lead > "+str(LValues[i])+"; p_{T} soft > "+str(SValues[i]), "l")	
		hist.Draw("same")	#non si distinguono molto bene i dati
	
	legend.Draw("")
	
	c1.Draw("")
	c1.SaveAs(d+"/norm_MassKK_qual.pdf")
	p.WriteObject(c1, "qual_MassKK_norm")
	os.system(f"xdg-open {d}/norm_MassKK_qual.pdf")
	
#	quadrature function
def quadrature (a,b):
	return pow( (pow(a,2) + pow(b,2)), 0.5 )	

#	FIT OF KK_PT CUT	#####################
#	take the cut pt_kL > 1.2 and pt_kS > 1.0

def fit (dataset):
	# variable

	kkmass = ROOT.RooRealVar("ditrack_mass", "m(KK) [GeV]", 1.01, 1.03)
	
	# dataset
	dfkk = dataset.AsNumpy(columns=["ditrack_mass"])

	kkroodata = ROOT.RooDataSet.from_numpy({"ditrack_mass": dfkk["ditrack_mass"]}, [kkmass])
	kkroohist = kkroodata.binnedClone()
	
		#	Number of entries
	entries = kkroohist.sumEntries()
	
		#	Create frame
	phiframe = kkmass.frame(Title="Dikaon Candidate Mass")

	
		#	Signal parameters
	mean = ROOT.RooRealVar("#mu_{#phi}", "mean of gaussian", 1.02, 1.019, 1.021) #mean value, min value, max value
	sigma = ROOT.RooRealVar("#sigma_{#phi}", "resolution", 0.00125, 0.001, 0.002) 
	width = ROOT.RooRealVar("#Gamma_{#phi}", "width", 0.00439, 0.001, 0.006)
	
	sigma.setConstant(1)

		#	Chebyshev coefficients
	f0 = ROOT.RooRealVar("f0", "f0", 0.1, 0., 1.)	#-5., -50.,50.
	f1 = ROOT.RooRealVar("f1", "f1", -0.0001, -0.01, 0.)	#5., -50.,50.
	#f2 = ROOT.RooRealVar("f2", "f2", -1.1, -50.,50.)

	
		#	Number of events
	Nbkg = ROOT.RooRealVar("N_{bkg}", "N bkg events", 50, 0., entries)
	Nsig = ROOT.RooRealVar("N_{sig}", "N sig events", 50, 0., entries)

		#	Model Functions
	sig = ROOT.RooVoigtian("signal", "signal", kkmass, mean, width, sigma)
	bkg = ROOT.RooChebychev("bkg", "Background", kkmass, [f0, f1]) 

		#	Total model
	model = ROOT.RooAddPdf("model", "voigt+cheb", [bkg, sig], [Nbkg, Nsig])
	
	
		# fit	
	model.fitTo(kkroohist)
	
		# ploton
	kkroohist.plotOn(phiframe)
	model.plotOn(phiframe, Components={sig}, LineStyle=":", LineColor="r")
	model.plotOn(phiframe, Components={bkg}, LineStyle=":", LineColor="g")
	model.plotOn(phiframe) # By default only fitted range is shown
	model.paramOn(phiframe, ROOT.RooFit.Parameters([mean, width, f0, f1]), ROOT.RooFit.Layout(0.65, 0.9, 0.9))

	xmin = mean.getVal() - 2*quadrature(width.getVal()/2, sigma.getVal())
	xmax = mean.getVal() + 2*quadrature(width.getVal()/2, sigma.getVal())

	
	kkmass.setRange("window", xmin,xmax)
	S = sig.createIntegral(kkmass, ROOT.RooFit.NormSet(kkmass), Range="window").getVal() * Nsig.getVal() 
	B = bkg.createIntegral(kkmass,  ROOT.RooFit.NormSet(kkmass), Range="window").getVal() * Nbkg.getVal()
	T = model.createIntegral(kkmass,  ROOT.RooFit.NormSet(kkmass), Range="window").getVal() * entries

	# PURITY of the signal
	P = S/T
	S2 = P*S
	
	print (f"\nthe evaluated bkg is {B}\n")
	print (f"the evaluated sig is {S}\n")
	print (f"the total integral are {T}\n")
	print (f"the entries are {entries}\n")
	print (f"\nPurity = {P}")
	print (f"\nSignificance = {S2}")
	phiframe.Draw()
	p.WriteObject(phiframe,"phi_mass_fit")
	pull = phiframe.pullHist()

	ymax = max( pull.GetY() )
	ymin = min( pull.GetY() )

	if ymax > 7 or ymin < -7: 
		os.system(f"rootbrowse {d}/cand.root")
		exit()
	
	return P
	

#______________________________________________________________________
def m_kk_fit(ptL = 1.2, ptS = 1.0, Q =11):	#ptL = 1.2, ptS = 1.0, Q =11
	dfkk = dataKKrs.Filter(f"trackL_pT > {ptL} & trackS_pT > {ptS}")\
	.Filter("ditrack_mass > 1.0 & ditrack_mass < 1.04")\
	.Filter(f"track1_pvAssocQ + track2_pvAssocQ > {Q}")\
	.AsNumpy(columns=["ditrack_mass"])
		
		#variable
	kkmass = ROOT.RooRealVar("ditrack_mass", "m(KK) [GeV]", 1.01, 1.03)

	kkroodata = ROOT.RooDataSet.from_numpy({"ditrack_mass": dfkk["ditrack_mass"]}, [kkmass])
	kkroohist = kkroodata.binnedClone()
	
		#	Number of entries
	entries = kkroohist.sumEntries()
	
		#	Create frame
	phiframe = kkmass.frame(Title="Dikaon Candidate Mass")

	
		#	Signal parameters
	mean = ROOT.RooRealVar("#mu_{#phi}", "mean of gaussian", 1.02, 1.019, 1.021) #mean value, min value, max value
	sigma = ROOT.RooRealVar("#sigma_{#phi}", "resolution", 0.00125, 0.001, 0.002) 
	width = ROOT.RooRealVar("#Gamma_{#phi}", "width", 0.00547, 0.002, 0.006)
	
	sigma.setConstant(1)		# fixed res from MC????????
	width.setConstant(1)
		#	Chebyshev coefficients
	f0 = ROOT.RooRealVar("f0", "f0", .75, 0., 5.)	#-5., -50.,50.
	f1 = ROOT.RooRealVar("f1", "f1", -0.0001, -0.05, 0.)	#5., -50.,50.
	#f2 = ROOT.RooRealVar("f2", "f2", 0.001, -0.1,0.1)

	
		#	Number of events
	Nbkg = ROOT.RooRealVar("N_{bkg}", "N bkg events", 50, 0., entries)
	Nsig = ROOT.RooRealVar("N_{sig}", "N sig events", 50, 0., entries)

		#	Model Functions
	sig = ROOT.RooVoigtian("signal", "signal", kkmass, mean, width, sigma)
	bkg = ROOT.RooChebychev("bkg", "Background", kkmass, [f0, f1]) 

		#	Total model
	model = ROOT.RooAddPdf("model", "voigt+cheb", [bkg, sig], [Nbkg, Nsig])
		
	model.fitTo(kkroohist)
		
		#print
	c0 = ROOT.TCanvas("canvas0", "canvas0", 1200, 600)

	kkroohist.plotOn(phiframe)
	model.plotOn(phiframe, Components={sig}, LineStyle=":", LineColor="r")
	model.plotOn(phiframe, Components={bkg}, LineStyle=":", LineColor="g")
	model.plotOn(phiframe) # By default only fitted range is shown
	model.paramOn(phiframe, ROOT.RooFit.Parameters([mean, width, Nbkg, Nsig, f0,f1]), ROOT.RooFit.Layout(0.65, 0.9, 0.9))

	xmin = mean.getVal() - 2*quadrature(width.getVal()/2, sigma.getVal())
	xmax = mean.getVal() + 2*quadrature(width.getVal()/2, sigma.getVal())

	
	kkmass.setRange("window", xmin,xmax)
	S = sig.createIntegral(kkmass, ROOT.RooFit.NormSet(kkmass), Range="window").getVal() * Nsig.getVal() 
	B = bkg.createIntegral(kkmass,  ROOT.RooFit.NormSet(kkmass), Range="window").getVal() * Nbkg.getVal()
	T = model.createIntegral(kkmass,  ROOT.RooFit.NormSet(kkmass), Range="window").getVal() * entries

	# PURITY of the signal
	P = S/T
	S2 = P*S
	
	print (f"\nthe evaluated bkg is {B}\n")
	print (f"the evaluated sig is {S}\n")
	print (f"the total integral are {T}\n")
	print (f"the entries are {entries}\n")
	print (f"\nPurity = {P}")
	print (f"\nSignificance = {S2}\n")
	line0 = ROOT.TLine(xmin, 0., xmin, 10000)
	line1 = ROOT.TLine(xmax, 0., xmax, 10000)
	line0.SetLineStyle(2)
	line0.SetLineColor(7)
	line0.SetLineWidth(4)
	line1.SetLineStyle(2)
	line1.SetLineColor(7)
	line1.SetLineWidth(4)

	c0.Divide(1,2)
	
	c0.cd(1)
	phiframe.Draw()
	line0.Draw("same")
	line1.Draw("same")
	
	c0.cd(2)
	pull = phiframe.pullHist()
	pull.GetYaxis().SetTitle("#chi^{2}")
	pull.SetTitle("")
	pull.Draw()

	

	#c0.Draw()
	
	p.WriteObject(phiframe,"phi_mass_fit")
	c0.SaveAs(d+"/PhiMassPlot.pdf")
	os.system(f"xdg-open {d}/PhiMassPlot.pdf")
	
	ymax = max( pull.GetY() )
	ymin = min( pull.GetY() )

	if ymax > 7 or ymin < -7: 
		os.system(f"rootbrowse {d}/cand.root")
		exit()
		# End timer
	end = time()

		# Calculate elapsed time
	elapsed = end - start
	print("\nTime for Phi mass plot: ", elapsed, "\n") 
	
	return P

#_______________________________________________________________________
# to test purity i have to do multiple fits with all compbinations of cuts:
# cuts in pt and cuts in vertex quality.

def purity ():

	c = ROOT.TCanvas()
	
	LValues = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]
	pvQ = [8., 9., 10., 11., 12., 13.]
	
	P1 = [0.08308301722436447,
		  0.11495195386832995,
		  0.12296124063941502,
		  0.1418172453034686,
		  0.15627704995988054,
		  0.168422194274125,
		  0.1768962788277655]
	P2 = [0.1553377081555507,
	 	  0.15546651659083266,
	 	  0.15562173896399706,
	 	  0.15636357553645472,
	 	  0.16132010661155846,
	 	  0.1629469241426943]
	
	
	pt_cuts = ROOT.TGraph()
	quality_cuts = ROOT.TGraph()
	#quality_cuts = ROOT.TH1D("purity vs quality","purity vs quality",len(pvQ),8.,14.)	
	

	 	  
#	plot = ROOT.TH2D("purity vs pt","purity vs pt",len(LValues),0.8,1.6, len(pvQ), 8.,14.)		
#	dataset = dataKKrs.Filter("ditrack_mass > 1.01 & ditrack_mass < 1.03")

	for i,P in zip(LValues,P1):
		
		print(f"Plot {i+1}/14 ....")
	#	data = dataset.Filter((f"trackL_pT > {pt} & trackS_pT > {pt-0.2}"))
		pt_cuts.AddPoint(i, P)
		
	for j,P in zip(pvQ,P2):
			
		print(f"Plot {j+9}/14 ....")
	#	print(f"Plot {6*i+j+1}/42 ....")
	#	data = dataset.Filter(f"track1_pvAssocQ + track2_pvAssocQ > {q}")
		quality_cuts.AddPoint(j, P) #SetBinContent(j, P2[j])fit(data)
	#	plot.SetBinContent(i,j, fit(data))
#		
	pt_cuts.SetStats(0)	
	quality_cuts.SetStats(0)
	#plot.SetStats(0)
	
	c.Divide(2,1)
	c.cd(1)	
	pt_cuts.Draw()
	
	c.cd(2)	
	quality_cuts.Draw()
#	plot.Draw("colz")
	c.Draw()
	
	p.WriteObject(c,"pt vs quality")
	c.SaveAs(d+"/pt_vs_qual.pdf")
	os.system(f"xdg-open {d}/pt_vs_qual.pdf")














#_______________________________________________________________________
ymumu_filter= "dimuon_mass > 9.881 & dimuon_mass < 10.147 & dimuon_pT > 12 & "
phiKKSelection = '''candidate_vProb > 0.1 &
ditrack_mass > 1.01438 &
ditrack_mass < 1.02449 &
trackL_pT > 1.2 & 
trackS_pT > 1.0'''
#candProb = " & candidate_vProb > 0.05"
cand_constr = "candidate_vMass > 11.11 & candidate_vMass < 11.3 & "
quality_filter = " & track1_pvAssocQ + track2_pvAssocQ > 11"


#binning = 250
#edge = [10.8, 13.7]

binning = 100
edge = [11.11, 11.3]

#	definitions of the selected candidates
candSelectionRS = dataKKrs.Filter(cand_constr +
								  ymumu_filter+
								  phiKKSelection +
								  quality_filter
								  
									)

candSelectionWS = dataKKws.Filter(ymumu_filter+phiKKSelection
									+ quality_filter)
#_____________________________________________________________________________
# ditrack mass: interval of phi rest mass
# vProb: vertex probability

#	MASS CONSTRAIN ON PI PI (identify if Y(3S)->Y(2S)pipi)
def mass_constr():

	hist0 = dataKKrs.Filter(ymumu_filter 
							+ "candidate_vProb > 0.1 &"
							+ "trackL_pT > 0.4 & trackS_pT > 0.4"
							).Histo1D(("MuMuKK cands", "Y(3S)#rightarrow Y(2S)(#rightarrow #mu^{+}#mu^{-})(#pi^{+}#pi^{-});m(#mu#mu#pi#pi) - m(#mu#mu) + m^{PDG}(Y(2S)) [GeV];Counts", binning, 10.1, 13.), "mumupipi_MC")
	
	#hist1 = dataKKrs.Filter(ymumu_filter + phiKKSelection).Histo1D(("MuMuKK cands", "Y(3S)#rightarrow Y(2S)(#rightarrow #mu^{+}#mu^{-})(#pi^{+}#pi^{-});m(#mu#mu#pi#pi) [GeV];Counts", binning, 10.0, 13.), "mumupipi_mass")
	
	c = ROOT.TCanvas()
	
	hist0.SetLineColor(1)
	#hist1.SetLineColor(2)
	
	hist0.Draw()
	#hist1.Draw("same")
	
		
	c.SaveAs("mumupipiMC.pdf")
	os.system("xdg-open mumupipiMC.pdf")


# PHI CANDIDATE PLOT	
#ws = wrong sign (take K+K+ and K-K-)
#compare the selection distribution with one random distribution


	#Γ φ = 0.00446 ± 0.00018
	#μφ  = 1.019445 ± 0.000036
	
	#	Y mass cuts: µ ± 2 sigma 
	#	φ mass cuts: µ ± 2 Sigma


# new dir with cuts written on file

tagli = ymumu_filter+phiKKSelection
tagli = tagli.replace(" & ", "\n")
os.system(f"echo \"{tagli}\nbinning: {binning}\" > {d}/tagli.txt")

#histograms
def mumukk(zoom = False):

	hist0 = candSelectionRS.Histo1D(("MuMuKK cands", "Y(2S)(#rightarrow #mu^{+}#mu^{-})#phi(#rightarrow K^{+}K^{-});m(#mu#muKK) - m(#mu#mu) + m^{PDG}(Y) [GeV];Counts RS", binning, edge[0], edge[1]), "candidate_vMass")
	
	hist1 = candSelectionWS.Histo1D(("MuMuKK cands", "Y(2S)(#rightarrow #mu^{+}#mu^{-})#phi(#rightarrow K^{+}K^{-});m(#mu#muKK) - m(#mu#mu) + m^{PDG}(Y) [GeV];Counts RS", binning, edge[0], edge[1]), "candidate_vMass")

	c0 = ROOT.TCanvas()

	hist0.SetStats(0)

	hist0.SetLineColor(1)
	hist1.SetLineColor(2)
	
	legend = ROOT.TLegend(0.7, 0.1, 0.89, 0.3) #(xmin,ymin,xmax,ymax)

	legend.AddEntry(hist0.GetPtr(), "RS kaons", "l")
	legend.AddEntry(hist1.GetPtr(), "WS kaons (norm)", "l")

	hist1.Scale(hist0.Integral()/hist1.Integral())
	        
	hist0.Draw("")
	hist1.Draw("same")
	legend.Draw("")
	c0.Draw("")
	c0.SaveAs(d+"/Candidate.pdf")
	os.system(f"xdg-open {d}/Candidate.pdf")
	
	p.WriteObject(c0,"candidate")
	
	if zoom:

		hist0.GetXaxis().SetRangeUser(11., 11.6)
		hist1.GetXaxis().SetRangeUser(11., 11.6)
		
		hist2 = hist0
		hist3 = hist1
		
		
		hist3.Scale(hist2.Integral()/hist3.Integral())
		
		hist2.SetTitle("Y(2S)(#rightarrow #mu^{+}#mu^{-})#phi(#rightarrow K^{+}K^{-}) (Zoom)")
		
		hist2.Draw("")
		hist3.Draw("same")
		
		legend.Draw("")
		c0.Draw("")
		c0.SaveAs(d+"/PhiCandidateZoom.pdf")
		
		p.WriteObject(c0,"phi_candidate(zoom)")
		os.system(f"xdg-open {d}/PhiCandidateZoom.pdf")
#_________________________________________________________________________

#	NUMBER OF CANDIDATES PER EVENT
#	arbitration of event multiplicity in Jupiter
def Ncand():
	
	candSelectionRS = dataKKrs.Filter(ymumu_filter+phiKKSelection + quality_filter)


	dfcands = candSelectionRS.AsNumpy(columns=["event", "run", "candidate_pT",
									 "candidate_vProb", "candidate_vMass"])
	
	candDF = pd.DataFrame(dfcands)
	
	
	candxEvent = candDF.groupby("event")["candidate_vProb"]
	NcandxEvent = candxEvent.count()
	
	
	cleanData = candDF.loc[candDF.groupby("event")["candidate_pT"].idxmax()]
	#hist = candxEvent.plot(kind="hist", logy=True)
	
	plt.hist(NcandxEvent, bins=range(1,12), align="left", edgecolor='black')
	plt.title("Candidate multiplicity")
	plt.xlabel("Multiplicity")
	plt.ylabel("Counts")
	plt.yscale('log')

	plt.show()
	
	multiplicity_ratio = (NcandxEvent > 1).sum()/(NcandxEvent > 0).sum()
	print(f"\nmultiplicity ratio is {round(multiplicity_ratio*100, 2)}%")
	
	
##	plt.hist(cleanData)
##	plt.show()
#	
	mumukkmass = ROOT.RooRealVar("candidate_vMass", "m(#mu#muKK) [GeV]", 11.0, 11.6)
	mumukkroodata = ROOT.RooDataSet.from_numpy({"candidate_vMass": numpy.array(cleanData["candidate_vMass"])}, [mumukkmass])
	
	c = ROOT.TCanvas()
	
	mumukkframe = mumukkmass.frame(Title="YPhi Candidate Mass")
	mumukkroodata.plotOn(mumukkframe)
	mumukkframe.Draw()


	c.Draw()
	c.SaveAs(d+"/clean_mumukk.pdf")
	
	p.WriteObject(c,"clean_mumukk")
	os.system(f"xdg-open {d}/clean_mumukk.pdf")

	
	

#	
########################################################################
#	MENU

compute = {	"1" : kaonL_pt,
			"2" : kaonS_pt,
			"3" : m_kk,
			"4" : m_kk_fit,
			"5"	: mass_constr,	
			"6" : mumukk,		
			"7" : mumukk,
			"8" : Ncand,
			"9" : purity,
			"0" : prob,
			"q" : exit
		  }

lang=[]
while("q" not in lang):

	lang = input("\nSelect plots (Separate by spacing):\n"+
				 "1. Leading Kaon pt\n"+
				 "2. Soft Kaon pt\n"+
				 "3. KK invariant mass (with K_pt cuts)\n"+
				 "4. Fit KK invariant mass\n"+
				 "5. Mass Constraint for Y3S into Y2S pi pi\n"+
				 "6. Phi Candidate plot\n"+
				 "7. Phi candidate plot with Zoom\n"+
				 "8. Candidate multiplicity\n"+
	#			 "ENTER: 1-5 Plots\n"+
				 "Press \"q\" to EXIT.\n").split()
	
	if "q" not in lang:
		print("Processing...") 

		# Start timer
	start = time()

	if not lang:	#print all plots		
		for func in compute.values() : func()

	else:			#print only selected-by-key plots
		for i in lang: 
		
			while i not in compute.keys(): 	#in case of multiple misdigit
				i = input(f"\"{i}\" is not valid. Please insert a valid key:\n")
			
			if i == "7" : compute[i](zoom = True)
			elif i == "4" : 
				ptL = float(input("ptL = "))
				ptS = ptL-0.2
				#quality = float(input("Q = "))
				compute[i](ptL, ptS)
			else : compute[i]()

		# End timer
	end = time()

		# Calculate elapsed time
	elapsed = end - start
	print("\nComputing time: ", elapsed, "\n") 
	print("What to do next?")
	
#p.Close()



