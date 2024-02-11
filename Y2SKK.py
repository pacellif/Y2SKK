print("\n^^^^^^^^\tY(2S)->µµ ANALYSER\t^^^^^^^^\n\nImporting modules...")


import ROOT 
import os
from time import time
from sys import argv
from definitions import UpsTreeDefinitions
import declarations

ROOT.gROOT.SetBatch(True)


#directory: /lustre/cms/store/user/vmastrap/MuMuKKRun2/Y2SKKfromSametV4
# comando scp per copiare macro in account remoto
# cambiare il path per aprire i file

#ROOT.ROOT.EnableImplicitMT(4)


#	FILE MANAGEMENT:	---------------------------------------------| 
#	- open .root dataset: it is possible to select a subset
#	- make custom directories in which store plots and .root files

with open('Y2SPhiRun2List.txt') as f:
    allFiles = f.readlines()
    f.close()

for i in range(len(allFiles)):			
    allFiles[i] = allFiles[i].replace("\n", "")


#		trees reading - using command line arguments

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

#	file root
fileroot = ROOT.TFile.Open(d+"/ups_mass.root","RECREATE")


#	plot template	-------------------------------------------------|
def cprint (hist, name, opt="same", stats=False, f = fileroot):
	title= "Y2S "+name
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
#------------------------------------------------------------------
def binCount (var, rangeName):
	varBinning = var.getBinning()
	a = var.getRange(rangeName)[0]
	b = var.getRange(rangeName)[1]
	nbins = varBinning.binNumber(b) - varBinning.binNumber(a)
	return nbins
#------------------------------------------------------------------	


print("Creating Dataset...")
dataY2S = UpsTreeDefinitions(ROOT.RDataFrame("rootuple/UpsTree",sample))



#	MUON TRANSVERSE MOMENTUM
def mu_pt():
	hist1 = dataY2S.Histo1D(("mu1_pt Distribution", "mu1_pt Distribution;p_{T}(#mu^{+});Counts", 500, 0., 40.), "mu1_pt")

	hist2 = dataY2S.Histo1D(("mu2_pt Distribution", "mu2_pt Distribution;p_{T}(#mu^{-});Counts", 500, 0., 40.), "mu2_pt")

		#other two variables - same results
	hist3 = dataY2S.Histo1D(("muonP_pT Distribution", "muonP_pT Distribution;p_{T}(#mu^{+});Counts", 500, 0., 40.), "muonP_pT")

	hist4 = dataY2S.Histo1D(("muonN_pT Distribution", "muonN_pT Distribution;p_{T}(#mu^{-});Counts", 500, 0., 40.), "muonN_pT")

	c = ROOT.TCanvas("Muon pT distribution", "Muon pT distribution")

	c.Divide(1,2)

	c.cd(1)
	hist1.Draw()


	c.cd(2)
	hist2.Draw()


	c.SaveAs(d+"/muonPTdistribution.pdf")
	fileroot.WriteObject(c,"MuonPT")
	os.system(f"xdg-open {d}/muonPTdistribution.pdf")

######	   MASS PLOT	#########

def m_Y2S():	#corretto

	upsilon = dataY2S.Filter("ups_vMass > 9.6 & ups_vMass < 10.3")\
	.Filter("ups_pT > 12")
	
	hist12 = upsilon\
	.Histo1D(("dimuon p_{T} > 12", "Y(2S) #rightarrow #mu^{+}#mu^{-};m(#mu^{+}#mu^{-}) [GeV];Counts", 500, 9.6, 10.3), "ups_vMass")
	
	hist15 = upsilon\
	.Filter("ups_pT > 15")\
	.Histo1D(("dimuon p_{T} > 15", "Y(2S) #rightarrow #mu^{+}#mu^{-};m(#mu^{+}#mu^{-}) [GeV];Counts", 500, 9.6, 10.3), "ups_vMass") 
	
	
	hist18 = upsilon\
			 .Filter("ups_pT > 18")\
	.Histo1D(("dimuon p_{T} > 18", "Y(2S) #rightarrow #mu^{+}#mu^{-};m(#mu^{+}#mu^{-}) [GeV];Counts", 500, 9.6, 10.3), "ups_vMass") 
	
	
	
	cprint([hist12,hist15,hist18], "YinvMass")
	
	


#	FIT		------------------------------------------------
	
def fit_Y2S():	
	
	upsmass = ROOT.RooRealVar("ups_vMass", "m(#mu#mu) [GeV]", 9.6, 10.3)
	upsmass.setBins(500)
	
	#alternatively
	#mumuroohist = dataY2S.Book(ROOT.std.move(ROOT.RooDataSetHelper("dataset", "Title of dataset", ROOT.RooArgSet(YMass))), ["ups_vMass"])
	
	massDF = dataY2S\
	.Filter("ups_vMass > 9.6 & ups_vMass < 10.3")\
	.Filter("ups_pT > 15")\
	.AsNumpy(columns=["ups_vMass"])
	
	mumuroodata = ROOT.RooDataSet\
	.from_numpy({"ups_vMass": massDF["ups_vMass"]},ROOT.RooArgSet(upsmass))
	mumuroohist = mumuroodata.binnedClone()
	
	xframe = upsmass.frame(Title="Y(2S) Mass")
	# mass Y2S PDG = 10.02326
	# mass Y3S PDG = 10.35520

		#signal mean
	mean2s = ROOT.RooRealVar("#mu_{Y(2S)}", "mean of gaussians", 
							  10., 9.8, 10.2)

		#sigmas
	sigma2s = ROOT.RooRealVar("#sigma_{Y(2S)}", "width", 0.063, 0.01, 5.) 
		
		#Crystal Ball parameters
	alpha = ROOT.RooRealVar ("#alpha","#alpha", 1.62, 0., 5.)
	n= ROOT.RooRealVar ("n","n", 2.5, 0., 5.)
	
	alpha.setConstant(1)	# from MC simulation
	n.setConstant(1)		# from MC simulation

		#chebychev coefficients
	f0 = ROOT.RooRealVar("f0", "f0", 5., 0., 10.)
	f1 = ROOT.RooRealVar("f1", "f1", 0., -20., 1.)
	f2 = ROOT.RooRealVar("f2", "f2", 2., 0., 8.)
		
		#fractions
	bkgfrac = ROOT.RooRealVar("f_{bkg}", "fraction of background", 0.5, 0.001, 1.)

	# 	MODELS FOR MASS PLOT

		#signals
	sig2s1 = ROOT.RooGaussian("signal2s_1", "signal2s_1", upsmass, mean2s, sigma2s)
	cb = ROOT.RooCBShape("Double CB", "#upsilon(2s) Pdf", upsmass, mean2s, sigma2s, alpha, n);
		
		#backgrounds
	bkg2 = ROOT.RooChebychev("bkg", "Background", upsmass, ROOT.RooArgList(f0,f1))
	bkg3 = ROOT.RooChebychev("bkg", "Background", upsmass, ROOT.RooArgList(f0,f1,f2))

		#MODELS
	model1 = ROOT.RooAddPdf("model1", "Cheb2+Gaus", [bkg2,sig2s1], [bkgfrac]) 
	model2 = ROOT.RooAddPdf("model2", "Cheb2+CB", [bkg2,cb], [bkgfrac]) 
	model3 = ROOT.RooAddPdf("model3", "Cheb3+Gaus", [bkg3,sig2s1], [bkgfrac]) 
	model4 = ROOT.RooAddPdf("model4", "Cheb3+CB", [bkg3,cb], [bkgfrac]) 

	allModels = {	1 : model1, 
					2 : model2, 
					5 : model3, 
					6 : model4
				}


		#line for pull plot
	zero = ROOT.TLine(9.7, 0., 10.3, 0.)
	zero.SetLineColor(ROOT.kRed)
	zero.SetLineWidth(2)
	zero.SetLineStyle(2)

	########	FIT RESULT AND CHI SQUARED COMPUTATION  ·················
	upsmass.setRange("range", 9.8 ,10.25)	 #set range before fitting
		
		
	c = ROOT.TCanvas("MassPlotY2S", "MassPlotY2S",1600,900)
	c.Divide(1,2)
	c.SetTitle("Different fits for Y(2S) InvMass")
		
		#PLOT DATA WITH DIFFERENT FITS
		
	for i, model in [(1,model3)]: #for i, model in allModels.items()
			
			#build frame with data
		xframe = upsmass.frame(Title = " ")	#Title="Y(2S) Mass"
		mumuroohist.plotOn(xframe, MarkerSize=0.3)
		
			#get components (bkg, sig)
		component = model.pdfList()
		
			#fit of all models
		model.fitTo(mumuroohist, Range="range")
		fitResult = model.fitTo(mumuroohist, Range="range", Save=True)
			
			#plotting data and fit
		
		model.plotOn(xframe, Components={component[0]}, LineStyle=":", LineColor="g")
		model.plotOn(xframe, Components={component[1]}, LineStyle=":", LineColor="r")
		model.plotOn(xframe,LineColor="b") #total curve
		
			#print all parameters + chi squared
		model.paramOn(xframe, ROOT.RooFit.Parameters([mean2s, sigma2s, bkgfrac]), ROOT.RooFit.Layout(0.1, 0.9, 0.9))	#choose parameters to print!!

			#add chisquare stats in a TPaveText
		text_box = ROOT.TPaveText(0.65, 0.75, 0.9, 0.9, "NDC")
		text_box.AddText( str(model.getTitle()) )	#type of fit
		text_box.SetFillColor(0)
		text_box.SetBorderSize(1)
		xframe.addObject(text_box)
		
			#add pull histos
		pull = xframe.pullHist()
		pull.GetYaxis().SetTitle("#chi^{2}")
		#pull.setXAxisLabel(" ")
		pull.GetXaxis().SetLimits(9.6, 10.3)
		pull.SetTitle(" ")
		pull.SetMarkerSize(0.3)

			#print plots in sub canvas
		c.cd(i)
		xframe.Draw("")
		
		c.cd(i+1)
		pull.Draw("")
		zero.Draw("same")	
###############################################################	
	"""
	fitResult = model.fitTo(mumuroohist, Range="range", Save=True)

	########	PLOTTING	···········································
	mumuroohist.plotOn(xframe, MarkerSize=0.3)
	model.plotOn(xframe,LineColor="r")

	component = model.pdfList()


	model.plotOn(xframe, Components={component[0]}, LineStyle=":", LineColor="g")
	model.plotOn(xframe, Components={component[1]}, LineStyle=":", LineColor="r")
	model.plotOn(xframe,LineColor="b")	#total
	
	model.paramOn(xframe, ROOT.RooFit.Layout(0.1, 0.9, 0.9)) #print all parameters

		#Type of fit
	text_box = ROOT.TPaveText(0.65, 0.75, 0.9, 0.9, "NDC")
	text_box.AddText( str(model.getTitle()) )	#type of fit
	text_box.SetFillColor(0)
	text_box.SetBorderSize(1)
	
		# pull plot
	pull = xframe.pullHist()
	pull.GetYaxis().SetTitle("#chi^{2}")
	#pull.setXAxisLabel(" ")
	pull.GetXaxis().SetLimits(9.6, 10.3)
	pull.SetTitle(" ")
	pull.SetMarkerSize(0.3)
	
	
	c.cd(1)
	xframe.Draw()
	text_box.Draw()
	
	c.cd(2)
	pull.Draw()
	zero.Draw("same")
	"""

	fileroot.WriteObject(c,"UpsInvMass_bestfit")
	c.SaveAs(d+"/MassPlotY2S_bestfit.pdf")
	os.system(f"xdg-open {d}/MassPlotY2S_bestfit.pdf")


#######		TRANSVERSAL MOMENTUM PLOTS	#############


def Y_pt():
	hist16 = dataY2S.Filter("run < 290000").Histo1D(("Y(2S) transverse momentum", "Y(2S) transverse momentum;p_{T}(#mu#mu) [GeV];Counts", 500, 0., 50), "ups_pT")
	hist17 = dataY2S.Filter("run < 310000").Histo1D(("Y(2S) transverse momentum", "Y(2S) transverse momentum;p_{T}(#mu#mu) [GeV];Counts", 500, 0., 50), "ups_pT")

	hist18 = dataY2S.Histo1D(("Y(2S) transverse momentum", "Y(2S) transverse momentum (stacked);p_{T}(#mu#mu) [GeV];Counts", 500, 0., 50.), "ups_pT")

	c_pt = ROOT.TCanvas("Y(2S) pT", "Y(2S) pT")


	hist17.SetFillColor(5)
	hist16.SetFillColor(3)


	hist18.Draw("")
	hist17.Draw("same")
	hist16.Draw("same")

	fileroot.WriteObject(c_pt,"UpsPT")
	c_pt.SaveAs(d+"/pt.pdf")
	os.system(f"xdg-open {d}/pt.pdf")


###########		PROBABILITY PLOT	###############

def Y_vProb():
	p_hist = dataY2S.Histo1D(("Y(2S) Probability", "Y(2S) Probability ;p;Counts", 500, 0., 1.), "ups_vProb")

	c_p = ROOT.TCanvas("Y2S prob", "Y2S prob")
	
	p_hist.Draw()
	
	fileroot.WriteObject(c_p,"ProbVertex")
	c_p.SaveAs(d+"/prob.pdf")
	os.system(f"xdg-open {d}/prob.pdf")


##########		RAPIDITY PLOT	###########

def Y_rap():
	hist16 = dataY2S.Filter("run < 290000").Histo1D(("Y(2S) Rapidity", "Y(2S) Rapidity;y(#mu#mu);Counts", 500,-2.5,2.5), "ups_rap")
	
	hist17 = dataY2S.Filter("run < 310000").Histo1D(("Y(2S) rapidity", "Y(2S) rapidity;y(#mu#mu);Counts", 500,-2.5,2.5), "ups_rap")
	
	hist19 = dataY2S.Histo1D(("Y(2S) rapidity", "Y(2S) rapidity (stacked);y(#mu#mu);Counts", 500, -2.5, 2.5), "ups_rap")

	hist17.SetFillColor(5)
	hist16.SetFillColor(3)


	c_rap = ROOT.TCanvas("Y2S Rapidity", "Y2S Rapidity")

	hist19.Draw("")
	hist17.Draw("same")
	hist16.Draw("same")

	fileroot.WriteObject(c_rap,"UpsRapidity")
	c_rap.SaveAs(d+"/rapidity.pdf")
	os.system(f"xdg-open {d}/rapidity.pdf")


#	PSEUDRAPIDITY PLOT	##############################
def Y_pseudorap():

	hist18 = dataY2S.Histo1D(("Y(2S) pseudorapidity", "Y(2S) Pseudorapidity (stacked);#eta(#mu#mu);Counts", 500, -2.0, 2.0), "ups_eta")
	
	hist16 = dataY2S.Filter("run < 290000").Histo1D(("Y(2S) pseudorapidity", "Y(2S) Pseudorapidity (stacked);#eta(#mu#mu);Counts", 500, -2.0, 2.0), "ups_eta")
	
	hist17 = dataY2S.Filter("run < 310000").Histo1D(("Y(2S) pseudorapidity", "Y(2S) Pseudorapidity (stacked);#eta(#mu#mu);Counts", 500, -2.0, 2.0), "ups_eta")

	hist17.SetFillColor(5)
	hist16.SetFillColor(3)


	c = ROOT.TCanvas("Y2S Pseudorapidity", "Y2S Pseudorapidity")

	hist18.Draw("")
	hist17.Draw("same")
	hist16.Draw("same")

	fileroot.WriteObject(c,"UpsPseudorap")
	c.SaveAs(d+"/Pseudorapidity.pdf")
	os.system(f"xdg-open {d}/Pseudorapidity.pdf")



#	DIMUON DECAY GEOMETRY		##############################
def dimuon_decay():

		#	angular phase plot
	hist1 = dataY2S.Histo2D(("dimuon decay", "Y(2S) #rightarrow #mu^{+}#mu^{-};#Delta#eta(#mu^{+}#mu^{-});#Delta#phi(#mu^{+}#mu^{-})", 50, 0., 2.5, 50, 0., 3.), "ups_deltaEta", "ups_deltaPhi")
	cprint([hist1], "Dimuon", "colz" )


		#	phase plot
	hist2 = dataY2S.Histo2D(("Dimuon decay", "Y(2S) #rightarrow #mu^{+}#mu^{-};#DeltaR(#mu^{+}#mu^{-});p_{T}(#mu^{+}#mu^{-}) [GeV]", 50, 0, 3, 50, 0, 150), "ups_deltaR", "ups_pT")
	cprint([hist2], "Dimuon2", "colz" )

		#	phase profile
	hist3 = dataY2S.Profile1D(("Dimuon decay", "Y(2S) #rightarrow #mu^{+}#mu^{-};p_{T}(#mu^{+}#mu^{-}) [GeV];#DeltaR(#mu^{+}#mu^{-})", 50, 0, 150), "ups_pT", "ups_deltaR")
	cprint([hist3],"profile",stats=True)
	


#	MENU

#	call functions by inserting the key from command line
compute = {	"1" : mu_pt,
			"2" : m_Y2S,
			"3" : fit_Y2S,
			"4" : Y_pt,		
			"5" : Y_vProb,		
			"6" : Y_rap,
			"7" : Y_pseudorap,
			"8" : dimuon_decay,
			"q" : exit
		  }

	
lang = input("\nSelect plots (Separate by spacing):\n1. Single Muon pT Plot\n2. Y(2S) Mass Plot\n3. Fit of Y(2S) Mass Plot\n4. Y(2S) pT\n5. Y(2S) Vertex Probability\n6. Y(2S) Rapidity\n7. Y(2S) Pseudorapidity\n8. Plots for Geometry of Di-Muon Decay\nENTER to print all plots.\nPress \"q\" to EXIT.\n").split()

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
		
		compute[i]()


	# End timer
end = time()

	# Calculate elapsed time
elapsed = end - start
print("\nComputing time: ", elapsed, "\n") 


fileroot.Close()

#From root to Pandas
# Y2S_tree = myfileY2SKK["rootuple/UpsTree;1"]
# print(Y2S_tree.keys())

# Y2SKK_tree = myfileY2SKK["rootuple/CandidateTree;1"]
# print(Y2SKK_tree.keys())

# data_Y2S = Y2S_tree.arrays(library="pd")
# print(data_Y2S.head())


#Upsilon:
['run', 'event', 'numPrimaryVertices', 'trigger', 'ups_p4', 'muonP_p4', 'muonN_p4', 'iPVwithmuons_ups', 'ups_vertexWeight', 'ups_vProb', 'ups_vMass', 'ups_vNChi2', 'ups_DCA', 'ups_ctauPV', 'ups_ctauErrPV', 'ups_lxyPV', 'ups_lxyErrPV', 'ups_cosAlpha', 'ups_ctauBS', 'ups_ctauErrBS', 'ups_lxyBS', 'ups_lxyErrBS', 'mu1_pt', 'mu1_ptErr', 'mu1_d0', 'mu1_d0Err', 'mu1_dz', 'mu1_dzErr', 'mu1_dxy', 'mu1_dxyErr', 'mu1_nvsh', 'mu1_nvph', 'mu1_charge', 'mu2_pt', 'mu2_ptErr', 'mu2_d0', 'mu2_d0Err', 'mu2_dz', 'mu2_dzErr', 'mu2_dxy', 'mu2_dxyErr', 'mu2_nvsh', 'mu2_nvph', 'mu2_charge']
#Upsilon KKcandidate_vMass
['run', 'event', 'nCandPerEvent', 'numPrimaryVertices', 'trigger', 'candidate_p4', 'track1_p4', 'track2_p4', 'ditrack_p4', 'dimuon_p4', 'muonp_p4', 
'muonn_p4', 'iPVwithmuons', 'dimuon_diMuIndx', 'dimuon_vertexWeight', 'dimuon_vProb', 'dimuon_vMass', 'dimuon_vNChi2', 'dimuon_DCA', 'dimuon_ctauPV', 
'dimuon_ctauErrPV', 'dimuon_lxyPV', 'dimuon_lxyErrPV', 'dimuon_cosAlpha', 'dimuon_ctauBS', 'dimuon_ctauErrBS', 'dimuon_lxyBS', 'dimuon_lxyErrBS', 
'candidate_vMass', 'candidate_vProb', 'candidate_vChi2', 'candidate_cosAlpha', 'candidate_ctauPV', 'candidate_ctauErrPV', 'candidate_charge', 
'candidate_lxy', 'candidate_lxyErr', 'candidate_lxyz', 'candidate_lxyzErr', 'thePrimaryV_X', 'thePrimaryV_Y', 'thePrimaryV_Z', 
'TheDecayVertex_X', 'TheDecayVertex_Y', 'TheDecayVertex_Z', 'thePrimaryV_2D_position', 'thePrimaryV_3D_position', 'TheDecayVertex_2D_position', 
'TheDecayVertex_3D_position', 'TheVertexDistance_2D', 'TheVertexDistance_3D', 'track1_d0', 'track1_d0Err', 'track1_dz', 'track1_dzErr', 'track1_dxy',
'track1_dxyErr', 'track1_nvsh', 'track1_nvph', 'track1_dRdimuon', 'track1_charge', 'track1_PV', 'track1_refVtx', 'track1_pvAssocQ', 'track1_dzAssocPV',
'track2_d0', 'track2_d0Err', 'track2_dz', 'track2_dzErr', 'track2_dxy', 'track2_dxyErr', 'track2_nvsh', 'track2_nvph', 'track2_dRdimuon', 
'track2_charge', 'track2_PV', 'track2_refVtx', 'track2_pvAssocQ', 'track2_dzAssocPV', 'ditrack_dRdimuon']
