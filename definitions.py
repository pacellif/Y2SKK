#	DEFINITIONS OF NEW COLUMNS

#Definitions for Upsilon dataset

def UpsTreeDefinitions (dataframe):
	return dataframe\
	.Define("ups_pT", "sqrt(ups_p4_fX*ups_p4_fX + ups_p4_fY*ups_p4_fY)")\
	.Define("ups_eta", "ComputePseudoRapidity(ups_p4_fX, ups_p4_fY, ups_p4_fZ, ups_p4_fE)")\
	.Define("ups_phi", "ComputePhi(ups_p4_fX, ups_p4_fY, ups_p4_fZ, ups_p4_fE)")\
	.Define("ups_rap", "ComputeRapidity(ups_p4_fX, ups_p4_fY, ups_p4_fZ, ups_p4_fE)")\
	.Define("muonP_pT", "sqrt(muonP_p4_fX*muonP_p4_fX + muonP_p4_fY*muonP_p4_fY)")\
	.Define("muonP_eta", "ComputePseudoRapidity(muonP_p4_fX, muonP_p4_fY, muonP_p4_fZ, muonP_p4_fE)")\
	.Define("muonP_phi", "ComputePhi(muonP_p4_fX, muonP_p4_fY, muonP_p4_fZ, muonP_p4_fE)")\
	.Define("muonP_rap", "ComputeRapidity(muonP_p4_fX, muonP_p4_fY, muonP_p4_fZ, muonP_p4_fE)")\
	.Define("muonN_pT", "sqrt(muonN_p4_fX*muonN_p4_fX + muonN_p4_fY*muonN_p4_fY)")\
	.Define("muonN_eta", "ComputePseudoRapidity(muonN_p4_fX, muonN_p4_fY, muonN_p4_fZ, muonN_p4_fE)")\
	.Define("muonN_phi", "ComputePhi(muonN_p4_fX, muonN_p4_fY, muonN_p4_fZ, muonN_p4_fE)")\
	.Define("muonN_rap", "ComputeRapidity(muonN_p4_fX, muonN_p4_fY, muonN_p4_fZ, muonN_p4_fE)")\
	.Define("ups_deltaEta", "abs(muonP_eta - muonN_eta)")\
	.Define("ups_deltaPhi", "Double_t angle = abs(muonP_phi - muonN_phi); if (angle > ROOT::Math::Pi()) angle = 2*ROOT::Math::Pi() - angle; return angle;")\
	.Define("ups_deltaR", "sqrt(ups_deltaEta*ups_deltaEta + ups_deltaPhi*ups_deltaPhi)")
	
def CandTreeDefinitions (dataframe):
	return dataframe\
	.Define("ditrack_mass", "sqrt(ditrack_p4_fE*ditrack_p4_fE - ditrack_p4_fX*ditrack_p4_fX - ditrack_p4_fY*ditrack_p4_fY - ditrack_p4_fZ*ditrack_p4_fZ)")\
	.Define("dimuon_mass", "sqrt(dimuon_p4_fE*dimuon_p4_fE - dimuon_p4_fX*dimuon_p4_fX - dimuon_p4_fY*dimuon_p4_fY - dimuon_p4_fZ*dimuon_p4_fZ)")\
	.Define("track1_pT", "sqrt(track1_p4_fX*track1_p4_fX + track1_p4_fY*track1_p4_fY)")\
	.Define("track1_eta", "ComputePseudoRapidity(track1_p4_fX, track1_p4_fY, track1_p4_fZ, track1_p4_fE)")\
	.Define("track2_pT", "sqrt(track2_p4_fX*track2_p4_fX + track2_p4_fY*track2_p4_fY)")\
	.Define("track2_eta", "ComputePseudoRapidity(track2_p4_fX, track2_p4_fY, track2_p4_fZ, track2_p4_fE)")\
	.Define("ditrack_pT", "sqrt(ditrack_p4_fX*ditrack_p4_fX + ditrack_p4_fY*ditrack_p4_fY)")\
	.Define("ditrack_eta", "ComputePseudoRapidity(ditrack_p4_fX, ditrack_p4_fY, ditrack_p4_fZ, ditrack_p4_fE)")\
	.Define("dimuon_pT", "sqrt(dimuon_p4_fX*dimuon_p4_fX + dimuon_p4_fY*dimuon_p4_fY)")\
	.Define("dimuon_eta", "ComputePseudoRapidity(dimuon_p4_fX, dimuon_p4_fY, dimuon_p4_fZ, dimuon_p4_fE)")\
	.Define("candidate_pT", "sqrt(candidate_p4_fX*candidate_p4_fX + candidate_p4_fY*candidate_p4_fY)")\
	.Define("candidate_eta", "ComputePseudoRapidity(candidate_p4_fX, candidate_p4_fY, candidate_p4_fZ, candidate_p4_fE)")\
	.Define("trackL_pT", "max(track1_pT, track2_pT)")\
	.Define("trackS_pT", "min(track1_pT, track2_pT)")\
	.Define("mumupipi_mass", "ComputeMuMuPiPiInvMass(muonp_p4_fX, muonp_p4_fY, muonp_p4_fZ, muonn_p4_fX, muonn_p4_fY, muonn_p4_fZ, track1_p4_fX, track1_p4_fY, track1_p4_fZ, track2_p4_fX, track2_p4_fY, track2_p4_fZ)")\
	.Define("pipi_mass", "ComputePiPiInvMass(track1_p4_fX, track1_p4_fY, track1_p4_fZ, track2_p4_fX, track2_p4_fY, track2_p4_fZ)")\
	.Define("mumupipi_MC", "ComputeMuMuPiPiInvMassMC(muonp_p4_fX, muonp_p4_fY, muonp_p4_fZ, muonn_p4_fX, muonn_p4_fY, muonn_p4_fZ, track1_p4_fX, track1_p4_fY, track1_p4_fZ, track2_p4_fX, track2_p4_fY, track2_p4_fZ)")	#Mass constraint
	
	
	
