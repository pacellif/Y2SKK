#	COMPUTING FUNCTIONS
from ROOT import gInterpreter

gInterpreter.Declare('''
Double_t ComputePseudoRapidity(Double_t px, Double_t py, Double_t pz, Double_t ene) {
    const ROOT::Math::PxPyPzE4D p(px, py, pz, ene);
    return p.Eta();
}
Double_t ComputePhi(Double_t px, Double_t py, Double_t pz, Double_t ene) {
    const ROOT::Math::PxPyPzE4D p(px, py, pz, ene);
    return p.Phi();
}
Double_t ComputeRapidity(Double_t px, Double_t py, Double_t pz, Double_t ene) {
    return 0.5 * log( (ene+pz)/(ene-pz) );
}
''')

gInterpreter.Declare('''
Double_t ComputeMuMuPiPiInvMass(Double_t px1, Double_t py1, Double_t pz1, 
                                    Double_t px2, Double_t py2, Double_t pz2,
                                    Double_t px3, Double_t py3, Double_t pz3,
                                    Double_t px4, Double_t py4, Double_t pz4) {
    const Double_t muonMass = 0.105658;
    const Double_t pionMass = 0.139570;
    const ROOT::Math::PxPyPzMVector mu1(px1, py1, pz1, muonMass);
    const ROOT::Math::PxPyPzMVector mu2(px2, py2, pz2, muonMass);
    const ROOT::Math::PxPyPzMVector pi1(px3, py3, pz3, pionMass);
    const ROOT::Math::PxPyPzMVector pi2(px4, py4, pz4, pionMass);
    return (mu1+mu2+pi1+pi2).mass();
}
Double_t ComputePiPiInvMass(Double_t px3, Double_t py3, Double_t pz3,
                     Double_t px4, Double_t py4, Double_t pz4) {
    const Double_t pionMass = 0.139570;
    const ROOT::Math::PxPyPzMVector pi1(px3, py3, pz3, pionMass);
    const ROOT::Math::PxPyPzMVector pi2(px4, py4, pz4, pionMass);
    return (pi1+pi2).mass();
}
''')

gInterpreter.Declare('''
Double_t ComputeMuMuPiPiInvMassMC(Double_t px1, Double_t py1, Double_t pz1,
                                    Double_t px2, Double_t py2, Double_t pz2,
                                    Double_t px3, Double_t py3, Double_t pz3,
                                    Double_t px4, Double_t py4, Double_t pz4) {
    const Double_t muonMass = 0.105658;
    const Double_t pionMass = 0.139570;
    const Double_t ups2SMass = 10.02326;
    const ROOT::Math::PxPyPzMVector mu1(px1, py1, pz1, muonMass);
    const ROOT::Math::PxPyPzMVector mu2(px2, py2, pz2, muonMass);
    const ROOT::Math::PxPyPzMVector pi1(px3, py3, pz3, pionMass);
    const ROOT::Math::PxPyPzMVector pi2(px4, py4, pz4, pionMass);
    return (mu1+mu2+pi1+pi2).mass() - (mu1+mu2).mass() + ups2SMass;
}
''')
