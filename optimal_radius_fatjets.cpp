/*
 * This code calculates the optimal R parameter to cluster w jets.
 * In the input root file, give a process that gives w jets which decay fully hadronically
 */

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <exception>
#include <cstdlib>
#include <cstdio>
#include <limits>
#include <numeric>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLorentzVector.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/CDFJetCluPlugin.hh"
#include "fastjet/internal/BasicRandom.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootUtilities.h"

#include "classes/DelphesClasses.h"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

int main(int argc, char* argv[]){
    gSystem->Load("libDelphes");
    TChain chain("Delphes");
    if (argc < 2){
        cerr << "usage: ./rsearch.cpp.ex <input ROOT file>\n";
        return EXIT_FAILURE;
    }
    char* filename = argv[1];

    chain.Add(filename);
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    TClonesArray* towerbr = treeReader->UseBranch("Tower");

    // store a list of all jet definitions for different R to make computation easier later
    map<double,JetDefinition> jetdefs;
    for ( double R = 0.4; R < 2.1; R += 0.1 )
        jetdefs.insert({R, JetDefinition(antikt_algorithm, R)});

    vector<map<double,long>> njet;
    for ( Long64_t entry = 0; entry < treeReader->GetEntries(); entry++ ) {
        // look at each event
        treeReader->ReadEntry(entry);
        // get a list of all particles at the tower
        vector<PseudoJet> particles;
        for ( int i = 0; i < towerbr->GetEntries(); i++ ){
            TLorentzVector p = ((Tower*)towerbr->At(i))->P4();
            particles.push_back(PseudoJet(p.Px(), p.Py(), p.Pz(), p.E()));
        }
        map<double,long> this_njet;
        for ( auto& jetdef : jetdefs ) {
            // for each jet definition, cluster all the particles
            ClusterSequence cs(particles, jetdef.second);
            // only consider the jets with pt > 300 to reduce noise
            vector<PseudoJet> jets = cs.inclusive_jets(300.0);
            long num_wjets = 0;
            // for all the jets,
            // look at the total invariant mass and note it if it is approx the mass of a w boson
            for ( auto& jet : jets )
                if ( 70 < jet.m() && jet.m() < 90 )
                    num_wjets++;
            this_njet[jetdef.first] = num_wjets;
        }
        // store all observations
        njet.push_back(this_njet);
    }

    // print data
    printf("R param   number of jets\n");
    for ( double R = 0.4; R < 2.1; R += 0.1 ) {
        long totjets = 0;
        for ( auto& entry : njet )
            totjets += entry[R];
        printf("  %1.1f %13ld\n", R, totjets);
    }

    return 0;
}
