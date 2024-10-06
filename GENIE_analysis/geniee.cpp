#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TObjArrayIter.h>
#include "Framework/Messenger/Messenger.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Conventions/Constants.h"

int main() {
    // Open the GHEP/ROOT file
    std::string filename = "/data/sample.ghep.root";
    TFile file(filename.c_str(), "READ");
    if (file.IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return 1;
    }

    // Get the tree header & print it
    NtpMCTreeHeader* header = dynamic_cast<NtpMCTreeHeader*>(file.Get("header"));
    if (!header) {
        std::cerr << "Error reading header from file: " << filename << std::endl;
        return 1;
    }
    LOG("test", pINFO) << *header;

    // Get the GENIE GHEP tree and set its branch address
    TTree* tree = dynamic_cast<TTree*>(file.Get("gtree"));
    if (!tree) {
        std::cerr << "Error reading tree from file: " << filename << std::endl;
        return 1;
    }

    NtpMCEventRecord* mcrec = nullptr;
    tree->SetBranchAddress("gmrec", &mcrec);

    // Event loop
    for (Long64_t i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);

        // Print-out the event
        EventRecord& event = *(mcrec->event);
        LOG("test", pINFO) << event;

        // Count final state muons
        int n_muons = 0;
        TObjArrayIter iter(&event);
        GHepParticle* p = nullptr;

        // Loop over event particles
        while ((p = dynamic_cast<GHepParticle*>(iter.Next()))) {
            int pdgc = p->Pdg();
            int status = p->Status();
            if (status != kIStStableFinalState) continue;
            bool is_muon = (pdgc == kPdgMuon || pdgc == -kPdgMuon);
            if (is_muon) n_muons++;
        }

        std::cout << "Event " << i << " has " << n_muons << " final state muons." << std::endl;
        mcrec->Clear();
    }

    file.Close();
    return 0;
}
