#include "TTree.h"

#include "hxx_tree.h"

#include <iostream>

hxx_tree::hxx_tree()
{
   elec_pt     = new std::vector<double>;
   elec_eta    = new std::vector<double>;
   elec_phi    = new std::vector<double>;
   elec_ch     = new std::vector<double>;
   dielec_minv = new std::vector<double>;
   //elec_P4     = new std::vector<TLorentzVector>;
   
   muon_pt     = new std::vector<double>;
   muon_eta    = new std::vector<double>;
   muon_phi    = new std::vector<double>;
   muon_ch     = new std::vector<double>;
   dimuon_minv = new std::vector<double>;   

   jet_pt     = new std::vector<double>;
   jet_eta    = new std::vector<double>;
   jet_phi    = new std::vector<double>;
   jet_btag   = new std::vector<int>;
   jet_tautag = new std::vector<int>;
}

void hxx_tree::WriteTree(TTree * tree){   
   
   tree->Branch("testvar", &testvar);
   tree->Branch("sample", &sample);
   tree->Branch("weight", &weight);
   tree->Branch("lepton_flavor", &lepton_flavor);

   // multiplicities
   tree->Branch("nelec", &nelec);
   tree->Branch("nmuon", &nmuon);
   tree->Branch("njet",  &njet);
   
   // dilepton variables:
   tree->Branch("l1_pt",  &l1_pt  );
   tree->Branch("l1_eta", &l1_eta );
   tree->Branch("l1_phi", &l1_phi );
   tree->Branch("l2_pt",  &l2_pt  );
   tree->Branch("l2_eta", &l2_eta );
   tree->Branch("l2_phi", &l2_phi );
   tree->Branch("mll",    &mll    );

   // generator-level variables:
   tree->Branch("gl1_pt",  &gl1_pt  );
   tree->Branch("gl1_eta", &gl1_eta );
   tree->Branch("gl1_phi", &gl1_phi );
   tree->Branch("gl2_pt",  &gl2_pt  );
   tree->Branch("gl2_eta", &gl2_eta );
   tree->Branch("gl2_phi", &gl2_phi );
   tree->Branch("gmll",    &gmll    );


   // met variables:
   tree->Branch("nopu_met",     &nopu_met);
   tree->Branch("nopu_met_phi", &nopu_met_phi);
   
   tree->Branch("ht",     &ht);

   // jets:
   tree->Branch("jet_pt",     jet_pt);
   tree->Branch("jet_eta",    jet_eta);
   tree->Branch("jet_phi",    jet_phi);
   tree->Branch("jet_btag",   jet_btag);
   tree->Branch("jet_tautag", jet_tautag);

   // elecs:
   tree->Branch("elec_pt",     elec_pt);
   tree->Branch("elec_eta",    elec_eta);
   tree->Branch("elec_phi",    elec_phi);
   tree->Branch("dielec_minv", dielec_minv);
   tree->Branch("elec_ch",     elec_ch);
   //tree->Branch("elec_P4",      elecP4); 

   // muons:
   tree->Branch("muon_pt",     muon_pt);
   tree->Branch("muon_eta",    muon_eta);
   tree->Branch("muon_phi",    muon_phi);
   tree->Branch("muon_ch",     muon_ch);
   tree->Branch("dimuon_minv", dimuon_minv);

   // dilepton variables
   tree->Branch("leadingm",    &leadingm);
   tree->Branch("subleadingm", &subleadingm);

   // four lepton variables
   tree->Branch("total_minv",  &total_minv);
}
      
void hxx_tree::ReadTree(TTree * tree){
   tree->SetBranchAddress("l1_pt", &l1_pt     );

   tree->SetBranchAddress("testvar", &testvar );
   tree->SetBranchAddress("sample", &sample   );
   tree->SetBranchAddress("weight", &weight   );
   tree->SetBranchAddress("lepton_flavor", &lepton_flavor );

   // multiplicities
   tree->SetBranchAddress("nelec", &nelec     );
   tree->SetBranchAddress("nmuon", &nmuon     );
   tree->SetBranchAddress("njet",  &njet      );
   
   // dilepton variables:
   tree->SetBranchAddress("l1_pt",  &l1_pt    );
   tree->SetBranchAddress("l1_eta", &l1_eta   );
   tree->SetBranchAddress("l1_phi", &l1_phi   );
   tree->SetBranchAddress("l2_pt",  &l2_pt    );
   tree->SetBranchAddress("l2_eta", &l2_eta   );
   tree->SetBranchAddress("l2_phi", &l2_phi   );
   tree->SetBranchAddress("mll",    &mll      );

   // generator-level variables:
   tree->SetBranchAddress("gl1_pt",  &gl1_pt    );
   tree->SetBranchAddress("gl1_eta", &gl1_eta   );
   tree->SetBranchAddress("gl1_phi", &gl1_phi   );
   tree->SetBranchAddress("gl2_pt",  &gl2_pt    );
   tree->SetBranchAddress("gl2_eta", &gl2_eta   );
   tree->SetBranchAddress("gl2_phi", &gl2_phi   );
   tree->SetBranchAddress("gmll",    &gmll      );

   // met variables:
   tree->SetBranchAddress("nopu_met",     &nopu_met     );
   tree->SetBranchAddress("nopu_met_phi", &nopu_met_phi );

   tree->SetBranchAddress("ht",           &ht);

   // ROOT makes this complicated...
   delete jet_pt;     jet_pt     = NULL;
   delete jet_phi;    jet_phi    = NULL;
   delete jet_eta;    jet_eta    = NULL;
   delete jet_btag;   jet_btag   = NULL;
   delete jet_tautag; jet_tautag = NULL;
   
   delete elec_pt;     elec_pt     = NULL;
   delete elec_phi;    elec_phi    = NULL;
   delete elec_eta;    elec_eta    = NULL;
   delete elec_ch;     elec_ch     = NULL;
   delete dielec_minv; dielec_minv = NULL;
   //delete elec_P4;     elec_P4     = NULL;

   delete muon_pt;     muon_pt     = NULL;
   delete muon_phi;    muon_phi    = NULL;
   delete muon_eta;    muon_eta    = NULL;
   delete muon_ch;     muon_ch     = NULL;
   delete dimuon_minv; dimuon_minv = NULL;
   
   tree->SetBranchAddress("jet_pt",     &jet_pt     );
   tree->SetBranchAddress("jet_phi",    &jet_phi    );
   tree->SetBranchAddress("jet_eta",    &jet_eta    );
   tree->SetBranchAddress("jet_btag",   &jet_btag   );
   tree->SetBranchAddress("jet_tautag", &jet_tautag );
   
   tree->SetBranchAddress("elec_pt",     &elec_pt     );
   tree->SetBranchAddress("elec_phi",    &elec_phi    );
   tree->SetBranchAddress("elec_eta",    &elec_eta    );
   tree->SetBranchAddress("dielec_minv", &dielec_minv );
   tree->SetBranchAddress("elec_ch",     &elec_ch     );
   //tree->SetBranchAddress("elec_P4",     &elec_P4     );

   tree->SetBranchAddress("muon_pt",     &muon_pt     );
   tree->SetBranchAddress("muon_phi",    &muon_phi    );
   tree->SetBranchAddress("muon_eta",    &muon_eta    );
   tree->SetBranchAddress("muon_ch",     &muon_ch     );
   tree->SetBranchAddress("dimuon_minv", &dimuon_minv );

   // dilepton variables
   tree->SetBranchAddress("leadingm",    &leadingm    );
   tree->SetBranchAddress("subleadingm", &subleadingm );

   // four lepton variables
   tree->SetBranchAddress("total_minv",  &total_minv  );

}

void hxx_tree::erase_jet(int i){
   jet_pt     ->erase(jet_pt->begin() + i);
   jet_eta    ->erase(jet_eta->begin() + i);
   jet_phi    ->erase(jet_phi->begin() + i);
   jet_btag   ->erase(jet_btag->begin() + i);
   jet_tautag ->erase(jet_tautag->begin() + i);
}

void hxx_tree::Clear(){
   testvar  = 0;
   sample   = 0;
   weight   = 1;
   weight_met = 0;
   lepton_flavor = 0;
   nelec    = 0;
   nmuon    = 0;
   njet     = 0;
   l1_pt    = 0;
   l1_eta   = 0;
   l1_phi   = 0;
   l2_pt    = 0;
   l2_eta   = 0;
   l2_phi   = 0;
   mll      = 0;
   gl1_pt   = 0;
   gl1_eta  = 0;
   gl1_phi  = 0;
   gl2_pt   = 0;
   gl2_eta  = 0;
   gl2_phi  = 0;
   gmll     = 0;
   nopu_met = 0;
   nopu_met_phi = 0;
   ht = 0;
   jet_pt ->clear();
   jet_eta->clear();
   jet_phi->clear();
   jet_btag->clear();
   jet_tautag->clear();
   elec_pt ->clear();
   elec_eta->clear();
   elec_phi->clear();
   elec_ch->clear();
   dielec_minv->clear();
   //elec_P4->clear();
   muon_pt-> clear();
   muon_eta->clear();
   muon_phi->clear();
   muon_ch->clear();
   dimuon_minv->clear();
   leadingm    = 0;
   subleadingm = 0;
   total_minv  = 0;
}
