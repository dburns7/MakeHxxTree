// -*- C++ -*-
//
// Package:    MakeHxxTree
// Class:      MakeHxxTree
// 
/**\class MakeHxxTree MakeHxxTree.cc hxx_4l/MakeHxxTree/src/MakeHxxTree.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  dburns@ucdavis.edu
//         Created:  Wed May  7 12:49:25 PDT 2014
// $Id$
//
//


// system include files
#include <memory>
#include <map>
#include <string>
#include <iostream>
using namespace std;

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TTree.h"
#include "hxx_tree.h"

using namespace edm;
//
// class declaration
//

class MakeHxxTree : public edm::EDAnalyzer {
   public:
      explicit MakeHxxTree(const edm::ParameterSet&);
      ~MakeHxxTree();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      edm::InputTag elecSrc_;
      edm::InputTag muonSrc_;
      edm::InputTag metSrc_;
      
      TTree * evt_tree;
      hxx_tree data;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MakeHxxTree::MakeHxxTree(const edm::ParameterSet& iConfig):
  elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc")),
  muonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
  metSrc_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc"))
{
   //now do what ever initialization is needed
  Service<TFileService> fs;
  evt_tree = fs->make<TTree>("hxxtree", "hxxtree");//("Ntuple",    "Ntuple");
  data.WriteTree(evt_tree);
  //data.AddBranches(evt_tree);
}


MakeHxxTree::~MakeHxxTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MakeHxxTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  int sample = 100;
  double lumi = 20;
  double Ngen = 1e4;
  double xsec = 1;
  double weight = lumi * xsec / Ngen;
  
  // get electron collection
  edm::Handle<edm::View<pat::Electron> > elecs;
  iEvent.getByLabel(elecSrc_,elecs);

  // get muon collection
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muonSrc_,muons);

  // get met collection
  edm::Handle<edm::View<pat::MET> > mets;
  iEvent.getByLabel(metSrc_,mets);

  
  int muonCnt = muons->size();
  int elecCnt = elecs->size();

  data.Clear();
  
  data.testvar = 2;
  data.sample = sample;
  data.weight = weight;
  data.nelec = elecCnt;
  data.nmuon = muonCnt;

  int es = 0;
  int ms = 0;

  if(muonCnt == 2 && elecCnt == 2)
  {
    for(edm::View<pat::Muon>::const_iterator muon1=muons->begin(); muon1!=muons->end(); ++muon1)
    {
      ms += 1;

      for(edm::View<pat::Muon>::const_iterator muon2=muons->begin(); muon2!=muons->end(); ++muon2)
      {
        if(muon2 < muon1) continue;
        //if(muon1->charge()*muon2->charge()>0) continue;
        ms += 1;
        
        for(edm::View<pat::Electron>::const_iterator elec1=elecs->begin(); elec1!=elecs->end(); ++elec1)
        {
          es += 1;
          
          for(edm::View<pat::Electron>::const_iterator elec2=elecs->begin(); elec2!=elecs->end(); ++elec2)
          {
            if(elec2<elec1) continue;
            //if(electron1->charge()*electron2->charge()>0) continue;
            es += 1;
            if(ms > 2 && es > 2) break;
            
            double dielec_minv, dimuon_minv, total_minv, leadingm, subleadingm;
            dielec_minv = (elec1->p4()+elec2->p4()).mass();
            dimuon_minv = (muon1->p4()+muon2->p4()).mass();
            total_minv = (muon1->p4()+muon2->p4()+elec1->p4()+elec2->p4()).mass();
            if(abs(dielec_minv-91) < abs(dimuon_minv-91))
            {
              leadingm = dielec_minv;
              subleadingm = dimuon_minv;
            }
            else if(abs(dimuon_minv-91) < abs(dielec_minv-91))
            {
              leadingm = dimuon_minv;
              subleadingm = dielec_minv;
            }
            
            data.muon_pt->push_back(muon1->pt());
            data.muon_eta->push_back(muon1->eta());
            data.muon_phi->push_back(muon1->phi());
            data.muon_ch->push_back(muon1->charge());
            data.muon_pt->push_back(muon2->pt());
            data.muon_eta->push_back(muon2->eta());
            data.muon_phi->push_back(muon2->phi());
            data.muon_ch->push_back(muon2->charge());
            data.elec_pt->push_back(elec1->pt());
            data.elec_eta->push_back(elec1->eta());
            data.elec_phi->push_back(elec1->phi());
            data.elec_ch->push_back(elec1->charge());  
            data.elec_pt->push_back(elec2->pt());
            data.elec_eta->push_back(elec2->eta());
            data.elec_phi->push_back(elec2->phi());
            data.elec_ch->push_back(elec2->charge());
            data.dielec_minv->push_back(dielec_minv);
            data.dimuon_minv->push_back(dimuon_minv);
            data.leadingm = leadingm;
            data.subleadingm = subleadingm;
            data.total_minv = total_minv;
          }
        }
      } 
    }
  }

  if(muonCnt == 4 && elecCnt == 0)
  {
    for(edm::View<pat::Muon>::const_iterator muon1=muons->begin(); muon1!=muons->end(); ++muon1)
    {
      ms += 1;

      for(edm::View<pat::Muon>::const_iterator muon2=muons->begin(); muon2!=muons->end(); ++muon2)
      {
        if(muon2 < muon1) continue;
        //if(muon1->charge()*muon2->charge()>0) continue;
        ms += 1;

        for(edm::View<pat::Muon>::const_iterator muon3=muons->begin(); muon3!=muons->end(); ++muon3)
        {
          if(muon3 < muon1 | muon3 < muon2) continue;
          ms+=1;

          for(edm::View<pat::Muon>::const_iterator muon4=muons->begin(); muon4!=muons->end(); ++muon4)
          {
            if(muon4 < muon1 | muon4 < muon2 | muon4 < muon3) continue;
            ms += 1;
            if(ms > 4) break;
            
            double dimuon_minv1, dimuon_minv2, total_minv, leadingm, subleadingm;
            dimuon_minv1 = (muon1->p4()+muon2->p4()).mass();
            dimuon_minv2 = (muon3->p4()+muon4->p4()).mass();
            total_minv = (muon1->p4()+muon2->p4()+muon3->p4()+muon4->p4()).mass();
            if(abs(dimuon_minv1-91) < abs(dimuon_minv2-91))
            {
              leadingm = dimuon_minv1;
              subleadingm = dimuon_minv2;
            }
            else if(abs(dimuon_minv2-91) < abs(dimuon_minv1-91))
            {
              leadingm = dimuon_minv2;
              subleadingm = dimuon_minv1;
            }

            data.muon_pt->push_back(muon1->pt());
            data.muon_eta->push_back(muon1->eta());
            data.muon_phi->push_back(muon1->phi());
            data.muon_ch->push_back(muon1->charge());
            data.muon_pt->push_back(muon2->pt());
            data.muon_eta->push_back(muon2->eta());
            data.muon_phi->push_back(muon2->phi());
            data.muon_ch->push_back(muon2->charge());
            data.muon_pt->push_back(muon3->pt());
            data.muon_eta->push_back(muon3->eta());
            data.muon_phi->push_back(muon3->phi());
            data.muon_ch->push_back(muon3->charge());
            data.muon_pt->push_back(muon4->pt());
            data.muon_eta->push_back(muon4->eta());
            data.muon_phi->push_back(muon4->phi());
            data.muon_ch->push_back(muon4->charge());
            data.dimuon_minv->push_back(dimuon_minv1);
            data.dimuon_minv->push_back(dimuon_minv2);
            data.leadingm = leadingm;
            data.subleadingm = subleadingm;
            data.total_minv = total_minv;

          }
        }
      }
    }  
  }

  if(muonCnt == 0 && elecCnt == 4)
  {
    for(edm::View<pat::Electron>::const_iterator elec1=elecs->begin(); elec1!=elecs->end(); ++elec1)
    {
      es += 1;
      
      for(edm::View<pat::Electron>::const_iterator elec2=elecs->begin(); elec2!=elecs->end(); ++elec2)
      {
        if(elec2 < elec1) continue;
        es += 1;
        
        for(edm::View<pat::Electron>::const_iterator elec3=elecs->begin(); elec3!=elecs->end(); ++elec3)
        {
          if(elec3 < elec1 | elec3 < elec2) continue;
          es+=1;
           
          for(edm::View<pat::Electron>::const_iterator elec4=elecs->begin(); elec4!=elecs->end(); ++elec4)
          {
            if(elec4 < elec1 | elec4 < elec2 | elec4 < elec3) continue;
            es+=1;
            if(es > 4) break;
            
            double dielec_minv1, dielec_minv2, total_minv, leadingm, subleadingm;
            dielec_minv1 = (elec1->p4()+elec2->p4()).mass();
            dielec_minv2 = (elec3->p4()+elec4->p4()).mass();
            total_minv = (elec1->p4()+elec2->p4()+elec3->p4()+elec4->p4()).mass();
            if(abs(dielec_minv1-91) < abs(dielec_minv2-91))
            {
              leadingm = dielec_minv1;
              subleadingm = dielec_minv2;
            }
            else if(abs(dielec_minv2-91) < abs(dielec_minv1-91))
            {
              leadingm = dielec_minv2;
              subleadingm = dielec_minv1;
            }

            data.elec_pt->push_back(elec1->pt());
            data.elec_eta->push_back(elec1->eta());
            data.elec_phi->push_back(elec1->phi());
            data.elec_ch->push_back(elec1->charge());
            data.elec_pt->push_back(elec2->pt());
            data.elec_eta->push_back(elec2->eta());
            data.elec_phi->push_back(elec2->phi());
            data.elec_ch->push_back(elec2->charge());
            data.elec_pt->push_back(elec3->pt());
            data.elec_eta->push_back(elec3->eta());
            data.elec_phi->push_back(elec3->phi());
            data.elec_ch->push_back(elec3->charge());
            data.elec_pt->push_back(elec4->pt());
            data.elec_eta->push_back(elec4->eta());
            data.elec_phi->push_back(elec4->phi());
            data.elec_ch->push_back(elec4->charge());
            data.dielec_minv->push_back(dielec_minv1);
            data.dielec_minv->push_back(dielec_minv2);
            data.leadingm = leadingm;
            data.subleadingm = subleadingm;
            data.total_minv = total_minv;
          }
        }
      }
    }    
  }



  evt_tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
MakeHxxTree::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MakeHxxTree::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MakeHxxTree::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MakeHxxTree::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MakeHxxTree::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MakeHxxTree::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MakeHxxTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MakeHxxTree);
