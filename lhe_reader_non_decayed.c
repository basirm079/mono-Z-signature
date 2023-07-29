//g++ -std=c++11 lhe_reader_non_decayed.c -o lhe_reader_non_decayed -I`root-config --incdir` `root-config --libs`
//&>/dev/null;x="${0%.*}";[ ! "$x" -ot "$0" ]||(rm -f "$x";g++ -o "$x" "$0" -I`root-config --incdir` `root-config --libs`);exit

// Build (for ubuntu 20... with g++ 9 and above): g++ lhe_reader_non_decayed.c -o lhe_reader_non_decayed -I`root-config --incdir` `root-config --libs`
// OR, g++ -std=c++11 lhe_reader_non_decayed.c -o lhe_reader_non_decayed -I`root-config --incdir` `root-config --libs`
// The above line will produce executable file name lhe_reader_non_decayed
// ./lhe_reader_non_decayed unweighted_events
#include <cmath>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <TH1F.h>
#include <TTree.h>
#include <TH2F.h>
#include <TFile.h>
//#include </home/muhammad/root6/build/include/TLorentzVector.h>
#include <TLorentzVector.h>
#include "TCanvas.h"
using namespace std;

// Pour une description du format leshouches
// hep-ph/0609017
//
// pour chaque evenement
// une ligne générale : NbPart idprocess poids scale alpha_em alpha_s
// pour chaque particule : id status mere1 mere2 couleur1 couleur2 px py pz E m lifetime spin  
int main(int argc, char **argv) {

  if (argc != 2) {
    cout << "Missing argument. Expected LHE filename without '.lhe'"<< endl;
    exit(-5);
  }
  float DeltaR(float eta1, float eta2, float phi1, float phi2);
  float DeltaPhi(float phi1, float phi2);

  string basename=argv[1];

  string lhefname = basename+".lhe";
  string rootfname = basename+".root";

  string tt;
  int event=0;
  int npart,idprocess;
  double weight,scale,alpha_em,alpha_s;

  TH1::SetDefaultSumw2(true);

  TH1F* hPt_xd_xd = new TH1F("Et_xd_xd","xd_xd Et",100,0,510);//momentum or energy
  TH1F* hphi_l_l = new TH1F("phi_l_l","l_l phi",100,-4,4);//phi btw lepton pair
  TH1F* hDR_l_l = new TH1F("DR_l_l","l_l DR",100,0,5) ;//phi-eta separation
  TH1F* hPt_l_l = new TH1F("Pt_l_l","l_l Pt",100,0,510);//momentum of lepton pair
  TH1F* hMt = new TH1F("Mt","MT",100,0,700);//transverse mass
  TH1F* hHt = new TH1F("Ht_all","all Ht",100,0,1850);//scalar sum
  TH1F* hTpv= new TH1F("new_all","all new",100,0,12);//Etmiss/scalar sum
  TH1F* hM_l_l = new TH1F("DM_l_l","l_l M",100,0,300);//invariant mass
  TH1F* hPt_lead_l = new TH1F("pt_lead_l","lead_l pt",100,0,370);//momentum of leading lepton
  TH1F* hPt_sublead_l = new TH1F("pt_sublead_l","sublead_l pt",100,0,250);//momentum of sub-leading lepton
  TH1F* hEta_lead_l = new TH1F("eta_lead_l","lead_l eta",100,-6,6);//eta of leading lepton
  TH1F* hEta_sublead_l = new TH1F("eta_sublead_l","sublead_l eta",100,-7,7);//eta of subleading lepton
  TH1F* hphi_l_l_ptmiss = new TH1F("phi_l_l_ptmissing","l_l_ptmissing phi",100,0,5);//phi-eta separation of lepton pair and missing energy
  TTree *bgtree5 = new TTree("bgtree5","bdt");
  double etmiss=0,del_phi=0,del_R=0,pt_ll=0,Ht=0,new_all=0,in_mass=0,lead_l=0,sublead_l=0,elead_l=0,esublead_l=0,l_l_ptmiss=0,MT;
  bgtree5->Branch("MT",&MT,"MT/D"); // info stored in trees for BDT
  bgtree5->Branch("etmiss",&etmiss,"etmiss/D");
  bgtree5->Branch("del_phi",&del_phi,"del_phi/D");
  bgtree5->Branch("del_R",&del_R,"del_R/D");
  bgtree5->Branch("pt_ll",&pt_ll,"pt_ll/D");
  bgtree5->Branch("Ht",&Ht,"Ht/D");
  bgtree5->Branch("new_all",&new_all,"new_all/D");
  bgtree5->Branch("in_mass",&in_mass,"in_mass/D");
  bgtree5->Branch("lead_l",&lead_l,"lead_l/D");
  bgtree5->Branch("sublead_l",&sublead_l,"sublead_l/D");
  bgtree5->Branch("elead_l",&elead_l,"elead_l/D");
  bgtree5->Branch("esublead_l",&esublead_l,"esublead_l/D");
  bgtree5->Branch("l_l_ptmiss",&l_l_ptmiss,"l_l_ptmiss/D");
  

  int nlept=0, nsemi=0, nhadr=0,z=0;
  float cut;   //for cut optimization
  //cout<<"cut"<<endl;
  //cin>>cut;
  ifstream ff(lhefname.c_str(),ios::in); //ouverture du fichier .lhe
  //ifstream ff("test.lhe",ios::in); //ouverture du fichier .lhe
  //ifstream ff("/home/cms/perries/madgraph-newphysics/s1/madevent/Events/zp4000_unweighted_events.lhe",ios::in); //ouverture du fichier .lhe
  //ifstream ff("/home/cms/perries/madgraph-newphysics/QCD/madevent/Events/qcd_unweighted_events.lhe",ios::in); //ouverture du fichier .lhe
  int negativeWeight = 0;
  long long line = 0;
  while(!ff.eof()) {
    std::stringstream buffer;
    ff>>tt;
    buffer << tt << std::endl;
    line++;
//    cout<<"this is line: "<<line<<endl;

    if(tt=="<event>") {
      ff>>npart>>idprocess>>weight>>scale>>alpha_em>>alpha_s; //event definition
      buffer << npart << " " << idprocess << " " << weight << " " << scale << " " << alpha_em << " " << alpha_s << std::endl;
      line++;
//    cout<<"this is line inside: "<<line<<endl;
      event++;
     // cout<<"weight="<<weight<<endl;
      if (weight < 0) {
        negativeWeight++;
        weight = -1;
      } else {
        weight = 1;
      }
      /*weight = 1.;*/
//      if (event==2)return false;
      if (event%100==0) cout << "reading event "<< event << endl;
      int lmin=-1, lmax=-1, bmin=-1, met1=-1, met2=-1, bmax=-1, muon=-1, met=-1, topLep,topbJet, tbarLep, tbarbJet, top_jet, tbar_jet;
      int xd = -1, xd_ = -1; //dark matter work
      int Part_Id, Moth_Id, Part_absId, Moth_absId;
int n_lep=0, n_alep=0, n_topjets=0, n_tbarjets=0, n_topbjets=0, n_bbar=0, n_tbar_nu=0, n_top_nu=0;
      int n_h =0;
      int c1=0,c2=0,c3=0,c4=0,l1=0,l2=0,l3=0,l4=0;
      int n_xd=0, n_xd_=0; //DM matter work
      int n_bj=0, n_bj_=0; //DM matter work
int  bj[2]={-1,-1};
int lp[2]={-1,-1};
      int q[4]={-1,-1,-1,-1};
      int top=-1,topbar=-1,zprime=-1;
      int *Id      = new int[npart+1];
      int *Status  = new int[npart+1];
      int *Mother1 = new int[npart+1];
      int *Mother2 = new int[npart+1];
      int *Color1  = new int[npart+1];
      int *Color2  = new int[npart+1];
      double *px = new double[npart+1];
      double *py = new double[npart+1];
      double *pz = new double[npart+1];
      double *E = new double[npart+1];
      double *m = new double[npart+1];
      double *lifetime = new double[npart+1];
      double *spin = new double[npart+1];
      TLorentzVector **v = new TLorentzVector*[npart+1];
      TLorentzVector v_top_lep, v_top_bJet, v_tbar_lep, v_tbar_bJet, v_top_jet, v_tbar_jet, v_tbar_nu, v_top_nu;
      TLorentzVector v_arb,v_xd, v_xd_,v_xd1,v_xd1_,v_xd2, v_xd2_,v_lep1,v_lep1_,v_lep2,v_lep2_,v_1,v_2; //DM lorentz vec
      TLorentzVector v_bj, v_bj_; //bj lorentz vec
      TLorentzVector v_higgs, v_ps_a, v_ps_A; //higgs lorentz vec
      // in lhe first line is number 1, so fill unused array [0] with a crazy value;
      Id[0]= -99999;
      Status[0]= -99999;
      Mother1[0]= -99999;
      Mother2[0]= -99999;
      Color1[0]= -99999;
      Color2[0]= -99999;
      px[0]= -99999;
      py[0]= -99999;
      pz[0]= -99999;
      E[0]= -99999;
      m[0]= -99999;
      lifetime[0]= -99999;
      float val[4]={0,0,0,0};
      float val1[4]={0,0,0,0};
      spin[0]= -99999;double sum1=0;Ht=0;
     for (int i=1 ; i<npart+1 ; i++) { //start at one
        ff >> Id[i] >> Status[i] >> Mother1[i] >> Mother2[i] >> Color1[i] >> Color2[i]
           >> px[i] >> py[i] >> pz[i] >> E[i] >> m[i] >> lifetime[i] >> spin[i] ;
        buffer << Id[i] << " " << Status[i] << " " << std::endl;
        line++;
        v[i] = new TLorentzVector(px[i], py[i], pz[i], E[i]);
        if (Status[i]==-1) continue; // status -1 = initial quark ==> skip
        int id = Id[i];  
        Part_absId = abs(Id[i]);
        Moth_absId =abs( Id[Mother1[i]]);
        Part_Id = Id[i];
        Moth_Id = Id[Mother1[i]];
        v_arb.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
        Ht = Ht + v_arb.Pt();
        if ( Id[i]==12 ) 
        {
          v_xd.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          c1=1;
        }
        if ( Id[i]==-12 ) 
        { 
          v_xd_.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          c2=1;
        }        
        if ( Id[i]==14 ) 
        {
          v_xd2.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          c3=1;
        }
        if ( Id[i]==-14 ) 
        { 
          v_xd2_.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          c4=1;
        }
        if ( Id[i]==-52 ) 
        { 
          v_xd1_.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          n_xd_++;
        }
        if ( Id[i]==52 ) 
        { 
          v_xd1.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          n_xd++;
        }
         if(Id[i]==11)
        {
        v_lep1.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
        val[0]=v_lep1.Pt();
        val1[0]=v_lep1.Eta();
        l1=1;
        }
        if(Id[i]==-11)
        {
        v_lep1_.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
        val[1]=v_lep1_.Pt();
        val1[1]=v_lep1_.Eta();
        l2=1;
        }
          if(Id[i]==13)
        {
        v_lep2.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
        val[2]=v_lep2.Pt();
        val1[2]=v_lep2.Eta();
        l3=1;
        }
        if(Id[i]==-13)
        {
        v_lep2_.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
        val[3]=v_lep2_.Pt();
        val1[3]=v_lep2_.Eta();
        l4=1;
        }
        if(l1==1&&l2==1)
        {
          v_1=v_lep1+v_lep1_;
          in_mass=(v_lep1+v_lep1_).M();
          hM_l_l->Fill(in_mass);
          del_phi = DeltaPhi(v_lep1.Phi(),v_lep1_.Phi());
          del_R = DeltaR(v_lep1.Eta(),v_lep1_.Eta(),v_lep1.Phi(),v_lep1_.Phi());
          pt_ll= v_1.Pt();
          l1=0;l2=0;
        }
        if(l3==1&&l4==1)
        { v_1 = v_lep2+v_lep2_;
          in_mass=(v_lep2+v_lep2_).M();
          hM_l_l->Fill(in_mass);
          del_phi = DeltaPhi(v_lep2.Phi(),v_lep2_.Phi());
          del_R = DeltaR(v_lep2.Eta(),v_lep2_.Eta(),v_lep2.Phi(),v_lep2_.Phi());
          pt_ll= v_1.Pt();
          l4=0;l3=0;
        }
        if(l3==1&&l2==1)
        {
          v_1 = v_lep2+v_lep1_;
          in_mass=(v_lep2+v_lep1_).M();
          hM_l_l->Fill(in_mass);
          del_phi= DeltaPhi(v_lep2.Phi(),v_lep1_.Phi());
          del_R = DeltaR(v_lep2.Eta(),v_lep1_.Eta(),v_lep2.Phi(),v_lep1_.Phi());
          pt_ll= v_1.Pt();
          l3=0;l2=0;
        }
        if(l1==1&&l4==1)
        {
          v_1 = v_lep1+v_lep2_;
          in_mass=(v_lep1+v_lep2_).M();
          hM_l_l->Fill(in_mass);
          del_phi= DeltaPhi(v_lep1.Phi(),v_lep2_.Phi());
          del_R = DeltaR(v_lep1.Eta(),v_lep2_.Eta(),v_lep1.Phi(),v_lep2_.Phi());
          pt_ll= v_1.Pt();
          l4=0;l1=0;
        }
if(c1==1 && c2==1)
{
          v_2 = v_xd_+v_xd;
          etmiss=v_2.Pt();
  c1=0;c2=0;
}
if(c1==1 && c3==1)
{

         v_2 = v_xd2+v_xd;
         etmiss=v_2.Pt();
c1=0;c3=0;
}
if(c1==1 && c4==1)
{
         v_2 = v_xd2_+v_xd;
         etmiss=v_2.Pt();
          c1=0;c4=0;
}
if(c2==1 && c3==1)
{
          v_2 = v_xd_+v_xd2;
          etmiss=v_2.Pt();
c2=0;c3=0;
}
if(c2==1 && c4==1)
{
          v_2 = v_xd_+v_xd2_;
          etmiss=v_2.Pt();
c2=0;c4=0;
}
if(c3==1 && c4==1)
{
         v_2 = v_xd2_+v_xd2;
         etmiss=v_2.Pt();
         c3=0;c4=0;
}     
       if ((n_xd + n_xd_) != 2)continue; // event conditions
          v_2 = v_xd1_+v_xd1;
          etmiss=v_2.Pt();
} //loop of i
l_l_ptmiss = DeltaPhi(v_1.Phi(),v_2.Phi());
float del = DeltaPhi(v_1.Phi(),v_2.Phi());
MT=sqrt(2*v_1.Pt()*v_2.Pt()*(1-cos(del)));
float val2[2]={0,0};int u=0;
for(int i=0;i<4;i++)
{
if(val1[i]==0) continue;
else 
{val2[u]=val1[i];
u++;
}
}
float val3[2]={0,0};int o=0;
for(int i=0;i<4;i++)
{
if(val[i]==0) continue;
else 
{val3[o]=val[i];
o++;
}
}
//for leading and subleading momentum
if(val3[0]>val3[1])
{
lead_l=val3[0];sublead_l=val3[1];
}
else
{
sublead_l=val3[0];lead_l=val3[1];
}
//for leading and subleading pseudorapidity
if(val2[0]>val2[1])
{
esublead_l=val2[0];elead_l=val2[1];
}
else
{
esublead_l=val2[1];elead_l=val2[0];
}
 
//if(etmiss >= 150 && new_all >= 5 && MT>=200 lead_l >= 120&& Ht >= 800) //&& MT>=200 && del_R <=1.5)//(etmiss >= 240 && new_all >= 8.5 && lead_l >= 140 && Ht >= 900)
//{
//z++;
//} //&& del_R <= 1.5 && fabs(del_phi) <=0.8 && new_all <= 8.5/* && in_mass <= 91 /*&& fabs(elead_l) <= 2.15  && lead_l <= 140  && Ht <= 900 && sublead_l <= 70 /*&& fabs(l_l_ptmiss) >= 3.0 )

if(etmiss >= 120 && new_all >= 5 && lead_l >= 90&& MT>=200&& Ht >= 500&& del_R <=1.5)//&& Ht >= 500)//&& new_all >= 5 && lead_l >= 120&& MT>=200&& Ht >= 800&& del_R <=1.5)
{z++;
}
//ws=0.000156,wbt=12.75,wbz=0.0038985,wbw=1.7805;
double wei = 0.0038985;
hPt_xd_xd->Fill(etmiss,wei);
hMt->Fill(MT,wei);
hDR_l_l->Fill(del_R,wei);
hPt_l_l->Fill(pt_ll,wei);
hphi_l_l->Fill(del_phi,wei);
hphi_l_l_ptmiss->Fill(l_l_ptmiss,wei);
hPt_sublead_l->Fill(sublead_l,wei);
hPt_lead_l->Fill(lead_l,wei);
hEta_sublead_l->Fill(esublead_l,wei);
hEta_lead_l->Fill(elead_l,wei);
hHt->Fill(Ht,wei);
new_all=etmiss/sqrt(Ht);
hTpv->Fill(new_all,wei);
bgtree5->Fill();

// --- end filling the histograms
      ff>>tt;
      line++;
      //if (event==100)  break;
      delete Id;
      delete Status;
      delete Mother1;
      delete Mother2;
      delete Color1;
      delete Color2;
      delete px;
      delete py;
      delete pz;
      delete E;
      delete m;
      delete lifetime;
      delete spin;
      for (int k=1 ; k<npart+1 ; delete v[k++]);
    }

}
  cout << " Total number of events --> " << event << endl;
  TFile *rootfile = new TFile(rootfname.c_str(),"recreate");
  hPt_xd_xd->Write();
  hphi_l_l->Write();
  hDR_l_l->Write();
  hPt_l_l->Write();
  hMt->Write();
  hHt->Write();
  hTpv->Write();
  hM_l_l->Write();
  hPt_sublead_l->Write();
  hPt_lead_l->Write();
  hEta_sublead_l->Write();
  hEta_lead_l->Write();
  hphi_l_l_ptmiss->Write();
  bgtree5->Print();
  bgtree5->Write();
  rootfile->Close();
//double cr;
//cr= z/lum;
cout<<z<<endl;
  cout << "Events with negative weight: " << negativeWeight << endl;
  cout << "lept decay = " << nlept*1.0/event << endl;
  cout << "hadr decay = " << nhadr*1.0/event << endl;
  cout << "semi decay = " << nsemi*1.0/event << endl;
  exit(0);
}
float DeltaPhi(float phi1, float phi2)
{
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  return dphi;
}
float DeltaR(float eta1, float eta2, float phi1, float phi2)
{
  float deta = eta2 - eta1;
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  float DELTAR = sqrt(pow(dphi,2)+pow(deta,2))*1.0;
  return DELTAR;
}
