#define ZMuSkimTree_cxx
#include "ZMuSkimTree2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TVector.h>
#include <iostream>
#include <TH2F.h>
#include <TGraph.h>
#include <vector>
#include <TFile.h>
#include <fstream>
#include <sstream>
#include <TBrowser.h>
#include <TSystem.h>
#include <iterator>
#include <TLegend.h>
#include <TDirectory.h>
#include <algorithm>

void ZMuSkimTree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ZMuSkimTree.C
//      root> ZMuSkimTree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   //if (maxentries>-1) nentries=maxentries;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<10;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      clog << "Entry\t" << jentry << endl;
      // if (Cut(ientry) < 0) continue;
      clog << ltSimTwinMuxOut_phi->size()<< endl;
   }
}


void ZMuSkimTree::Count()
{
  
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();
	clog << nentries<< "\n";
	//if (maxentries>-1) nentries=maxentries;
	Long64_t nbytes = 0, nb = 0;
	Long64_t phicount = 0, simphicount = 0;
  TH2F* h = new TH2F ("wheel0", "wheel0", 28, -0.5, 12.5,10,0,5);
  TH2F* h2 = new TH2F ("wheel1", "wheel1", 28, -0.5, 12.5,10,0,5);
  TH2F* h3 = new TH2F ("wheel2", "wheel2", 28, -0.5, 12.5,10,0,5);
  TH2F* h4 = new TH2F ("wheel-1", "wheel-1", 28, -0.5, 12.5,10,0,5);
  TH2F* h5 = new TH2F ("wheel-2", "wheel-2", 28, -0.5, 12.5,10,0,5);
  TH2F* Lh = new TH2F ("wheel0", "wheel0", 28, -0.5, 12.5,10,0,5);
  TH2F* Lh2 = new TH2F ("wheel1", "wheel1", 28, -0.5, 12.5,10,0,5);
  TH2F* Lh3 = new TH2F ("wheel2", "wheel2", 28, -0.5, 12.5,10,0,5);
  TH2F* Lh4 = new TH2F ("wheel-1", "wheel-1", 28, -0.5, 12.5,10,0,5);
  TH2F* Lh5 = new TH2F ("wheel-2", "wheel-2", 28, -0.5, 12.5,10,0,5);

  long int j = 0;
	Long64_t diff = 0;
  vector<int> fillingplus[10];
  vector<int> fillingless[10];
  TCanvas* c1 = new TCanvas("c1","c1");
  TCanvas* c2 = new TCanvas("c2","c2");
  //TCanvas* c3 = new TCanvas("c3","c3");
  //TCanvas* c4 = new TCanvas("c4","c4");
  //TCanvas* c5 = new TCanvas("c5","c5");
  c1->Divide(3,2);
  c2->Divide(3,2);
  //TH1F* hwh0 = new TH1F("hwh0","number of hits wheel 0",);
	for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
    	Long64_t ientry = LoadTree(jentry);
    	if (ientry < 0) break;
    	fChain->SetBranchStatus("*",0); 
    	fChain->SetBranchStatus("ltSimTwinMuxOut_phi",1);
      fChain->SetBranchStatus("ltTwinMuxOut_phi",1);
      fChain->SetBranchStatus("ltTwinMuxOut_wheel",1);
      fChain->SetBranchStatus("ltTwinMuxOut_sector",1);
      fChain->SetBranchStatus("ltTwinMuxOut_station",1);
      fChain->SetBranchStatus("ltSimTwinMuxOut_wheel",1);
      fChain->SetBranchStatus("ltSimTwinMuxOut_sector",1);
      fChain->SetBranchStatus("ltSimTwinMuxOut_station",1);

    	nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (ltSimTwinMuxOut_phi->size()!= ltTwinMuxOut_phi->size())
      {
        if (ltSimTwinMuxOut_phi->size() > ltTwinMuxOut_phi->size())
        {
          for (unsigned int i = (ltSimTwinMuxOut_phi->size()-ltTwinMuxOut_phi->size()); i <ltSimTwinMuxOut_phi->size() ; ++i)
          {
            int Swh=ltSimTwinMuxOut_wheel->at(i);
            int Ssec=ltSimTwinMuxOut_sector->at(i);
            int Ssta=ltSimTwinMuxOut_station->at(i);
            clog << i << endl;

            switch(Swh)
            {
              case 0:
              fillingplus[0].push_back(Ssec);
              fillingplus[1].push_back(Ssta);
              break;
              case 1:
              fillingplus[2].push_back(Ssec);
              fillingplus[3].push_back(Ssta);
              break;
              case 2:
              fillingplus[4].push_back(Ssec);
              fillingplus[5].push_back(Ssta);
              break;
              case -1:
              fillingplus[6].push_back(Ssec);
              fillingplus[7].push_back(Ssta);
              break;
              case -2:
              fillingplus[8].push_back(Ssec);
              fillingplus[9].push_back(Ssta);
              break;
              default: break;
            }
          }

          std::vector<int>::iterator j=fillingplus[1].begin();
          for (std::vector<int>::iterator i = fillingplus[0].begin(); i != fillingplus[0].end(); ++i)
          {
            h->Fill(*i,*j);
            ++j;
          }
          std::vector<int>::iterator k=fillingplus[3].begin();
          for (std::vector<int>::iterator kk = fillingplus[2].begin(); kk != fillingplus[2].end(); ++kk)
          {
            h2->Fill(*kk,*k);
            ++k;
          }
          std::vector<int>::iterator l=fillingplus[5].begin();
          for (std::vector<int>::iterator ll = fillingplus[4].begin(); ll != fillingplus[4].end(); ++ll)
          {
            h3->Fill(*ll,*l);
            ++l;
          }
          std::vector<int>::iterator m=fillingplus[7].begin();
          for (std::vector<int>::iterator mm = fillingplus[6].begin(); mm != fillingplus[6].end(); ++mm)
          {
            h4->Fill(*mm,*m);
            ++m;
          }
          std::vector<int>::iterator n=fillingplus[9].begin();
          for (std::vector<int>::iterator nn = fillingplus[8].begin(); nn != fillingplus[8].end(); ++nn)
          {
            h5->Fill(*nn,*n);
            ++n; 
          }
        }
        else
        {
          for (unsigned int i = (ltTwinMuxOut_phi->size()-ltSimTwinMuxOut_phi->size()); i <ltTwinMuxOut_phi->size() ; ++i)
          {
            int wh=ltTwinMuxOut_wheel->at(i);
            int sec=ltTwinMuxOut_sector->at(i);
            int sta=ltTwinMuxOut_station->at(i);
            switch(wh)
            {
              case 0:
              fillingless[0].push_back(sec);
              fillingless[1].push_back(sta);
              break;
              case 1:
              fillingless[2].push_back(sec);
              fillingless[3].push_back(sta);
              break;
              case 2:
              fillingless[4].push_back(sec);
              fillingless[5].push_back(sta);
              break;
              case -1:
              fillingless[6].push_back(sec);
              fillingless[7].push_back(sta);
              break;
              case -2:
              fillingless[8].push_back(sec);
              fillingless[9].push_back(sta);
              break;
              default: break;
            }
          }

          std::vector<int>::iterator j=fillingless[1].begin();
          for (std::vector<int>::iterator i = fillingless[0].begin(); i != fillingless[0].end(); ++i)
          {
            clog << "*j " <<*j << '\t';
            clog << "fillingh\n";
            Lh->Fill(*i,*j);
            ++j;          
          }
          std::vector<int>::iterator k=fillingless[3].begin();
          for (std::vector<int>::iterator kk = fillingless[2].begin(); kk != fillingless[2].end(); ++kk)
          {
            clog << "*k " <<*k << '\t';
            clog << "fillingh2\n";
            Lh2->Fill(*kk,*k);
            ++k;          
          }
          std::vector<int>::iterator l=fillingless[5].begin();
          for (std::vector<int>::iterator ll = fillingless[4].begin(); ll != fillingless[4].end(); ++ll)
          {
            clog << "*l " <<*l << '\t';
            clog << "fillingh3\n";
            Lh3->Fill(*ll,*l);
            ++l;         
          }
          std::vector<int>::iterator m=fillingless[7].begin();
          for (std::vector<int>::iterator mm = fillingless[6].begin(); mm != fillingless[6].end(); ++mm)
          {
            clog << "*m " <<*m << '\t';
            clog << "fillingh4\n";
            Lh4->Fill(*mm,*m);
            ++m;    
          }
          std::vector<int>::iterator n=fillingless[9].begin();
          for (std::vector<int>::iterator nn = fillingless[8].begin(); nn != fillingless[8].end(); ++nn)
          {
            clog << "*n " <<*n << '\t';
            clog << "fillingh5\n";
            Lh5->Fill(*nn,*n);
            ++n;  
          }

        }
        fillingplus[0].clear();
        fillingplus[1].clear();
        fillingplus[2].clear();
        fillingplus[3].clear();
        fillingplus[4].clear();
        fillingplus[5].clear();
        fillingplus[6].clear();
        fillingplus[7].clear();
        fillingplus[8].clear();
        fillingplus[9].clear();
        fillingless[0].clear();
        fillingless[1].clear();
        fillingless[2].clear();
        fillingless[3].clear();
        fillingless[4].clear();
        fillingless[5].clear();
        fillingless[6].clear();
        fillingless[7].clear();
        fillingless[8].clear();
        fillingless[9].clear();
      }
    }	
    
    

    c1->cd(1);
    h->Draw(); 
    c1->cd(2);
    h2->Draw();    
    c1->cd(3);
    h3->Draw();   
    c1->cd(4);
    h4->Draw();
    c1->cd(5);
    h5->Draw();
    c1->Draw();

    c2->cd(1);
    Lh->Draw(); 
    c2->cd(2);
    Lh2->Draw();    
    c2->cd(3);
    Lh3->Draw();   
    c2->cd(4);
    Lh4->Draw();
    c2->cd(5);
    Lh5->Draw();
    c2->Draw();    
}


void ZMuSkimTree::PhiWheelDist()
{
  if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
   //if (maxentries>-1) nentries=maxentries;
    vector<float> phi[10];
    TCanvas* c3 = new TCanvas("c3","c3");
    //TCanvas* c2 = new TCanvas("c2","c2");
    c3->Divide(3,2);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      clog <<"jentry" << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->SetBranchStatus("*",0); 
      fChain->SetBranchStatus("ltSimTwinMuxOut_phi",1);
      fChain->SetBranchStatus("ltTwinMuxOut_phi",1);
      fChain->SetBranchStatus("ltTwinMuxOut_wheel",1);
      fChain->SetBranchStatus("ltTwinMuxOut_sector",1);
      fChain->SetBranchStatus("ltTwinMuxOut_station",1);
      fChain->SetBranchStatus("ltSimTwinMuxOut_wheel",1);
      fChain->SetBranchStatus("ltSimTwinMuxOut_sector",1);
      fChain->SetBranchStatus("ltSimTwinMuxOut_station",1);

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      if (ltSimTwinMuxOut_phi->size()== ltTwinMuxOut_phi->size())
      {
        clog<< "primo if\n";
        for (unsigned int i = 0; i < ltSimTwinMuxOut_phi->size(); ++i)
        {
          clog<<"for\n";
          switch(ltSimTwinMuxOut_wheel->at(i))
          {
            case 0:
            clog<<"caso0\n";
            phi[0].push_back(ltTwinMuxOut_phi->at(i));
            phi[1].push_back(ltSimTwinMuxOut_phi->at(i));
            break;
            case 1:
            clog<<"caso1\n";
            phi[2].push_back(ltTwinMuxOut_phi->at(i));
            phi[3].push_back(ltSimTwinMuxOut_phi->at(i));
            break;
            case 2:
            clog<<"caso2\n";
            phi[4].push_back(ltTwinMuxOut_phi->at(i));
            phi[5].push_back(ltSimTwinMuxOut_phi->at(i));
            break;
            case -1:
            clog<<"caso-1\n";
            phi[6].push_back(ltTwinMuxOut_phi->at(i));
            phi[7].push_back(ltSimTwinMuxOut_phi->at(i));
            break;
            case -2:
            clog<<"caso-2\n";
            phi[8].push_back(ltTwinMuxOut_phi->at(i));
            phi[9].push_back(ltSimTwinMuxOut_phi->at(i));
            break;
            default:
            break;
          }
        }
      }
      else 
      continue;
   }
    TGraph* g1 = new TGraph(phi[0].size(),phi[0].data(),phi[1].data());
    TGraph* g2 = new TGraph(phi[2].size(),phi[2].data(),phi[3].data());
    TGraph* g3 = new TGraph(phi[4].size(),phi[4].data(),phi[5].data());
    TGraph* g4 = new TGraph(phi[6].size(),phi[6].data(),phi[7].data());
    TGraph* g5 = new TGraph(phi[8].size(),phi[8].data(),phi[9].data());
    c3->cd(1);
    g1->Draw("AP");
    c3->cd(2);
    g2->Draw("AP");
    c3->cd(3);
    g3->Draw("AP");
    c3->cd(4);
    g4->Draw("AP");
    c3->cd(5);
    g5->Draw("AP");
}

void ZMuSkimTree::PhiWheelDist(short w,short sc) 
{
  if (abs(w)>2) return;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //if (maxentries>-1) nentries=maxentries;
  vector<float> phi[2];
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      clog << jentry <<'\t';
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (ltSimTwinMuxOut_phi->size()== ltTwinMuxOut_phi->size())
      {
        for (unsigned int i = 0; i < ltSimTwinMuxOut_phi->size(); ++i)
        {

          if (ltSimTwinMuxOut_wheel->at(i) == w && ltTwinMuxOut_wheel->at(i) == w && ltSimTwinMuxOut_sector->at(i) == sc && ltTwinMuxOut_sector->at(i) == sc)
          {
            phi[0].push_back(ltTwinMuxOut_phi->at(i));
            phi[1].push_back(ltSimTwinMuxOut_phi->at(i));
          }
          else
          continue;
        }
      }
      else
        continue;
    }
    gSystem->Exec("echo  \'\a\' \'\a\' \'\a\' \'\a\' \'\a\' \'\a\' \'\a\' " );
   

    TGraph* g1 = new TGraph(phi[0].size(),phi[0].data(),phi[1].data());
      g1->Draw("AP");

}

void ZMuSkimTree::PhiWheelDist(short w, short sc, short st) 
{
  if (abs(w)>2 || abs(sc)>12 ) return;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //if (maxentries>-1) nentries=maxentries;
  vector<float> phi[7000];
  
  TGraph* g[3500];
   unsigned int iter=100;
   Long64_t nbytes = 0, nb = 0;
   Long64_t nbytes2 = 0, nb2 = 0;
   int* bxs = new int[iter];
   for (Long64_t jentry2=0; jentry2<iter;jentry2++) 
    {
      Long64_t ientry2 = LoadTree(jentry2);
      if (ientry2 < 0) break;
      //clog << "/////////////////////////////////////////////////" <<jentry2 << endl;
      nb2 = fChain->GetEntry(jentry2);   nbytes2 += nb2;
      fChain->SetBranchStatus("*",0); 
      fChain->SetBranchStatus("bunchXing",1);
      bxs[jentry2] = bunchXing;
    }

    for (unsigned int j=0; j<iter; j++)
    {
   for (Long64_t jentry=0; jentry<iter;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      fChain->SetBranchStatus("*",1); 
      
      if (ltSimTwinMuxOut_phi->size()== ltTwinMuxOut_phi->size())
      {
        for (unsigned int i = 0; i < ltSimTwinMuxOut_phi->size(); ++i)
        {
          if (ltTwinMuxOut_wheel->at(i) == w && ltTwinMuxOut_sector->at(i) == sc && ltTwinMuxOut_station->at(i) == st && (bunchXing + ltTwinMuxOut_bx->at(i)) == bxs[j])
          {
            phi[bxs[j]].push_back(ltTwinMuxOut_phi->at(i));
          }
          if (ltSimTwinMuxOut_wheel->at(i) == w && ltSimTwinMuxOut_sector->at(i) == sc && ltSimTwinMuxOut_station->at(i) == st && (bunchXing + ltSimTwinMuxOut_bx->at(i)) == bxs[j])
          {
            phi[2*bxs[j] - bxs[j]].push_back(ltSimTwinMuxOut_phi->at(i));
          }
          else
          continue;
        }
      
      }
      else
      continue;
    
   }
 }
 
 string f = "bx";
for (unsigned int bxsa = 0; bxsa < iter; ++bxsa)
 {
   
  
    if (phi[bxs[bxsa]].size())
    {
      string path = "/home/marco/root/macros/";
      string nam = path +"bxprscatter_w"+ to_string(w)+ "sc"+ to_string(sc)+ "st" + to_string(st) +".root";
      char const *addr = nam.c_str();
      TFile* out = new TFile(addr, "UPDATE");
      if ( out->IsOpen() ) printf("File opened successfully\n");
      g[bxs[bxsa]] = new TGraph(phi[bxs[bxsa]].size(),phi[bxs[bxsa]].data(),phi[2*(bxs[bxsa]) - bxs[bxsa]].data());
      clog << "real\t" << phi[bxs[bxsa]].data()[0] << '\t' << "sim\t" << phi[2*(bxs[bxsa]) - bxs[bxsa]].data()[0] << endl;
      string s = f + to_string(bxs[bxsa]);
      clog << s << endl; 
      char const *pchar = s.c_str();
      g[bxs[bxsa]]->SetMarkerStyle(20);
      g[bxs[bxsa]]->SetLineWidth(0);
      g[bxs[bxsa]]->Write(pchar);
      out->Write();
      out->Close();
    }
 }
      
    





    

}
void ZMuSkimTree::Scanner() 
{
  /////////////////Aggiungi bx e altro dati
  //if (abs(w)>2 || abs(se)>12 ) return;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //if (maxentries>-1) nentries=maxentries;
  //vector<float> phi[2];
   Long64_t nbytes = 0, nb = 0;
   int second = 0;
   ofstream listoutput("outscanner.csv");
   ostringstream bufflist;
   listoutput <<"Entry\t "  << "phi \t"  << "Wheel "<< "\tSector " << "\tStation "  << "\tNprimitive"  << "\tis2nd"<< "\tBx(diffbx)"  <<  endl;

   for (Long64_t jentry=0; jentry<10000;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (ltTwinMuxOut_wheel->size())
      {
      for (int w = -2; w < 3; ++w)
      {
        for (int sc = 1; sc < 13; ++sc)
        {
          
          for (int st = 1; st < 6; ++st)
          {
            int d=0;
        
        

            
              for (unsigned int i = 0; i < ltTwinMuxOut_wheel->size(); ++i)
              {
             
                if (ltTwinMuxOut_wheel->at(i) == w && ltTwinMuxOut_sector->at(i) == sc && ltTwinMuxOut_station->at(i) == st && this->Finder(w,sc,st) == 0)
                {
                  ++d;
                  listoutput << jentry  << '\t' <<ltTwinMuxOut_phi->data()[i]  << '\t'<< w   << '\t'<< sc  << '\t'<< st  << '\t'<< d  << '\t'<< ltTwinMuxOut_is2nd->data()[i]  << '\t' << (bunchXing + ltTwinMuxOut_bx->data()[i]) << '('<< ltTwinMuxOut_bx->data()[i] << ')' << "\t\t" << ltTwinMuxOut_rpcbit->data()[i] <<  endl;
                  
                  /*for (unsigned int k = 0; k <  ltSimTwinMuxOut_phi->size(); ++k)
                  {
                  if (ltSimTwinMuxOut_wheel->at(k) == w && ltSimTwinMuxOut_sector->at(k) == sc && ltSimTwinMuxOut_station->at(k) == st && d==1) 
                    {
                      listoutput << '\t' <<ltSimTwinMuxOut_phi->data()[k]  << '\t'<< w   << '\t'<< sc  << '\t'<< st  << '\t'<< d  << '\t'<< ltSimTwinMuxOut_is2nd->data()[k]  << '\t' << (bunchXing + ltSimTwinMuxOut_bx->data()[k]) << '('<< ltSimTwinMuxOut_bx->data()[k] << ')' << "\tSIMULATED" << '\t' << ltSimTwinMuxOut_rpcbit->data()[k] <<  endl;
                      ++second;
                    }
                    else
                      continue;
                  }*/
            
                }

              else
              continue;


              }
              
              
            }

          }
        }
      }
      else
        continue;

    }
    listoutput.close();
    clog << "is2nd\t" << second << endl;
}

void ZMuSkimTree::ScannerInv() 
{
  /////////////////Aggiungi bx e altro dati
  //if (abs(w)>2 || abs(se)>12 ) return;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //if (maxentries>-1) nentries=maxentries;
  //vector<float> phi[2];
   Long64_t nbytes = 0, nb = 0;
   int second = 0;
   ofstream listoutput("outputrpcbitinv.csv");
   ostringstream bufflist;
   listoutput <<"Entry\t "  << "phi \t"  << "Wheel "<< "\tSector " << "\tStation "  << "\tNprimitive"  << "\tis2nd"<< "\tBx(diffbx)"  <<  endl;

   for (Long64_t jentry=0; jentry<100000;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (ltTwinMuxOut_phi->size())
      {
      for (int w = -2; w < 3; ++w)
      {
        for (int sc = 1; sc < 13; ++sc)
        {
          
          for (int st = 1; st < 6; ++st)
          {
            int d=0;
        
        

            for (unsigned int k = 0; k <  ltSimTwinMuxOut_phi->size(); ++k)
                  {
                  if (ltSimTwinMuxOut_wheel->at(k) == w && ltSimTwinMuxOut_sector->at(k) == sc && ltSimTwinMuxOut_station->at(k) == st ) 
                    {
                      listoutput << '\t' <<ltSimTwinMuxOut_phi->data()[k]  << '\t'<< w   << '\t'<< sc  << '\t'<< st  << '\t'<< d  << '\t'<< ltSimTwinMuxOut_is2nd->data()[k]  << '\t' << (bunchXing + ltSimTwinMuxOut_bx->data()[k]) << '('<< ltSimTwinMuxOut_bx->data()[k] << ')' << "\tSIMULATED" << '\t' << ltSimTwinMuxOut_rpcbit->data()[k] <<  endl;
                      ++second;
              for (unsigned int i = 0; i < ltTwinMuxOut_phi->size(); ++i)
              {
             
                if (ltTwinMuxOut_wheel->at(i) == w && ltTwinMuxOut_sector->at(i) == sc && ltTwinMuxOut_station->at(i) == st)
                {
                  
                  listoutput << jentry  << '\t' <<ltTwinMuxOut_phi->data()[i]  << '\t'<< w   << '\t'<< sc  << '\t'<< st  << '\t'<< d  << '\t'<< ltTwinMuxOut_is2nd->data()[i]  << '\t' << (bunchXing + ltTwinMuxOut_bx->data()[i]) << '('<< ltTwinMuxOut_bx->data()[i] << ')' << "\t\t" << ltTwinMuxOut_rpcbit->data()[i] <<  endl;
                  
                  
                    }
                    else
                      continue;
                  }
            
                }

              else
              continue;


              }
              
              
            }

          }
        }
      }
      else
        continue;

    }
    listoutput.close();
    clog << "is2nd\t" << second << endl;
}
void ZMuSkimTree::PhiDiff()
{
  //DISTRIBUZIONE DIFFERENZE NELLE PHI
  if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
   //if (maxentries>-1) nentries=maxentries;
    vector<float> phi[2];
    TCanvas* c3 = new TCanvas("c3","c3");
    TH1F* h = new TH1F ("h", "h", 28, -0.5, 12.5);

    

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<10000;jentry++) 
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
    if (ltSimTwinMuxOut_phi->size() == ltTwinMuxOut_phi->size())
    {
      for (unsigned int i = 0; i < ltSimTwinMuxOut_phi->size(); ++i)
      {
        phi[0].push_back(abs(ltSimTwinMuxOut_phi->at(i) - ltTwinMuxOut_phi->at(i)));
      }

    }

    for (std::vector<float>::iterator i = phi[0].begin(); i != phi[0].end(); ++i)
      {
        h->Fill(*i);
      }

  }
  h->Draw();

}



void ZMuSkimTree::asd()
{
 for (int w = -2; w < 3; ++w)
      {
        for (int sc = 1; sc < 13; ++sc)
        {
          
          for (int st = 1; st < 6; ++st)
          {
            clog << "iteration\t" << w << '\t'<< sc << '\t'<<st << endl;
            this->PhiWheelDist(w,sc,st);
          }
        }
      }
}


void ZMuSkimTree::SkimmerNoDT()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //if (maxentries>-1) nentries=maxentries;
  //vector<float> phi[2];
   Long64_t nbytes = 0, nb = 0;
   //ofstream listoutput("70kreal.csv");
   //ofstream listoutput2("70kemulover2.csv");

   //listoutput <<"Entry\t "  << "phi \t"  << "Wheel "<< "\tSector " << "\tStation "  << "\tBx" << endl;
   //listoutput2 <<"Entry\t "  << "phi \t"  << "Wheel "<< "\tSector " << "\tStation "  << "\tBx" << endl;
   //////////////////////Vector declarations
   vector<int> skimbx;
   vector<int> simskimbx;
   vector<int> skimentry;
   vector<int> simskimentry;
   vector<float> skimphi;
   vector<float> simskimphi;
   vector<float> skimphiB;
   vector<float> simskimphiB;
   vector<int> wheel;
   vector<int> sector;
   vector<int> station;
   //vector<float> pos;
   //vector<float> dir;
   vector<int> simwheel;
   vector<int> simsector;
   vector<int> simstation;
   //vector<float> simpos;
   //vector<float> simdir;
   vector<int> sample;
   vector<int> asd;
   //vector<float> segmpos;
   //vector<float> segmdirx;
   //vector<float> segmdirz;



  ////////////////Vectors overcounts
  vector<int> pskimbx;
  vector<int> pwheel;
  vector<int> psector;
  vector<int> pstation;
  vector<int> pskimentry;



  vector<int> qsimskimbx;
  vector<int> qsimwheel;
  vector<int> qsimsector;
  vector<int> qsimstation;
  vector<int> qsimskimentry;
  

  ////////////////Vectors equalcounts

  vector<int> tskimbx;
  vector<int> twheel;
  vector<int> tsector;
  vector<int> tstation;
  vector<int> tskimentry;
  vector<float> tskimphi;
  vector<float> tskimphiB;
  vector<float> tsimskimphi;
  vector<float> tsimskimphiB;
  vector<int> tsimskimbx;
  //vector<float> tpos;
  //vector<float> tdir;
  //vector<float> tsimpos;
  //vector<float> tsimdir;



  /////////Real and emulated data reading



   
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

     
      if (ltTwinMuxOut_phi->size())
      {
        for (int w = -2; w < 3; ++w)
        {
          for (int sc = 1; sc < 13; ++sc)
          {
          
            for (int st = 1; st < 5; ++st)
            {
              

              for (unsigned int i = 0; i < ltTwinMuxOut_phi->size(); ++i)
                {
             
                  if (ltTwinMuxOut_wheel->at(i) == w && ltTwinMuxOut_sector->at(i) == sc && ltTwinMuxOut_station->at(i) == st && ltTwinMuxOut_phi->at(i) != 2044)
                  {
                    if(!ltTwinMuxIn_phi || !(this->Finder(w,sc,st)))
                    {
                    //listoutput << jentry  << '\t' <<ltTwinMuxOut_phi->at(i)  << '\t'<< w   << '\t'<< sc  << '\t'<< st   << '\t' << ltTwinMuxOut_bx->at(i) << endl;
                    skimbx.push_back(ltTwinMuxOut_bx->at(i));
                    skimphi.push_back(ltTwinMuxOut_phi->at(i));
                    skimphiB.push_back(ltTwinMuxOut_phiB->at(i));
                    skimentry.push_back(jentry);
                    wheel.push_back(w);
                    sector.push_back(sc);
                    station.push_back(st);
                    ////pos.push_back(ltTwinMuxOut_pos->at(i));
                    ////dir.push_back(ltTwinMuxOut_dir->at(i));
                   
                    }
                  
                  }
                  else
                    continue;
                }
             
            }

          }
        }
      }
     if (ltSimTwinMuxOut_phi->size())
      {
        for (int w = -2; w < 3; ++w)
        {
          for (int sc = 1; sc < 13; ++sc)
          {
          
            for (int st = 1; st < 5; ++st)
            {
              

              for (unsigned int i = 0; i < ltSimTwinMuxOut_phi->size(); ++i)
                {
             
                  if (ltSimTwinMuxOut_wheel->at(i) == w && ltSimTwinMuxOut_sector->at(i) == sc && ltSimTwinMuxOut_station->at(i) == st)
                  {
                    if(!ltTwinMuxIn_phi || !(this->Finder(w,sc,st)))
                    {
                    //listoutput2  << jentry  << '\t' <<ltSimTwinMuxOut_phi->at(i)  << '\t'<< w   << '\t'<< sc  << '\t'<< st   << '\t' << ltSimTwinMuxOut_bx->at(i) << endl;
                    simskimbx.push_back(ltSimTwinMuxOut_bx->at(i));
                    simskimphi.push_back(ltSimTwinMuxOut_phi->at(i));
                    simskimphiB.push_back(ltSimTwinMuxOut_phiB->at(i));
                    simskimentry.push_back(jentry);
                    simwheel.push_back(w);
                    simsector.push_back(sc);
                    simstation.push_back(st);
                    //simpos.push_back(ltSimTwinMuxOut_pos->at(i));
                    //simdir.push_back(ltSimTwinMuxOut_dir->at(i));
                    
                    }
                  
                  }
                  else
                    continue;
                }
              }
            }
          }
        } 
    
    }

    /*vector<int>::iterator j=simskimbx.begin();
    for (std::vector<int>::iterator i = skimbx.begin(); i != skimbx.end(); ++i)
    {
     listoutput << *i << '\t' << *j << endl;
     ++j; 
    }*/
    
      /*
      Pensa a farlo coi vector ma non sono sicuro sia il modo migliore e funzionante
      Valuta se è meglio fare questo confronto tra i tripletti rilevati con tutti quelli che non hanno dato input
      all'interno della entry (molto probabilmente sì poichè sono storie diverse una entry dall'altra) 
      Una volta trovate le tracce che non hanno input (da fare separatamente con reali e simulati) confrontare il bx.
      Poi devo trovare tutti gli input con un solo DT input e vedere come il Twin ha corretto e come il simulatore.
      Ricordati la storia degli rpcbit (se serve chiedi)
      MUOVITI E NON PERDERE TEMPO CAZZO

      */
  
////////////////////////Overcount algorithm
clog << "starting overcounts\n";

////Real
int ov;
for (std::vector<int>::iterator i = skimentry.begin(); i != skimentry.end(); ++i)
{
  ov = 0;
  sample.push_back(*i);
  sample.push_back(wheel.at(distance(skimentry.begin(),i)));
  sample.push_back(sector.at(distance(skimentry.begin(),i)));
  sample.push_back(station.at(distance(skimentry.begin(),i)));
  sample.push_back(skimbx.at(distance(skimentry.begin(),i)));

  for (std::vector<int>::iterator j = simskimentry.begin(); j != simskimentry.end(); ++j)
  {
    asd.push_back(*j);
    asd.push_back(simwheel.at(distance(simskimentry.begin(),j)));
    asd.push_back(simsector.at(distance(simskimentry.begin(),j)));
    asd.push_back(simstation.at(distance(simskimentry.begin(),j)));
    asd.push_back(simskimbx.at(distance(simskimentry.begin(),j)));
    if (sample == asd)
    {
      ++ov;
      asd.clear();
      break;
    }
    else asd.clear();

  }
  if(ov==0)
  {
    pskimbx.push_back(skimbx.at(distance(skimentry.begin(),i)));
    pwheel.push_back(wheel.at(distance(skimentry.begin(),i)));
    psector.push_back(sector.at(distance(skimentry.begin(),i)));
    pstation.push_back(station.at(distance(skimentry.begin(),i)));
    pskimentry.push_back(*i);
  }
  sample.clear();
}

//Emul

for (std::vector<int>::iterator i = simskimentry.begin(); i != simskimentry.end(); ++i)
{
  ov = 0;
  sample.push_back(*i);
  sample.push_back(simwheel.at(distance(simskimentry.begin(),i)));
  sample.push_back(simsector.at(distance(simskimentry.begin(),i)));
  sample.push_back(simstation.at(distance(simskimentry.begin(),i)));
  sample.push_back(simskimbx.at(distance(simskimentry.begin(),i)));
  for (std::vector<int>::iterator j = skimentry.begin(); j != skimentry.end(); ++j)
  {
    asd.push_back(*j);
    asd.push_back(wheel.at(distance(skimentry.begin(),j)));
    asd.push_back(sector.at(distance(skimentry.begin(),j)));
    asd.push_back(station.at(distance(skimentry.begin(),j)));
    asd.push_back(skimbx.at(distance(skimentry.begin(),j)));
    if (sample == asd)
    {
      ++ov;
      asd.clear();
      break;
    }
    else asd.clear();

  }
  if(ov==0)
  {
    qsimskimbx.push_back(simskimbx.at(distance(simskimentry.begin(),i)));
    qsimwheel.push_back(simwheel.at(distance(simskimentry.begin(),i)));
    qsimsector.push_back(simsector.at(distance(simskimentry.begin(),i)));
    qsimstation.push_back(simstation.at(distance(simskimentry.begin(),i)));
    qsimskimentry.push_back(*i);
  }
  sample.clear();
}



///////////////////////Equalcount algorithm
clog << "starting equalcount\n";

for (std::vector<int>::iterator i = skimentry.begin(); i != skimentry.end(); ++i)
{
  
  sample.push_back(*i);
  sample.push_back(wheel.at(distance(skimentry.begin(),i)));
  sample.push_back(sector.at(distance(skimentry.begin(),i)));
  sample.push_back(station.at(distance(skimentry.begin(),i)));
  sample.push_back(skimbx.at(distance(skimentry.begin(),i)));
  for (std::vector<int>::iterator j = simskimentry.begin(); j != simskimentry.end(); ++j)
  {
    asd.push_back(*j);
    asd.push_back(simwheel.at(distance(simskimentry.begin(),j)));
    asd.push_back(simsector.at(distance(simskimentry.begin(),j)));
    asd.push_back(simstation.at(distance(simskimentry.begin(),j)));
    asd.push_back(simskimbx.at(distance(simskimentry.begin(),j)));
    if (sample == asd)
    {
      tskimbx.push_back(skimbx.at(distance(skimentry.begin(),i)));
      twheel.push_back(wheel.at(distance(skimentry.begin(),i)));
      tsector.push_back(sector.at(distance(skimentry.begin(),i)));
      tstation.push_back(station.at(distance(skimentry.begin(),i)));
      tskimphi.push_back(skimphi.at(distance(skimentry.begin(),i)));
      tskimphiB.push_back(skimphiB.at(distance(skimentry.begin(),i)));
      //tpos.push_back(pos.at(distance(skimentry.begin(),i)));
      //tdir.push_back(dir.at(distance(skimentry.begin(),i)));
      tsimskimphi.push_back(simskimphi.at(distance(simskimentry.begin(),j)));
      tsimskimphiB.push_back(simskimphiB.at(distance(simskimentry.begin(),j)));
      tsimskimbx.push_back(simskimbx.at(distance(simskimentry.begin(),j)));
      //tsimpos.push_back(simpos.at(distance(simskimentry.begin(),j)));
      //tsimdir.push_back(simdir.at(distance(simskimentry.begin(),j)));
      tskimentry.push_back(*i);

      asd.clear();
      break;
    }
    else asd.clear();

  }
  sample.clear();
}
////////////////////////////LAMBDA FUNCTION to see which data is nearer to the reconstructed segments
/*
[this,nentries,&tskimentry, &twheel,&tsector,&tstation,&segmdirx,&segmdirz,&segmpos]()
{



   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = this->LoadTree(jentry);
      if (ientry < 0) break;
      nb = this->fChain->GetEntry(jentry);   nbytes += nb;

      for (std::vector<int>::iterator i = tskimentry.begin(); *i <= jentry; ++i)
      {
        if (*i == jentry)
        {
          clog <<"YATTA!\n";
          for (unsigned int j = 0; j < dtsegm4D_wheel->size(); ++j)
          {

            if (dtsegm4D_wheel->at(j) == twheel.at(distance(tskimentry.begin(),i)) && dtsegm4D_sector->at(j) == tsector.at(distance(tskimentry.begin(),i)) && dtsegm4D_station->at(j) == tstation.at(distance(tskimentry.begin(),i)) && dtsegm4D_phinhits->at(j)>5)
            {
              segmpos.push_back(dtsegm4D_x_pos_loc->at(j));
              segmdirx.push_back(dtsegm4D_x_dir_loc->at(j));
              segmdirz.push_back(dtsegm4D_z_dir_loc->at(j));
clog << this->Finder(dtsegm4D_wheel->at(j),dtsegm4D_sector->at(j),dtsegm4D_station->at(j)) << endl;
            }
          }
        }
      }
     clog <<"LAMBDA FUNCTION FTW!! " << jentry << endl; 

   }

}();
     clog <<"tpos " << tpos.size() << endl;
     clog << "segmpos " << segmpos.size() << endl;
TH1F* hh = new TH1F("hh","hh",500,-100,100);
TH1F* hhsim = new TH1F("hhsim","hhsim",500,-100,100);
TH1F* hhh = new TH1F("hhh","hhh",500,-100,100);
TH1F* hhhsim = new TH1F("hhh","hhh",500,-100,100);
TCanvas* ccc = new TCanvas("ccc","ccc"); 
for (std::vector<float>::iterator i = segmpos.begin(); i != segmpos.end(); ++i)
{
  hh->Fill(tpos.at(distance(segmpos.begin(),i))-*i);
}

for (std::vector<float>::iterator i = segmpos.begin(); i != segmpos.end(); ++i)
{
  hhsim->Fill(tsimpos.at(distance(segmpos.begin(),i))-*i);
}

for (std::vector<float>::iterator i = segmdirx.begin(); i != segmdirx.end(); ++i)
{
  hhh->Fill(tdir.at(distance(segmdirx.begin(),i))-180*atan(*i/segmdirz.at(distance(segmdirx.begin(),i)))/3.1415);
}
for (std::vector<float>::iterator i = segmdirx.begin(); i != segmdirx.end(); ++i)
{
  hhhsim->Fill(tsimdir.at(distance(segmdirx.begin(),i))-180*atan(*i/segmdirz.at(distance(segmdirx.begin(),i)))/3.1415);
}
//ltTwinMuxIn_dir-180.0*atan(dtsegm4D_x_dir_loc/dtsegm4D_z_dir_loc)/3.1415>>h(200,-10,10)

ccc->Divide(2,2);
ccc->cd(1);
hh->Draw();
ccc->cd(2);
hhh->Draw();
ccc->cd(3);
hhsim->Draw();
ccc->cd(4);
hhhsim->Draw();
ccc->Draw();
*/


//listoutput.close();
//listoutput2.close();

////////////////Delta Phis
float u;
vector<float> phid;
vector<float> phibd;
TH1F* hphid = new TH1F("","",100,-0.5,100.5);
TH1F* hphibd = new TH1F("","",100,-0.5,99.5);
for (std::vector<float>::iterator i = tskimphi.begin(); i != tskimphi.end(); ++i)
{
  u = sqrt((*i - tsimskimphi.at(distance(tskimphi.begin(),i)))*(*i - tsimskimphi.at(distance(tskimphi.begin(),i))));
  phid.push_back(u);
  hphid->Fill(u);
}
u = 0;
for (std::vector<float>::iterator i = tskimphiB.begin(); i != tskimphiB.end(); ++i)
{
  u = sqrt((*i - tsimskimphiB.at(distance(tskimphiB.begin(),i)))*(*i - tsimskimphiB.at(distance(tskimphiB.begin(),i))));
  phibd.push_back(u);
  hphibd->Fill(u);
}

/////////////Drawing

this->Histofiller(station,sector,wheel, "real");
this->Histofiller(simstation,simsector,simwheel,  "emulated");
this->Histofiller(pstation,psector,pwheel, "real_over");
this->Histofiller(qsimstation,qsimsector,qsimwheel,  "emulated_over");
//this->Histofiller(tsimstation,tsimsector,tsimwheel,  "emulated_equal");
this->Histofiller(tstation,tsector,twheel,  "real_equal");

TFile* outg = new TFile("/home/marco/root/macros/graph.root", "UPDATE");
  if ( outg->IsOpen() ) printf("File opened successfully\n");


hphid->Write("phidiff");
hphibd->Write("phibdiff");



TGraph* gg = new TGraph(tskimphi.size(),tskimphi.data(),tsimskimphi.data());
TGraph* gg2 = new TGraph(tskimphiB.size(),tskimphiB.data(),tsimskimphiB.data());
//TGraph* gg3 = new TGraph(tskimbx.size(),tskimbx.data(),tsimskimbx.data());
TH2F* h1 = new TH2F ("BunchXing","BunchXing",7,-3.5,3.5,7,-3.5,3.5);
for (std::vector<int>::iterator i = tskimbx.begin(); i != tskimbx.end(); ++i)
{
  h1->Fill(*i,tsimskimbx.at(i - tskimbx.begin()));
}

gg->SetDrawOption("ap");
gg2->SetDrawOption("ap");
h1->SetOption("colztext");

gg->Write("Phi");
gg2->Write("PhiBending");
h1->Write("BunchXing");
outg->Write();
outg->Close();


}

int ZMuSkimTree::Finder(short w, short sc, short st)
{
  vector<int> wss;
  vector<int> sample;
  wss.push_back(w);
  wss.push_back(sc);
  wss.push_back(st);
  int d= 0;
  for (unsigned int i = 0; i < ltTwinMuxIn_phi->size(); ++i)
  {
    sample.push_back(ltTwinMuxIn_wheel->at(i));
    sample.push_back(ltTwinMuxIn_sector->at(i));
    sample.push_back(ltTwinMuxIn_station->at(i));
    if(sample == wss) ++d;
    sample.clear();
  }
  
  wss.clear();
  if(d == 0)
  {
    return 0;
    
  }
  
  else
  {
    return d;
    
  }

}
void ZMuSkimTree::InputScanner()
{
if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //if (maxentries>-1) nentries=maxentries;
  //vector<float> phi[2];
   Long64_t nbytes = 0, nb = 0;
   ofstream listoutput("inputscanner.csv");
   ostringstream bufflist;
   listoutput <<"Entry\t "  << "Wheel "<< "\tSector " << "\tStation " <<  endl;
   vector<int> noinp;
   for (Long64_t jentry=0; jentry<100;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


        for (int w = -2; w < 3; ++w)
        {
          for (int sc = 1; sc < 13; ++sc)
          {
          
            for (int st = 1; st < 6; ++st)
            {
              if(!ltTwinMuxIn_phi || !(this->Finder(w,sc,st)))
              {
                listoutput << jentry  << '\t'<< w   << '\t'<< sc  << '\t'<< st   <<  endl;
              }
              else
                continue;
            }
          }
        }
      }
}



void ZMuSkimTree::Skimmer1DT()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //if (maxentries>-1) nentries=maxentries;

  Long64_t nbytes = 0, nb = 0;
  ofstream listoutput("Skimmer1DT.csv");
  ofstream listoutput2("SimSkimmer1DT.csv");
  //listoutput <<"Entry\t "  << "phi \t"  << "Wheel "<< "\tSector " << "\tStation "  << "\tBx" << endl;
  //listoutput2 <<"Entry\t "  << "phi \t"  << "Wheel "<< "\tSector " << "\tStation "  << "\tBx" << endl;
  vector<int> skimbx;
  vector<int> simskimbx;
  vector<int> skimentry;
  vector<int> simskimentry;
  vector<float> skimphi;
  vector<float> simskimphi;
  vector<float> skimphiB;
  vector<float> simskimphiB;
  vector<int> wheel;
  vector<int> sector;
  vector<int> station;
  vector<int> simwheel;
  vector<int> simsector;
  vector<int> simstation;
  vector<int> sample;
  vector<int> asd;

  TCanvas* c = new TCanvas("c","c");
  TH1F* h = new TH1F("TwinMuxOut","TwinMuxOut",8,-4,4);
  TH1F* h2 = new TH1F("SimTwinMuxOut","SimTwinMuxOut",8,-4,4);
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

  
    if (ltTwinMuxOut_phi->size())
    {
      for (int w = -2; w < 3; ++w)
      {
        for (int sc = 1; sc < 13; ++sc)
        {
        
          for (int st = 1; st < 5; ++st)
          {
              

            for (unsigned int i = 0; i < ltTwinMuxOut_phi->size(); ++i)
            {
           
              if (ltTwinMuxOut_wheel->at(i) == w && ltTwinMuxOut_sector->at(i) == sc && ltTwinMuxOut_station->at(i) == st && ltTwinMuxOut_phi->at(i) != 2044)
              {
                if(this->Finder(w,sc,st) == 1)
                {
                  skimbx.push_back(ltTwinMuxOut_bx->at(i));
                  skimphi.push_back(ltTwinMuxOut_phi->at(i));
                  skimphiB.push_back(ltTwinMuxOut_phiB->at(i));
                  skimentry.push_back(jentry);
                  wheel.push_back(w);
                  sector.push_back(sc);
                  station.push_back(st);                
                }
              }
              else
                continue;
            }
          }
        }
      }
    }

  
    if (ltSimTwinMuxOut_phi->size())
      {
        for (int w = -2; w < 3; ++w)
        {
          for (int sc = 1; sc < 13; ++sc)
          {
          
            for (int st = 1; st < 5; ++st)
            {
              

              for (unsigned int i = 0; i < ltSimTwinMuxOut_phi->size(); ++i)
                {
             
                  if (ltSimTwinMuxOut_wheel->at(i) == w && ltSimTwinMuxOut_sector->at(i) == sc && ltSimTwinMuxOut_station->at(i) == st)
                  {
                    if(this->Finder(w,sc,st)==1)
                    {

                      simskimbx.push_back(ltSimTwinMuxOut_bx->at(i));
                      simskimphi.push_back(ltSimTwinMuxOut_phi->at(i));
                      simskimphiB.push_back(ltSimTwinMuxOut_phiB->at(i));
                      simskimentry.push_back(jentry);
                      simwheel.push_back(w);
                      simsector.push_back(sc);
                      simstation.push_back(st);                    
                    }
                  
                  }
                  else
                    continue;
                }
              }
            }
          }
        }






  }
  clog << "starting overcounts\n";
/*for (std::vector<int>::iterator i = skimentry.end(); i != skimentry.begin(); --i)
{
  clog << *i << endl;
}*/
vector<int> pskimbx;
vector<int> pwheel;
vector<int> psector;
vector<int> pstation;
vector<int> pskimentry;
int ov;
for (std::vector<int>::iterator i = skimentry.begin(); i != skimentry.end(); ++i)
{
  clog << distance(skimentry.begin(),i) << '\t';
  ov = 0;
  sample.push_back(*i);
  sample.push_back(wheel.at(distance(skimentry.begin(),i)));
  sample.push_back(sector.at(distance(skimentry.begin(),i)));
  sample.push_back(station.at(distance(skimentry.begin(),i)));
  sample.push_back(skimbx.at(distance(skimentry.begin(),i)));

  for (std::vector<int>::iterator j = simskimentry.begin(); j != simskimentry.end(); ++j)
  {
    asd.push_back(*j);
    asd.push_back(simwheel.at(distance(simskimentry.begin(),j)));
    asd.push_back(simsector.at(distance(simskimentry.begin(),j)));
    asd.push_back(simstation.at(distance(simskimentry.begin(),j)));
    asd.push_back(simskimbx.at(distance(simskimentry.begin(),j)));
    if (sample == asd)
    {
      ++ov;
      asd.clear();
      break;
    }
    else asd.clear();

  }
  if(ov==0)
  {
    pskimbx.push_back(skimbx.at(distance(skimentry.begin(),i)));
    pwheel.push_back(wheel.at(distance(skimentry.begin(),i)));
    psector.push_back(sector.at(distance(skimentry.begin(),i)));
    pstation.push_back(station.at(distance(skimentry.begin(),i)));
    pskimentry.push_back(*i);
  }
  sample.clear();
}
clog << endl;
/*for (std::vector<int>::iterator i = skimentry.end(); i != skimentry.begin(); --i)
{ 
  for (std::vector<int>::iterator j = simskimentry.end(); j != simskimentry.begin(); --j)
  { 

    if (*i == *j) 
    {
      
      pskimbx.erase(pskimbx.begin() + distance(skimentry.begin(),i));
      pwheel.erase(pwheel.begin() + distance(skimentry.begin(),i));
      psector.erase(psector.begin() + distance(skimentry.begin(),i));
      pstation.erase(pstation.begin() + distance(skimentry.begin(),i));
      break;
      
    }
    
  }

}*/

vector<int> qsimskimbx;
vector<int> qsimwheel;
vector<int> qsimsector;
vector<int> qsimstation;
vector<int> qsimskimentry;


for (std::vector<int>::iterator i = simskimentry.begin(); i != simskimentry.end(); ++i)
{
  clog << "sim\t" << distance(simskimentry.begin(),i) << '\t';
  ov = 0;
  sample.push_back(*i);
  sample.push_back(simwheel.at(distance(simskimentry.begin(),i)));
  sample.push_back(simsector.at(distance(simskimentry.begin(),i)));
  sample.push_back(simstation.at(distance(simskimentry.begin(),i)));
  sample.push_back(simskimbx.at(distance(simskimentry.begin(),i)));
  for (std::vector<int>::iterator j = skimentry.begin(); j != skimentry.end(); ++j)
  {
    asd.push_back(*j);
    asd.push_back(wheel.at(distance(skimentry.begin(),j)));
    asd.push_back(sector.at(distance(skimentry.begin(),j)));
    asd.push_back(station.at(distance(skimentry.begin(),j)));
    asd.push_back(skimbx.at(distance(skimentry.begin(),j)));
    if (sample == asd)
    {
      ++ov;
      asd.clear();
      break;
    }
    else asd.clear();

  }
  if(ov==0)
  {
    qsimskimbx.push_back(simskimbx.at(distance(simskimentry.begin(),i)));
    qsimwheel.push_back(simwheel.at(distance(simskimentry.begin(),i)));
    qsimsector.push_back(simsector.at(distance(simskimentry.begin(),i)));
    qsimstation.push_back(simstation.at(distance(simskimentry.begin(),i)));
    qsimskimentry.push_back(*i);
  }
  sample.clear();
}
clog << endl;










/*for (std::vector<int>::iterator i = simskimentry.end(); i != simskimentry.begin(); --i)
{

  for (std::vector<int>::iterator j = skimentry.end(); j != skimentry.begin(); --j)
  {
    if (*i == *j) 
    {
      qsimskimbx.erase(qsimskimbx.begin() + distance(simskimentry.begin(),i));
      qsimwheel.erase(qsimwheel.begin() + distance(simskimentry.begin(),i));
      qsimsector.erase(qsimsector.begin() + distance(simskimentry.begin(),i));
      qsimstation.erase(qsimstation.begin() + distance(simskimentry.begin(),i));
      break;
    }
  }
}*/
clog << "starting equalcount\n";
vector<int> tskimbx;
vector<int> twheel;
vector<int> tsector;
vector<int> tstation;
vector<int> tskimentry;
vector<float> tskimphi;
vector<float> tskimphiB;
vector<float> tsimskimphi;
vector<float> tsimskimphiB;
vector<int> tsimskimbx;
for (std::vector<int>::iterator i = skimentry.begin(); i != skimentry.end(); ++i)
{
  
  sample.push_back(*i);
  sample.push_back(wheel.at(distance(skimentry.begin(),i)));
  sample.push_back(sector.at(distance(skimentry.begin(),i)));
  sample.push_back(station.at(distance(skimentry.begin(),i)));
  sample.push_back(skimbx.at(distance(skimentry.begin(),i)));
  for (std::vector<int>::iterator j = simskimentry.begin(); j != simskimentry.end(); ++j)
  {
    asd.push_back(*j);
    asd.push_back(simwheel.at(distance(simskimentry.begin(),j)));
    asd.push_back(simsector.at(distance(simskimentry.begin(),j)));
    asd.push_back(simstation.at(distance(simskimentry.begin(),j)));
    asd.push_back(simskimbx.at(distance(simskimentry.begin(),j)));
    if (sample == asd)
    {
      tskimbx.push_back(skimbx.at(distance(skimentry.begin(),i)));
      twheel.push_back(wheel.at(distance(skimentry.begin(),i)));
      tsector.push_back(sector.at(distance(skimentry.begin(),i)));
      tstation.push_back(station.at(distance(skimentry.begin(),i)));
      tskimphi.push_back(skimphi.at(distance(skimentry.begin(),i)));
      tskimphiB.push_back(skimphiB.at(distance(skimentry.begin(),i)));
      tsimskimphi.push_back(simskimphi.at(distance(simskimentry.begin(),j)));
      tsimskimphiB.push_back(simskimphiB.at(distance(simskimentry.begin(),j)));
      tsimskimbx.push_back(simskimbx.at(distance(simskimentry.begin(),j)));
      tskimentry.push_back(*i);

      asd.clear();
      break;
    }
    else asd.clear();

  }
  sample.clear();
}



      
float u;
vector<float> phid;
vector<float> phibd;
TH1F* hphid = new TH1F("","",501,-0.5,500.5);
TH1F* hphibd = new TH1F("","",501,-0.5,500.5);
for (std::vector<float>::iterator i = tskimphi.begin(); i != tskimphi.end(); ++i)
{

  u = sqrt((*i - tsimskimphi.at(distance(tskimphi.begin(),i)))*(*i - tsimskimphi.at(distance(tskimphi.begin(),i))));
  clog << u << endl;
  phid.push_back(u);
  hphid->Fill(u);
}
u = 0;
clog << "phibd\n";
for (std::vector<float>::iterator i = tskimphiB.begin(); i != tskimphiB.end(); ++i)
{
  clog << u << endl;
  u = sqrt((*i - tsimskimphiB.at(distance(tskimphiB.begin(),i)))*(*i - tsimskimphiB.at(distance(tskimphiB.begin(),i))));
  phibd.push_back(u);
  hphibd->Fill(u);
}


this->Histofiller(station,sector,wheel, "real1DT");
this->Histofiller(simstation,simsector,simwheel,  "emulated1DT");
this->Histofiller(pstation,psector,pwheel, "real_over1DT");
this->Histofiller(qsimstation,qsimsector,qsimwheel,  "emulated_over1DT");
//this->Histofiller(tsimstation,tsimsector,tsimwheel,  "emulated_equal");
this->Histofiller(tstation,tsector,twheel,  "real_equal1DT");

TFile* outg = new TFile("/home/marco/root/macros/graph1DT.root", "UPDATE");
  if ( outg->IsOpen() ) printf("File opened successfully\n");
hphid->Write("phidiff");
hphibd->Write("phibdiff");

TGraph* gg = new TGraph(tskimphi.size(),tskimphi.data(),tsimskimphi.data());
TGraph* gg2 = new TGraph(tskimphiB.size(),tskimphiB.data(),tsimskimphiB.data());
//TGraph* gg3 = new TGraph(tskimbx.size(),tskimbx.data(),tsimskimbx.data());
TH2F* h1 = new TH2F ("BunchXing","BunchXing",7,-3.5,3.5,7,-3.5,3.5);
for (std::vector<int>::iterator i = tskimbx.begin(); i != tskimbx.end(); ++i)
{
  h1->Fill(*i,tsimskimbx.at(i - tskimbx.begin()));
}

gg->SetDrawOption("ap");
gg2->SetDrawOption("ap");
h1->SetOption("colztext");

gg->Write("Phi");
gg2->Write("PhiBending");
h1->Write("BunchXing");
outg->Write();
outg->Close();





}



void ZMuSkimTree::Histofiller(vector<int> station,vector<int> sector, vector<int> wheel, string type)
{
  


  TH2F* h1 = new TH2F ("Station 1","Station 1",5,-2.5,2.5,12,0.5,12.5);
  TH2F* h2 = new TH2F ("Station 2","Station 2",5,-2.5,2.5,12,0.5,12.5);
  TH2F* h3 = new TH2F ("Station 3","Station 3",5,-2.5,2.5,12,0.5,12.5);
  TH2F* h4 = new TH2F ("Station 4","Station 4",5,-2.5,2.5,12,0.5,12.5);



  for (std::vector<int>::iterator i = station.begin(); i != station.end(); ++i)
  {
    switch(*i)
    {
      case 1: 
      h1->Fill(wheel.at(i - station.begin()),sector.at(i - station.begin()));
      break;
      case 2:
      h2->Fill(wheel.at(i - station.begin()),sector.at(i - station.begin()));
      break;
      case 3:
      h3->Fill(wheel.at(i - station.begin()),sector.at(i - station.begin()));
      break;
      case 4:
      h4->Fill(wheel.at(i - station.begin()),sector.at(i - station.begin()));
      break;
      default:
      break;
    }
  }

  //c->Divide(2,2);
  //c->cd(1);
  //h1->Draw();
  //c->cd(2);
  //h2->Draw();
  //c->cd(3);
  //h3->Draw();
  //c->cd(4);
  //h4->Draw();

  h1->Scale((1/(h1->Integral()))*100);
  h2->Scale((1/(h2->Integral()))*100);
  h3->Scale((1/(h3->Integral()))*100);
  h4->Scale((1/(h4->Integral()))*100);










  TFile* out = new TFile("/home/marco/root/macros/hist.root", "UPDATE");
  if ( out->IsOpen() ) printf("File opened successfully\n");
  TDirectory* fold;
  if(!(out->GetDirectory(type.c_str()))) fold = out->mkdir(type.c_str());
  else fold = out->GetDirectory(type.c_str());
  fold->cd();
  TCanvas* c = new TCanvas("c","c");
  c->Divide(2,2);
  h1->SetXTitle("Wheel");
  h1->SetYTitle("Sector");
  h2->SetXTitle("Wheel");
  h2->SetYTitle("Sector");
  h3->SetXTitle("Wheel");
  h3->SetYTitle("Sector");
  h4->SetXTitle("Wheel");
  h4->SetYTitle("Sector");
  h1->SetOption("colztext");
  h2->SetOption("colztext");
  h3->SetOption("colztext");
  h4->SetOption("colztext");
  //c->cd(1);
  //h1->Draw();
  //c->cd(2);
  //h2->Draw();
  //c->cd(3);
  //h3->Draw();
  //c->cd(4);
  //h4->Draw();
  //c->Write();
  h1->Write("Section_1");
  h2->Write("Section_2");
  h3->Write("Section_3");
  h4->Write("Section_4");
  fold->cd();
  out->Write();
  out->Close();

  delete h1;
  delete h2;
  delete h3;
  delete h4;


}







////////////////////////////////*****Bisogna controllare l'analisi delle tracce 1DT
////////////////////////////////Per il confronto tra segmenti ricostruiti e dati twinmux per vedere se hanno ragione i dati o emu devo vedere se è un problema di numero di tracce per la posizione perchè l'istogramma fa proprio cagare mentre per la direzione è discreto ma comunque necessita di più entries
////////////////////////////////a proposito del modo per controllare velocemente se le tracce del twinmux sono associate a un muone veramente rilevato dovrei utilizzare i dati Mu... che sono tclonesarray e dovrei chiedere come sono organizzati perché ogni entry ce n'è almeno uno (di numero variabile) di array "vuoti" con size 1 e dato 999.
////////////////////////////////*****Ho già normalizzato le mappe e bisogna stamparle per bene.
////////////////////////////////Per la distribuzione del deltaphi opterei per una scala logaritmica poichè lo 0 (che volendo si potrebbe togliere ma direi di no) e il 4 sono belli pieni mentre altri valori hanno occorrenze minime
////////////////////////////////Mentre per i deltaphibend c'è una distribuzione più regolare e la scala logaritmica rovinerebbe tutto
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////