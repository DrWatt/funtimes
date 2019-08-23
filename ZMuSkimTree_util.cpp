#include "ZMuSkimTree2posdir.h"
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
#include <TExec.h>


void ZMuSkimTree::Loop()
{
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


int ZMuSkimTree::Finder(short w, short sc, short st,short a, short bx)
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
    if(sample == wss)
    {
      ++d;
      if(a == 1 && ltTwinMuxIn_bx->at(i) != bx) d = 20;
    }
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


void ZMuSkimTree::Histofiller(vector<int> station,vector<int> sector, vector<int> wheel, string type, short normalization, vector<int> pstation,vector<int> psector, vector<int> pwheel, string ptype)
{
  
  TCanvas* c = new TCanvas("c","c");
  TCanvas* c2 = new TCanvas("c2","c2");
  TCanvas* c3 = new TCanvas("c3","c3");
  TCanvas* c4 = new TCanvas("c4","c4");
  //c->Divide(2,1);
  //c2->Divide(2,1);
  //c3->Divide(2,2);
  TLegend* leg = new TLegend(0.1,0.1);
  TH2F* h1 = new TH2F (type.c_str(),type.c_str(),5,-2.5,2.5,12,0.5,12.5);
  TH2F* h2 = new TH2F (type.c_str(),type.c_str(),5,-2.5,2.5,12,0.5,12.5);
  TH2F* h3 = new TH2F (type.c_str(),type.c_str(),5,-2.5,2.5,12,0.5,12.5);
  TH2F* h4 = new TH2F (type.c_str(),type.c_str(),5,-2.5,2.5,12,0.5,12.5);


  TH2F* hh1 = new TH2F (ptype.c_str(),ptype.c_str(),5,-2.5,2.5,12,0.5,12.5);
  TH2F* hh2 = new TH2F (ptype.c_str(),ptype.c_str(),5,-2.5,2.5,12,0.5,12.5);
  TH2F* hh3 = new TH2F (ptype.c_str(),ptype.c_str(),5,-2.5,2.5,12,0.5,12.5);
  TH2F* hh4 = new TH2F (ptype.c_str(),ptype.c_str(),5,-2.5,2.5,12,0.5,12.5);

  //h1->SetStats(kFALSE);
  //h2->SetStats(kFALSE);
  //h3->SetStats(kFALSE);
  //h4->SetStats(kFALSE);
//
//  //hh1->SetStats(kFALSE);
//  //hh2->SetStats(kFALSE);
//  //hh3->SetStats(kFALSE);
  //hh4->SetStats(kFALSE);





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



 if(normalization)
 {


  for (std::vector<int>::iterator i = pstation.begin(); i != pstation.end(); ++i)
  {
    switch(*i)
    {
      case 1: 
      hh1->Fill(pwheel.at(i - pstation.begin()),psector.at(i - pstation.begin()));
      break;
      case 2:
      hh2->Fill(pwheel.at(i - pstation.begin()),psector.at(i - pstation.begin()));
      break;
      case 3:
      hh3->Fill(pwheel.at(i - pstation.begin()),psector.at(i - pstation.begin()));
      break;
      case 4:
      hh4->Fill(pwheel.at(i - pstation.begin()),psector.at(i - pstation.begin()));
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
  clog << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\t" << ptype.c_str() << endl;
  clog << "MB1\t" << hh1->GetEntries() << endl;
  clog << "MB2\t" << hh2->GetEntries() << endl;
  clog << "MB3\t" << hh3->GetEntries() << endl;
  clog << "MB4\t" << hh4->GetEntries() << endl;
  clog << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"<< endl;
  hh1->Divide(h1);
  hh2->Divide(h2);
  hh3->Divide(h3);
  hh4->Divide(h4);



  hh1->Scale(100);
  hh2->Scale(100);
  hh3->Scale(100);
  hh4->Scale(100);






  TFile* out = new TFile("/home/marco/root/macros/hist.root", "UPDATE");
  if ( out->IsOpen() ) printf("File opened successfully\n");
  TDirectory* fold;
  if(!(out->GetDirectory(ptype.c_str()))) fold = out->mkdir(ptype.c_str());
  else fold = out->GetDirectory(ptype.c_str());
  fold->cd();

  hh1->SetXTitle("Wheel");
  hh1->SetYTitle("Sector");
  hh2->SetXTitle("Wheel");
  hh2->SetYTitle("Sector");
  hh3->SetXTitle("Wheel");
  hh3->SetYTitle("Sector");
  hh4->SetXTitle("Wheel");
  hh4->SetYTitle("Sector");
  hh1->GetXaxis()->SetNdivisions(5);
  hh2->GetXaxis()->SetNdivisions(5);
  hh3->GetXaxis()->SetNdivisions(5);
  hh4->GetXaxis()->SetNdivisions(5);
  hh1->GetYaxis()->SetNdivisions(12);
  hh2->GetYaxis()->SetNdivisions(12);
  hh3->GetYaxis()->SetNdivisions(12);
  hh4->GetYaxis()->SetNdivisions(12);
  hh1->SetOption("colztext");
  hh2->SetOption("colztext");
  hh3->SetOption("colztext");
  hh4->SetOption("colztext");
  hh1->GetZaxis()->SetRangeUser(0,11);
  hh2->GetZaxis()->SetRangeUser(0,11);
  hh3->GetZaxis()->SetRangeUser(0,11);
  hh4->GetZaxis()->SetRangeUser(0,11);
  hh1->SetMarkerStyle(3);
  hh2->SetMarkerStyle(3);
  hh3->SetMarkerStyle(3);
  hh4->SetMarkerStyle(3);
  hh1->SetMarkerSize(1.5);
  hh2->SetMarkerSize(1.5);
  hh3->SetMarkerSize(1.5);
  hh4->SetMarkerSize(1.5);
  gStyle->SetPaintTextFormat("4.3f");
  //hh1->SetPaintTextFormat("4.3f");
  //hh2->SetPaintTextFormat("4.3f");
  //hh3->SetPaintTextFormat("4.3f");
  //hh4->SetPaintTextFormat("4.3f");
  //c->cd(1);
  //h1->Draw();
  c3->cd();
  hh1->Draw();
  ////leg->Draw();
  c4->cd();
  hh2->Draw();
  //c3->cd(2);
  //hh3->Draw();
  //c3->cd(4);
  //hh4->Draw();
  ////leg->Draw();
  //c->cd(4);
  //h4->Draw();
  //c->Write();
  //hh1->Write("Stationo_1");
  //hh2->Write("Stationo_2");
  //hh3->Write("Stationo_3");
  //hh4->Write("Stationo_4");
  fold->cd();
  out->Write();
  out->Close();
  clog << "integral1\t" << hh1->Integral() << endl;
  clog << "integral2\t" << hh2->Integral() << endl;

}
  //h1->Scale((1/(h1->Integral()))*100);
  //h2->Scale((1/(h2->Integral()))*100);
  //h3->Scale((1/(h3->Integral()))*100);
  //h4->Scale((1/(h4->Integral()))*100);







  TFile* out = new TFile("/home/marco/root/macros/hist.root", "UPDATE");
  if ( out->IsOpen() ) printf("File opened successfully\n");
  TDirectory* fold;
  if(!(out->GetDirectory(type.c_str()))) fold = out->mkdir(type.c_str());
  else fold = out->GetDirectory(type.c_str());
  fold->cd();
  //TCanvas* c = new TCanvas("c","c");
  //c->Divide(2,2);
  h1->SetXTitle("Wheel");
  h1->SetYTitle("Sector");
  h2->SetXTitle("Wheel");
  h2->SetYTitle("Sector");
  h3->SetXTitle("Wheel");
  h3->SetYTitle("Sector");
  h4->SetXTitle("Wheel");
  h4->SetYTitle("Sector");
  h1->GetXaxis()->SetNdivisions(5);
  h2->GetXaxis()->SetNdivisions(5);
  h3->GetXaxis()->SetNdivisions(5);
  h4->GetXaxis()->SetNdivisions(5);
  h1->GetYaxis()->SetNdivisions(12);
  h2->GetYaxis()->SetNdivisions(12);
  h3->GetYaxis()->SetNdivisions(12);
  h4->GetYaxis()->SetNdivisions(12);
  h1->SetOption("colztext");
  h2->SetOption("colztext");
  h3->SetOption("colztext");
  h4->SetOption("colztext");
  h1->SetMarkerStyle(3);
  h2->SetMarkerStyle(3);
  h3->SetMarkerStyle(3);
  h4->SetMarkerStyle(3);
  h1->SetMarkerSize(1.5);
  h2->SetMarkerSize(1.5);
  h3->SetMarkerSize(1.5);
  h4->SetMarkerSize(1.5);
  //TExec *ex2 = new TExec("ex2","gStyle->SetPaintTextFormat(\"4.0f\") ");

  //h1->SetPaintTextFormat("4.0f");
  //h2->SetPaintTextFormat("4.0f");
  //h3->SetPaintTextFormat("4.0f");
  //h4->SetPaintTextFormat("4.0f");
  c->cd();
  
  h1->Draw();
  //leg->Draw();
  c2->cd();
  h2->Draw();
  //leg->Draw();
  //c->cd(3);
  //h3->Draw();
  //c->cd(4);
  //h4->Draw();
  c->Write();
  c2->Write();
  //c3->cd(1);
  //h3->Draw();
  //c3->cd(3);
  //h4->Draw();
  //h1->Write("Station_1");
  //h2->Write("Station_2");
  //h3->Write("Station_3");
  //h4->Write("Station_4");
  c3->Write();
  c4->Write();
  fold->cd();

  out->Write();
  out->Close();
  

  delete hh1;
  delete hh2;
  delete hh3;
  delete hh4;

  delete h1;
  delete h2;
  delete h3;
  delete h4;


}
