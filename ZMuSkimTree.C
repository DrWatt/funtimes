#define ZMuSkimTree_cxx
#include "ZMuSkimTree_util.cpp"

void ZMuSkimTree::SkimmerNoDT()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //if (maxentries>-1) nentries=maxentries;
  Long64_t nbytes = 0, nb = 0;
  //ofstream listoutput("70kreal.csv");
  //ofstream listoutput2("70kemul.csv");

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
  vector<float> pos;
  vector<float> dir;
  vector<int> simwheel;
  vector<int> simsector;
  vector<int> simstation;
  vector<float> simpos;
  vector<float> simdir;
  vector<int> sample;
  vector<int> asd;
  vector<float> segmpos;
  vector<float> segmdirx;
  vector<float> segmdirz;

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
  vector<float> tpos;
  vector<float> tdir;
  vector<float> tsimpos;
  vector<float> tsimdir;

  ///////////////////Segment stuff
  //TCanvas* ccc = new TCanvas("ccc","ccc"); 
  TH1F* hh = new TH1F("Track x_pos","Track x_pos",60,-10,10);
  TH1F* hhsim = new TH1F("Emulator Track x_pos","Emulator Track x_pos",60,-10,10);
  TH1F* hhh = new TH1F("Track Direction","Track Direction",120,-20,20);
  TH1F* hhhsim = new TH1F("Emulator Track Direction","Emulator Track Direction",120,-20,20);
  vector<int> segmentry;
  vector<int> segwheel;
  vector<int> segsector;
  vector<int> segstation;

  /////////Real and emulated data reading
  clog << "Extracting data from ROOT tree\n";
   
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    float t = jentry/1000.;
    if (t == (int)t) clog <<'\r' << "Entry " << jentry;

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
              if (ltTwinMuxOut_wheel->at(i) == w && ltTwinMuxOut_sector->at(i) == sc && ltTwinMuxOut_station->at(i) == st && ltTwinMuxOut_phi->at(i) != 2044 && abs(ltTwinMuxOut_bx->at(i)) < 3)
              {
                if(!(this->Finder(w,sc,st)))
                {
                  skimbx.push_back(ltTwinMuxOut_bx->at(i));
                  skimphi.push_back(ltTwinMuxOut_phi->at(i));
                  skimphiB.push_back(ltTwinMuxOut_phiB->at(i));
                  skimentry.push_back(jentry);
                  wheel.push_back(w);
                  sector.push_back(sc);
                  station.push_back(st);
                  pos.push_back(ltTwinMuxOut_pos->at(i));
                  dir.push_back(ltTwinMuxOut_dir->at(i));   
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
            if (ltSimTwinMuxOut_wheel->at(i) == w && ltSimTwinMuxOut_sector->at(i) == sc && ltSimTwinMuxOut_station->at(i) == st && abs(ltSimTwinMuxOut_bx->at(i)) < 3)
            {
              if(!(this->Finder(w,sc,st)))
              {
                simskimbx.push_back(ltSimTwinMuxOut_bx->at(i));
                simskimphi.push_back(ltSimTwinMuxOut_phi->at(i));
                simskimphiB.push_back(ltSimTwinMuxOut_phiB->at(i));
                simskimentry.push_back(jentry);
                simwheel.push_back(w);
                simsector.push_back(sc);
                simstation.push_back(st);
                simpos.push_back(ltSimTwinMuxOut_pos->at(i));
                simdir.push_back(ltSimTwinMuxOut_dir->at(i));                    
              }      
            }
            else
              continue;
            }
          }
        }
      }
    } 
    if (dtsegm4D_wheel->size())
    {
      for (unsigned int i = 0; i < dtsegm4D_wheel->size(); ++i)
      {
        for (int w = -2; w < 3; ++w)
        {
          for (int sc = 1; sc < 13; ++sc)
          {
          
            for (int st = 1; st < 5; ++st)
            {
              
              if (dtsegm4D_wheel->at(i) == w && dtsegm4D_sector->at(i) == sc && dtsegm4D_station->at(i) == st && this->Finder(w,sc,st) == 0 && dtsegm4D_phinhits->at(i)>5 )
              {
                segmentry.push_back(jentry);
                segwheel.push_back(w);
                segsector.push_back(sc);
                segstation.push_back(st);
                segmpos.push_back(dtsegm4D_x_pos_loc->at(i));
                segmdirx.push_back(dtsegm4D_x_dir_loc->at(i));
                segmdirz.push_back(dtsegm4D_z_dir_loc->at(i));
              }
            }
          }
        }
      }
    } 
  }
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
      if(*j == *i)
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
      if(*j == *i)
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
    //sample.push_back(skimbx.at(distance(skimentry.begin(),i)));
    for (std::vector<int>::iterator j = simskimentry.begin(); j != simskimentry.end(); ++j)
    {
      if(*j == *i)
      {
        asd.push_back(*j);
        asd.push_back(simwheel.at(distance(simskimentry.begin(),j)));
        asd.push_back(simsector.at(distance(simskimentry.begin(),j)));
        asd.push_back(simstation.at(distance(simskimentry.begin(),j)));
        //asd.push_back(simskimbx.at(distance(simskimentry.begin(),j)));
        if (sample == asd)
        {
          tskimbx.push_back(skimbx.at(distance(skimentry.begin(),i)));
          twheel.push_back(wheel.at(distance(skimentry.begin(),i)));
          tsector.push_back(sector.at(distance(skimentry.begin(),i)));
          tstation.push_back(station.at(distance(skimentry.begin(),i)));
          tskimphi.push_back(skimphi.at(distance(skimentry.begin(),i)));
          tskimphiB.push_back(skimphiB.at(distance(skimentry.begin(),i)));
          tpos.push_back(pos.at(distance(skimentry.begin(),i)));
          tdir.push_back(dir.at(distance(skimentry.begin(),i)));
          tsimskimphi.push_back(simskimphi.at(distance(simskimentry.begin(),j)));
          tsimskimphiB.push_back(simskimphiB.at(distance(simskimentry.begin(),j)));
          tsimskimbx.push_back(simskimbx.at(distance(simskimentry.begin(),j)));
          tsimpos.push_back(simpos.at(distance(simskimentry.begin(),j)));
          tsimdir.push_back(simdir.at(distance(simskimentry.begin(),j)));
          tskimentry.push_back(*i);
  
          asd.clear();
          break;
        }
        else asd.clear();
      }
  
    }
    sample.clear();
  }
////////////////////////////IN MEMORIAM: LAMBDA FUNCTION to see which data is nearer to the reconstructed segments

  for (std::vector<int>::iterator i = tskimentry.begin(); i != tskimentry.end(); ++i)
  {
    
    sample.push_back(*i);
    sample.push_back(twheel.at(distance(tskimentry.begin(),i)));
    sample.push_back(tsector.at(distance(tskimentry.begin(),i)));
    sample.push_back(tstation.at(distance(tskimentry.begin(),i)));
    for (std::vector<int>::iterator j = segmentry.begin(); j != segmentry.end(); ++j)
    {
      if(*j == *i)
      {
        asd.push_back(*j);
        asd.push_back(segwheel.at(distance(segmentry.begin(),j)));
        asd.push_back(segsector.at(distance(segmentry.begin(),j)));
        asd.push_back(segstation.at(distance(segmentry.begin(),j)));
        if (sample == asd)
        {
          hh->Fill(tpos.at(distance(tskimentry.begin(),i))-segmpos.at(distance(segmentry.begin(),j)));
          hhsim->Fill(tsimpos.at(distance(tskimentry.begin(),i))-segmpos.at(distance(segmentry.begin(),j)));
          hhh->Fill(tdir.at(distance(tskimentry.begin(),i))-180*atan(segmdirx.at(distance(segmentry.begin(),j))/segmdirz.at(distance(segmentry.begin(),j)))/3.1415);
          hhhsim->Fill(tsimdir.at(distance(tskimentry.begin(),i))-180*atan(segmdirx.at(distance(segmentry.begin(),j))/segmdirz.at(distance(segmentry.begin(),j)))/3.1415);

          asd.clear();
          break;
        }
        else asd.clear();
      }
  
    }
    sample.clear();
  }

 
//////dir formula:  ltTwinMuxIn_dir-180.0*atan(dtsegm4D_x_dir_loc/dtsegm4D_z_dir_loc)/3.1415>>h(200,-10,10)

////////////////Delta Phis
  float u;
  vector<float> phid;
  vector<float> phibd;
  TH1F* hphid = new TH1F("","",80,-39.5,39.5);
  TH1F* hphibd = new TH1F("","",80,-39.5,39.5);
  for (std::vector<float>::iterator i = tskimphi.begin(); i != tskimphi.end(); ++i)
  {
    u = *i - tsimskimphi.at(distance(tskimphi.begin(),i));
    phid.push_back(u);
    hphid->Fill(u);
  }
  u = 0;
  for (std::vector<float>::iterator i = tskimphiB.begin(); i != tskimphiB.end(); ++i)
  {
    u =*i - tsimskimphiB.at(distance(tskimphiB.begin(),i));
    phibd.push_back(u);
    hphibd->Fill(u);
  }
  gStyle->SetPaintTextFormat("4.3f");
  /////////////Drawing
  hh->SetXTitle("#Deltax (cm)");
  hh->SetYTitle("Counts");
  hhsim->SetXTitle("#Deltax (cm)");
  hhsim->SetYTitle("Counts");
  hhh->SetXTitle("#DeltaDir (deg)");
  hhh->SetYTitle("Counts");
  hhhsim->SetXTitle("#DeltaDir (deg)");
  hhhsim->SetYTitle("Counts");
  hh->Fit("gaus");
  
  hhsim->Fit("gaus");
  hhh->Fit("gaus");//
  hhhsim->Fit("gaus","","",-5,5);

  TFile* outs = new TFile("/home/marco/root/macros/Segm.root", "UPDATE");
    if ( outs->IsOpen() ) printf("Segm.root opened successfully\n");
  //ccc->Divide(2,2);
  //ccc->cd(1);
    gStyle->SetOptStat("e");
    gStyle->SetOptFit(0001);
  hh->Write("Track x_pos");
  //ccc->cd(2);
  hhh->Write("Track Direction");
  //ccc->cd(3);
  hhsim->Write("Emulator Track x_pos");
  //ccc->cd(4);
  hhhsim->Write("Emulator Track Direction");
  //ccc->Draw();

  outs->Write();
  outs->Close();
  delete hh;
  delete hhh;
  delete hhsim;
  delete hhhsim;
  this->Histofiller(station,sector,wheel, "TwinMux Data",1,pstation,psector,pwheel, "TwinMux Overcounts");
  this->Histofiller(simstation,simsector,simwheel,  "Emulator Data",1,qsimstation,qsimsector,qsimwheel,  "Emulator Overcounts");
  this->Histofiller(tstation,tsector,twheel,  "real_equal");
  
  TFile* outg = new TFile("/home/marco/root/macros/graph.root", "UPDATE");
    if ( outg->IsOpen() ) printf("graph.root opened successfully\n");

  hphid->Write("phidiff");
  hphibd->Write("phibdiff");

  TGraph* gg = new TGraph(tskimphi.size(),tskimphi.data(),tsimskimphi.data());
  TGraph* gg2 = new TGraph(tskimphiB.size(),tskimphiB.data(),tsimskimphiB.data());
  TH2F* h1 = new TH2F ("BunchXing","BunchXing",7,-3.5,3.5,7,-3.5,3.5);
  for (std::vector<int>::iterator i = tskimbx.begin(); i != tskimbx.end(); ++i)
  {
    h1->Fill(*i,tsimskimbx.at(i - tskimbx.begin()));
  }
  h1->SetXTitle("TwinMux BX");
  h1->SetYTitle("Emulator BX");

  gg->SetDrawOption("ap");
  gg2->SetDrawOption("ap");
  h1->SetOption("colztext");
  h1->Scale((1/(h1->Integral()))*100);
  gg->Write("Phi");
  gg2->Write("PhiBending");
  h1->Write("BunchXing");
  outg->Write();
  outg->Close();


}





void ZMuSkimTree::Skimmer1DT()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //if (maxentries>-1) nentries=maxentries;

  Long64_t nbytes = 0, nb = 0;
  //ofstream listoutput("Skimmer1DT.csv");
  //ofstream listoutput2("SimSkimmer1DT.csv");
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
  vector<int> simwheel;
  vector<int> simsector;
  vector<int> simstation;
  vector<short> rpcbit;
  vector<short> quality;
  vector<short> simquality;
  vector<int> sample;
  vector<int> asd;

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
  vector<int> buffer = z;

  /////////Real and emulated data reading
  clog << "Extracting data from ROOT tree\n";
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    float t = jentry/1000.;
    if (t == (int)t) clog <<'\r' << "Entry " << jentry;
  
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
           
              if (ltTwinMuxOut_wheel->at(i) == w && ltTwinMuxOut_sector->at(i) == sc && ltTwinMuxOut_station->at(i) == st && ltTwinMuxOut_phi->at(i) != 2044 && abs(ltTwinMuxOut_bx->at(i)) < 3 &&  ltTwinMuxOut_rpcbit->at(i) == 1 )
              {
                if(this->Finder(w,sc,st,1,ltTwinMuxOut_bx->at(i))>= 20 && this->Finder(w,sc,st) == 1 )
                {
                  skimbx.push_back(ltTwinMuxOut_bx->at(i));
                  skimphi.push_back(ltTwinMuxOut_phi->at(i));
                  skimphiB.push_back(ltTwinMuxOut_phiB->at(i));
                  skimentry.push_back(jentry);
                  wheel.push_back(w);
                  sector.push_back(sc);
                  station.push_back(st);  
                  rpcbit.push_back(ltTwinMuxOut_rpcbit->at(i));
                  quality.push_back(ltTwinMuxOut_quality->at(i));              
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
             
              if (ltSimTwinMuxOut_wheel->at(i) == w && ltSimTwinMuxOut_sector->at(i) == sc && ltSimTwinMuxOut_station->at(i) == st && abs(ltSimTwinMuxOut_bx->at(i)) < 3)
              {
                if(this->Finder(w,sc,st,1,ltSimTwinMuxOut_bx->at(i))>= 20 &&this->Finder(w,sc,st) == 1)
                {

                  simskimbx.push_back(ltSimTwinMuxOut_bx->at(i));
                  simskimphi.push_back(ltSimTwinMuxOut_phi->at(i));
                  simskimphiB.push_back(ltSimTwinMuxOut_phiB->at(i));
                  simskimentry.push_back(jentry);
                  simwheel.push_back(w);
                  simsector.push_back(sc);
                  simstation.push_back(st);  
                  simquality.push_back(ltSimTwinMuxOut_quality->at(i));
                  
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
  clog << endl;
  ////////////////////////Overcount algorithm
  clog << "starting overcounts\n";

  ///Real
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
      if(*j == *i)
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

  ///Emul
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
      if(*j == *i)
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

  ///////////////////////Equalcount algorithm
  clog << "starting equalcount\n";
  int pindex = 0;
  short f;
  bool qualityflag = 0;
  float phibuff;
  for (std::vector<int>::iterator i = skimentry.begin(); i != skimentry.end(); ++i)
  {
    f=0;
    qualityflag = 0;
    buffer = z;
    sample.push_back(*i);
    sample.push_back(wheel.at(distance(skimentry.begin(),i)));
    sample.push_back(sector.at(distance(skimentry.begin(),i)));
    sample.push_back(station.at(distance(skimentry.begin(),i)));
    //sample.push_back(skimbx.at(distance(skimentry.begin(),i)));
    phibuff = skimphi.at(distance(skimentry.begin(),i));
    for (std::vector<int>::iterator j = simskimentry.begin(); j != simskimentry.end(); ++j)
    {
      if(*j == *i)
      {

        asd.push_back(*j);
        asd.push_back(simwheel.at(distance(simskimentry.begin(),j)));
        asd.push_back(simsector.at(distance(simskimentry.begin(),j)));
        asd.push_back(simstation.at(distance(simskimentry.begin(),j)));

        if(qualityflag) ++f;
        if(qualityflag && asd == buffer)
        {
          if (abs(simskimphi.at(distance(simskimentry.begin(),j))-phibuff) < abs(simskimphi.at(pindex))-phibuff)
          {
            tsimskimphi.pop_back();
            tsimskimphiB.pop_back();
            tsimskimbx.pop_back();
  
            tsimskimphi.push_back(simskimphi.at(distance(simskimentry.begin(),j)));
            tsimskimphiB.push_back(simskimphiB.at(distance(simskimentry.begin(),j)));
            tsimskimbx.push_back(simskimbx.at(distance(simskimentry.begin(),j)));
            pindex = distance(simskimentry.begin(),j);
          }
          
          
        }
        //asd.push_back(simskimbx.at(distance(simskimentry.begin(),j)));
        if (sample == asd && !qualityflag)
        {
          qualityflag = 1;
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

          pindex = distance(simskimentry.begin(),j);
          buffer = asd;
          asd.clear();
        }
        else asd.clear();
      }
      if (f>10) break;
    }
    sample.clear();
  }
      
  ////////////////Delta Phis
  float u;
  vector<float> phid;
  vector<float> phibd;
  TH1F* hphid = new TH1F("","",20,-42,42);
  TH1F* hphibd = new TH1F("","",20,-42,42);
  for (std::vector<float>::iterator i = tskimphi.begin(); i != tskimphi.end(); ++i)
  {
  
    u = *i - tsimskimphi.at(distance(tskimphi.begin(),i));
    //clog << u << endl;
    phid.push_back(u);
    hphid->Fill(u);
  }
  u = 0;
  //clog << "phibd\n";
  for (std::vector<float>::iterator i = tskimphiB.begin(); i != tskimphiB.end(); ++i)
  {
    //clog << u << endl;
    u = *i - tsimskimphiB.at(distance(tskimphiB.begin(),i));
    phibd.push_back(u);
    hphibd->Fill(u);
  }
  gStyle->SetPaintTextFormat("4.2f");
  /////////////Drawing
  this->Histofiller(station,sector,wheel, "TwinMux Data",1,pstation,psector,pwheel, "TwinMux Overcounts");
  this->Histofiller(simstation,simsector,simwheel,  "Emulator Data",1,qsimstation,qsimsector,qsimwheel,  "Emulator Overcounts");
  this->Histofiller(tstation,tsector,twheel,  "real_equal1DTbit");
  
  
  TFile* outg = new TFile("/home/marco/root/macros/graph1DT.root", "UPDATE");
    if ( outg->IsOpen() ) printf("File opened successfully\n");
  
  hphid->Write("phidiff");
  hphibd->Write("phibdiff");
  
  TGraph* gg = new TGraph(tskimphi.size(),tskimphi.data(),tsimskimphi.data());
  TGraph* gg2 = new TGraph(tskimphiB.size(),tskimphiB.data(),tsimskimphiB.data());

  TH2F* h1 = new TH2F ("BunchXing","BunchXing",7,-3.5,3.5,7,-3.5,3.5);
  for (std::vector<int>::iterator i = tskimbx.begin(); i != tskimbx.end(); ++i)
  {
    h1->Fill(*i,tsimskimbx.at(i - tskimbx.begin()));
  }
  h1->SetXTitle("TwinMux BX");
  h1->SetYTitle("Emulator BX");
  gg->SetDrawOption("ap");
  gg2->SetDrawOption("ap");
  h1->SetOption("colztext");
  h1->Scale((1/(h1->Integral()))*100);
  
  gg->Write("Phi");
  gg2->Write("PhiBending");
  h1->Write("BunchXing");
  outg->Write();
  outg->Close();


  /*TH2F* h2 = new TH2F ("Station 1","Station 1",5,-2.5,2.5,12,0.5,12.5);
  for (std::vector<short>::iterator i = rpcbit.begin(); i != rpcbit.end(); ++i)
  {
    if(*i == 1)
    {
      h2->Fill(wheel.at(i - rpcbit.begin()),sector.at(i - rpcbit.begin()));
    }
  }
  
  h2->SetOption("colztext");
  h2->Draw();
  */

}



void ZMuSkimTree::QualityCheck()
{
  TH1F* bx = new TH1F("BX Distribution","BX Distribution", 10,-4.5,5.5);
  TH1F* bxsim = new TH1F("","", 10,-4.5,5.5);
  TH1F* hq = new TH1F ("TwinMux Quality","TwinMux Quality", 9,-0.5,8.5);
  TH1F* hqb = new TH1F ("","", 9,-0.5,8.5);
  TH1F* shq = new TH1F ("Emulator Quality","Emulator Quality", 9,-0.5,8.5);
  TH1F* shqb = new TH1F ("","", 9,-0.5,8.5);
  TCanvas* c = new TCanvas("c","c");
    TCanvas* c2 = new TCanvas("c2","c2");


  Long64_t nentries = fChain->GetEntriesFast();
  //if (maxentries>-1) nentries=maxentries;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    float t = jentry/1000.;
    if (t == (int)t) 
      {

        clog << '\r' << "Entry " << jentry;
        }
    if (ltTwinMuxOut_phi->size())
    {
      for (unsigned int i = 0; i < ltTwinMuxOut_phi->size(); ++i)
      {
        if (ltTwinMuxOut_phi->at(i) != 2044 && abs(ltTwinMuxOut_bx->at(i)) < 3)
        {
          bx->Fill(ltTwinMuxOut_bx->at(i));
          //clog << jentry << '\t';
          hq->Fill(ltTwinMuxOut_quality->at(i));
          if(this->Finder(ltTwinMuxOut_wheel->at(i),ltTwinMuxOut_sector->at(i),ltTwinMuxOut_station->at(i),1,ltTwinMuxOut_bx->at(i))>= 20 && ltTwinMuxOut_rpcbit->at(i)==1 && this->Finder(ltTwinMuxOut_wheel->at(i),ltTwinMuxOut_sector->at(i),ltTwinMuxOut_station->at(i)) == 1) hqb->Fill(ltTwinMuxOut_quality->at(i));    
        }
        else
          continue;
      }
    }

    if (ltSimTwinMuxOut_phi->size())
    {
      for (unsigned int i = 0; i < ltSimTwinMuxOut_phi->size(); ++i)
      {
        if (ltSimTwinMuxOut_phi->at(i) != 2044 && abs(ltSimTwinMuxOut_bx->at(i)) < 3)
        {
          bxsim->Fill(ltSimTwinMuxOut_bx->at(i));
          //clog << jentry << '\t';
          shq->Fill(ltSimTwinMuxOut_quality->at(i));
          if(this->Finder(ltSimTwinMuxOut_wheel->at(i),ltSimTwinMuxOut_sector->at(i),ltSimTwinMuxOut_station->at(i),1,ltSimTwinMuxOut_bx->at(i))>= 20 && this->Finder(ltSimTwinMuxOut_wheel->at(i),ltSimTwinMuxOut_sector->at(i),ltSimTwinMuxOut_station->at(i)) == 1) shqb->Fill(ltSimTwinMuxOut_quality->at(i));         
        }
        else
          continue;
      }
    }
  }
  clog << endl;
  TH1F* bxdif = new TH1F(*bxsim);
  bxdif->SetTitle("bxdif");
  bxdif->Add(bx, -1);
  hq->GetYaxis()->SetRangeUser(1,1000000);
  shq->GetYaxis()->SetRangeUser(1,1000000);
  hqb->SetLineColor(kRed);
  shqb->SetLineColor(kRed);
  bxsim->SetLineColor(kRed);
  
  c->Divide(2,1);

  c->cd(1);
  hq->Draw();
  hqb->Draw("SAME");
  c->cd(2);
  shq->Draw();
  shqb->Draw("SAME");
  c->Draw();
  c2->Divide(2,1);
  c2->cd(1);
  bx->Draw();
  bxsim->Draw("SAME");
  c2->cd(2);
  bxdif->Draw();

  c2->Draw();
  
}









////////////////////////////////*****Bisogna controllare l'analisi delle tracce 1DT
////////////////////////////////Per il confronto tra segmenti ricostruiti e dati twinmux per vedere se hanno ragione i dati o emu devo vedere se è un problema di numero di tracce per la posizione perchè l'istogramma fa proprio cagare mentre per la direzione è discreto ma comunque necessita di più entries
////////////////////////////////a proposito del modo per controllare velocemente se le tracce del twinmux sono associate a un muone veramente rilevato dovrei utilizzare i dati Mu... che sono tclonesarray e dovrei chiedere come sono organizzati perché ogni entry ce n'è almeno uno (di numero variabile) di array "vuoti" con size 1 e dato 999.
////////////////////////////////*****Ho già normalizzato le mappe e bisogna stamparle per bene.
////////////////////////////////NODT: Per la distribuzione del deltaphi opterei per una scala logaritmica poichè lo 0 (che volendo si potrebbe togliere ma direi di no) e il 4 sono belli pieni mentre altri valori hanno occorrenze minime
////////////////////////////////Mentre per i deltaphibend c'è una distribuzione più regolare e la scala logaritmica rovinerebbe tutto
////////////////////////////////
////////////////////////////////****AGGIUNTA normalizzazione negli overcounts rispetto ai dati totali reali o emulati
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////TAglia BX nel plot
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////
////////////////////////////////1DT: bitrpc quanti dati hanno il bit 1.
////////////////////////////////Quanti dati emulati di questi hanno bx cabiato
////////////////////////////////
////////////////////////////////
////////////////////////////////Bit 1 dati c'è stato match
////////////////////////////////bit 1 emul solo se c'è stato cambio bx
////////////////////////////////
////////////////////////////////Problema molteplicità mappe*****Dall'equalcount non dovrebbero spuntare valori doppi su cui loopa più volte, se c'è una traccia nell'emulatore due volte, viene presa quella
////////////////////////////////
////////////////////////////////
////////////////////////////////controllo anche quality
////////////////////////////////qualità tracce twinmuxout tutte sovrapposte a qulla di dat con cambio bx
////////////////////////////////Lo stesso con emulati
////////////////////////////////
////////////////////////////////Panoramica bx con tutti i dati tra input e output per vedere chi migliora il picco in 0
////////////////////////////////
////////////////////////////////Migliorare il plot di quality introducendo le richieste sul bit e su un solo DT
////////////////////////////////
////////////////////////////////








//////////////////////MERDACCE
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

}
void ZMuSkimTree::Segm()
{
 if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //if (maxentries>-1) nentries=maxentries;
  //vector<float> phi[2];
   Long64_t nbytes = 0, nb = 0;
   TH1F* h = new TH1F("j","j",400,-200,200);
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
    float t = jentry/1000.;
    if (t == (int)t) clog << "Entry " << jentry << endl;

     
      if (dtsegm4D_wheel->size())
      {
        for (unsigned int i = 0; i < dtsegm4D_wheel->size(); ++i)
        {
        for (int w = -2; w < 3; ++w)
        {
          for (int sc = 1; sc < 13; ++sc)
          {
          
            for (int st = 1; st < 5; ++st)
            {
              
              if (dtsegm4D_wheel->at(i) == w && dtsegm4D_sector->at(i) == sc && dtsegm4D_station->at(i) == st && this->Finder(w,sc,st) == 0 && dtsegm4D_phinhits->at(i)>5 )
              h->Fill(dtsegm4D_x_pos_loc->at(i));



            }
          }
        }
      }
      }
    }

h->Draw();

}*/