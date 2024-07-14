#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TText.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include <math.h>
#include <iostream>
#include <fstream>
using namespace Pythia8;
using namespace std;

int main() { 
  gStyle->SetOptStat(1111111); //Opção da caixa de estatística

  Pythia pythia;
  // Setando flags do Pythia: energia de CM, partículas incidentes, processos requeridos...
  // Colisão pp numa energia de centro de massa de 7 TeV
  pythia.readString("Beams:eCM = 7000."); // energia do CM
  pythia.readString("Beams:idA = 2212");  // próton incidente no beam A
  pythia.readString("Beams:idB = 2212");  // próton incidente no beam B

  //pythia.readString("Onia:all(3S1)=on"); // esta é a única flag de processo ligada, ativa a produção de todas as ressonâncias de bottomonium
 // pythia.readString("charmonium:all");
  pythia.readString("Charmonium:all=on");
  pythia.readString("443:onMode = off");
  pythia.readString("443:onIfMatch = 13 -13");
  
  // Inicializando o Pythia
  pythia.init();
  Event *event = &pythia.event;

  // Declarando os histrogramas:

  // Variáveis cinemáticas do Upsilon (eta, pT e phi) sem e com corte no pT e eta dos múons (pT>3,5 GeV e |eta|<2,1).
  // Histogramas serão preenchidos com eta e pT do Upsilon em 2 condições:
  // Com corte: ambos múons devem ter pT>3,5 GeV e |eta|<2,1
  // Sem corte: não há corte na variável dos múons
  // Para calcular a aceptância bin a bin, 
  // basta dividir o histogramas com corte pelos 
  // histogramas sem corte (isto é feito ao final do programa)
  TH1D* UpsilonPt = new TH1D("UpsilonPt","Upsilon p_{T}",30,0,30);
  TH1D* UpsilonCutPt = new TH1D("UpsilonCutPt","Upsilon withCuts p_{T}",30,0,30);
  TH1D* UpsilonEta = new TH1D("UpsilonEta","Upsilon Eta",20,-4.,4.);
  TH1D* UpsilonCutEta = new TH1D("UpsilonCutEta","Upsilon withCut Eta",20,-4.,4.);
  TH1D* UpsilonPhi = new TH1D("UpsilonPhi","Upsilon ",64,-3.2,3.2);
  TH1D* UpsilonCutPhi = new TH1D("UpsilonCutPhi","UpsilonCutPt",64,-3.2,3.2);

  // Também vamos considerar a polarização do Upsilon: I ~ 1 +/- cos^2(theta)  
  // PLUS: alpha=1 e MINUS = -1
  TH1D* UpsilonPt_PLUS = new TH1D("UpsilonPt_PLUS","Upsilon p_{T} I ~ 1+cos^{2}#theta",30,0,30);
  TH1D* UpsilonCutPt_PLUS = new TH1D("UpsilonCutPt_PLUS","Upsilon withCut p_{T} I ~ 1+cos^{2}#theta",30,0,30);
  TH1D* UpsilonPt_MINUS = new TH1D("UpsilonPt_MINUS","Upsilon p_{T} I ~ 1-cos^{2}#theta",30,0,30);
  TH1D* UpsilonCutPt_MINUS = new TH1D("UpsilonCutPt_MINUS","Upsilon withCut p_{T} I ~ 1-cos^{2}#theta",30,0,30);

  // Variáveis cinemáticas (pT, eta e phi) dos múons positivos (mup) e negativos (mun)       
  TH1D* munPt = new TH1D("munPt","mu- p_{T}",100,0,30);
  TH1D* mupPt = new TH1D("mupPt","mu+ p_{T}",100,0,30);
  TH1D* munEta = new TH1D("munEta","mu- Eta",100,-4.,4.);
  TH1D* mupEta = new TH1D("mupEta","mu+ Eta",100,-4.,4.);
  TH1D* munPhi = new TH1D("munPhi","mu- Phi",100,-3.2,3.2);
  TH1D* mupPhi = new TH1D("mupPhi","mu+ Phi",100,-3.2,3.2);

  // Histograma 2D para mostrar a região de aceptância do CMS
  TH2D* mu_pt_eta = new TH2D("mu_pt_eta", "p_{T} x #eta",100,0,30,100,-4.,4.);
  gStyle->SetOptStat(0);

  double acept(0.),aceptPLUS(0),aceptMINUS(0); // Cálculo da aceptância para alpha =0 (acept), =1 (aceptPLUS) e =-1 (aceptMINUS).
  int ncut(0),ncutPLUS(0),ncutMINUS(0); // Contador de eventos com corte na variável dos múons
  int ntotal(0); // Contador do número total de eventos
  int nev(100000); // Número total de eventos gerados
  int UpsilonDaughter1(0), UpsilonDaughter2(0); // Variáveis que vão receber as partículas filhas do Upsilon.
  int indexUpsilon(0); // Índice do Upsilon no evento gerado
  int munIndex(0), mupIndex(0); // Índice do par de múons no evento gerado
  double thetastar(0), IPLUS(0), IMINUS(0); // theta* e os pesos correspondentes a cada polarização (alpha=1 (PLUS) e =-1 (MINUS))

  // Começa o loop de eventos
  for (int iEvent = 0; iEvent < nev; ++iEvent) {
    if (!pythia.next()) continue;
    if (iEvent < 1) {pythia.info.list(); pythia.event.list();} // Imprime o primeiro evento

    indexUpsilon = -1; // Índice do Upsilon recebe inicialmente uma flag -1
    // Começa o loop de partículas
    for (int i = 0; i < pythia.event.size(); ++i){

      Particle& theParticle = pythia.event[i];

      // Procura por um Upsilon (id=553)
      if (abs(theParticle.id()) == 443) {
        indexUpsilon = i; // Pega o índice do Upsilon
        break; // Se acha o Upsilon, sai do loop de partículas
      }     
    } // Fim do loop de partículas.
    if (indexUpsilon == -1) continue; // Se não houve identificação de Upsilon no evento, segue para o próximo evento, continua o loop.

    //Encontra as filhas do Upsilon
    UpsilonDaughter1 = pythia.event[indexUpsilon].daughter1();
    UpsilonDaughter2 = pythia.event[indexUpsilon].daughter2();

    // Zerando o índice dos múons
    munIndex=0;
    mupIndex=0;

    // Procurando por um par de múon-antimúon (id=+/-13) entre as filhas do Upsilon
    if (UpsilonDaughter1<UpsilonDaughter2) {
      // Varredura sobre todas as filhas do Upsilon
      for (int i=UpsilonDaughter1; i<=UpsilonDaughter2; ++i) {
        if (pythia.event[i].id()==13)  munIndex=i;
        if (pythia.event[i].id()==-13) mupIndex=i;
      }
    }
    // Checando se encontrou um par muon/antimuon entre as filhas do Upsilon
    if (munIndex!=0 && mupIndex!=0) {
      // Associando a cada partícula (muon+, muon- e Upsilon) um TLorentzVector
      TLorentzVector MuonP(pythia.event[mupIndex].px(),pythia.event[mupIndex].py(),pythia.event[mupIndex].pz(),pythia.event[mupIndex].e());
      TLorentzVector MuonN(pythia.event[munIndex].px(),pythia.event[munIndex].py(),pythia.event[munIndex].pz(),pythia.event[munIndex].e());
      TLorentzVector MuonP_CM(pythia.event[mupIndex].px(),pythia.event[mupIndex].py(),pythia.event[mupIndex].pz(),pythia.event[mupIndex].e()); // A seguir faremos um boost
      TLorentzVector MuonN_CM(pythia.event[munIndex].px(),pythia.event[munIndex].py(),pythia.event[munIndex].pz(),pythia.event[munIndex].e()); // A seguir faremos um boost
      TLorentzVector Upsilon(pythia.event[munIndex].px(),pythia.event[munIndex].py(),pythia.event[munIndex].pz(),pythia.event[munIndex].e());

      // Obter o centro de massa do par de múons.
      TVector3 Upsilon_CM = -(MuonP+MuonN).BoostVector(); 

      // Realizando um boost, ou seja, obtendo as variáveis cinemáticas dos múons no sistema de referência do seu CM (onde o  Upsilon está em repouso).
      MuonP_CM.Boost(Upsilon_CM);
      MuonN_CM.Boost(Upsilon_CM);

      // Imprime na tela informações sobre o evento:
      cout << "Event number " << iEvent << endl;
      cout << "Found an event " << pythia.event[indexUpsilon].name() << " -> " << pythia.event[munIndex].name() << " " << pythia.event[mupIndex].name() << endl;
      cout << "Mu+ 4-mom = " << pythia.event[munIndex].p() << endl;
      cout << "Mu- 4-mom = " << pythia.event[mupIndex].p() << endl;
      
      // Preenchendo os histogramas com as variáveis cinemática dos múons
      mupPt->Fill(pythia.event[mupIndex].pT());
      munPt->Fill(pythia.event[munIndex].pT());   
      mupEta->Fill(pythia.event[mupIndex].phi());
      munEta->Fill(pythia.event[munIndex].phi()); 
      mupPhi->Fill(pythia.event[mupIndex].eta());
      munPhi->Fill(pythia.event[munIndex].eta());  

      // ângulo theta*, que é o ângulo entre o momentum do múon no sistema de repouso do Upsilon e o momentum do Upsilon no sistema de laboratório.
      thetastar = abs(MuonP_CM.Theta() - pythia.event[indexUpsilon].theta());

      // cálculo do peso devido a polarização
      IPLUS = (3./4.)*(1 + pow(cos(thetastar),2)); // alpha = 1
      IMINUS = (3./2.)*(1 - pow(cos(thetastar),2)); // alpha = -1

      // Preenchendo os histogramas com as variáveis cinemáticas do Upsilon
      UpsilonPt->Fill(pythia.event[indexUpsilon].pT());
      UpsilonPt_PLUS->Fill(IPLUS*pythia.event[indexUpsilon].pT());
      UpsilonPt_MINUS->Fill(IMINUS*pythia.event[indexUpsilon].pT());            
      UpsilonEta->Fill(pythia.event[indexUpsilon].eta());
      UpsilonPhi->Fill(pythia.event[indexUpsilon].phi());

      mu_pt_eta->Fill(pythia.event[munIndex].pT(),pythia.event[munIndex].eta());
      mu_pt_eta->Fill(pythia.event[mupIndex].pT(),pythia.event[mupIndex].eta());

      // Verificar se o Upsilon não polarizado está na região de aceptância
      if (pythia.event[munIndex].pT() > 1.0      && pythia.event[mupIndex].pT() > 1.0 &&
        abs(pythia.event[munIndex].eta()) < 2.4 && abs(pythia.event[mupIndex].eta()) < 2.4) {
        ncut++;
        UpsilonCutPt->Fill(pythia.event[indexUpsilon].pT());                             
        UpsilonCutEta->Fill(pythia.event[indexUpsilon].eta());
        UpsilonCutPhi->Fill(pythia.event[indexUpsilon].phi());
      }
      // Verificar se o Upsilon polarizado com alpha=1 está na região de aceptância
      if (IPLUS*pythia.event[munIndex].pT() > 1.0      && IPLUS*pythia.event[mupIndex].pT() > 1.0 &&
        abs(IPLUS*pythia.event[munIndex].eta()) < 2.4 && abs(IPLUS*pythia.event[mupIndex].eta()) < 2.4) {
        ncutPLUS++;
        UpsilonCutPt_PLUS->Fill(IPLUS*pythia.event[indexUpsilon].pT());          
      }      
      // Verificar se o Upsilon polarizado com alpha=-1 está na região de aceptância
      if (IMINUS*pythia.event[munIndex].pT() > 1.0      && IMINUS*pythia.event[mupIndex].pT() > 1.0 &&
        abs(IMINUS*pythia.event[munIndex].eta()) < 2.4 && abs(IMINUS*pythia.event[mupIndex].eta()) < 2.4) {
        ncutMINUS++;
        UpsilonCutPt_MINUS->Fill(IMINUS*pythia.event[indexUpsilon].pT());          
      }       
      ntotal++;
    }
  } // Fim do loop de eventos

  // Cálculo da aceptância para as diferentes polarizações assumidas
  acept = double(ncut)/double(ntotal); // não-polarizado
  aceptPLUS = double(ncutPLUS)/double(ntotal); // alpha=1
  aceptMINUS = double(ncutMINUS)/double(ntotal); // alpha=-1

  ofstream myfile;
  myfile.open("acept.txt");
  myfile << "Aceptancia para diferentes polarizacoes.\n";
  myfile << "Nominal: " << ncut << "/" << ntotal << " = " << acept << "\n";
  myfile << "Transversal: " << ncutPLUS << "/" << ntotal << " = "  << aceptPLUS << "\n";
  myfile << "Longitudinal: " << ncutMINUS << "/" << ntotal << " = " << aceptMINUS << "\n";
  myfile << "Numero Total de J/psi que decairam em mu+mu- " << ntotal << "\n";
  myfile << "Numero Total de eventos em que os dois muons do J/psi estavam dentro da area de Aceptancia " << ncut << "\n";
  myfile.close();

  // Imprimir o resultado na tela
  if (ntotal!=0) cout << "Acceptance = #cut/#total = " << ncut << "/" << ntotal << " = " << acept << endl;
  if (ntotal!=0) cout << "AcceptancePLUS = #cut/#total = " << ncutPLUS << "/" << ntotal << " = " << aceptPLUS << endl;
  if (ntotal!=0) cout << "AcceptanceMINUS = #cut/#total = " << ncutMINUS << "/" << ntotal << " = " << aceptMINUS << endl;

  // Informação sobre a estatística da geração dos eventos
  pythia.stat();
  // Criar o arquivo de output onde os histogramas serão armazenados
  TFile* outFile = new TFile("gen.root","RECREATE");
  UpsilonPhi->Write(); UpsilonCutPhi->Write();
  UpsilonEta->Write(); UpsilonCutEta->Write();
  UpsilonPt->Write();  UpsilonCutPt->Write();
  UpsilonPt_PLUS->Write(); UpsilonCutPt_PLUS->Write();    
  UpsilonPt_MINUS->Write(); UpsilonCutPt_MINUS->Write();     
  munPt->Write();  mupPt->Write();
  munEta->Write(); mupEta->Write();
  munPhi->Write(); mupPhi->Write();

  mu_pt_eta->GetXaxis()->SetRangeUser(0, 15);

  mu_pt_eta->Write();

  // Calculando a aceptância em função do pT para diferentes polarizações
  TH1 *Accept = (TH1*)UpsilonCutPt->Clone("Accept"); 
  TH1 *AcceptPLUS = (TH1*)UpsilonCutPt_PLUS->Clone("AcceptPLUS"); 
  TH1 *AcceptMINUS = (TH1*)UpsilonCutPt_MINUS->Clone("AcceptMINUS"); 

  Accept->Divide(UpsilonCutPt,UpsilonPt);
  Accept->Write();
  AcceptPLUS->Divide(UpsilonCutPt_PLUS,UpsilonPt_PLUS);
  AcceptPLUS->Write();
  AcceptMINUS->Divide(UpsilonCutPt_MINUS,UpsilonPt_MINUS);
  AcceptMINUS->Write();

  // Definir os limites dos eixos X e Y para o histograma AcceptMINUS
  Accept->GetXaxis()->SetRangeUser(0, 10);
  Accept->GetYaxis()->SetRangeUser(0, 1.);

  AcceptPLUS->GetXaxis()->SetRangeUser(0, 10);
  AcceptPLUS->GetYaxis()->SetRangeUser(0, 1.);

  AcceptMINUS->GetXaxis()->SetRangeUser(0, 10);
  AcceptMINUS->GetYaxis()->SetRangeUser(0, 1.);

  // Desenhar e salvar gráfico da acceptância
  TCanvas* c1 = new TCanvas("c1","",800,800);
  AcceptPLUS->SetTitle("Acceptance J/#psi(1S)");
  AcceptPLUS->GetXaxis()->SetTitle("p_{T}^{J/#psi} (GeV)");  
  AcceptPLUS->GetYaxis()->SetTitle("Acceptance");  

  AcceptPLUS->SetMarkerStyle(21);
  Accept->SetMarkerStyle(21);  
  AcceptMINUS->SetMarkerStyle(21);

  AcceptPLUS->SetMarkerColor(kBlue);
  Accept->SetMarkerColor(kBlack);
  AcceptMINUS->SetMarkerColor(kRed);

  auto legend = new TLegend(0.1,0.75,0.3,0.9);
  legend->AddEntry(Accept,"Nominal","p");
  legend->AddEntry(AcceptPLUS,"Transverse","p");
  legend->AddEntry(AcceptMINUS,"Longitudinal","p");

  AcceptPLUS->Draw("P");
  Accept->Draw("P""SAME");
  AcceptMINUS->Draw("P""SAME");  
  legend->Draw("SAME");

  c1->SetGrid();
  c1->SaveAs("CMS_acceptance.pdf","pdf");
  c1->SaveAs("CMS_acceptance.root","root");
  c1->Close();
  gStyle->SetOptStat(0); 

  // Desenhar e salvar gráfico da regiao de acceptância
  TCanvas* c2 = new TCanvas("c2","",800,800);
  mu_pt_eta->GetXaxis()->SetTitle("p_{T} (GeV)");
  mu_pt_eta->GetYaxis()->SetTitle("#eta");
  mu_pt_eta->Draw("colz");
  TLine *line1 = new TLine(0.,-2.4,15.,-2.4);
  TLine *line2 = new TLine(0.,2.4,15.,2.4);
  TLine *line3 = new TLine(1.0,-2.4,1.0,2.4);
  line1->SetLineColor(kRed);
  line2->SetLineColor(kRed);
  line3->SetLineColor(kRed);
  line1->SetLineWidth(4); line1->SetLineStyle(9);
  line2->SetLineWidth(4); line2->SetLineStyle(9);
  line3->SetLineWidth(4); line3->SetLineStyle(9);

  TText *t = new TText(10.,-2.3,"Acceptance region");
  t->SetTextColor(kRed);
  t->SetTextFont(43);
  t->SetTextSize(20);

  line1->Draw();
  line2->Draw();
  line3->Draw();

  t->Draw();  
  c2->SaveAs("CMS_acceptance2.pdf","pdf");
  c2->SaveAs("CMS_acceptance2.root","root");
  c2->Close();
  outFile->Close();
  return 0;
}
