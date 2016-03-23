
/*

Leitor do arquivo .txt criado pelo gerador principal, onde retornará o os histogramas de momento transversal e rapidez
de proton, pion, e kaon junto a uma curva da função da TACNE.

*/
void Leito_dados(){

//   massas = 9.46, 5.275, 5.275, 5.271, 5.271, 3.097, 2.981, 2.281, 2.010, 2.010, 2.010, 2.010, 1.971, 1.971, 1.869, 1.869, 1.865, 1.865, 1.672, 1.533, 1.533, 1.385, 1.385, 1.385, 1.3213, 1.3149, 1.232, 1.232, 1.232, 1.232, 1.1973, 1.1925, 1.1894, 1.1156, 1.02, 0.9576, 0.939573, 0.938280, 0.892, 0.892, 0.892, 0.892, 0.783, 0.77,0.77, 0.77, 0.5488, 0.49772, 0.49772, 0.49367, 0.49367, 0.139569, 0.139569, 0.134964 
  
  
        
        const int num = 100000; 
        double xmaxpt = 3;
        double xminpt = 0;
        double xmin = -10;
	double xmax = 10;
	double massa[num];
        double carga[num]; 
        double m_t[num];      
        double y[num];
        int div = 12;
	double yr = xmin;
	double lim = xmax-xmin;
        double pion0 = 0.134964;
        double pionpn = 0.139569;
        double pion = (pion0+pionpn)/2;
        double kaon0  = 0.49772;
        double kaonpn = 0.49367;
        double kaon   = (kaon0+kaonpn)/2;
        double proton = 0.938280;
	
	
	
     

        TString NomeDoArquivo = "Particulas.txt";// Implantar caminho dos dados // 
        ifstream ifs(NomeDoArquivo.Data());

	
	
	TH1F *hpion = new TH1F("pion (pt 7TeV)","TRANSVERSE MOMENTUM",100,xminpt,xmaxpt);
	TH1F *hpiony = new TH1F("pion (y 7TeV)","RAPIDITY",100,xmin,xmax);
	TH1F *hkaon = new TH1F("kaon (pt 7TeV)","TRANSVERSE MOMENTUM",100,xminpt,xmaxpt);
	TH1F *hkaony = new TH1F("kaon (y 7TeV)","RAPIDITY",100,xmin,xmax);
	TH1F *hproton = new TH1F("proton (pt 7TeV)","TRANSVERSE MOMENTUM",100,xminpt,xmaxpt);
	TH1F *hprotony = new TH1F("proton (y 7TeV)","RAPIDITY",100,xmin,xmax);
	
	TGraph *gr = new TGraph(div+1); 
	
	
	
        for (int l=0; l < num; l++){
 		ifs >>massa[l] >> carga[l] >>  m_t[l] >> y[l];
    
	}
	
	for(int t=0; t < num; t++){
	  
	  if(massa[t] == pion0 ){
	    
	    hpion->Fill(m_t[t]);
	    hpiony-> Fill(y[t]);
	    
	  }
	  
	  if(massa[t] == pionpn ){
	    
	    hpion->Fill(m_t[t]);
	    hpiony-> Fill(y[t]);
	    
	  }
	 
	  if(massa[t] == kaon0 ){
	    
	    hkaon->Fill(m_t[t]);
	    hkaony-> Fill(y[t]);
	    
	  }
	  
	  if(massa[t] == kaonpn ){
	    
	    hkaon->Fill(m_t[t]);
	    hkaony-> Fill(y[t]);
	    
	  }
	  
	  if(massa[t] == proton ){
	    
	    hproton->Fill(m_t[t]);
	    hprotony-> Fill(y[t]);
	    
	  }
	  
	  
}

 
 
 
 
 for(int t=0; t <= div ; t++){
  
                  TF1 *fr = new TF1("Function", "((1/(sqrt(2*pi)*1.8))*(pow((1-(1.146-1)*(-(pow((x-2.3),2)/(2*pow(1.8,2))))),(-(1/(1.146-1)))))+(1/(sqrt(2*pi)*1.8))*(pow((1-(1.146-1)*(-(pow((x+2.3),2)/(2*pow(1.8,2))))),(-1/(1.146-1)))))*((pow(0.068,3)*pow((1/cosh([1]-x)),2))*(pow(((1.146-1)*cosh([1]-x))/0.068,(2*1.146-3)/(1.146-1)))*(([2]*(1.146-1))+(0.068*(1/cosh([1]-x))))*(((pow(([2]+((0.068*(1/cosh([1]-x)))/(1.146-1))),(-1.146/(1.146-1))))*((-pow([2],2)*(1.146-2)+(2*[2]*0.068*(1/cosh([1]-x))+(2*pow(0.068,2)*pow((1/cosh([1]-x)),2))))))/((4*pow(pi,2))*(1.146-2)*(pow((1.146-1),3))*(2*1.146-3))))", xmin, xmax);  //Função para a integração.
                  fr->SetParameter(2,pion);
                  fr->SetParameter(1,yr);

                  const int np = 1000;
                  double *xf=new double[np];
                  double *wf=new double[np];

                  fr->CalcGaussLegendreSamplingPoints(np,xf,wf,1e-15);

                  gr->GetX()[t] =  yr;
                  gr->GetY()[t] = fr->IntegralFast(np,xf,wf,xmin,xmax);
                  // Os valores são guardados no gráfico gr. 
                  
                  
                  yr = yr + lim/div;
                  // y representa o valor do range da integral.
                  }
    
            yr = xmin;

    int gMaxpiony = hpiony->GetMaximum();
    TSpline3 *sp = new TSpline3("sp",gr);
    sp->SaveAs("spf.C");
    gROOT->LoadMacro("spf.C++");
    TF1 *sFunction = new TF1("Function","([1]*121000)*spf(x)",xmin,xmax);
    sFunction->SetParameter(1,gMaxpiony);
    
    TCanvas *hx = new TCanvas("hx","Graph",200,10,700,500);
    hpiony->SetXTitle("y");
    hpiony->Draw();
    sFunction->Draw("same");
    hx->Print("histograma_y_particula_PION.pdf");
    delete hx;
     
    
    
    
    for(int ti=0; ti <= div ; ti++){
  
                  TF1 *frkaon = new TF1("Function", "((1/(sqrt(2*pi)*1.8))*(pow((1-(1.146-1)*(-(pow((x-2.3),2)/(2*pow(1.8,2))))),(-(1/(1.146-1)))))+(1/(sqrt(2*pi)*1.8))*(pow((1-(1.146-1)*(-(pow((x+2.3),2)/(2*pow(1.8,2))))),(-1/(1.146-1)))))*((pow(0.068,3)*pow((1/cosh([1]-x)),2))*(pow(((1.146-1)*cosh([1]-x))/0.068,(2*1.146-3)/(1.146-1)))*(([2]*(1.146-1))+(0.068*(1/cosh([1]-x))))*(((pow(([2]+((0.068*(1/cosh([1]-x)))/(1.146-1))),(-1.146/(1.146-1))))*((-pow([2],2)*(1.146-2)+(2*[2]*0.068*(1/cosh([1]-x))+(2*pow(0.068,2)*pow((1/cosh([1]-x)),2))))))/((4*pow(pi,2))*(1.146-2)*(pow((1.146-1),3))*(2*1.146-3))))", xmin, xmax);  //Função para a integração.
                  frkaon->SetParameter(2,kaon);
                  frkaon->SetParameter(1,yr);

                  const int npkaon = 1000;
                  double *xfkaon=new double[npkaon];
                  double *wfkaon=new double[npkaon];

                  frkaon->CalcGaussLegendreSamplingPoints(npkaon,xfkaon,wfkaon,1e-15);

                  gr->GetX()[ti] =  yr;
                  gr->GetY()[ti] = frkaon->IntegralFast(npkaon,xfkaon,wfkaon,xmin,xmax);
                  // Os valores são guardados no gráfico gr. 
                  
                  
                  yr = yr + lim/div;
                  // y representa o valor do range da integral.
                  }
    
            yr = xmin;

    int gMaxkaony = hkaony->GetMaximum();
    TSpline3 *spkaon = new TSpline3("spkaon",gr);
    sp->SaveAs("spfkaon.C");
    gROOT->LoadMacro("spfkaon.C++");
    TF1 *sFunctionkaon = new TF1("Function","([1]*120000)*spfkaon(x)",xmin,xmax);
    sFunctionkaon->SetParameter(1,gMaxkaony);
    
    TCanvas *hx = new TCanvas("hx","Graph",200,10,700,500);
    hkaony->SetXTitle("y");
    hkaony->Draw();
    sFunctionkaon->Draw("same");
    hx->Print("histograma_y_particulaKAON.pdf");
    delete hx;
    
  
    
    
    
    
    for(int tii=0; tii <= div ; tii++){
  
                  TF1 *frproton = new TF1("Function", "((1/(sqrt(2*pi)*1.8))*(pow((1-(1.146-1)*(-(pow((x-2.3),2)/(2*pow(1.8,2))))),(-(1/(1.146-1)))))+(1/(sqrt(2*pi)*1.8))*(pow((1-(1.146-1)*(-(pow((x+2.3),2)/(2*pow(1.8,2))))),(-1/(1.146-1)))))*((pow(0.068,3)*pow((1/cosh([1]-x)),2))*(pow(((1.146-1)*cosh([1]-x))/0.068,(2*1.146-3)/(1.146-1)))*(([2]*(1.146-1))+(0.068*(1/cosh([1]-x))))*(((pow(([2]+((0.068*(1/cosh([1]-x)))/(1.146-1))),(-1.146/(1.146-1))))*((-pow([2],2)*(1.146-2)+(2*[2]*0.068*(1/cosh([1]-x))+(2*pow(0.068,2)*pow((1/cosh([1]-x)),2))))))/((4*pow(pi,2))*(1.146-2)*(pow((1.146-1),3))*(2*1.146-3))))", xmin, xmax);  //Função para a integração.
                  frproton->SetParameter(2,proton);
                  frproton->SetParameter(1,yr);

                  const int npproton = 1000;
                  double *xfproton=new double[npproton];
                  double *wfproton=new double[npproton];

                  frproton->CalcGaussLegendreSamplingPoints(npproton,xfproton,wfproton,1e-15);

                  gr->GetX()[ti] =  yr;
                  gr->GetY()[ti] = frproton->IntegralFast(npkaon,xfproton,wfproton,xmin,xmax);
                  // Os valores são guardados no gráfico gr. 
                  
                  
                  yr = yr + lim/div;
                  // y representa o valor do range da integral.
                  }
    
            yr = xmin;

    int gMaxprotony = hprotony->GetMaximum();
    TSpline3 *spproton = new TSpline3("spproton",gr);
    sp->SaveAs("spfproton.C");
    gROOT->LoadMacro("spfproton.C++");
    TF1 *sFunctionproton = new TF1("Function","([1]*118000)*spfproton(x)",xmin,xmax);
    sFunctionproton->SetParameter(1,gMaxprotony);
    
    TCanvas *hx = new TCanvas("hx","Graph",200,10,700,500);
    hprotony->SetXTitle("y");
    hprotony->Draw();
    sFunctionproton->Draw("same");
    hx->Print("histograma_y_particulaproton.pdf");
    delete hx;
   
    
    
    
    
    
    
    
    int gMaxpion = hpion->GetMaximum();
    TF1 *fGraphPion = new TF1("Function", "([1]/2.15)*(((x*sqrt(pow(x,2)+pow([2],2)))/0.068)*1*((2-1.146)*(3-2*1.146))/((2-1.146)*pow([2],2)+2*[2]*0.068+2*pow(0.068,2)))*(pow((1+(1.146-1)*([2]/0.068)),(1/(1.146-1))))*(pow((1+(1.146-1)*((sqrt(pow(x,2)+pow([2],2)))/0.068)),(-1.146/(1.146-1))))", xminpt, xmaxpt);  //Função para a integração
    fGraphPion->SetParameter(1,gMaxpion);
    fGraphPion->SetParameter(2,pion);
 
    TCanvas *hxx = new TCanvas("hxx","Graph",200,10,700,500);
    hpion->SetXTitle("pt (GeV/c)");
    hpion->Draw(); 
    fGraphPion->Draw("same");
    hxx->Print("Momento_Trasversal_pion_fct.pdf");
    delete hxx;

    int gMaxkaon = hkaon->GetMaximum();
    TF1 *fGraphkaon = new TF1("Function", "([1]/1.35)*(((x*sqrt(pow(x,2)+pow([2],2)))/0.068)*1*((2-1.146)*(3-2*1.146))/((2-1.146)*pow([2],2)+2*[2]*0.068+2*pow(0.068,2)))*(pow((1+(1.146-1)*([2]/0.068)),(1/(1.146-1))))*(pow((1+(1.146-1)*((sqrt(pow(x,2)+pow([2],2)))/0.068)),(-1.146/(1.146-1))))", xminpt, xmaxpt);  //Função para a integração
    fGraphkaon->SetParameter(1,gMaxkaon);
    fGraphkaon->SetParameter(2,kaon);
 
    TCanvas *hxx = new TCanvas("hxx","Graph",200,10,700,500);
    hkaon->SetXTitle("pt (GeV/c)");
    hkaon->Draw(); 
    fGraphkaon->Draw("same");
    hxx->Print("Momento_Trasversal_kaon_fct.pdf");
    delete hxx;
  
    int gMaxproton = hproton->GetMaximum();
    TF1 *fGraphproton = new TF1("Function", "([1]*1.1)*(((x*sqrt(pow(x,2)+pow([2],2)))/0.068)*1*((2-1.146)*(3-2*1.146))/((2-1.146)*pow([2],2)+2*[2]*0.068+2*pow(0.068,2)))*(pow((1+(1.146-1)*([2]/0.068)),(1/(1.146-1))))*(pow((1+(1.146-1)*((sqrt(pow(x,2)+pow([2],2)))/0.068)),(-1.146/(1.146-1))))", xminpt, xmaxpt);  //Função para a integração
    fGraphproton->SetParameter(1,gMaxproton);
    fGraphproton->SetParameter(2,proton);
 
    TCanvas *hxx = new TCanvas("hxx","Graph",200,10,700,500);
    hproton->SetXTitle("pt (GeV/c)");
    hproton->Draw(); 
    fGraphproton->Draw("same");
    hxx->Print("Momento_Trasversal_proton_fct.pdf");
    delete hxx;
 
   }
