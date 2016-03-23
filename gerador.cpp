/* 

Gerador principal em que são gerados todos os arquivos para criação do .txt
Os arquivos gerados são massa, carga, momento transversal, rapidez.

*/


void gerador(){
  
  // Entradas do usuário.
  
  int eventos = 1;
  double energia = 2300;
  
  // Eventos: é a quantidade de eventos em que o código vai rodar. 
  // Energia: é a energia do centro de massa da particula.
  
/*****************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************/
  
 //Variáveis
  

  double xminpt = 0;                          // Valor do limite mínimo para pt
  double xmaxpt = 3;                          // Valor do limite máximo para pt
  double xmin = -10;                          // X min para a integral de rapidez  
  double xmax = 10;                           // X max para a integral de rapidez
  const Double_t div =12;                     // Número de pontos para plotar no gráfico da rapidez
  double lim = xmax - xmin;                   // valor do limite de x
  double y = xmin;                            // valores para o eixo x 
  const int pint = 100;                       // pontos para a interpolação
  double xi = xmin;                           // valor para o eixo x da interpolação 
  const int size = 54;                        // Número de particulas  
  double n[size]={0.000000000178527, 0.000000002432434666636607, 0.000000002432434666636607, 0.000000002440493057336401, 0.000000002440493057336401, 0.0000000230402, 0.0000000268926, 0.000000077143, 0.0000001243885392283889, 0.0000001243885392283889, 0.0000001243885392283889, 0.0000001243885392283889, 0.000000133774146395552, 0.000000133774146395552, 0.0000001626575153468772, 0.0000001626575153468772, 0.0000001639355950507426, 0.0000001639355950507426, 0.000000242905, 0.0000003291610694918138, 0.0000003291610694918138, 0.0000004649454399855365, 0.0000004649454399855365, 0.0000004649454399855365, 0.0000005436225477345856, 0.0000005523795903395713, 0.0000006826355989218414, 0.0000006826355989218414, 0.0000006826355989218414, 0.0000006826355989218414, 0.0000007479215986494419, 0.0000007575305024323641, 0.0000007638148820053583, 0.000000933644, 0.00000122603, 0.0000014766141547977551, 0.000001560097463554025, 0.000001566298958394324, 0.0000018089421064631, 0.0000018089421064631, 0.0000018089421064631, 0.0000018089421064631, 0.00000258166, 0.000002697942688005066, 0.000002697942688005066, 0.000002697942688005066, 0.000006059602317495987, 0.000007431177664583774, 0.000007431177664583774, 0.000007554558250517257, 0.000007554558250517257, 0.0000356944, 0.0000356944, 0.0000363897};     // Resultado de N para cada massa, vindo de uma equação da TACNE. Resultados obtidos pelo programa Mathematica.
  double massa[size]={9.46, 5.275, 5.275, 5.271, 5.271, 3.097, 2.981, 2.281, 2.010, 2.010, 2.010, 2.010, 1.971, 1.971, 1.869, 1.869, 1.865, 1.865, 1.672, 1.533, 1.533, 1.385, 1.385, 1.385, 1.3213, 1.3149, 1.232, 1.232, 1.232, 1.232, 1.1973, 1.1925, 1.1894, 1.1156, 1.02, 0.9576, 0.939573, 0.938280, 0.892, 0.892, 0.892, 0.892, 0.783, 0.77,0.77, 0.77, 0.5488, 0.49772, 0.49772, 0.49367, 0.49367, 0.139569, 0.139569, 0.134964};        // Valor da massa de cada particulda.
  int carga[size]={0,0,0,1,-1,0,0,1,1,-1,0,0,1,-1,1,-1,0,0,-1,-1,0,1,-1,0,-1,0,1,-1,0,2,-1,0,1,0,0,0,0,1,1,-1,0,0,0,1,-1,0,0,0,0,1,-1,1,-1,0};   // Valor da carga das particulas. Obs: particulas e cargas na ordem.
  double m = 0;                               // Variável utilizada para guardar os valores de N.
  double p[size];                             // Vetor para guardar os valores de N normalizada.
  double r = 0;                               // Variável utilizada para guardar os valores de N na normalização.
  double x[size *10000];                      // Vetor utilizado para guardar os valores do sorteio da particula.
  double resultado_integral[size];            // Vetor utilizado para guardar os valores das integrais de rapidez.
  double mx[pint+1][size+1];                  // Matriz que guarda os valores de x do gráfico após a interpolação. 
  double my[pint+1][size+1];                  // Matriz que guarda os valores de y do gráfico após a interpolação.
  const double divpt = 200;                   // Número de pontos para o gráfico pt
  double mxpt[divpt+1][size+1];               // Matriz que guarda os valores de x do gráfico de pt
  double mypt[divpt+1][size+1];               // Matriz que guarda os valores de y do gráfico de pt
  double xipt = xminpt;                       // Valor para o eixo x da interpolação de pt
  double limpt = xmaxpt-xminpt;               // valor do limite de x pt


  
  /***********************************************************************************************************************************************************************************************************************************************************************************************************************************************************/

  
  TGraph *gr = new TGraph(div+1);                                                    // Classe para guardar os pontos do gráfico de rapidez.
  TGraph *gInt = new TGraph(pint+1);                                                 // Classe para guardar os pontos do gráfico interpolado, para depois integrar. 
  TGraph *gIntegral = new TGraph(pint+1);                                            // Classe para guardar os pontos do gráfico para fazer a normalização.     
  TGraph *gIntegralpt = new TGraph(divpt+1);                                         // Classe para guardar os pontos do gráfico para fazer a normalização.
  gRandom = new TRandom3(0);
  gRandom->SetSeed(0);                                                   // Semente usada para gerar números aleatórios
  
  ofstream myfile;                                                                   // Objeto para a criação de um arquivo externo.
  TString EnderecoSalva = "Particulas.txt";
  myfile.open(EnderecoSalva);

 /***********************************************************************************************************************************************************************************************************************************************************************************************************************************************************/
  
  // Rapidez   
 
 
 // Integral feita para todas as massas, guardando os valores em vetores.  
  
 // São atríbuidos valores a y e feito a integral para cada valor, formando assim um gráfico para rapidez   
  
  for(int e=0; e<size; e++){                 // Inicio do "for" com a variável (e) "for(e)" modifica o valor da massa[] nas funções. 
  
          for(int t=0; t <= div ; t++){      // Início do "for" com a variável (t) "for(t)" modifica o valor de y para rapidez.
  
                  TF1 *fr = new TF1("Function", "((1/(sqrt(2*pi)*1.8))*(pow((1-(1.146-1)*(-(pow((x-2.3),2)/(2*pow(1.8,2))))),(-(1/(1.146-1)))))+(1/(sqrt(2*pi)*1.8))*(pow((1-(1.146-1)*(-(pow((x+2.3),2)/(2*pow(1.8,2))))),(-1/(1.146-1)))))*((pow(0.068,3)*pow((1/cosh([1]-x)),2))*(pow(((1.146-1)*cosh([1]-x))/0.068,(2*1.146-3)/(1.146-1)))*(([2]*(1.146-1))+(0.068*(1/cosh([1]-x))))*(((pow(([2]+((0.068*(1/cosh([1]-x)))/(1.146-1))),(-1.146/(1.146-1))))*((-pow([2],2)*(1.146-2)+(2*[2]*0.068*(1/cosh([1]-x))+(2*pow(0.068,2)*pow((1/cosh([1]-x)),2))))))/((4*pow(pi,2))*(1.146-2)*(pow((1.146-1),3))*(2*1.146-3))))", xmin, xmax);  //Função para a integração.
                  fr->SetParameter(2,massa[e]);    // valor da massa atribuido a cada vez que "for(e)" roda.
                  fr->SetParameter(1,y);           // valor de y atribuido a cada vez que o "for(t)" roda

                  const int np = 1000;
                  double *xf=new double[np];
                  double *wf=new double[np];

                  fr->CalcGaussLegendreSamplingPoints(np,xf,wf,1e-15);                          // Classe chamada para fazer a integral de da rapidez

                  gr->GetX()[t] =  y;                                                           // y = ao valor da abscissa do gráfico
                  gr->GetY()[t] = fr->IntegralFast(np,xf,wf,xmin,xmax);                         // integral da rapidez = a ordenada do gráfico
                  // Os valores são guardados no gráfico gr. 
                  
                  
                  y = y + lim/div;
                  // y representa o valor do range da integral.
                  }
    
            y = xmin;



            // faz a interpolação com os valores do gráfico "gr" e guarda nas matrizes mx e my.   
            for(int v=0 ; v <= pint ; v++){                         // Inicio do "for" com a variáve (pint) "for(pint)" para mudar os valores de xi
    
                    mx[v][e] = xi;                                  // valores atríbuidos a abscissa do gráfico na interpolação 
                    my[v][e] = gr->Eval(xi,0,"S");                  // dado um valor da abcissa retorna o ponto da ordenada do gráfico usando a interpolação TSpline3
                    gInt->GetX()[v] = mx[v][e];                     // guardando os valores no gráfico gInt para fazer a integral da interpolação 
                    gInt->GetY()[v] = my[v][e];                     // guardando os valores no gráfico gInt para fazer a integral da interpolação 
                                        
                    xi = xi + lim/pint;                             // Atribuindo os valores a abscissa do gráfico

                    }                                               // fim do for(pint)

             xi = xmin;
             resultado_integral[e]=gInt->Integral();                // faz a integral de gInt

     
     
     

/***********************************************************************************************************************************************************************************************************************************************************************************************************************************************************/
 
   // Momento transversal

   
   // Com a função já normalizada, basta atribuir valores na abscissa para se obter valores da ordenada com a função Eval() onde os pontos estão sendo guardados nas matrize mxpt e mypt.
   
      for(int ptP=0; ptP <= divpt ;ptP++ ){           // Inicio do "for" com variável (ptP) "for(ptP)" muda o valor de xipt para obter o valor da ordenada 
              
              TF1 *f = new TF1("Function", "(((x*sqrt(pow(x,2)+pow([2],2)))/0.068)*1*((2-1.146)*(3-2*1.146))/((2-1.146)*pow([2],2)+2*[2]*0.068+2*pow(0.068,2)))*(pow((1+(1.146-1)*([2]/0.068)),(1/(1.146-1))))*(pow((1+(1.146-1)*((sqrt(pow(x,2)+pow([2],2)))/0.068)),(-1.146/(1.146-1))))", xminpt, xmaxpt);  //Função para a integração
              f->SetParameter(2,massa[e]);            // atribui o valor de massa[] do que "for(e)" determina. 
                      
              mxpt[ptP][e] = xipt;                    // atribui valor para a abscissa 
              mypt[ptP][e] = f->Eval(xipt);           // obtem o valor da ordenada.
       
    
       xipt = xipt + limpt/divpt;                     // Atribuindo os valores a abscissa do gráfico
    }                                                 // Fim de "for(ptP)"        
    
    xipt = xminpt;
    
    
    
    
  }                       // fim de "for(e)"
  
   
/***********************************************************************************************************************************************************************************************************************************************************************************************************************************************************/
   
 // somatória de todos os valores de n obtendo o valor de m, e depois dividindo todos os valores de n por m
 
  for(int i=0 ; i < size;i++){
  
          m = m + n[i];            //somando todos os valores de N[].

          }
    
  for(int j=0;j<size;j++){
 
          r = r+(n[j]/m);          // dividindo cada valor de n por m
          p[j] = r;                // Vetor onde esta sendo guardado os valores de n normalizado
      
          }
  m = 0;
  r = 0;

   
  for(int g=0;g<eventos;g++){                    // Inicio do "for" com a variavel (g) "for(g)" determina o número de eventos 
 
/*********************************************************************************************************************************************************************************************************************************************************************************************************************************************************/  

  // Usando a multiplicidade para calcular o numero de partículas dado a energia 
  
  double num_part_log , num_part;

  num_part_log =(7.8*log(energia))-19;   // Equação logarítimica              

  num_part = (0.004 * energia)+31;       // Equação linear

  cout<< num_part_log <<endl;
  
  // Obs: para os calculos pode se utilizar de apenas uma equação que será determinada na variável n_p. 
  
  const int n_p = num_part_log;          // Número de particulas geradas. "neste caso está sendo usado a equação logarítmica para calcular o número de partículas" 
  double result[n_p][4];                 // Matriz criada para guardar os resultados.
     
  for (int l=0; l < n_p; l++){
  
           x[l] = gRandom->Uniform(0,1);         // Sorteio das partículas.                                   
                  
           for(int k=0; k < size; k++){
 
                   if(x[l] < p[k]){              // se o numero sorteado "x[]" for menor que "n" normalizado "p[]" utilizar o valor de massa  "massa[]"
 
  
/*********************************************************************************************************************************************************************************************************************************************************************************************************************************************************/  

                   // Momento transversal
               
                   // Carrega os valores do gráfico para a massa[] determinada 
		   for(int valpt=0 ; valpt <= divpt ; valpt++){
   
                           gIntegralpt->GetX()[valpt] = mxpt[valpt][k];                    // Obtem os valores do gráfico para determinada massa m0. 
                           gIntegralpt->GetY()[valpt] = mypt[valpt][k];
                    
                   }
               
               
                   for(int mt=0;mt< n_p*10;mt++){
 
                           double Ypt = gRandom->Uniform(xminpt,xmaxpt);                         // Sorteio do valor da abscissa
                           double Rpt = gIntegralpt->Eval(Ypt,0,"s");                            // Obtem o valor da ordenada do gráfico.
                           double temppt = gRandom->Uniform(0,2.5);                              // sorteio de uma valor para a ordenada 
 
                           if(temppt <= Rpt){                       // Se o valor da ordenada "temppt" for menor ou igual a o valor da ordenada "Rpt" atribuir o valor da abscissa "Ypt"
    
                                     result[l][2] = Ypt;              // Lugar na matriz onde é guardado o valor da ordenada de pt.
                                 
                                     break;
                           }
 
                   }
 

/*********************************************************************************************************************************************************************************************************************************************************************************************************************************************************/  
                   
                   
                   // Rapidez 
  
           
           
              // Carrega os valores do gráfico para a massa[] determinada 
                   for(int val=0 ; val <= pint ; val++){
   
                           gIntegral->GetX()[val] = mx[val][k];                    // Obtem os valores do gráfico para determinada massa. 
                           gIntegral->GetY()[val] = my[val][k];
                    
                           }
 
                  
               
 
                      for(int rp=0;rp< n_p*100;rp++){    
                  
                                double Y = gRandom->Uniform(xmin,xmax);                         // Sorteio do valor de x
                                double R = (gIntegral->Eval(Y,0,"s"))/(resultado_integral[k]);  // Normaliza o valor de y do gráfico para dado x.
                                double temp = gRandom->Uniform(0,1);                          // sorteio de uma valor de y 

                                if(temp <= R){
    
                                      result[l][3] = Y;              // Lugar na matriz onde é guardado o valor de Y da rapidez.
                              
                                      break;
                               }
		      }
                              
                          
                   result[l][0] = massa[k];               // Valor da massa da partícula
                   result[l][1] = carga[k];               // Valor da carga da partícula

                   break;
 
                   }
           }
  } 
      

// salvando os valores em um arquivo .txt.
      
  for(const int c=0; c < n_p; c++){

                myfile << result[c][0] << "\t" << result[c][1] << "\t" << result[c][2] << "\t" << result[c][3] <<endl;       // Salvando os valores em um arquivo externo       
    
                }
  }                   // Fim de "for(g)"
   
  myfile.close();
     
  }
