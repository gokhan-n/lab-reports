{
  float del_r = 0.002;
  float del_v = 0.1;
  float del_a = 0.01;
  float constant = 0.0006925657;
  float constant2 = constant * constant;
  float K = 2*0.2*.02*5*5*5/(1.257e-6*1.257e-6*154*154*4*4*4);
  std::vector<float> results;
  std::vector<float> resultSigmas;
  float meanupper;
  float meanlower;




  std::vector<std::string> filenames = {
        "qmr2.txt",
        "qmr3.txt",
        "qmr4.txt",
        "qmr5.txt",
        "qmallr.txt"


     };

     std::vector<std::string> titles = {
        "2V/r^2 vs B^2 for r = 2 cm",
        "2V/r^2 vs B^2 for r = 3 cm",
        "2V/r^2 vs B^2 for r = 4 cm",
        "2V/r^2 vs B^2 for r = 5 cm",
        "2V/r^2 vs B^2 for all r values taken"
    };
    int ndata;  



  for (int k=0; k<5;k++){
    if (k != 4){
      ndata = 4;
    }
    else {
      ndata =16;
    }
    float amp[ndata], volt[ndata], radius[ndata];
    float x[ndata], y[ndata], sx[ndata], sy[ndata];

    
    std::ifstream file(filenames[k].c_str());

    for (int l = 0; l < ndata; ++l) {
        file >> amp[l] >> volt[l] >> radius[l];
    }
  
    for (int i=0; i<ndata; ++i){ // Error Propagation
      float r2 = radius[i]*radius[i];
      float a2 = amp[i]*amp[i];
      float del_r2 = r2 * 2 * del_r / radius[i]; 
      x[i] = constant2 * a2;
      y[i] = 2 * volt[i] / r2;
      sx[i] = x[i] * sqrt((del_v / volt[i])*(del_v / volt[i]) + (del_r2 / r2)*(del_r2 / r2));
      sy[i] = y[i] * 2 * del_a / a2;
  }
  

    float weight = 0;
    float totw = 0; // Total weight
    float xybar = 0, xbar = 0, ybar = 0, x2bar = 0, y2bar = 0; // weighted averages

    for (int j=0; j<ndata; ++j) {
      weight = 1./((sy[j]*sy[j]) + (sx[j]*sx[j]));
      totw += weight;
      xybar += (x[j]*y[j]*weight);
      xbar += (x[j]*weight);
      ybar += (y[j]*weight);
      x2bar += (x[j]*x[j]*weight);
      y2bar += (y[j]*y[j]*weight);
    } 


    float sy2bar = ndata / totw; // weighted average error squared
    float slope = (xybar - xbar*ybar) / (x2bar - xbar*xbar);
    float itcpt = ybar - slope * xbar;
    float slopeerr = sqrt ( sy2bar / (ndata * (x2bar - xbar*xbar) ) );
    float itcpterr = sqrt ( x2bar ) * slopeerr;
  

    TGraphErrors *mygraph = new TGraphErrors(ndata,x,y,sx,sy);
    TCanvas *c = new TCanvas();
    mygraph->SetTitle(titles[k].c_str());
    mygraph->Draw("A*");

    gStyle->SetOptFit(1111);
 


    TF1 *fnew = new TF1("fnew","[0]*x+[1]",0,900e3);
    fnew->SetParameters(slope, itcpt); 
    mygraph->Fit(fnew);
    results.push_back(fnew->GetParameter(0));
    resultSigmas.push_back(fnew->GetParError(1));
  
  }

  
  for (int u = 0; u < 5; ++u) {
       
        meanupper += results[u]/(pow(resultSigmas[u],2));
        meanlower += 1/(pow(resultSigmas[u],2));

  }
  
  float result = meanupper/meanlower;
  float resultSigma = sqrt(1/meanlower);

  cout << "the charge to mass ratio off the electron is found to be " << result <<"C/kg"<< "+- "<<resultSigma<<"C/kg"<<"\n";
  
  }
