{
    float del_time = 0.5;

    std::vector<std::string> filenames = {
        "3k 2s.txt",
        "3k 3s.txt",
        "3k 4s.txt",
        "3k 5s.txt",
        "35k 3s.txt",
        "35k 4s.txt",
        "35k 5s.txt",
        "4k 2s.txt",
        "4k 3s.txt",
        "4k 4s.txt",
        "4k 5s.txt",
        "45k 2s.txt",
        "45k 3s.txt",
        "45k 4s.txt",
        "45k 5s.txt",
     };

     std::vector<std::string> titles = {
        "Charging frequency(1/s) vs average time(s) for 3000 volt, 2 squeezes",
        "Charging frequency(1/s) vs average time(s) for 3000 volt, 3 squeezes",
        "Charging frequency(1/s) vs average time(s) for 3000 volt, 4 squeezes",
        "Charging frequency(1/s) vs average time(s) for 3000 volt, 5 squeezes",
        "Charging frequency(1/s) vs average time(s) for 3500 volt, 3 squeezes",
        "Charging frequency(1/s) vs average time(s) for 3500 volt, 4 squeezes",
        "Charging frequency(1/s) vs average time(s) for 3500 volt, 5 squeezes",
        "Charging frequency(1/s) vs average time(s) for 4000 volt, 2 squeezes",
        "Charging frequency(1/s) vs average time(s) for 4000 volt, 3 squeezes",
        "Charging frequency(1/s) vs average time(s) for 4000 volt, 4 squeezes",
        "Charging frequency(1/s) vs average time(s) for 4000 volt, 5 squeezes",
        "Charging frequency(1/s) vs average time(s) for 4500 volt, 2 squeezes",
        "Charging frequency(1/s) vs average time(s) for 4500 volt, 3 squeezes",
        "Charging frequency(1/s) vs average time(s) for 4500 volt, 4 squeezes",
        "Charging frequency(1/s) vs average time(s) for 4500 volt, 5 squeezes",
    };

      std::vector<float> halfLifes;
      std::vector<float> uncertainties;
      float weightedupper = 0.0;
      float weightedlower = 0.0;




    for(int j =0;j<filenames.size();j++){
 
        int ndata;

        std::ifstream inputFile(filenames[j].c_str());
        std::string line;
        int linecount = 0; // Initialize a counter
        while (std::getline(inputFile, line)) {  // Loop through each line in the file
            linecount++;  // Incrementing line count for each line read
    }
     

        ndata = linecount;
    
        double times[ndata]; 
        double x[ndata], y[ndata], sx[ndata], sy[ndata];


        std::ifstream file(filenames[j].c_str());
        for (int l = 0; l < ndata; ++l) {
            file >> times[l];
    }
     

       

        for (int i=0; i<ndata - 1; ++i){
            if(((times[i+1]+times[i])/2)<0){
                x[i]=-1*(times[i+1]+times[i])/2;
            }
            else if(0<((times[i+1]+times[i])/2)<0.1){
                x[i] = 0.01;
            }
            else{
                x[i]=(times[i+1]+times[i])/2;
            }
            y[i] = 1/(times[i+1]-times[i]);
            sx[i] = sqrt(2*((1/4)*(1/4)));
            sy[i] = y[i]*y[i]*sqrt(2*(del_time*del_time));
    }

    
        TGraphErrors *mygraph = new TGraphErrors(ndata-1,x,y,sx,sy);
        TCanvas *c = new TCanvas();
        mygraph->Draw("A*");
        mygraph->SetTitle(titles[j].c_str());


        mygraph->GetXaxis()->SetTitle("Average Time(s)");
        mygraph->GetYaxis()->SetTitle("Charging Frequency(1/s)");
    
    
        TF1 *expo_fit = new TF1("expo_fit","[0]*exp([1]*x)",4,400);
        float par_0 = 0.2;
        float par_1 = -TMath::Log(2)/55;
        expo_fit->SetParameters(par_0,par_1);
        expo_fit->SetParNames("Scale","Exponential Decay");
        expo_fit->SetLineColor(kRed);
        expo_fit->SetLineWidth(2);
        mygraph->Fit("expo_fit","R");
        double lambda_error = expo_fit->GetParError(1);
        double scale = expo_fit->GetParameter(0);
        double lambda = expo_fit->GetParameter(1);
        gStyle->SetOptFit(1);
        cout << "Scale = " << scale  << "\n";
        cout << "Lambda = "<< lambda <<"+-" << lambda_error <<"\n";

        float result = -log(2)/lambda;
        float uncertainty = (log(2)/(lambda*lambda))*lambda_error;
        halfLifes.push_back(result);
        uncertainties.push_back(uncertainty);

   
    }

    for(int p=0;p<15;p++){
        //cout<<"half lifes are "<<halfLifes[p]<< "uncertainties are"<<uncertainties[p]<<endl;
        weightedupper += halfLifes[p] /(pow(uncertainties[p], 2));
        weightedlower += 1/(pow(uncertainties[p], 2));
        
    }
    

    float halfLife = weightedupper/weightedlower;
    float finalUncertainty = sqrt(1/weightedlower);
    cout<<"so half life of Radon is found to be "<<halfLife<<"s +-"<<finalUncertainty<<"s"<<endl;
    


}