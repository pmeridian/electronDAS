{
  std::cout << "Set defaults" << std::endl;
  gSystem->Load("libPhysics");
  gSystem->Load("libRooFit");
  gSystem->SetIncludePath("-I${ROOFITSYS}/include/");                                                                              
  //                                                                                                                                                                               
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);  // Show overflow, underflow + SumOfWeights	\
    
    //  gStyle->SetStatStyle(0); // for a completely transparent stat box \
    
    gStyle->SetOptFit(0);
    gStyle->SetOptFile(1);
    gStyle->SetOptTitle(0);
}


