{

  int num_channels = 16;

  stringstream infile_name_stream;
  infile_name_stream << "/farmshare/user_data/blenardo/struck_data/"
              << "tier1_SIS3316Raw_20181003013401_Ambe_CH1NaI1700V_CH0PSD1100V_22inchSTD_PSDTrigg150mVTh_PSDX10_125MHz__1-ngm.root";

  stringstream outfile_name_template_stream;
  outfile_name_template_stream << "/farmshare/user_data/blenardo/struck_data/"
            << "tier1_SIS3316Raw_20181003013401_Ambe_CH1NaI1700V_CH0PSD1100V_22inchSTD_PSDTrigg150mVTh_PSDX10_125MHz__1-ngm_";

  string infile_name = infile_name_stream.str();
  string outfile_name_template = outfile_name_template_stream.str();

  TFile * f = new TFile( infile_name.c_str() );

  
  TTree * hit_tree = (TTree *)f->Get("HitTree");  
  NGMHitv8 * hit = new NGMHitv8();

  hit_tree->SetBranchAddress("HitTree",&hit);


  printf("Number of entries: %d\n",hit_tree->GetEntries());
  printf("This means there are %d events\n",hit_tree->GetEntries()/num_channels);

  n_entries = hit_tree->GetEntries();

  int i_entry = 0;
  int i_outfile = 0;
  char outfilename[100];

  while( i_entry < n_entries ) {
          
       sprintf(outfilename,"%s%04d.root",outfile_name_template.c_str(),i_outfile);
       TFile * out_file = new TFile(outfilename,"RECREATE");
       TTree * out_hit_tree = new TTree();

       out_hit_tree->Branch("HitTree",&hit);
       out_hit_tree->SetName("HitTree");     
 
       int j_entry_out = 0;
       printf("============================= OUTFILE %d ================================\n",i_outfile);
       printf("%s\n",outfilename); 
       while( j_entry_out < 10000 && i_entry + j_entry_out < n_entries ) {
          hit_tree->GetEntry(i_entry);
          out_hit_tree->Fill(); 
          j_entry_out += 1;
          i_entry += 1;
       }
       out_file->cd();
       out_hit_tree->Write();
       out_file->Close();
       delete out_hit_tree;
       delete out_file;
       i_entry += 1;
       i_outfile += 1;
  }

}
