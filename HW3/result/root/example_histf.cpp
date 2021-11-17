/// @file example_hist.cpp
/// @brief example root scripts to read binary file and 
/// make a histogram (single precision version)
/// @author pistack (Junho Lee)
/// @date 2021. 11. 17.

{
TFile *file = new TFile("example_histf.root", "recreate");
TTree *tree = new TTree("example_histf", "tree title");
// canvas for zeta
TCanvas *plot_monitor = new TCanvas("monitor", "Monitor", 
800, 600);
// canvas for theta
TCanvas *plot_dist = new TCanvas("distribution", "Distribution", 
800, 600);

int burn = 10000; // discard first 10000 markov chain
std::size_t n=1; // n th move
Float_t action; // action
Float_t e; // estimated error

tree -> Branch("n", &n);
tree -> Branch("S", &action);
tree -> Branch("Error", &e);

std::ifstream infile("../monitor_bin/monitor6f_b.bin", ios::in | ios::binary);
std::size_t size; // size of file
infile.seekg(0, ios::end);
size = infile.tellg();
infile.seekg(0, ios::beg);
char * buffer = new char [size];
infile.read(buffer, size);
infile.close();
for(std::size_t j = 2*sizeof(Float_t)*burn; j < size; j = j + 2*sizeof(Float_t))
{
    char * tmp = new char [sizeof(Float_t)];
    for(int k = 0; k < sizeof(Float_t); ++k)
    tmp[k] = buffer[j+k];
    action = *reinterpret_cast<Float_t*>(tmp);
    for(int k = 0; k < sizeof(Float_t); ++k)
    tmp[k] = buffer[j+k+sizeof(Float_t)];
    e = *reinterpret_cast<Float_t*>(tmp);
    tree -> Fill();
    ++n;
}
cout << n << endl;

tree -> Write(); // write tree to file
plot_monitor -> cd();
tree -> Draw("S:n", ""); // monitor plot
// histogram
plot_dist -> cd();
tree -> Draw("S", "");
file -> Close(); // close file 
}

