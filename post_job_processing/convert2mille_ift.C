#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <cmath> 

using namespace std;

class Mille 
{
public:
  Mille(const char *outFileName, bool asBinary = true, bool writeZero = true);
  ~Mille();

  void mille(int NLC, const float *derLc, int NGL, const float *derGl,
	     const int *label, float rMeas, float sigma);
  void special(int nSpecial, const float *floatings, const int *integers);
  void kill();
  void end();

private:
  void newSet();
  bool checkBufferSize(int nLocal, int nGlobal);

  std::ofstream myOutFile; ///< C-binary for output
  bool myAsBinary;         ///< if false output as text
  bool myWriteZero;        ///< if true also write out derivatives/labels ==0
  /// buffer size for ints and floats
  enum {myBufferSize = 5000};  ///< buffer size for ints and floats
  int   myBufferInt[myBufferSize];   ///< to collect labels etc.
  float myBufferFloat[myBufferSize]; ///< to collect derivatives etc.
  int   myBufferPos; ///< position in buffer
  bool  myHasSpecial; ///< if true, special(..) already called for this record
  /// largest label allowed: 2^31 - 1
  enum {myMaxLabel = (0xFFFFFFFF - (1 << 31))};
};

Mille::Mille(const char *outFileName, bool asBinary, bool writeZero) : 
  myOutFile(outFileName, (asBinary ? (std::ios::binary | std::ios::out) : std::ios::out)),
  myAsBinary(asBinary), myWriteZero(writeZero), myBufferPos(-1), myHasSpecial(false)
{
  // Instead myBufferPos(-1), myHasSpecial(false) and the following two lines
  // we could call newSet() and kill()...
  myBufferInt[0]   = 0;
  myBufferFloat[0] = 0.;

  if (!myOutFile.is_open()) {
    std::cerr << "Mille::Mille: Could not open " << outFileName 
	      << " as output file." << std::endl;
  }
}

//___________________________________________________________________________
/// Closes file.
Mille::~Mille()
{
  myOutFile.close();
}

//___________________________________________________________________________
/// Add measurement to buffer.
/**
 * \param[in]    NLC    number of local derivatives
 * \param[in]    derLc  local derivatives
 * \param[in]    NGL    number of global derivatives
 * \param[in]    derGl  global derivatives
 * \param[in]    label  global labels
 * \param[in]    rMeas  measurement (residuum)
 * \param[in]    sigma  error
 */
void Mille::mille(int NLC, const float *derLc,
		  int NGL, const float *derGl, const int *label,
		  float rMeas, float sigma)
{
  if (sigma <= 0.) return;
  if (myBufferPos == -1) this->newSet(); // start, e.g. new track
  if (!this->checkBufferSize(NLC, NGL)) return;

  // first store measurement
  ++myBufferPos;
  myBufferFloat[myBufferPos] = rMeas;
  myBufferInt  [myBufferPos] = 0;

  // store local derivatives and local 'lables' 1,...,NLC
  for (int i = 0; i < NLC; ++i) {
    if (derLc[i] || myWriteZero) { // by default store only non-zero derivatives
      ++myBufferPos;
      myBufferFloat[myBufferPos] = derLc[i]; // local derivatives
      myBufferInt  [myBufferPos] = i+1;      // index of local parameter
    }
  }

  // store uncertainty of measurement in between locals and globals
  ++myBufferPos;
  myBufferFloat[myBufferPos] = sigma;
  myBufferInt  [myBufferPos] = 0;

  // store global derivatives and their labels
  for (int i = 0; i < NGL; ++i) {
    if (derGl[i] || myWriteZero) { // by default store only non-zero derivatives
      if ((label[i] > 0 || myWriteZero) && label[i] <= myMaxLabel) { // and for valid labels
	++myBufferPos;
	myBufferFloat[myBufferPos] = derGl[i]; // global derivatives
	myBufferInt  [myBufferPos] = label[i]; // index of global parameter
      } else {
	std::cerr << "Mille::mille: Invalid label " << label[i] 
		  << " <= 0 or > " << myMaxLabel << std::endl; 
      }
    }
  }
}

//___________________________________________________________________________
/// Add special data to buffer.
/**
 * \param[in]    nSpecial   number of floats/ints
 * \param[in]    floatings  floats
 * \param[in]    integers   ints
 */
void Mille::special(int nSpecial, const float *floatings, const int *integers)
{
  if (nSpecial == 0) return;
  if (myBufferPos == -1) this->newSet(); // start, e.g. new track
  if (myHasSpecial) {
    std::cerr << "Mille::special: Special values already stored for this record."
	      << std::endl; 
    return;
  }
  if (!this->checkBufferSize(nSpecial, 0)) return;
  myHasSpecial = true; // after newSet() (Note: MILLSP sets to buffer position...)

  //  myBufferFloat[.]  | myBufferInt[.]
  // ------------------------------------
  //      0.0           |      0
  //  -float(nSpecial)  |      0
  //  The above indicates special data, following are nSpecial floating and nSpecial integer data.

  ++myBufferPos; // zero pair
  myBufferFloat[myBufferPos] = 0.;
  myBufferInt  [myBufferPos] = 0;

  ++myBufferPos; // nSpecial and zero
  myBufferFloat[myBufferPos] = -nSpecial; // automatic conversion to float
  myBufferInt  [myBufferPos] = 0;

  for (int i = 0; i < nSpecial; ++i) {
    ++myBufferPos;
    myBufferFloat[myBufferPos] = floatings[i];
    myBufferInt  [myBufferPos] = integers[i];
  }
}

//___________________________________________________________________________
/// Reset buffers, i.e. kill derivatives accumulated for current set.
void Mille::kill()
{
  myBufferPos = -1;
}

//___________________________________________________________________________
/// Write buffer (set of derivatives with same local parameters) to file.
void Mille::end()
{
  if (myBufferPos > 0) { // only if anything stored...
    const int numWordsToWrite = (myBufferPos + 1)*2;

    if (myAsBinary) {
      myOutFile.write(reinterpret_cast<const char*>(&numWordsToWrite), 
		      sizeof(numWordsToWrite));
      myOutFile.write(reinterpret_cast<char*>(myBufferFloat), 
		      (myBufferPos+1) * sizeof(myBufferFloat[0]));
      myOutFile.write(reinterpret_cast<char*>(myBufferInt), 
		      (myBufferPos+1) * sizeof(myBufferInt[0]));
    } else {
      myOutFile << numWordsToWrite << "\n";
      for (int i = 0; i < myBufferPos+1; ++i) {
	myOutFile << myBufferFloat[i] << " ";
      }
      myOutFile << "\n";
      
      for (int i = 0; i < myBufferPos+1; ++i) {
	myOutFile << myBufferInt[i] << " ";
      }
      myOutFile << "\n";
    }
  }
  myBufferPos = -1; // reset buffer for next set of derivatives
}

//___________________________________________________________________________
/// Initialize for new set of locals, e.g. new track.
void Mille::newSet()
{
  myBufferPos = 0;
  myHasSpecial = false;
  myBufferFloat[0] = 0.0;
  myBufferInt  [0] = 0;   // position 0 used as error counter
}

//___________________________________________________________________________
/// Enough space for next nLocal + nGlobal derivatives incl. measurement?
/**
 * \param[in]   nLocal  number of local derivatives
 * \param[in]   nGlobal number of global derivatives
 * \return      true if sufficient space available (else false)
 */
bool Mille::checkBufferSize(int nLocal, int nGlobal)
{
  if (myBufferPos + nLocal + nGlobal + 2 >= myBufferSize) {
    ++(myBufferInt[0]); // increase error count
    std::cerr << "Mille::checkBufferSize: Buffer too short (" 
	      << myBufferSize << "),"
	      << "\n need space for nLocal (" << nLocal<< ")"
	      << "/nGlobal (" << nGlobal << ") local/global derivatives, " 
	      << myBufferPos + 1 << " already stored!"
	      << std::endl;
    return false;
  } else {
    return true;
  }
}

void convert2mille_v2_withlocaly_ift(){
  //data23
  // Mille mille_file("/afs/cern.ch/user/a/agarabag/EOS/millepede/ift_align_data/mp2input_mc_data_misalgn_5var_new.bin");
  Mille mille_file("/home/agarabag/millepede/ift_mc_align/mp2input_ke_ift_spoint_ini.bin");

  //data22
  //TFile* f1=new TFile("/afs/cern.ch/user/k/keli/eos/Faser/alignment/global/8023_8025_8115_8301_8730_9073/kfalignment_data_iter5_noIFT_noZ.root");
  //Mille mille_file("/afs/cern.ch/user/k/keli/eos/Faser/alignment/global/8023_8025_8115_8301_8730_9073/mp2input.bin");
  //MC
  //TFile* f1=new TFile("/afs/cern.ch/user/k/keli/work/FASER/test/kfalignment_mc_0.root");
  //TFile* f1=new TFile("/afs/cern.ch/user/k/keli/eos/Faser/alignment/global/misalign_MC/kfalignment_mc_iter1_2ndf.root");
  //TFile* f1=new TFile("/afs/cern.ch/user/k/keli/eos/Faser/alignment/global/misalign_MC/kfalignment_mc_iter3.root");
  //TFile* f1=new TFile("/afs/cern.ch/user/k/keli/eos/Faser/alignment/global/misalign_MC/inputformp2.root");
  //TFile* f1=new TFile("/afs/cern.ch/user/k/keli/eos/Faser/alignment/global/misalign_MC/inputformp2_iter0.root");
  //Mille mille_file("/afs/cern.ch/user/k/keli/eos/Faser/alignment/global/misalign_MC/mp2input.bin");
  int Nfile=48;

  // for(int fileId=0;fileId<Nfile;fileId++){
  // string InputFileName="/eos/user/t/tarai/software/result/kfalignment_010738_"+std::to_string(fileId)+".root";
  // TFile* f1=new TFile(InputFileName.c_str());
  //  TFile* f1=new TFile("/afs/cern.ch/user/a/agarabag/EOS/ift_alignment/run/data_3_station_align/kfalignment_data_skimmed.root");
  TFile* f1=new TFile("/data/agarabag/alignment_output/new_alg_spoint_ift_ini/kfalignment_mc.root");


  TTree* t1=(TTree*)f1->Get("trackParam");
  bool dump6ndf_modules=false;
  bool dumpz_modules=false;
  bool dump3ndf_modules=false;

  bool dumplayers=true;
  bool dump6ndf_layers=true;
  bool dumpz_layers=false;
  
  bool dumpstations=false;
  bool dump5ndf_stations=false;
  bool dumpz_station = false;

  double m_fitParam_x=0;
  double m_fitParam_y=0;
  double m_fitParam_pull_x=0;
  double m_fitParam_pull_y=0;
  double m_fitParam_chi2=0;
  double m_fitParam_ndf=0;
  double m_fitParam_nMeasurements=0;
  double m_fitParam_px=0;
  double m_fitParam_py=0;
  double m_fitParam_pz=0;
  double m_fitParam_charge=0;
  std::vector<double>* m_ift_id=0;
  std::vector<double>* m_ift_local_x=0;
  std::vector<double>* m_ift_local_y=0;
  std::vector<double>* m_ift_local_xe=0;
  std::vector<double>* m_ift_local_ye=0;
  std::vector<double>* m_ift_local_residual_x=0;
  std::vector<double>* m_ift_local_residual_y=0;
  std::vector<double>* m_ift_align_derivative_x_x=0;
  std::vector<double>* m_ift_align_derivative_x_y=0;
  std::vector<double>* m_ift_align_derivative_x_z=0;
  std::vector<double>* m_ift_align_derivative_x_rx=0;
  std::vector<double>* m_ift_align_derivative_x_ry=0;
  std::vector<double>* m_ift_align_derivative_x_rz=0;
  std::vector<double>* m_ift_align_derivative_x_par_x=0;
  std::vector<double>* m_ift_align_derivative_x_par_y=0;
  std::vector<double>* m_ift_align_derivative_x_par_theta=0;
  std::vector<double>* m_ift_align_derivative_x_par_phi=0;
  std::vector<double>* m_ift_align_derivative_x_par_qop=0;
  std::vector<double>* m_ift_align_derivative_global_y_x=0;
  std::vector<double>* m_ift_align_derivative_global_y_y=0;
  std::vector<double>* m_ift_align_derivative_global_y_z=0;
  std::vector<double>* m_ift_align_derivative_global_y_rx=0;
  std::vector<double>* m_ift_align_derivative_global_y_ry=0;
  std::vector<double>* m_ift_align_derivative_global_y_rz=0;
  std::vector<double>* m_ift_align_derivative_y_x=0;
  std::vector<double>* m_ift_align_derivative_y_y=0;
  std::vector<double>* m_ift_align_derivative_y_z=0;
  std::vector<double>* m_ift_align_derivative_y_rx=0;
  std::vector<double>* m_ift_align_derivative_y_ry=0;
  std::vector<double>* m_ift_align_derivative_y_rz=0;
  std::vector<double>* m_ift_align_derivative_y_par_x=0;
  std::vector<double>* m_ift_align_derivative_y_par_y=0;
  std::vector<double>* m_ift_align_derivative_y_par_theta=0;
  std::vector<double>* m_ift_align_derivative_y_par_phi=0;
  std::vector<double>* m_ift_align_derivative_y_par_qop=0;
  std::vector<double>* m_ift_align_derivative_global_x_x=0;
  std::vector<double>* m_ift_align_derivative_global_x_y=0;
  std::vector<double>* m_ift_align_derivative_global_x_z=0;
  std::vector<double>* m_ift_align_derivative_global_x_rx=0;
  std::vector<double>* m_ift_align_derivative_global_x_ry=0;
  std::vector<double>* m_ift_align_derivative_global_x_rz=0;

  t1->SetBranchAddress("fitParam_x",&m_fitParam_x);
  t1->SetBranchAddress("fitParam_y",&m_fitParam_y);
  t1->SetBranchAddress("fitParam_chi2",&m_fitParam_chi2);
  t1->SetBranchAddress("fitParam_pz",&m_fitParam_pz);
  t1->SetBranchAddress("fitParam_px",&m_fitParam_px);
  t1->SetBranchAddress("fitParam_py",&m_fitParam_py);
  t1->SetBranchAddress("fitParam_align_ift_id", &m_ift_id);
  t1->SetBranchAddress("fitParam_align_ift_local_measured_x", &m_ift_local_x);
  t1->SetBranchAddress("fitParam_align_ift_local_measured_y", &m_ift_local_y);
  t1->SetBranchAddress("fitParam_align_ift_local_measured_xe", &m_ift_local_xe);
  t1->SetBranchAddress("fitParam_align_ift_local_measured_ye", &m_ift_local_ye);
  t1->SetBranchAddress("fitParam_align_ift_local_residual_x", &m_ift_local_residual_x);
  t1->SetBranchAddress("fitParam_align_ift_local_residual_y", &m_ift_local_residual_y);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_x_par_x",     &m_ift_align_derivative_x_par_x);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_x_par_y",     &m_ift_align_derivative_x_par_y);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_x_par_theta", &m_ift_align_derivative_x_par_theta);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_x_par_phi",   &m_ift_align_derivative_x_par_phi);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_x_par_qop",   &m_ift_align_derivative_x_par_qop);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_x_x",         &m_ift_align_derivative_x_x);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_x_y",         &m_ift_align_derivative_x_y);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_x_z",         &m_ift_align_derivative_x_z);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_x_rx",        &m_ift_align_derivative_x_rx);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_x_ry",        &m_ift_align_derivative_x_ry);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_x_rz",        &m_ift_align_derivative_x_rz);
  t1->SetBranchAddress("fitParam_align_ift_global_derivation_y_x",        &m_ift_align_derivative_global_y_x);
  t1->SetBranchAddress("fitParam_align_ift_global_derivation_y_y",        &m_ift_align_derivative_global_y_y);
  t1->SetBranchAddress("fitParam_align_ift_global_derivation_y_z",        &m_ift_align_derivative_global_y_z);
  t1->SetBranchAddress("fitParam_align_ift_global_derivation_y_rx",       &m_ift_align_derivative_global_y_rx);
  t1->SetBranchAddress("fitParam_align_ift_global_derivation_y_ry",       &m_ift_align_derivative_global_y_ry);
  t1->SetBranchAddress("fitParam_align_ift_global_derivation_y_rz",       &m_ift_align_derivative_global_y_rz);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_y_par_x",     &m_ift_align_derivative_y_par_x);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_y_par_y",     &m_ift_align_derivative_y_par_y);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_y_par_theta", &m_ift_align_derivative_y_par_theta);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_y_par_phi",   &m_ift_align_derivative_y_par_phi);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_y_par_qop",   &m_ift_align_derivative_y_par_qop);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_y_x",         &m_ift_align_derivative_y_x);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_y_y",         &m_ift_align_derivative_y_y);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_y_z",         &m_ift_align_derivative_y_z);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_y_rx",        &m_ift_align_derivative_y_rx);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_y_ry",        &m_ift_align_derivative_y_ry);
  t1->SetBranchAddress("fitParam_align_ift_local_derivation_y_rz",        &m_ift_align_derivative_y_rz);
  t1->SetBranchAddress("fitParam_align_ift_global_derivation_x_x",        &m_ift_align_derivative_global_x_x);
  t1->SetBranchAddress("fitParam_align_ift_global_derivation_x_y",        &m_ift_align_derivative_global_x_y);
  t1->SetBranchAddress("fitParam_align_ift_global_derivation_x_z",        &m_ift_align_derivative_global_x_z);
  t1->SetBranchAddress("fitParam_align_ift_global_derivation_x_rx",       &m_ift_align_derivative_global_x_rx);
  t1->SetBranchAddress("fitParam_align_ift_global_derivation_x_ry",       &m_ift_align_derivative_global_x_ry);
  t1->SetBranchAddress("fitParam_align_ift_global_derivation_x_rz",       &m_ift_align_derivative_global_x_rz);


  std::vector<double>* m_fitParam_align_global_derivative_y_x=0;
  std::vector<double>* m_fitParam_align_global_derivative_y_y=0;
  std::vector<double>* m_fitParam_align_global_derivative_y_z=0;
  std::vector<double>* m_fitParam_align_global_derivative_y_rx=0;
  std::vector<double>* m_fitParam_align_global_derivative_y_ry=0;
  std::vector<double>* m_fitParam_align_global_derivative_y_rz=0;
  std::vector<double>* m_fitParam_align_local_derivative_x_x=0;
  std::vector<double>* m_fitParam_align_local_derivative_x_y=0;
  std::vector<double>* m_fitParam_align_local_derivative_x_z=0;
  std::vector<double>* m_fitParam_align_local_derivative_x_rx=0;
  std::vector<double>* m_fitParam_align_local_derivative_x_ry=0;
  std::vector<double>* m_fitParam_align_local_derivative_x_rz=0;
  std::vector<double>* m_fitParam_align_local_residual_x=0;
  std::vector<double>* m_fitParam_align_local_measured_x=0;
  std::vector<double>* m_fitParam_align_local_measured_xe=0;
  std::vector<double>* m_fitParam_align_id=0;
  std::vector<double>* m_fitParam_align_local_derivative_x_par_x=0;
  std::vector<double>* m_fitParam_align_local_derivative_x_par_y=0;
  std::vector<double>* m_fitParam_align_local_derivative_x_par_theta=0;
  std::vector<double>* m_fitParam_align_local_derivative_x_par_phi=0;
  std::vector<double>* m_fitParam_align_local_derivative_x_par_qop=0;
  t1->SetBranchAddress("fitParam_align_global_derivation_y_x",&m_fitParam_align_global_derivative_y_x);
  t1->SetBranchAddress("fitParam_align_global_derivation_y_y",&m_fitParam_align_global_derivative_y_y);
  t1->SetBranchAddress("fitParam_align_global_derivation_y_z",&m_fitParam_align_global_derivative_y_z);
  t1->SetBranchAddress("fitParam_align_global_derivation_y_rx",&m_fitParam_align_global_derivative_y_rx);
  t1->SetBranchAddress("fitParam_align_global_derivation_y_ry",&m_fitParam_align_global_derivative_y_ry);
  t1->SetBranchAddress("fitParam_align_global_derivation_y_rz",&m_fitParam_align_global_derivative_y_rz);
  t1->SetBranchAddress("fitParam_align_local_derivation_x_x",&m_fitParam_align_local_derivative_x_x);
  t1->SetBranchAddress("fitParam_align_local_derivation_x_y",&m_fitParam_align_local_derivative_x_y);
  t1->SetBranchAddress("fitParam_align_local_derivation_x_z",&m_fitParam_align_local_derivative_x_z);
  t1->SetBranchAddress("fitParam_align_local_derivation_x_rx",&m_fitParam_align_local_derivative_x_rx);
  t1->SetBranchAddress("fitParam_align_local_derivation_x_ry",&m_fitParam_align_local_derivative_x_ry);
  t1->SetBranchAddress("fitParam_align_local_derivation_x_rz",&m_fitParam_align_local_derivative_x_rz);
  t1->SetBranchAddress("fitParam_align_local_residual_x",&m_fitParam_align_local_residual_x);
  t1->SetBranchAddress("fitParam_align_local_measured_x",&m_fitParam_align_local_measured_x);
  t1->SetBranchAddress("fitParam_align_local_measured_xe",&m_fitParam_align_local_measured_xe);
  t1->SetBranchAddress("fitParam_align_id",&m_fitParam_align_id);
  t1->SetBranchAddress("fitParam_align_local_derivation_x_par_x",&m_fitParam_align_local_derivative_x_par_x);
  t1->SetBranchAddress("fitParam_align_local_derivation_x_par_y",&m_fitParam_align_local_derivative_x_par_y);
  t1->SetBranchAddress("fitParam_align_local_derivation_x_par_theta",&m_fitParam_align_local_derivative_x_par_theta);
  t1->SetBranchAddress("fitParam_align_local_derivation_x_par_phi",&m_fitParam_align_local_derivative_x_par_phi);
  t1->SetBranchAddress("fitParam_align_local_derivation_x_par_qop",&m_fitParam_align_local_derivative_x_par_qop);

  int nevt= t1->GetEntries();

  //loop over all the events
  std::vector<int> labels;
  std::vector<float> glo_der;
  std::vector<float> loc_der;
  float resi=0.;
  float resi_e=0.;

  std::vector<int> labels_y;
  std::vector<float> glo_der_y;
  std::vector<float> loc_der_y;
  float resi_y=0.;
  float resi_e_y=0.;

  std::vector<int> labels_og;
  std::vector<float> glo_der_og;
  std::vector<float> loc_der_og;
  float resi_og=0.;
  float resi_e_og=0.;

  // int ioutput=0;
  int nhits[8]={0};

  for(int ievt = 0;ievt< nevt; ++ievt){

    t1->GetEntry(ievt);
    //std::cout<<"like "<<ievt<<" "<<m_fitParam_chi2/m_fitParam_ndf<<" "<<m_fitParam_pz<<" "<<m_ift_id->size()<<std::endl;
    if(m_fitParam_chi2>2000||m_fitParam_pz<100||m_fitParam_pz>5000||m_ift_id->size()<3)continue;
    //if(m_fitParam_chi2>500||m_fitParam_pz<100||m_fitParam_pz>5000||m_ift_id->size()<15)continue;
    // ++ioutput;
    //    if(ioutput>1000)continue;
    for(unsigned int ihit = 0; ihit<m_ift_id->size();++ihit){
      labels.clear();
      glo_der.clear();
      loc_der.clear();
      labels_y.clear();
      glo_der_y.clear();
      loc_der_y.clear();

      if(fabs(m_ift_local_residual_x->at(ihit))>0.5)continue;
      if(m_ift_align_derivative_x_x->at(ihit)<-9000||m_ift_align_derivative_x_rz->at(ihit)<-9000||m_ift_align_derivative_global_y_x->at(ihit)<-9000||m_ift_align_derivative_global_y_y->at(ihit)<-9000||m_ift_align_derivative_global_y_z->at(ihit)<-9000||m_ift_align_derivative_global_y_rx->at(ihit)<-9000||m_ift_align_derivative_global_y_ry->at(ihit)<-9000||m_ift_align_derivative_global_y_rz->at(ihit)<-9000)continue;
      if(m_ift_align_derivative_y_x->at(ihit)<-9000||m_ift_align_derivative_y_rz->at(ihit)<-9000||m_ift_align_derivative_global_x_x->at(ihit)<-9000||m_ift_align_derivative_global_x_y->at(ihit)<-9000||m_ift_align_derivative_global_x_z->at(ihit)<-9000||m_ift_align_derivative_global_x_rx->at(ihit)<-9000||m_ift_align_derivative_global_x_ry->at(ihit)<-9000||m_ift_align_derivative_global_x_rz->at(ihit)<-9000)continue;

      int moduleid = m_ift_id->at(ihit);
      if(moduleid%10==1){
        --moduleid;
      }
      moduleid+=1000;//station from 1 not 0
      //	int layerid = moduleid/100;
      //	if(moduleid/100==10){
      //	  if(nhits[moduleid%100/10]>1000)continue;
      //	  nhits[moduleid%100/10]+=1;
      //	}


      // if(dump6ndf_modules){
	// labels.push_back(moduleid*10+0+1);//millepede can not have label at 0
	// labels.push_back(moduleid*10+1+1);
  //      	labels.push_back(moduleid*10+2+1);
	// labels.push_back(moduleid*10+3+1);
	// labels.push_back(moduleid*10+4+1);
	// glo_der.push_back(m_ift_align_derivative_x_x->at(ihit));
	// glo_der.push_back(m_ift_align_derivative_x_y->at(ihit));
	// glo_der.push_back(m_ift_align_derivative_x_rx->at(ihit));
	// glo_der.push_back(m_ift_align_derivative_x_ry->at(ihit));
	// glo_der.push_back(m_ift_align_derivative_x_rz->at(ihit));

	// if(dumpz_modules){
	//   labels.push_back(moduleid*10+5+1);
	//   glo_der.push_back(m_ift_align_derivative_x_z->at(ihit));
	//   labels_y.push_back(moduleid*10+5+1);
	//   glo_der_y.push_back(m_ift_align_derivative_y_z->at(ihit));
	// }

	// labels_y.push_back(moduleid*10+0+1);
	// labels_y.push_back(moduleid*10+1+1);
	// labels_y.push_back(moduleid*10+2+1);
	// labels_y.push_back(moduleid*10+3+1);
	// labels_y.push_back(moduleid*10+4+1);
	// glo_der_y.push_back(m_ift_align_derivative_y_x->at(ihit));
	// glo_der_y.push_back(m_ift_align_derivative_y_y->at(ihit));
	// glo_der_y.push_back(m_ift_align_derivative_y_rx->at(ihit));
	// glo_der_y.push_back(m_ift_align_derivative_y_ry->at(ihit));
	// glo_der_y.push_back(m_ift_align_derivative_y_rz->at(ihit));
  //     }
  //     else if(dump3ndf_modules){
	// labels.push_back(moduleid*10+0+1);//millepede can not have label at 0
	// labels.push_back(moduleid*10+1+1);
	// labels.push_back(moduleid*10+2+1);
	// glo_der.push_back(m_ift_align_derivative_x_x->at(ihit));
	// glo_der.push_back(m_ift_align_derivative_x_y->at(ihit));
	// glo_der.push_back(m_ift_align_derivative_x_rz->at(ihit));

	// labels_y.push_back(moduleid*10+0+1);
	// labels_y.push_back(moduleid*10+1+1);
	// labels_y.push_back(moduleid*10+2+1);
	// glo_der_y.push_back(m_ift_align_derivative_y_x->at(ihit));
	// glo_der_y.push_back(m_ift_align_derivative_y_y->at(ihit));
	// glo_der_y.push_back(m_ift_align_derivative_y_rz->at(ihit));
  //     }
  //     else{
	// labels.push_back(moduleid*10+0+1);//millepede can not have label at 0
	// labels.push_back(moduleid*10+1+1);
	// glo_der.push_back(m_ift_align_derivative_x_x->at(ihit));
	// glo_der.push_back(m_ift_align_derivative_x_rz->at(ihit));    

	// labels_y.push_back(moduleid*10+0+1);
	// labels_y.push_back(moduleid*10+1+1);
	// glo_der_y.push_back(m_ift_align_derivative_y_x->at(ihit));
	// glo_der_y.push_back(m_ift_align_derivative_y_rz->at(ihit));     
  //     }

      if(dumplayers){
	int layerid = moduleid/100;
	//	std::cout<<"id "<<moduleid<<" "<<layerid<<std::endl;
	if(dump6ndf_layers){
	  if(fabs(m_ift_align_derivative_global_y_rx->at(ihit))>2||fabs(m_ift_align_derivative_global_y_ry->at(ihit))>2)continue;
	  labels.push_back(layerid*10+0+1);//millepede can not have label at 0
	  labels.push_back(layerid*10+1+1);
	  labels.push_back(layerid*10+2+1);
	  labels.push_back(layerid*10+3+1);
	  labels.push_back(layerid*10+4+1);
	  glo_der.push_back(m_ift_align_derivative_global_y_x->at(ihit));
	  glo_der.push_back(m_ift_align_derivative_global_y_y->at(ihit));

	  labels_y.push_back(layerid*10+0+1);
	  labels_y.push_back(layerid*10+1+1);
	  labels_y.push_back(layerid*10+2+1);
	  labels_y.push_back(layerid*10+3+1);
	  labels_y.push_back(layerid*10+4+1);
	  glo_der_y.push_back(m_ift_align_derivative_global_x_x->at(ihit));
	  glo_der_y.push_back(m_ift_align_derivative_global_x_y->at(ihit));
	  if(dumpz_layers){
	    labels.push_back(layerid*10+5+1);
	    glo_der.push_back(m_ift_align_derivative_global_y_z->at(ihit));
	    labels_y.push_back(layerid*10+5+1);
	    glo_der_y.push_back(m_ift_align_derivative_global_x_z->at(ihit));
	  }
	  glo_der.push_back(m_ift_align_derivative_global_y_rx->at(ihit));
	  glo_der.push_back(m_ift_align_derivative_global_y_ry->at(ihit));
	  glo_der.push_back(m_ift_align_derivative_global_y_rz->at(ihit));

	  glo_der_y.push_back(m_ift_align_derivative_global_x_rx->at(ihit));
	  glo_der_y.push_back(m_ift_align_derivative_global_x_ry->at(ihit));
	  glo_der_y.push_back(m_ift_align_derivative_global_x_rz->at(ihit));
	} 
	else{
	  labels.push_back(layerid*10+0+1);//millepede can not have label at 0
	  labels.push_back(layerid*10+1+1);
	  glo_der.push_back(m_ift_align_derivative_global_y_y->at(ihit));
	  glo_der.push_back(m_ift_align_derivative_global_y_rz->at(ihit));
    
	  labels_y.push_back(layerid*10+0+1);
	  labels_y.push_back(layerid*10+1+1);
	  glo_der_y.push_back(m_ift_align_derivative_global_x_y->at(ihit));
	  glo_der_y.push_back(m_ift_align_derivative_global_x_rz->at(ihit));
	}
      }
      if(dumpstations){
	int stationid = moduleid/1000;
	if(dump5ndf_stations){
	  labels.push_back(stationid*10+0+1);//millepede can not have label at 0
	  labels.push_back(stationid*10+1+1);
	  labels.push_back(stationid*10+2+1);
	  labels.push_back(stationid*10+3+1);
	  labels.push_back(stationid*10+4+1);
	  glo_der.push_back(m_ift_align_derivative_global_y_x->at(ihit));
	  glo_der.push_back(m_ift_align_derivative_global_y_y->at(ihit));

	  labels_y.push_back(stationid*10+0+1);
	  labels_y.push_back(stationid*10+1+1);
	  labels_y.push_back(stationid*10+2+1);
	  labels_y.push_back(stationid*10+3+1);
	  labels_y.push_back(stationid*10+4+1);
	  glo_der_y.push_back(m_ift_align_derivative_global_x_x->at(ihit));
	  glo_der_y.push_back(m_ift_align_derivative_global_x_y->at(ihit));

	  if(dumpz_station){
	    labels.push_back(stationid*10+5+1);
	    glo_der.push_back(m_ift_align_derivative_global_y_z->at(ihit));
	    labels_y.push_back(stationid*10+5+1);
	    glo_der_y.push_back(m_ift_align_derivative_global_x_z->at(ihit));
	  }

    glo_der.push_back(m_ift_align_derivative_global_y_rx->at(ihit));
	  glo_der.push_back(m_ift_align_derivative_global_y_ry->at(ihit));
	  glo_der.push_back(m_ift_align_derivative_global_y_rz->at(ihit));

	  glo_der_y.push_back(m_ift_align_derivative_global_x_rx->at(ihit));
	  glo_der_y.push_back(m_ift_align_derivative_global_x_ry->at(ihit));
	  glo_der_y.push_back(m_ift_align_derivative_global_x_rz->at(ihit));
	}
      }

      loc_der.push_back(m_ift_align_derivative_x_par_x->at(ihit));
      loc_der.push_back(m_ift_align_derivative_x_par_y->at(ihit));
      loc_der.push_back(m_ift_align_derivative_x_par_theta->at(ihit));
      loc_der.push_back(m_ift_align_derivative_x_par_phi->at(ihit));
      loc_der.push_back(m_ift_align_derivative_x_par_qop->at(ihit));
      loc_der_y.push_back(m_ift_align_derivative_y_par_x->at(ihit));
      loc_der_y.push_back(m_ift_align_derivative_y_par_y->at(ihit));
      loc_der_y.push_back(m_ift_align_derivative_y_par_theta->at(ihit));
      loc_der_y.push_back(m_ift_align_derivative_y_par_phi->at(ihit));
      loc_der_y.push_back(m_ift_align_derivative_y_par_qop->at(ihit));
      //      std::cout<<"like "<<ievt<<" "<<ihit<<" "<<loc_der.size()<<std::endl;
      float* lcder = loc_der.data();
      float* glder = glo_der.data();
      int* label = labels.data();
      float* lcder_y = loc_der_y.data();
      float* glder_y = glo_der_y.data();
      int* label_y = labels_y.data();
      resi=m_ift_local_residual_x->at(ihit);
      resi_e=m_ift_local_xe->at(ihit);
      resi_y=m_ift_local_residual_y->at(ihit);
      resi_e_y=m_ift_local_ye->at(ihit);
      mille_file.mille(loc_der.size(), lcder, glo_der.size(), glder, label, resi, resi_e);
      mille_file.mille(loc_der_y.size(), lcder_y, glo_der_y.size(), glder_y, label_y, resi_y, resi_e_y);
    }



    for(int ihit = 0; ihit<m_fitParam_align_id->size();++ihit){
      labels_og.clear();
      glo_der_og.clear();
      loc_der_og.clear();
      if(fabs(m_fitParam_align_local_residual_x->at(ihit))>0.05)continue;
      if(m_fitParam_align_local_derivative_x_x->at(ihit)<-9000||m_fitParam_align_local_derivative_x_rz->at(ihit)<-9000||m_fitParam_align_global_derivative_y_x->at(ihit)<-9000||m_fitParam_align_global_derivative_y_y->at(ihit)<-9000||m_fitParam_align_global_derivative_y_z->at(ihit)<-9000||m_fitParam_align_global_derivative_y_rx->at(ihit)<-9000||m_fitParam_align_global_derivative_y_ry->at(ihit)<-9000||m_fitParam_align_global_derivative_y_rz->at(ihit)<-9000)continue;

      int moduleid = m_fitParam_align_id->at(ihit);
      //      std::cout<<"like moduleid "<<moduleid;
      if(moduleid%10==1)--moduleid;
      //      std::cout<<" "<<moduleid<<std::endl;
      moduleid+=1000;//station from 1 not 0
      //	int layerid = moduleid/100;
      //	if(moduleid/100==10){
      //	  if(nhits[moduleid%100/10]>1000)continue;
      //	  nhits[moduleid%100/10]+=1;
      //	}
  //     if(dump6ndf_modules){
	// labels_og.push_back(moduleid*10+0+1);//millepede can not have label at 0
	// labels_og.push_back(moduleid*10+1+1);
	// labels_og.push_back(moduleid*10+2+1);
	// labels_og.push_back(moduleid*10+3+1);
	// labels_og.push_back(moduleid*10+4+1);
	// if(dumpz_modules){
	//   labels_og.push_back(moduleid*10+5+1);
	//   glo_der_og.push_back(m_fitParam_align_local_derivative_x_z->at(ihit));
	// }
	// glo_der_og.push_back(m_fitParam_align_local_derivative_x_x->at(ihit));
	// glo_der_og.push_back(m_fitParam_align_local_derivative_x_y->at(ihit));
	// glo_der_og.push_back(m_fitParam_align_local_derivative_x_rx->at(ihit));
	// glo_der_og.push_back(m_fitParam_align_local_derivative_x_ry->at(ihit));
	// glo_der_og.push_back(m_fitParam_align_local_derivative_x_rz->at(ihit));
  //     }
  //     else if(dump3ndf_modules){
	// labels_og.push_back(moduleid*10+0+1);//millepede can not have label at 0
	// labels_og.push_back(moduleid*10+1+1);
	// labels_og.push_back(moduleid*10+2+1);
	// glo_der_og.push_back(m_fitParam_align_local_derivative_x_x->at(ihit));
	// glo_der_og.push_back(m_fitParam_align_local_derivative_x_y->at(ihit));
	// glo_der_og.push_back(m_fitParam_align_local_derivative_x_rz->at(ihit));
  //     }
  //     else{
	// labels_og.push_back(moduleid*10+0+1);//millepede can not have label at 0
	// labels_og.push_back(moduleid*10+1+1);
	// glo_der_og.push_back(m_fitParam_align_local_derivative_x_x->at(ihit));
	// glo_der_og.push_back(m_fitParam_align_local_derivative_x_rz->at(ihit));         
  //     }
      if(dumplayers){
	int layerid = moduleid/100;
	//	std::cout<<"id "<<moduleid<<" "<<layerid<<std::endl;
	if(dump6ndf_layers){
	  if(fabs(m_fitParam_align_global_derivative_y_rx->at(ihit))>2||fabs(m_fitParam_align_global_derivative_y_ry->at(ihit))>2)continue;
	  labels_og.push_back(layerid*10+0+1);//millepede can not have label at 0
	  labels_og.push_back(layerid*10+1+1);
	  labels_og.push_back(layerid*10+2+1);
	  labels_og.push_back(layerid*10+3+1);
	  labels_og.push_back(layerid*10+4+1);
	  glo_der_og.push_back(m_fitParam_align_global_derivative_y_x->at(ihit));
	  glo_der_og.push_back(m_fitParam_align_global_derivative_y_y->at(ihit));
	  if(dumpz_layers)
	    {
	      labels_og.push_back(layerid*10+5+1);
	      glo_der_og.push_back(m_fitParam_align_global_derivative_y_z->at(ihit));
	    }
	  glo_der_og.push_back(m_fitParam_align_global_derivative_y_rx->at(ihit));
	  glo_der_og.push_back(m_fitParam_align_global_derivative_y_ry->at(ihit));
	  glo_der_og.push_back(m_fitParam_align_global_derivative_y_rz->at(ihit));
	} 
	else{
	  labels_og.push_back(layerid*10+0+1);//millepede can not have label at 0
	  labels_og.push_back(layerid*10+1+1);
	  glo_der_og.push_back(m_fitParam_align_global_derivative_y_y->at(ihit));
	  glo_der_og.push_back(m_fitParam_align_global_derivative_y_rz->at(ihit));
	}
      }
      if(dumpstations){
	int stationid = moduleid/1000;
	if(dump5ndf_stations){
	  labels_og.push_back(stationid*10+0+1);//millepede can not have label at 0
	  labels_og.push_back(stationid*10+1+1);
	  labels_og.push_back(stationid*10+2+1);
	  labels_og.push_back(stationid*10+3+1);
	  labels_og.push_back(stationid*10+4+1);
	  glo_der_og.push_back(m_fitParam_align_global_derivative_y_x->at(ihit));
	  glo_der_og.push_back(m_fitParam_align_global_derivative_y_y->at(ihit));
    if(dumpz_station){
      labels_og.push_back(stationid*10+5+1);
      glo_der_og.push_back(m_fitParam_align_global_derivative_y_z->at(ihit));
	  }
	  glo_der_og.push_back(m_fitParam_align_global_derivative_y_rx->at(ihit));
	  glo_der_og.push_back(m_fitParam_align_global_derivative_y_ry->at(ihit));
	  glo_der_og.push_back(m_fitParam_align_global_derivative_y_rz->at(ihit));
	}
      }

      loc_der_og.push_back(m_fitParam_align_local_derivative_x_par_x->at(ihit));
      loc_der_og.push_back(m_fitParam_align_local_derivative_x_par_y->at(ihit));
      loc_der_og.push_back(m_fitParam_align_local_derivative_x_par_theta->at(ihit));
      loc_der_og.push_back(m_fitParam_align_local_derivative_x_par_phi->at(ihit));
      loc_der_og.push_back(m_fitParam_align_local_derivative_x_par_qop->at(ihit));
      //      std::cout<<"like "<<ievt<<" "<<ihit<<" "<<loc_der_og.size()<<std::endl;
      float* lcder = loc_der_og.data();
      float* glder = glo_der_og.data();
      int* label = labels_og.data();
      resi_og=m_fitParam_align_local_residual_x->at(ihit);
      resi_e_og=m_fitParam_align_local_measured_xe->at(ihit);
      mille_file.mille(loc_der_og.size(), lcder, glo_der_og.size(), glder, label, resi_og, resi_e_og);
    }


    mille_file.end();
    //    mille_file.flushTrack();
  }
  // } // file loop
  mille_file.kill();


}
