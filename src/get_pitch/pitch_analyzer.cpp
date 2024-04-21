/// @file

#include <iostream>
#include <math.h>
#include "pitch_analyzer.h"

using namespace std;

/// Name space of UPC
namespace upc {
  void PitchAnalyzer::autocorrelation(const vector<float> &x, vector<float> &r) const {

    for (unsigned int l = 0; l < r.size(); ++l) {
  		/// \TODO Compute the autocorrelation r[l] 
      /// \FET Autocorrelació ***computada***  Per afegir missatges al make doc
      
      r[l] = 0;
      // Formula Autocorrelación -> 1/N sum(x[n]*x[n-l])
      for (unsigned int n = l; n < x.size(); n++){
        r[l] += x[n]*x[n-l];
      }      
      r[l] /= x.size();  
    }

    if (r[0] == 0.0F) //to avoid log() and divide zero 
      r[0] = 1e-10; 
  }

  void PitchAnalyzer::set_window(Window win_type) {
    if (frameLen == 0)
      return;

    window.resize(frameLen);

    switch (win_type) {
    case HAMMING:
      /// \TODO Implement the Hamming window
      window.resize(frameLen);
      for (int i = 0; i < frameLen; ++i) {
        window[i] = 0.54 - 0.46 * cos(2 * M_PI * i / (frameLen - 1));
      }
      /// \FET Finestra de Hamming implementada. Utilitzarem la rectangular ja que funciona millor
      break;
    case RECT:
        window.assign(frameLen, 1);
    default:
      window.assign(frameLen, 1);
    }
  }

  void PitchAnalyzer::set_f0_range(float min_F0, float max_F0) {
    npitch_min = (unsigned int) samplingFreq/max_F0;
    if (npitch_min < 2)
      npitch_min = 2;  // samplingFreq/2

    npitch_max = 1 + (unsigned int) samplingFreq/min_F0;

    //frameLen should include at least 2*T0
    if (npitch_max > frameLen/2)
      npitch_max = frameLen/2;
  }

  bool PitchAnalyzer::unvoiced(float pot, float r1norm, float rmaxnorm,float ZCR) const {
    /// \TODO Implement a rule to decide whether the sound is voiced or not.
    /// * You can use the standard features (pot, r1norm, rmaxnorm),
    ///   or compute and use other ones.

    
float score = 0;
    static float power_first_window = 0;
    static int window = 0;
    const float potvalue = -46, r1value = 0.5, rmaxvalue = 0.41, zcrvalue= 0.1;

    if (pot < potvalue)
      score += 0.5;
    else if (r1norm < r1value)
      score += 0.5;
    else if (rmaxnorm < rmaxvalue)
      score += 0.5;
    
    if (ZCR > zcrvalue)
      score += 0.5;

    if (score >= 1)
      return true;
    else
      return false;
  }

//calul del ZCR
  float PitchAnalyzer::compute_zcr(const vector<float> &x) const {
    float suma=0;
    unsigned int N= x.size();
    
     for(int i=1; i<N; i++){
        if((x[i-1]>=0 && x[i]<=0)||(x[i-1]<=0 && x[i]>=0)){
        suma=suma+1;
        }
        
    }
    return (float) (suma)/(2*(N));
}

  float PitchAnalyzer::compute_pitch(vector<float> & x) const {
    if (x.size() != frameLen)
      return -1.0F;

    //Window input frame
    for (unsigned int i=0; i<x.size(); ++i)
      x[i] *= window[i];

    vector<float> r(npitch_max);

    //Compute correlation
    autocorrelation(x, r);
    //Compute ZCR
    float ZCR= compute_zcr(x);

    vector<float>::const_iterator iR = r.begin(), iRMax = iR;

    /// \TODO 
	/// Find the lag of the maximum value of the autocorrelation away from the origin.<br>
	/// Choices to set the minimum value of the lag are:
	///    - The first negative value of the autocorrelation.
	///    - The lag corresponding to the maximum value of the pitch.
    ///	   .
	/// In either case, the lag should not exceed that of the minimum value of the pitch.

    for (iRMax = iR = r.begin() + npitch_min; iR < r.begin() + npitch_max; iR++){
      if(*iR > *iRMax){
        iRMax = iR;
      }
    }

    unsigned int lag = iRMax - r.begin();

    float pot = 10 * log10(r[0]);

    //You can print these (and other) features, look at them using wavesurfer
    //Based on that, implement a rule for unvoiced
    //change to #if 1 and compile
#if 1
    if (r[0] > 0.0F)
      cout << roundf(pot* 100) / 100.0  << '\t' << roundf(r[1]/r[0]* 100) / 100.0 << '\t' << roundf(r[lag]/r[0]  * 100) / 100.0<< endl;
#endif 
    
    if (unvoiced(pot, r[1]/r[0], r[lag]/r[0],ZCR))
      return 0;
    else
      return (float) samplingFreq/(float) lag;
  }
}
