/* 	autotuned~.c

	An auto-tuning PD External, based on autotalent, an auto-tuning LADSPA plugin.

	Free software by Thomas A. Baran.
	http://web.mit.edu/tbaran/www/autotalent.html
	VERSION 0.2
	March 20, 2010
   
	Pd port by megrimm 2011

	This program is free software; you can redistribute it and/or modify        
	it under the terms of the GNU General Public License as published by        
	the Free Software Foundation; either version 2 of the License, or           
	(at your option) any later version.                                         
                                                                                
	This program is distributed in the hope that it will be useful,             
	but WITHOUT ANY WARRANTY; without even the implied warranty of              
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               
	GNU General Public License for more details.                                
                                                                                
	You should have received a copy of the GNU General Public License           
	along with this program; if not, write to the Free Software                 
	Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  

 */

/*****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#define PI (float)3.14159265358979323846
#define L2SC (float)3.32192809488736218171

/*****************************************************************************/

#include "m_pd.h"    // Required for all PD external objects
//#include "mayer_fft.h"

#ifndef CLIP
#define CLIP(a, lo, hi) ( (a)>(lo)?( (a)<(hi)?(a):(hi) ):(lo) )
#endif

#ifndef MAYER_H
#define MAYER_H

#define REAL float

void mayer_realfft(int n, REAL *real);
void mayer_realifft(int n, REAL *real);

#endif

// Variables for FFT routine
typedef struct
{
  int nfft;        // size of FFT
  int numfreqs;    // number of frequencies represented (nfft/2 + 1)
  float* fft_data; // array for writing/reading to/from FFT function
} fft_vars;

// Constructor for FFT routine
fft_vars* fft_con(int nfft)
{
  fft_vars* membvars = (fft_vars*) malloc(sizeof(fft_vars));

  membvars->nfft = nfft;
  membvars->numfreqs = nfft/2 + 1;

  membvars->fft_data = (float*) calloc(nfft, sizeof(float));

  return membvars;
}

// Destructor for FFT routine
void fft_des(fft_vars* membvars)
{
  free(membvars->fft_data);

  free(membvars);
}

// Perform forward FFT of real data
// Accepts:
//   membvars - pointer to struct of FFT variables
//   input - pointer to an array of (real) input values, size nfft
//   output_re - pointer to an array of the real part of the output,
//     size nfft/2 + 1
//   output_im - pointer to an array of the imaginary part of the output,
//     size nfft/2 + 1
void fft_forward(fft_vars* membvars, float* input, float* output_re, float* output_im)
{
  int ti;
  int nfft;
  int hnfft;
  int numfreqs;

  nfft = membvars->nfft;
  hnfft = nfft/2;
  numfreqs = membvars->numfreqs;

  for (ti=0; ti<nfft; ti++) {
    membvars->fft_data[ti] = input[ti];
  }

  mayer_realfft(nfft, membvars->fft_data);

  output_im[0] = 0;
  for (ti=0; ti<hnfft; ti++) {
    output_re[ti] = membvars->fft_data[ti];
    output_im[ti+1] = membvars->fft_data[nfft-1-ti];
  }
  output_re[hnfft] = membvars->fft_data[hnfft];
  output_im[hnfft] = 0;
}

// Perform inverse FFT, returning real data
// Accepts:
//   membvars - pointer to struct of FFT variables
//   input_re - pointer to an array of the real part of the output,
//     size nfft/2 + 1
//   input_im - pointer to an array of the imaginary part of the output,
//     size nfft/2 + 1
//   output - pointer to an array of (real) input values, size nfft
void fft_inverse(fft_vars* membvars, float* input_re, float* input_im, float* output)
{
  int ti;
  int nfft;
  int hnfft;
  int numfreqs;

  nfft = membvars->nfft;
  hnfft = nfft/2;
  numfreqs = membvars->numfreqs;

  for (ti=0; ti<hnfft; ti++) {
    membvars->fft_data[ti] = input_re[ti];
    membvars->fft_data[nfft-1-ti] = input_im[ti+1];
  }
  membvars->fft_data[hnfft] = input_re[hnfft];

  mayer_realifft(nfft, membvars->fft_data);

  for (ti=0; ti<nfft; ti++) {
    output[ti] = membvars->fft_data[ti];
  }
}

// DONE WITH FFT CODE



void *autotuned_class;

typedef struct _autotuned 	// Data structure for this object PD
{
    t_object x_obj;
    float x_f;
		
	// parameters
	
	float fMix;
	float fShift;
	float fTune;
	float fA;
	float fBb;
	float fB;
	float fC;
	float fDb;
	float fD;
	float fEb;
	float fE;
	float fF;
	float fGb;
	float fG;
	float fAb;
	float fAmount;
	float fGlide;

	fft_vars* fmembvars; 		// member variables for fft routine
	unsigned long fs; 			// Sample rate
	unsigned long cbsize; 		// size of circular buffer
	unsigned long corrsize; 	// cbsize/2 + 1
	unsigned long cbiwr;
	unsigned long cbord;
	float* cbi; 				// circular input buffer
//	float* cbf; 				// circular formant correction buffer
	float* cbo; 				// circular output buffer
	float* cbonorm;				// circular output buffer used to normalize signal
	float* cbwindow; 			// hann of length N/2, zeros for the rest
	float* acwinv; 				// inverse of autocorrelation of window
	float* hannwindow; 			// length-N hann
	int noverlap;
	float* ffttime;
	float* fftfreqre;
	float* fftfreqim;
	
	// VARIABLES FOR LOW-RATE SECTION
	float aref; 				// A tuning reference (Hz)
	float pperiod; 				// Pitch period (seconds)
	float pitch; 				// Pitch (semitones)
	float pitchpers; 			// Pitch persist (semitones)
	float conf; 				// Confidence of pitch period estimate (between 0 and 1)
	float vthresh; 				// Voiced speech threshold
	float pmax; 				// Maximum allowable pitch period (seconds)
	float pmin; 				// Minimum allowable pitch period (seconds)
	unsigned long nmax; 		// Maximum period index for pitch prd est
	unsigned long nmin; 		// Minimum period index for pitch prd est
	float lrshift; 				// Shift prescribed by low-rate section
	int ptarget; 				// Pitch target, between 0 and 11
	float sptarget; 			// Smoothed pitch target
	int wasvoiced; 				// 1 if previous frame was voiced
	float persistamt; 			// Proportion of previous pitch considered during next voiced period
	float glidepersist;
	
	// VARIABLES FOR PITCH SHIFTER
	float phprd; 				// phase period
	float phprdd;				// default (unvoiced) phase period
	float phinc; 				// input phase increment
	float phincfact; 			// factor determining output phase increment
	float phasein;
	float phaseout;
	float* frag; 				// windowed fragment of speech
	unsigned long fragsize; 	// size of fragment in samples	
	float clockinterval;
	void *periodout, *confout; 	// floatout for pitch
	void *clock;
	
} t_autotuned;	



//prototypes for methods
void *autotuned_new();
void autotuned_free(t_autotuned *x);	
void autotuned_init(t_autotuned *x, unsigned long sr);
t_int *autotuned_perform(t_int *w);
void autotuned_dsp(t_autotuned *x, t_signal **sp, short *count);
void autotuned_assist(t_autotuned *x, void *b, long m, long a, char *s);
//void autotuned_setparam(t_autotuned *x, t_symbol *m, short argc, t_atom *argv);
void autotuned_list(t_autotuned *x, t_symbol *m, short argc, t_atom *argv);
void autotuned_mix(t_autotuned *x, t_floatarg f);
void autotuned_shift(t_autotuned *x, t_floatarg f);
void autotuned_tune(t_autotuned *x, t_floatarg f);
void autotuned_correction(t_autotuned *x, t_floatarg f);
void autotuned_glide(t_autotuned *x, t_floatarg f);
void autotuned_processclock(t_autotuned *x);




void autotuned_tilde_setup(void)
{
	autotuned_class = class_new(gensym("autotuned~"), (t_newmethod)autotuned_new, 
							   (t_method)autotuned_free ,sizeof(t_autotuned), 0,A_GIMME,0);
	CLASS_MAINSIGNALIN(autotuned_class, t_autotuned, x_f );
	class_addmethod(autotuned_class,(t_method)autotuned_dsp, gensym("dsp"), 0);
	class_addmethod(autotuned_class,(t_method)autotuned_assist, gensym("assist"), 0);
	class_addmethod(autotuned_class,(t_method)autotuned_list,gensym("list"),A_GIMME,0);
	class_addmethod(autotuned_class,(t_method)autotuned_mix,gensym("mix"),A_FLOAT,0);
	class_addmethod(autotuned_class,(t_method)autotuned_shift,gensym("shift"),A_FLOAT,0);
	class_addmethod(autotuned_class,(t_method)autotuned_tune,gensym("tune"),A_FLOAT,0);
	class_addmethod(autotuned_class,(t_method)autotuned_correction,gensym("correction"),A_FLOAT,0);
	class_addmethod(autotuned_class,(t_method)autotuned_glide,gensym("glide"),A_FLOAT,0);
	post("autotuned~ v.0.2");
	post("megrimm 2011");
}



// Create - Contruction of signal inlets and outlets
void *autotuned_new()
{
    unsigned long sr;
    t_autotuned *x = (t_autotuned *)pd_new(autotuned_class);
	
	if(sys_getsr()) sr = sys_getsr();
	else sr = 44100;
	
	autotuned_init(x , sr);
	
	// second argument = number of signal inlets
    //dsp_setup((t_object *)x, 1);
	//inlet_new(&x->x_obj, &x->x_obj.ob_pd,gensym("signal"));
	//inlet_new(&x->x_obj, &x->x_obj.ob_pd,gensym("signal"), gensym("signal"));
	floatinlet_new (&x->x_obj, &x->fAmount);
	floatinlet_new (&x->x_obj, &x->fGlide);
	floatinlet_new (&x->x_obj, &x->fMix);
	floatinlet_new (&x->x_obj, &x->fShift);
	floatinlet_new (&x->x_obj, &x->fTune);
	//floatinlet_new (&x->x_obj, &x->fPersist);
	//symbolinlet_new (&x->x_obj, &x->fAmount);
	
	outlet_new(&x->x_obj, gensym("signal"));
	//x->f_out = outlet_new(&x->x_obj, &s_float);
	x->confout = outlet_new(&x->x_obj, &s_float);
	x->periodout = outlet_new(&x->x_obj, &s_float);
	//x->confout = outlet_new(x, "float");
	//x->periodout = outlet_new(x, "float");
    x->clock = clock_new(x,(t_method)autotuned_processclock);
	x->clockinterval = 10.;
    
    return (x);									// Return the pointer
}

// Destroy
void autotuned_free(t_autotuned *x)
{
	//dsp_free((t_object *)x);		// Always call dsp_free first in this routine
	clock_unset(x->clock);
	//object_free(x->clock,0); 
	fft_des(x->fmembvars);
	freebytes(x->cbi,0);
	freebytes(x->cbo,0);
	freebytes(x->cbonorm,0);
	freebytes(x->cbwindow,0);
	freebytes(x->hannwindow,0);
	freebytes(x->acwinv,0);
	freebytes(x->frag,0);
	freebytes(x->ffttime,0);
	freebytes(x->fftfreqre,0);
	freebytes(x->fftfreqim,0);

}

void autotuned_init(t_autotuned *x,unsigned long sr)
{
	unsigned long ti;

	x->fs = sr;
	x->aref = 440;
	
	if (x->fs >=88200) {
		x->cbsize = 4096;
	}
	else {
		x->cbsize = 2048;
	}
	x->corrsize = x->cbsize / 2 + 1;
	
	x->pmax = 1/(float)70;  // max and min periods (ms)
	x->pmin = 1/(float)700; // eventually may want to bring these out as sliders
	
	x->pperiod = x->pmax;
	
	x->nmax = (unsigned long)(x->fs * x->pmax);
	if (x->nmax > x->corrsize) {
		x->nmax = x->corrsize;
	}
	x->nmin = (unsigned long)(x->fs * x->pmin);
	
	x->cbi = (float*) calloc(x->cbsize, sizeof(float));
	x->cbo = (float*) calloc(x->cbsize, sizeof(float));
	x->cbonorm = (float*) calloc(x->cbsize, sizeof(float));
	
	x->cbiwr = 0;
	x->cbord = 0;
	
	// Standard raised cosine window, max height at N/2
	x->hannwindow = (float*) calloc(x->cbsize, sizeof(float));
	for (ti=0; ti<x->cbsize; ti++) {
		x->hannwindow[ti] = -0.5*cos(2*PI*ti/(x->cbsize - 1)) + 0.5;
	}
	
	// Generate a window with a single raised cosine from N/4 to 3N/4
	x->cbwindow = (float*) calloc(x->cbsize, sizeof(float));
	for (ti=0; ti<(x->cbsize / 2); ti++) {
		x->cbwindow[ti+x->cbsize/4] = -0.5*cos(4*PI*ti/(x->cbsize - 1)) + 0.5;
	}
	
	x->noverlap = 4;
	
	x->fmembvars = fft_con(x->cbsize);
	
	x->ffttime = (float*) calloc(x->cbsize, sizeof(float));
	x->fftfreqre = (float*) calloc(x->corrsize, sizeof(float));
	x->fftfreqim = (float*) calloc(x->corrsize, sizeof(float));
	
	
	// ---- Calculate autocorrelation of window ----
	x->acwinv = (float*) calloc(x->cbsize, sizeof(float));
	for (ti=0; ti<x->cbsize; ti++) {
		x->ffttime[ti] = x->cbwindow[ti];
	}
	fft_forward(x->fmembvars, x->cbwindow, x->fftfreqre, x->fftfreqim);
	for (ti=0; ti<x->corrsize; ti++) {
		x->fftfreqre[ti] = (x->fftfreqre[ti])*(x->fftfreqre[ti]) + (x->fftfreqim[ti])*(x->fftfreqim[ti]);
		x->fftfreqim[ti] = 0;
	}
	fft_inverse(x->fmembvars, x->fftfreqre, x->fftfreqim, x->ffttime);
	for (ti=1; ti<x->cbsize; ti++) {
		x->acwinv[ti] = x->ffttime[ti]/x->ffttime[0];
		if (x->acwinv[ti] > 0.000001) {
			x->acwinv[ti] = (float)1/x->acwinv[ti];
		}
		else {
			x->acwinv[ti] = 0;
		}
	}
	x->acwinv[0] = 1;
	// ---- END Calculate autocorrelation of window ----	
	
	x->lrshift = 0;
	x->ptarget = 0;
	x->sptarget = 0;
	x->wasvoiced = 0;
	x->persistamt = 0;
	
	x->glidepersist = 100; // 100 ms glide persist
	
	x->vthresh = 0.8;  //  The voiced confidence (unbiased peak) threshold level
	
	// Pitch shifter initialization
	x->phprdd = 0.01; // Default period
	x->phprd = x->phprdd;
	x->phinc = (float)1/(x->phprd * x->fs);
	x->phincfact = 1;
	x->phasein = 0;
	x->phaseout = 0;
	x->frag = (float*) calloc(x->cbsize, sizeof(float));
	x->fragsize = 0;
}

//*********************//
// DSP Methods PD //
//*********************//

void autotuned_dsp(t_autotuned *x, t_signal **sp, short *count)
{
	clock_delay(x->clock, 0.);
	
	if(x->fs != sp[0]->s_sr)  autotuned_init(x, sp[0]->s_sr);
	
	dsp_add(autotuned_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}


//*************************//
// Perform Routine PD//
//*************************//
t_int *autotuned_perform(t_int *w)
{
	t_autotuned *x = (t_autotuned *)(w[1]); // object is first arg 
	t_float *in = (t_float *)(w[2]);
	t_float *out = (t_float *)(w[3]);
	unsigned long SampleCount = (unsigned long)(w[4]);
	
	// copy struct variables to local
	
	float fMix = x->fMix;
	float fShift = x->fShift;
	float fTune = x->fTune;
	float fA = x->fA;
	float fBb = x->fBb;
	float fB = x->fB;
	float fC = x->fC;
	float fDb = x->fDb;
	float fD = x->fD;
	float fEb = x->fEb;
	float fE = x->fE;
	float fF = x->fF;
	float fGb = x->fGb;
	float fG = x->fG;
	float fAb = x->fAb;
	float fGlide = x->fGlide;
	float fPersist = x->glidepersist;
	float fAmount = x->fAmount;
	
	x->aref = (float)440*pow(2,fTune/12);
	
	unsigned long int lSampleIndex;
	
	unsigned long N = x->cbsize;
	unsigned long Nf = x->corrsize;
	unsigned long fs = x->fs;

	float pmax = x->pmax;
	float pmin = x->pmin;
	unsigned long nmax = x->nmax;
	unsigned long nmin = x->nmin;
	
	float pperiod = x->pperiod;
	float pitch = x->pitch;
	float conf = x->conf;
	float aref = x->aref;
	
		//
	
	long int ti;
	long int ti2;
	long int ti3;
	float tf;
	float tf2;
	float tf3;
	
	// Variables for cubic spline interpolator
	float indd;
  int ind0;
  int ind1;
  int ind2;
  int ind3;
  float vald;
  float val0;
  float val1;
  float val2;
  float val3;
	
	
	//******************//
	//  MAIN DSP LOOP   //
	//******************//
	for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++)  
	{
		
		// load data into circular buffer
		tf = (float) *(in++);
		x->cbi[x->cbiwr] = tf;
		x->cbiwr++;
		if (x->cbiwr >= N) {
			x->cbiwr = 0;
		}
		
		
		// ********************//
		// * Low-rate section *//
		// ********************//
		
		// Every N/noverlap samples, run pitch estimation / correction code
		if ((x->cbiwr)%(N/x->noverlap) == 0) {
			
			
			// ---- Obtain autocovariance ---- //
			
			// Window and fill FFT buffer
			ti2 = (long) x->cbiwr;
			for (ti=0; ti<(long)N; ti++) {
				x->ffttime[ti] = (float)(x->cbi[(ti2-ti)%N]*x->cbwindow[ti]);
			}
			
			// Calculate FFT
			fft_forward(x->fmembvars, x->ffttime, x->fftfreqre, x->fftfreqim);
			
			// Remove DC
			x->fftfreqre[0] = 0;
			x->fftfreqim[0] = 0;
			
			// Take magnitude squared
			for (ti=1; ti< (long) Nf; ti++) {
				x->fftfreqre[ti] = (x->fftfreqre[ti])*(x->fftfreqre[ti]) + (x->fftfreqim[ti])*(x->fftfreqim[ti]);
				x->fftfreqim[ti] = 0;
			}
			
			// Calculate IFFT
			fft_inverse(x->fmembvars, x->fftfreqre, x->fftfreqim, x->ffttime);
			
			// Normalize
			for (ti=1; ti<(long)N; ti++) {
				x->ffttime[ti] = x->ffttime[ti] / x->ffttime[0];
			}
			x->ffttime[0] = 1;
			
			//  ---- END Obtain autocovariance ----
			
			
			//  ---- Calculate pitch and confidence ----
			
			// Calculate pitch period
			//   Pitch period is determined by the location of the max (biased)
			//   peak within a given range
			//   Confidence is determined by the corresponding unbiased height
			tf2 = 0;
			pperiod = pmin;
			for (ti=nmin; ti<(long)nmax; ti++) {
				ti2 = ti-1;
				ti3 = ti+1;
				if (ti2<0) {
					ti2 = 0;
				}
				if (ti3>(long)Nf) {
					ti3 = Nf;
				}
				tf = x->ffttime[ti];
				
				if (tf>x->ffttime[ti2] && tf>=x->ffttime[ti3] && tf>tf2) {
					tf2 = tf;
					conf = tf*x->acwinv[ti];
					pperiod = (float)ti/fs;
				}
			}
			
			// Convert to semitones
			pitch = (float) -12*log10((float)aref*pperiod)*L2SC;
			x->pitch = pitch;
			x->pperiod = pperiod;
			x->conf = conf;
			
			//  ---- END Calculate pitch and confidence ----
			
			
			//  ---- Determine pitch target ----
			
			// If voiced
			if (conf>=x->vthresh) {
				// TODO: Scale sliders
				// Determine pitch target
				tf = -1;
				tf2 = 0;
				tf3 = 0;
				for (ti=0; ti<12; ti++) {
					switch (ti) {
						case 0:
							tf2 = fA;
							break;
						case 1:
							tf2 = fBb;
							break;
						case 2:
							tf2 = fB;
							break;
						case 3:
							tf2 = fC;
							break;
						case 4:
							tf2 = fDb;
							break;
						case 5:
							tf2 = fD;
							break;
						case 6:
							tf2 = fEb;
							break;
						case 7:
							tf2 = fE;
							break;
						case 8:
							tf2 = fF;
							break;
						case 9:
							tf2 = fGb;
							break;
						case 10:
							tf2 = fG;
							break;
						case 11:
							tf2 = fAb;
							break;
					}
					/* 	  if (ti==x->ptarget) { */
					/* 	    tf2 = tf2 + 0.01; // add a little hysteresis */
					/* 	  } */
					tf2 = tf2 - (float)fabs( (pitch-(float)ti)/6 - 2*floorf(((pitch-(float)ti)/12 + 0.5)) ); // like a likelihood function
					if (tf2>=tf) {                                                                           // that we're maximizing
						tf3 = (float)ti;                                                                       // to find the target pitch
						tf = tf2;
					}
				}
				x->ptarget = tf3;
				
				// Glide persist
				if (x->wasvoiced == 0) {
					x->wasvoiced = 1;
					tf = x->persistamt;
					x->sptarget = (1-tf)*x->ptarget + tf*x->sptarget;
					x->persistamt = 1;
				}
				
				// Glide on circular scale
				tf3 = (float)x->ptarget - x->sptarget;
				tf3 = tf3 - (float)12*floorf(tf3/12 + 0.5);
				if (fGlide>0) {
					tf2 = (float)1-pow((float)1/24, (float)N * 1000/ (x->noverlap*fs*fGlide));
				}
				else {
					tf2 = 1;
				}
				x->sptarget = x->sptarget + tf3*tf2;
			}
			// If not voiced
			else {
				x->wasvoiced = 0;
				
				// Keep track of persist amount
				if (fPersist>0) {
					tf = pow((float)1/2, (float)N * 1000/ (x->noverlap*fs*fPersist));
				}
				else {
					tf = 0;
				}
				x->persistamt = x->persistamt * tf; // Persist amount decays exponentially
			}
			// END If voiced
			
			//  ---- END Determine pitch target ----
			
			
			// ---- Determine correction to feed to the pitch shifter ----
			tf = x->sptarget - pitch; // Correction amount
			tf = tf - (float)12*floorf(tf/12 + 0.5); // Never do more than +- 6 semitones of correction
			if (conf<x->vthresh) {
				tf = 0;
			}
			x->lrshift = fShift + fAmount*tf;  // Add in pitch shift slider
			
			
			// ---- Compute variables for pitch shifter that depend on pitch ---
			x->phincfact = (float)pow(2, x->lrshift/12);
			if (conf>=x->vthresh) {  // Keep old period when unvoiced
				x->phinc = (float)1/(pperiod*fs);
				x->phprd = pperiod*2;
			}
		}
		// ************************
		// * END Low-Rate Section *
		// ************************
		
		
		// *****************
		// * Pitch Shifter *
		// *****************
		
    	// Pitch shifter (kind of like a pitch-synchronous version of Fairbanks' technique)
    	//   Note: pitch estimate is naturally N/2 samples old
		x->phasein = x->phasein + x->phinc;
		x->phaseout = x->phaseout + x->phinc*x->phincfact;
		
		//   When input phase resets, take a snippet from N/2 samples in the past
		if (x->phasein >= 1) {
			x->phasein = x->phasein - 1;
			ti2 = x->cbiwr - (long int)N/2;
			for (ti=-((long int)N)/2; ti<(long int)N/2; ti++) {
				x->frag[ti%N] = x->cbi[(ti + ti2)%N];
			}
		}
		
		//   When output phase resets, put a snippet N/2 samples in the future
		if (x->phaseout >= 1) {
			x->fragsize = x->fragsize*2;
			if (x->fragsize >= N) {
				x->fragsize = N;
			}
			x->phaseout = x->phaseout - 1;
			ti2 = x->cbord + N/2;
			ti3 = (long int)(((float)x->fragsize) / x->phincfact);
			for (ti=-ti3/2; ti<(ti3/2); ti++) {
				tf = x->hannwindow[(long int)N/2 + ti*(long int)N/ti3];
				x->cbo[(ti + ti2)%N] = x->cbo[(ti + ti2)%N] + x->frag[((int)(x->phincfact*ti))%N]*tf;
				x->cbonorm[(ti + ti2)%N] = x->cbonorm[(ti + ti2)%N] + tf;
			}
			x->fragsize = 0;
		}
		x->fragsize++;
		
		//   Get output signal from buffer
		tf = x->cbonorm[x->cbord];
		//   Normalize
		if (tf>0.5) {
			tf = (float)1/tf;
		}
		else {
			tf = 1;
		}
		tf = tf*x->cbo[x->cbord]; // read buffer
		tf = x->cbo[x->cbord];
		x->cbo[x->cbord] = 0; // erase for next cycle
		x->cbonorm[x->cbord] = 0;
		x->cbord++; // increment read pointer
		if (x->cbord >= N) {
			x->cbord = 0;
		}
		

/*		
		    if (x->phaseout >= 1) {
      x->fragsize = x->fragsize*2;
      if (x->fragsize > N) {
	x->fragsize = N;
      }
      x->phaseout = x->phaseout - 1;
      ti2 = x->cbord + N/2;
      ti3 = (long int)(((float)x->fragsize) / x->phincfact);
      if (ti3>=N/2) {
	ti3 = N/2 - 1;
      }
      for (ti=-ti3/2; ti<(ti3/2); ti++) {
	tf = x->hannwindow[(long int)N/2 + ti*(long int)N/ti3];
	// 3rd degree polynomial interpolator - based on eqns from Hal Chamberlin's book
	indd = x->phincfact*ti;
	ind1 = (int)indd;
	ind2 = ind1+1;
	ind3 = ind1+2;
	ind0 = ind1-1;
	val0 = x->frag[(ind0+N)%N];
	val1 = x->frag[(ind1+N)%N];
	val2 = x->frag[(ind2+N)%N];
	val3 = x->frag[(ind3+N)%N];
	vald = 0;
	vald = vald - (float)0.166666666667 * val0 * (indd - ind1) * (indd - ind2) * (indd - ind3);
	vald = vald + (float)0.5 * val1 * (indd - ind0) * (indd - ind2) * (indd - ind3);
	vald = vald - (float)0.5 * val2 * (indd - ind0) * (indd - ind1) * (indd - ind3);
	vald = vald + (float)0.166666666667 * val3 * (indd - ind0) * (indd - ind1) * (indd - ind2);
	x->cbo[(ti + ti2 + N)%N] = x->cbo[(ti + ti2 + N)%N] + vald*tf;
      }
      x->fragsize = 0;
    }
    x->fragsize++;

    //   Get output signal from buffer
    tf = x->cbo[x->cbord]; // read buffer

    x->cbo[x->cbord] = 0; // erase for next cycle
    x->cbord++; // increment read pointer
    if (x->cbord >= N) {
      x->cbord = 0;
    }
*/
		
		// *********************
		// * END Pitch Shifter *
		// *********************
		
		
		// Write audio to output
		// Mix (blend between original (delayed) =0 and shifted/corrected =1)
		*(out++) = fMix*tf + (1-fMix)*x->cbi[(x->cbiwr - N + 1)%N];
		
	}

    return (w + 5); // always add one more than the 2nd argument in dsp_add()
}

void autotuned_mix(t_autotuned *x, t_floatarg f)
{
	x->fMix   = CLIP(f,0.,1.);
}

void autotuned_shift(t_autotuned *x, t_floatarg f)
{
	x->fShift   = CLIP(f,-12.,12.);
}

void autotuned_tune(t_autotuned *x, t_floatarg f)
{
	x->fTune  = CLIP(f,-1.,1.);
}

void autotuned_correction(t_autotuned *x, t_floatarg f)
{
	x->fAmount   = CLIP(f,0.,1.);
}

void autotuned_glide(t_autotuned *x, t_floatarg f)
{
	x->fGlide   = CLIP(f,0.,500.);
}


void autotuned_list(t_autotuned *x, t_symbol *m, short argc, t_atom *argv)
{
	if(argc== 12)
	{
		x->fC     = CLIP(atom_getfloat(argv), 0., 1.);  
		x->fDb    = CLIP(atom_getfloat(argv+1), 0., 1.)  ;  
		x->fD     = CLIP(atom_getfloat(argv+2), 0., 1.)  ;  
		x->fEb    = CLIP(atom_getfloat(argv+3), 0., 1.)  ;  
		x->fE     = CLIP(atom_getfloat(argv+4), 0., 1.)  ;  
		x->fF     = CLIP(atom_getfloat(argv+5), 0., 1.)  ;  
		x->fGb    = CLIP(atom_getfloat(argv+6), 0., 1.)  ;  
		x->fG     = CLIP(atom_getfloat(argv+7), 0., 1.)  ;  
		x->fAb    = CLIP(atom_getfloat(argv+8), 0., 1.)  ;  
		x->fA     = CLIP(atom_getfloat(argv+9), 0., 1.)  ;  
		x->fBb    = CLIP(atom_getfloat(argv+10), 0., 1.)  ;  
		x->fB     = CLIP(atom_getfloat(argv+11), 0., 1.)  ;  
	}
	else {post("bad list"); };
}

void autotuned_processclock(t_autotuned *x)
{
	clock_delay(x->clock, x->clockinterval); // schedule the next clock
	
	outlet_float(x->confout, x->conf);
	outlet_float(x->periodout, 1.f / x->pperiod);
}

void autotuned_assist(t_autotuned *x, void *b, long m, long a, char *s)
{
	if (m == 1) //input
	{
		sprintf(s,"Signal Input, Messages");
	}
	else
	{
		switch (a) {	
		case 0:
			sprintf(s,"Signal Output");
			break;
		case 1:
				sprintf(s,"confidence (0-1)(float)");
				break;
		case 2:
				
				sprintf(s,"frequency (hz) (float)");
				break;
				
		}
	}
	
}
