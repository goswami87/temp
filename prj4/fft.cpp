/*
This is traditional 2-radix DIT FFT algorithm implementation.
It is based on conventional 3-loop structure. 
INPUT:
	In_R, In_I[]: Real and Imag parts of Complex signal

OUTPUT:
	In_R, In_I[]: Real and Imag parts of Complex signal
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fft.h"

#if 1 //too much time
unsigned int reverse_bits(unsigned int input) {
	int i, rev = 0;
	for (i = 0; i < M; i++) {
		#pragma HLS unroll
		rev = (rev << 1) | (input & 1);
		input = input >> 1;
	}
	return rev;
}

void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE],
                 DTYPE X_R_OUT[SIZE], DTYPE X_I_OUT[SIZE]) {
			
    unsigned int i, rev;

    for (i = 0; i < SIZE; i++) {
		//#pragma HLS unroll
        rev = reverse_bits(i);

        if (rev > i) {
            // swap
            X_R_OUT[i]   = X_R[rev];
            X_R_OUT[rev] = X_R[i];

            X_I_OUT[i]   = X_I[rev];
            X_I_OUT[rev] = X_I[i];
        } else if (rev < i) {
            // already swapped earlier
        } else {
            // rev == i, copy directly
            X_R_OUT[i] = X_R[i];
            X_I_OUT[i] = X_I[i];
        }
    }
}
void fft_stage(int stage, DTYPE X_R[SIZE], DTYPE X_I[SIZE],DTYPE Out_R[SIZE], DTYPE Out_I[SIZE]) {

	int DFTpts = 1 << stage; // DFT = 2^stage = points in sub DFT
	int numBF = DFTpts / 2; // Butterfly WIDTHS in sub−DFT
	int step = SIZE >> stage;
	DTYPE k = 0;
	DTYPE e = -6.283185307178 / DFTpts;
	DTYPE a = 0.0;
	// Perform butteries for j−th stage
	butterfly_loop:
	for (int j = 0; j < numBF; j++) {
		//#pragma HLS PIPELINE II=1
		//#pragma HLS UNROLL factor=4
		int tw_idx = j * step;

        // Read twiddle from ROM
        DTYPE c = W_real[tw_idx];
        DTYPE s = W_imag[tw_idx];

		a = a + e;
		// Compute butterflies that use same W**k
		dft_loop:
		for (int i = j; i < SIZE; i += DFTpts) {
			#pragma HLS PIPELINE II=1
			int i_lower = i + numBF; // index of lower point in butterfly
			DTYPE temp_R = X_R[i_lower] * c - X_I[i_lower] * s;
			DTYPE temp_I = X_I[i_lower] * c + X_R[i_lower] * s;
			Out_R[i_lower] = X_R[i] - temp_R;
			Out_I[i_lower] = X_I[i] - temp_I;
			Out_R[i] = X_R[i] + temp_R;
			Out_I[i] = X_I[i] + temp_I;
		}
		k += step;
	}
}

//void fft(DTYPE X_R[SIZE], DTYPE X_I[SIZE])
void fft(DTYPE X_R[SIZE], DTYPE X_I[SIZE],DTYPE OUT_X_R[SIZE], DTYPE OUT_X_I[SIZE])
{
#pragma HLS dataflow
	DTYPE Stage_R[M][SIZE], Stage_I[M][SIZE];

	bit_reverse(X_R, X_I, Stage_R[0], Stage_I[0]);
	stage_loop:
	//for (int stage = 1; stage < M; stage++) { // Do M−1 stages of butterflies
	for (int stage = 1; stage < M; stage++){
		#pragma HLS UNROLL
		fft_stage(stage, Stage_R[stage-1], Stage_I[stage-1], Stage_R[stage], Stage_I[stage]);
	}
	fft_stage(M, Stage_R[M-1], Stage_I[M-1], OUT_X_R, OUT_X_I);
}

#endif //too much time