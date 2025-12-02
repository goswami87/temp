/*
This is traditional 2-radix DIT FFT algorithm implementation.
INPUT:
	In_R, In_I[]: Real and Imag parts of Complex signal

OUTPUT:
	Out_R, Out_I[]: Real and Imag parts of Complex signal
*/

#include "fft.h"

void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE],DTYPE O_X_R[SIZE], DTYPE O_X_I[SIZE]);
void fft_stages(DTYPE X_R[SIZE], DTYPE X_I[SIZE], int STAGES, DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]);

void fft(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE])
{
#pragma HLS dataflow	

	//Call fft
	DTYPE Stage0_R[SIZE], Stage0_I[SIZE];	
	DTYPE Stage1_R[SIZE], Stage1_I[SIZE];
	DTYPE Stage2_R[SIZE], Stage2_I[SIZE];
	DTYPE Stage3_R[SIZE], Stage3_I[SIZE];
	DTYPE Stage4_R[SIZE], Stage4_I[SIZE];
	DTYPE Stage5_R[SIZE], Stage5_I[SIZE];
	DTYPE Stage6_R[SIZE], Stage6_I[SIZE];
	DTYPE Stage7_R[SIZE], Stage7_I[SIZE];
	DTYPE Stage8_R[SIZE], Stage8_I[SIZE];
	DTYPE Stage9_R[SIZE], Stage9_I[SIZE];

	bit_reverse(X_R, X_I, Stage0_R, Stage0_I);

	fft_stages(Stage0_R, Stage0_I, 1, Stage1_R, Stage1_I);
	fft_stages(Stage1_R, Stage1_I, 2, Stage2_R, Stage2_I);
	fft_stages(Stage2_R, Stage2_I, 3, Stage3_R, Stage3_I);
	fft_stages(Stage3_R, Stage3_I, 4, Stage4_R, Stage4_I);
	fft_stages(Stage4_R, Stage4_I, 5, Stage5_R, Stage5_I);
	fft_stages(Stage5_R, Stage5_I, 6, Stage6_R, Stage6_I);
	fft_stages(Stage6_R, Stage6_I, 7, Stage7_R, Stage7_I);
	fft_stages(Stage7_R, Stage7_I, 8, Stage8_R, Stage8_I);
	fft_stages(Stage8_R, Stage8_I, 9, Stage9_R, Stage9_I);
	fft_stages(Stage9_R, Stage9_I, 10, OUT_R, OUT_I);

}

void bit_reverse(
    float X_R_in[SIZE], float X_I_in[SIZE],
    float X_R_out[SIZE], float X_I_out[SIZE]
) {
    for (int i = 0; i < SIZE; i++) {
        #pragma HLS PIPELINE II=1

        int rev = 0;
        int temp = i;

        // Reverse M bits
        for (int j = 0; j < M; j++) {
            rev = (rev << 1) | (temp & 1);
            temp >>= 1;
        }

        // Out-of-place bit reverse:
        // No swap needed, simply place input[rev] into output[i]
        X_R_out[i] = X_R_in[rev];
        X_I_out[i] = X_I_in[rev];
    }
}

/*=======================BEGIN: FFT=========================*/

//stages
void fft_stages(DTYPE X_R[SIZE], DTYPE X_I[SIZE], int stage, DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {
#pragma HLS ARRAY_PARTITION variable=X_R block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=X_I block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=OUT_R block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=OUT_I block factor=2 dim=1

#pragma HLS ARRAY_PARTITION variable=W_real complete dim=1
#pragma HLS ARRAY_PARTITION variable=W_imag complete dim=1
//Insert your code here
	int DFTpts = 1 << stage; // DFT = 2^stage = points in sub DFT
	int numBF = DFTpts / 2; // Butterfly WIDTHS in sub−DFT
	int step = SIZE >> stage;
	DTYPE e = -6.283185307178 / DFTpts;
	DTYPE a = 0.0;
	// Perform butteries for j−th stage
	butterfly_loop:
	for (int j = 0; j < numBF; j++) {
		//#pragma HLS PIPELINE II=2
		//#pragma HLS UNROLL  //report time explodes but works
		int tw_idx = j * step;

        // Read twiddle from ROM
        DTYPE c = W_real[tw_idx];
        DTYPE s = W_imag[tw_idx];

		a = a + e;
		// Compute butterflies that use same W**k
		dft_loop:
		for (int i = j; i < SIZE; i += DFTpts) {
			#pragma HLS PIPELINE II=2
			//#pragma HLS UNROLL factor=2
			int i_lower = i + numBF; // index of lower point in butterfly
			DTYPE temp_R = X_R[i_lower] * c - X_I[i_lower] * s;
			DTYPE temp_I = X_I[i_lower] * c + X_R[i_lower] * s;
			OUT_R[i_lower] = X_R[i] - temp_R;
			OUT_I[i_lower] = X_I[i] - temp_I;
			OUT_R[i] = X_R[i] + temp_R;
			OUT_I[i] = X_I[i] + temp_I;
		}
	}
}
/*=======================END: FFT=========================*/
