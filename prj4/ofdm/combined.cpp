/*
 This is traditional 2-radix DIT FFT algorithm implementation.
 INPUT:
 In_R, In_I[]: Real and Imag parts of Complex signal

 OUTPUT:
 Out_R, Out_I[]: Real and Imag parts of Complex signal
 */

#include "fft.h"

void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE]);
void fft_stage_first(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE],
		DTYPE OUT_I[SIZE]);
void fft_stages(DTYPE X_R[SIZE], DTYPE X_I[SIZE], int STAGES, DTYPE OUT_R[SIZE],
		DTYPE OUT_I[SIZE]);
void fft_stage_last(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE],
		DTYPE OUT_I[SIZE]);
void qpsk_decode(DTYPE R[SIZE], DTYPE I[SIZE], int D[SIZE]);

void demod(DTYPE X_R[SIZE], DTYPE X_I[SIZE], int D[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE])
{

	fft(X_R, X_I, OUT_R, OUT_I);
	qpsk_decode(OUT_R, OUT_I, D);

}
#if 0
void fft(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE],
		DTYPE OUT_I[SIZE]) {

	bit_reverse(X_R, X_I);

	//Call fft
	DTYPE Stage1_R[SIZE], Stage1_I[SIZE];
	DTYPE Stage2_R[SIZE], Stage2_I[SIZE];
	DTYPE Stage3_R[SIZE], Stage3_I[SIZE];
	DTYPE Stage4_R[SIZE], Stage4_I[SIZE];
	DTYPE Stage5_R[SIZE], Stage5_I[SIZE];
	DTYPE Stage6_R[SIZE], Stage6_I[SIZE];
	DTYPE Stage7_R[SIZE], Stage7_I[SIZE];
	DTYPE Stage8_R[SIZE], Stage8_I[SIZE];
	DTYPE Stage9_R[SIZE], Stage9_I[SIZE];

	fft_stage_first(X_R, X_I, Stage1_R, Stage1_I);
	fft_stages(Stage1_R, Stage1_I, 2, Stage2_R, Stage2_I);
	fft_stages(Stage2_R, Stage2_I, 3, Stage3_R, Stage3_I);
	fft_stages(Stage3_R, Stage3_I, 4, Stage4_R, Stage4_I);
	fft_stages(Stage4_R, Stage4_I, 5, Stage5_R, Stage5_I);
	fft_stages(Stage5_R, Stage5_I, 6, Stage6_R, Stage6_I);
	fft_stages(Stage6_R, Stage6_I, 7, Stage7_R, Stage7_I);
	fft_stages(Stage7_R, Stage7_I, 8, Stage8_R, Stage8_I);
	fft_stages(Stage8_R, Stage8_I, 9, Stage9_R, Stage9_I);
	fft_stage_last(Stage9_R, Stage9_I, OUT_R, OUT_I);

	bit_reverse(OUT_R, OUT_I);

}

/*=======================BEGIN: FFT=========================*/
//stage 1
void fft_stage_first(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {

	//Insert your code here

}

//stages
void fft_stages(DTYPE X_R[SIZE], DTYPE X_I[SIZE], int stage, DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {

	//Insert your code here
}

//last stage
void fft_stage_last(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {

	//Insert your code here

}

void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE X_R_Copy[SIZE], DTYPE X_I_Copy[SIZE])
{

	//Insert your code here
}
#endif
/*=======================END: FFT=========================*/
/*
This is traditional 2-radix DIT FFT algorithm implementation.
INPUT:
	In_R, In_I[]: Real and Imag parts of Complex signal

OUTPUT:
	Out_R, Out_I[]: Real and Imag parts of Complex signal
*/

//#include "fft.h"

DTYPE e_LUT[M] = {
-3.141592653589,
-1.570796326795,
-0.785398163397,
-0.392699081699,
-0.196349540849,
-0.098174770425,
-0.049087385212,
-0.024543692606,
-0.012271846303,
-0.006135914655};

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
#pragma HLS ARRAY_PARTITION variable=Stage0_R block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage0_I block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage1_R block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage1_I block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage2_R block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage2_I block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage3_R block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage3_I block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage4_R block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage4_I block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage5_R block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage5_I block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage6_R block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage6_I block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage7_R block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage7_I block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage8_R block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage8_I block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage9_R block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage9_I block factor=2 dim=1

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
#pragma HLS ARRAY_PARTITION variable=X_R_in block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=X_I_in block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=X_R_out block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=X_I_out block factor=2 dim=1	
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
/*#pragma HLS ARRAY_PARTITION variable=X_R   complete dim=1
#pragma HLS ARRAY_PARTITION variable=X_I   complete dim=1
#pragma HLS ARRAY_PARTITION variable=OUT_R complete dim=1
#pragma HLS ARRAY_PARTITION variable=OUT_I complete dim=1*/

#pragma HLS ARRAY_PARTITION variable=X_R block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=X_I block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=OUT_R block factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=OUT_I block factor=2 dim=1

#pragma HLS ARRAY_PARTITION variable=W_real complete dim=1
#pragma HLS ARRAY_PARTITION variable=W_imag complete dim=1
//Insert your code here
	int DFTpts = 1 << stage; // DFT = 2^stage = points in sub DFT
	//int numBF = DFTpts / 2; // Butterfly WIDTHS in sub−DFT
	int numBF = DFTpts >> 1; // Butterfly WIDTHS in sub−DFT
	int step = SIZE >> stage;
	//DTYPE k = 0;
	//DTYPE e = -6.283185307178 / DFTpts;
	DTYPE e = e_LUT[stage-1];
	DTYPE a = 0.0;
	// Perform butteries for j−th stage
	butterfly_loop:
	for (int j = 0; j < numBF; j++) {
		//#pragma HLS PIPELINE II=2
		//#pragma HLS UNROLL factor=8
		int tw_idx = j * step;

        // Read twiddle from ROM
        DTYPE c = W_real[tw_idx];
        DTYPE s = W_imag[tw_idx];

		a = a + e;
		// Compute butterflies that use same W**k
		dft_loop:
		for (int i = j; i < SIZE; i += DFTpts) {
			#pragma HLS PIPELINE II=2 //value 1 cause 4% dsp overshoot
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

#include "fft.h"

static unsigned short count;
static DTYPE xr[ SIZE ];
static DTYPE xi[ SIZE ];
static DTYPE xr_out[ SIZE ];
static DTYPE xi_out[ SIZE ];
static int   dout[ SIZE ];
#if 0
void ofdm_receiver( volatile DTYPE *inptr, volatile uint32_t *outptr )
{
#pragma AP interface ap_fifo port=inptr
#pragma AP interface ap_fifo port=outptr
#pragma AP interface ap_ctrl_none port=return

	*outptr++ = dout[ count ];

	xr[ count ] = *inptr++;
	xi[ count ] = *inptr++;
	count++;
	if( count == 1024 ){
		count = 0;
		demod( xr, xi, dout, xr_out, xi_out );
	}
}
#else
void ofdm_receiver(volatile DTYPE *inptr, volatile uint32_t *outptr)
{
#pragma AP interface ap_fifo port=inptr
#pragma AP interface ap_fifo port=outptr
#pragma AP interface ap_ctrl_none port=return

    // --- READ INPUT SAMPLE ---
    DTYPE r = *inptr++;
    DTYPE i = *inptr++;

    xr[count] = r;
    xi[count] = i;

    // --- PROCESS BLOCK WHEN FULL ---
    if(count == 1023){
        demod(xr, xi, dout, xr_out, xi_out);
    }

    // --- OUTPUT SYMBOL ---
    *outptr = dout[count];

    // Increment count AFTER using dout[count]
    count++;
    if(count == 1024)
        count = 0;
}
#endif
#include "fft.h"
#include <stdio.h>

void qpsk_decode(DTYPE R[SIZE], DTYPE I[SIZE], int D[SIZE]) {

	//Write your code here
#pragma HLS PIPELINE II=1
#pragma HLS ARRAY_PARTITION variable=R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=I cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=D cyclic factor=4 dim=1

decode_loop:
    for (int n = 0; n < SIZE; n++) {

        // Hard-decision:
        int bit_I = (R[n] >= 0) ? 0 : 1;  // MSB
        int bit_Q = (I[n] >= 0) ? 0 : 1;  // LSB

        // Pack 2 bits into one integer
        int symbol = (bit_I << 1) | bit_Q;

        D[n] = symbol;
    }
}

