/*
This is traditional 2-radix DIT FFT algorithm implementation.
INPUT:
	In_R, In_I[]: Real and Imag parts of Complex signal

OUTPUT:
	Out_R, Out_I[]: Real and Imag parts of Complex signal
*/

#include "fft.h"

//void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE]);
void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE],DTYPE O_X_R[SIZE], DTYPE O_X_I[SIZE]);
void fft_stage_first(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]);
void fft_stages(DTYPE X_R[SIZE], DTYPE X_I[SIZE], int STAGES, DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]);
void fft_stage_last(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]);

#if 0
void fft(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE])
{
#pragma HLS dataflow	
	//bit_reverse(X_R, X_I);
	//Call fft
	static DTYPE Stage0_R[SIZE], Stage0_I[SIZE];	
	static DTYPE Stage1_R[SIZE], Stage1_I[SIZE];
	static DTYPE Stage2_R[SIZE], Stage2_I[SIZE];
	static DTYPE Stage3_R[SIZE], Stage3_I[SIZE];
	static DTYPE Stage4_R[SIZE], Stage4_I[SIZE];
	static DTYPE Stage5_R[SIZE], Stage5_I[SIZE];
	static DTYPE Stage6_R[SIZE], Stage6_I[SIZE];
	static DTYPE Stage7_R[SIZE], Stage7_I[SIZE];
	static DTYPE Stage8_R[SIZE], Stage8_I[SIZE];
	static DTYPE Stage9_R[SIZE], Stage9_I[SIZE];
#pragma HLS ARRAY_PARTITION variable=Stage0_R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage0_I cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage1_R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage1_I cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage2_R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage2_I cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage3_R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage3_I cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage4_R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage4_I cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage5_R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage5_I cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage6_R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage6_I cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage7_R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage7_I cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage8_R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage8_I cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage9_R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=Stage9_I cyclic factor=4 dim=1

	bit_reverse(X_R, X_I, Stage0_R, Stage0_I);

	//fft_stage_first(Stage0_R, Stage0_I, Stage1_R, Stage1_I);
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
	//fft_stage_last(Stage9_R, Stage9_I, OUT_R, OUT_I);

}
#else
void fft(DTYPE X_R[SIZE], DTYPE X_I[SIZE],
         DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE])
{
#pragma HLS DATAFLOW

    // Ping-pong buffers for dataflow compatibility
    static DTYPE buf0_R[SIZE], buf0_I[SIZE];
    static DTYPE buf1_R[SIZE], buf1_I[SIZE];
    static DTYPE buf2_R[SIZE], buf2_I[SIZE];
    static DTYPE buf3_R[SIZE], buf3_I[SIZE];
    static DTYPE buf4_R[SIZE], buf4_I[SIZE];
    static DTYPE buf5_R[SIZE], buf5_I[SIZE];
    static DTYPE buf6_R[SIZE], buf6_I[SIZE];
    static DTYPE buf7_R[SIZE], buf7_I[SIZE];
    static DTYPE buf8_R[SIZE], buf8_I[SIZE];
    static DTYPE buf9_R[SIZE], buf9_I[SIZE];
/*#pragma HLS ARRAY_PARTITION variable=X_R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=X_I cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=OUT_R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=OUT_I cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=buf0_R complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf0_I complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf1_R complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf1_I complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf2_R complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf2_I complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf3_R complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf3_I complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf4_R complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf4_I complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf5_R complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf5_I complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf6_R complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf6_I complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf7_R complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf7_I complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf8_R complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf8_I complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf9_R complete dim=1
#pragma HLS ARRAY_PARTITION variable=buf9_I complete dim=1*/
// (repeat for all buffers – this allows pipelining)

    // Stage 0: bit reverse
    bit_reverse(X_R, X_I, buf0_R, buf0_I);

    // FFT stages
    fft_stages(buf0_R, buf0_I, 1, buf1_R, buf1_I);
    fft_stages(buf1_R, buf1_I, 2, buf2_R, buf2_I);
    fft_stages(buf2_R, buf2_I, 3, buf3_R, buf3_I);
    fft_stages(buf3_R, buf3_I, 4, buf4_R, buf4_I);
    fft_stages(buf4_R, buf4_I, 5, buf5_R, buf5_I);
    fft_stages(buf5_R, buf5_I, 6, buf6_R, buf6_I);
    fft_stages(buf6_R, buf6_I, 7, buf7_R, buf7_I);
    fft_stages(buf7_R, buf7_I, 8, buf8_R, buf8_I);
    fft_stages(buf8_R, buf8_I, 9, buf9_R, buf9_I);

    // Final stage
    fft_stages(buf9_R, buf9_I, 10, OUT_R, OUT_I);
}
#endif
#if 0
void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE]){
	//Insert your code here
	for (int i = 0; i < SIZE; i++) {
		int rev = 0;
		int temp = i;

		// Reverse M bits
		for (int j = 0; j < M; j++) {
			rev = (rev << 1) | (temp & 1);
			temp >>= 1;
		}

		// Swap only if rev > i (avoid double swap)
		if (rev > i) {
			int temp_R = X_R[i];
			int temp_I = X_I[i];
			X_R[i] = X_R[rev];
			X_I[i] = X_I[rev];
			X_R[rev] = temp_R;
			X_I[rev] = temp_I;
		}
	}

}
#else
void bit_reverse(
    DTYPE X_R_in[SIZE], DTYPE X_I_in[SIZE],
    DTYPE X_R_out[SIZE], DTYPE X_I_out[SIZE])
{
bitrev_loop:
    for (int i = 0; i < SIZE; i++) {
        #pragma HLS PIPELINE II=1

        // Unrolled bit reversal
        int rev = 0;
        int tmp = i;

    bitcalc:
        for (int b = 0; b < M; b++) {
            #pragma HLS UNROLL
            rev = (rev << 1) | (tmp & 1);
            tmp >>= 1;
        }

        // Out-of-place write
        X_R_out[i] = X_R_in[rev];
        X_I_out[i] = X_I_in[rev];
    }
}
#endif
#if 0
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

/*
unsigned int reverse_bits(unsigned int input) {
//#pragma HLS INLINE
    int i, rev = 0;
    for (i = 0; i < M; i++) {
        #pragma HLS UNROLL
		//#pragma HLS PIPELINE II=1
        rev = (rev << 1) | (input & 1);
        input = input >> 1;
    }
    return rev;
}

void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE],
                 DTYPE X_R_OUT[SIZE], DTYPE X_I_OUT[SIZE]) {
    // DIRECT MAPPING OPTIMIZATION
    // Instead of swapping, we simply wire input[rev] to output[i].
    // With full unroll, this becomes pure routing (0 cycles logic, only wire delay).
    loop_bit_reverse:
    for (int i = 0; i < SIZE; i++) {
        #pragma HLS UNROLL
		//#pragma HLS PIPELINE II=1
        unsigned int rev = reverse_bits(i);
        X_R_OUT[i] = X_R[rev];
        X_I_OUT[i] = X_I[rev];
    }
}
*/
#endif
/*=======================BEGIN: FFT=========================*/
//stage 1
void fft_stage_first(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {

//Insert your code here
	fft_stages(X_R,  X_I,1, OUT_R, OUT_I);
}
#if 0
//stages
void fft_stages(DTYPE X_R[SIZE], DTYPE X_I[SIZE], int stage, DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {

//Insert your code here
	int DFTpts = 1 << stage; // DFT = 2^stage = points in sub DFT
	int numBF = DFTpts / 2; // Butterfly WIDTHS in sub−DFT
	int step = SIZE >> stage;
	DTYPE k = 0;
	DTYPE e = -6.283185307178 / DFTpts;
	DTYPE a = 0.0;
	// Perform butteries for j−th stage
	butterfly_loop:
	for (int j = 0; j < numBF; j++) {
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
			OUT_R[i_lower] = X_R[i] - temp_R;
			OUT_I[i_lower] = X_I[i] - temp_I;
			OUT_R[i] = X_R[i] + temp_R;
			OUT_I[i] = X_I[i] + temp_I;
		}
		k += step;
	}
}
#else
void fft_stages(
    DTYPE X_R[SIZE], DTYPE X_I[SIZE],
    int stage,
    DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE])
{/*
#pragma HLS ARRAY_PARTITION variable=X_R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=X_I cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=OUT_R cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=OUT_I cyclic factor=4 dim=1	
    #pragma HLS ARRAY_PARTITION variable=X_R   complete dim=1
    #pragma HLS ARRAY_PARTITION variable=X_I   complete dim=1
    #pragma HLS ARRAY_PARTITION variable=OUT_R complete dim=1
    #pragma HLS ARRAY_PARTITION variable=OUT_I complete dim=1 */
    int DFTpts = 1 << stage; 
    int half = DFTpts >> 1;
    int step = SIZE >> stage;

j_loop:
    for (int j = 0; j < half; j++) {
        #pragma HLS PIPELINE II=1

        int tw = j * step;
        DTYPE wr = W_real[tw];
        DTYPE wi = W_imag[tw];

    butterfly_loop:
        for (int i = j; i < SIZE; i += DFTpts) {
            //#pragma HLS UNROLL factor=2  // improves ILP without exploding DSP usage

            int lo = i + half;

            // Complex multiply
            DTYPE xr = X_R[lo];
            DTYPE xi = X_I[lo];

            DTYPE tr = xr * wr - xi * wi;
            DTYPE ti = xi * wr + xr * wi;

            // Store upper and lower
            OUT_R[i] = X_R[i] + tr;
            OUT_I[i] = X_I[i] + ti;

            OUT_R[lo] = X_R[i] - tr;
            OUT_I[lo] = X_I[i] - ti;
        }
    }
}

#endif

//last stage
void fft_stage_last(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {

//Insert your code here
	fft_stages(X_R,  X_I,10, OUT_R, OUT_I);
}
/*=======================END: FFT=========================*/




