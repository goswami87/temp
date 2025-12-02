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
#if 0
void fft(DTYPE X_R[SIZE], DTYPE X_I[SIZE])
{

	DTYPE temp_R;		/*temporary storage complex variable*/
	DTYPE temp_I;		/*temporary storage complex variable*/


	int i,j,k;			/* loop indexes */
	int i_lower;		/* Index of lower point in butterfly */
	int step;

	int stage;
	int DFTpts;
	int numBF;			/*Butterfly Width*/

	//int N2 = SIZE2;	/* N2=N>>1 */

	/*=====================BEGIN BIT REBERSAL===========================*/
	// write your code here
	for (i = 0; i < SIZE; i++) {
		int rev = 0;
		int temp = i;

		// Reverse M bits
		for (j = 0; j < M; j++) {
			rev = (rev << 1) | (temp & 1);
			temp >>= 1;
		}

		// Swap only if rev > i (avoid double swap)
		if (rev > i) {
			temp_R = X_R[i];
			temp_I = X_I[i];
			X_R[i] = X_R[rev];
			X_I[i] = X_I[rev];
			X_R[rev] = temp_R;
			X_I[rev] = temp_I;
		}
	}
	/*++++++++++++++++++++++END OF BIT REVERSAL++++++++++++++++++++++++++*/

	/*=======================BEGIN: FFT=========================*/
	// Do M stages of butterflies
	//step=N2;
	DTYPE a, e, c, s;

	stages:
	for(stage=1; stage<= M; stage++)
	{
		DFTpts = 1 << stage;		// DFT = 2^stage = points in sub DFT
		numBF = DFTpts/2; 			// Butterfly WIDTHS in sub-DFT
		//k=0;

		e = -6.283185307178/DFTpts;

		a = 0.0;
		// Perform butterflies for j-th stage
		butterfly:
		for(j=0; j<numBF; j++)
		{

			c = cos(a);
			s = sin(a);
			a = a + e;

			// Compute butterflies that use same W**k
			DFT:
			for(i=j; i<SIZE; i += DFTpts)
			{

				i_lower = i + numBF;			//index of lower point in butterfly
				temp_R = X_R[i_lower]*c- X_I[i_lower]*s;
				temp_I = X_I[i_lower]*c+ X_R[i_lower]*s;

				X_R[i_lower] = X_R[i] - temp_R;
				X_I[i_lower] = X_I[i] - temp_I;
				X_R[i] = X_R[i] + temp_R;
				X_I[i] = X_I[i] + temp_I;
			}
			//k+=step;
		}
		//step=step/2;
	}
}
#endif
#if 0
/*=======================END: FFT=========================*/
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
		#pragma HLS unroll
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
		int tw_idx = j * step;

        // Read twiddle from ROM
        DTYPE c = W_real[tw_idx];
        DTYPE s = W_imag[tw_idx];

		//DTYPE c = cos(a);
		//DTYPE s = sin(a);
		a = a + e;
		// Compute butterflies that use same W**k
		dft_loop:
		for (int i = j; i < SIZE; i += DFTpts) {
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

void fft(DTYPE X_R[SIZE], DTYPE X_I[SIZE])
{
//#pragma HLS dataflow
	DTYPE Stage_R[M][SIZE], Stage_I[M][SIZE];
#pragma HLS array partition variable=Stage_R dim=1 complete
#pragma HLS array partition variable=Stage_I dim=1 complete
	//bit_reverse(X_R, X_I, Stage_R[0], Stage_I[0]);
	bit_reverse(X_R, X_I, Stage_R[0], Stage_I[0]);
	stage_loop:
	//for (int stage = 1; stage < M; stage++) { // Do M−1 stages of butterflies
	for (int stage = 1; stage < M; stage++){
		//#pragma HLS unroll
		#pragma HLS pipeline
		fft_stage(stage, Stage_R[stage-1], Stage_I[stage-1], Stage_R[stage], Stage_I[stage]);
	}
	fft_stage(M, Stage_R[M-1], Stage_I[M-1], X_R, X_I);
}
#endif
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
/*					
#pragma HLS ARRAY_PARTITION variable=X_R complete dim=1
#pragma HLS ARRAY_PARTITION variable=X_I complete dim=1	
#pragma HLS ARRAY_PARTITION variable=X_R_OUT complete dim=1
#pragma HLS ARRAY_PARTITION variable=X_I_OUT complete dim=1		
*/				
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
/*	
#pragma HLS ARRAY_PARTITION variable=X_R complete dim=1
#pragma HLS ARRAY_PARTITION variable=X_I complete dim=1	
#pragma HLS ARRAY_PARTITION variable=Out_R complete dim=1
#pragma HLS ARRAY_PARTITION variable=Out_I complete dim=1	
*/
	int DFTpts = 1 << stage; // DFT = 2^stage = points in sub DFT
	int numBF = DFTpts / 2; // Butterfly WIDTHS in sub−DFT
	int step = SIZE >> stage;
	DTYPE k = 0;
	DTYPE e = -6.283185307178 / DFTpts;
	DTYPE a = 0.0;
	// Perform butteries for j−th stage
	butterfly_loop:
	for (int j = 0; j < numBF; j++) {
		#pragma HLS PIPELINE II=1
		#pragma HLS UNROLL factor=4
		int tw_idx = j * step;

        // Read twiddle from ROM
        DTYPE c = W_real[tw_idx];
        DTYPE s = W_imag[tw_idx];

		//DTYPE c = cos(a);
		//DTYPE s = sin(a);
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

void fft(DTYPE X_R[SIZE], DTYPE X_I[SIZE])
{
#pragma HLS dataflow
	DTYPE Stage_R[M][SIZE], Stage_I[M][SIZE];
//#pragma HLS array_partition variable=Stage_R dim=1 complete
//#pragma HLS array_partition variable=Stage_I dim=1 complete
	//bit_reverse(X_R, X_I, Stage_R[0], Stage_I[0]);
	bit_reverse(X_R, X_I, Stage_R[0], Stage_I[0]);
	stage_loop:
	//for (int stage = 1; stage < M; stage++) { // Do M−1 stages of butterflies
	for (int stage = 1; stage < M; stage++){
		//#pragma HLS unroll
		//#pragma HLS pipeline
		fft_stage(stage, Stage_R[stage-1], Stage_I[stage-1], Stage_R[stage], Stage_I[stage]);
	}
	fft_stage(M, Stage_R[M-1], Stage_I[M-1], X_R, X_I);
}

#endif //too much time

#if 0
#include <math.h>
#include <complex>
// Assuming these are defined in your header
// #define SIZE 1024
// #define M 10 
// typedef float DTYPE;

// Precompute reverse bits to avoid logic overhead during runtime
// or rely on HLS constant propagation if M is static.
unsigned int reverse_bits(unsigned int input) {
#pragma HLS INLINE
    int i, rev = 0;
    for (i = 0; i < M; i++) {
        //#pragma HLS UNROLL
		#pragma HLS PIPELINE II=1
        rev = (rev << 1) | (input & 1);
        input = input >> 1;
    }
    return rev;
}

void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE],
                 DTYPE X_R_OUT[SIZE], DTYPE X_I_OUT[SIZE]) {
#pragma HLS INLINE off
//#pragma HLS RESOURCE variable=X_R core=RAM_2P_BRAM
//#pragma HLS RESOURCE variable=X_I core=RAM_2P_BRAM
    // Partitioning allows parallel access to all elements (Wires, not RAM)
    //#pragma HLS ARRAY_PARTITION variable=X_R complete dim=1
    //#pragma HLS ARRAY_PARTITION variable=X_I complete dim=1 
    //#pragma HLS ARRAY_PARTITION variable=X_R_OUT complete dim=1
    //#pragma HLS ARRAY_PARTITION variable=X_I_OUT complete dim=1                     

    // DIRECT MAPPING OPTIMIZATION
    // Instead of swapping, we simply wire input[rev] to output[i].
    // With full unroll, this becomes pure routing (0 cycles logic, only wire delay).
    loop_bit_reverse:
    for (int i = 0; i < SIZE; i++) {
        //#pragma HLS UNROLL
		#pragma HLS PIPELINE II=1
        unsigned int rev = reverse_bits(i);
        X_R_OUT[i] = X_R[rev];
        X_I_OUT[i] = X_I[rev];
    }
}

void fft_stage(int stage, DTYPE X_R[SIZE], DTYPE X_I[SIZE], 
               DTYPE Out_R[SIZE], DTYPE Out_I[SIZE]) {
#pragma HLS INLINE off
    // Partitioning ensures we can read/write multiple values per clock
    //#pragma HLS ARRAY_PARTITION variable=X_R complete dim=1
    //#pragma HLS ARRAY_PARTITION variable=X_I complete dim=1 
    //#pragma HLS ARRAY_PARTITION variable=Out_R complete dim=1
    //#pragma HLS ARRAY_PARTITION variable=Out_I complete dim=1   
//#pragma HLS RESOURCE variable=X_R core=RAM_2P_BRAM
//#pragma HLS RESOURCE variable=X_I core=RAM_2P_BRAM

    int DFTpts = 1 << stage;       // Points in sub-DFT
    int numBF = DFTpts / 2;        // Butterfly width
    int step = SIZE >> stage;      // Stride for twiddle factors

    // LOOP FLATTENING OPTIMIZATION
    // We flatten the nested loop to ensure continuous pipelining (II=1)
    // Total Butterflies per stage = SIZE / 2
    
    // Note: HLS can automatically flatten nested loops if bounds are constant.
    // Since 'stage' is constant (due to unrolling in top function), 
    // DFTpts and numBF are constants.
    
    butterfly_loop:
    for (int j = 0; j < numBF; j++) {
        // Compute Twiddle Index once per row
        int tw_idx = j * step;
        DTYPE c = W_real[tw_idx];
        DTYPE s = W_imag[tw_idx];

        dft_loop:
        for (int i = j; i < SIZE; i += DFTpts) {
            #pragma HLS LOOP_FLATTEN
            #pragma HLS PIPELINE II=1
            
            int i_lower = i + numBF; 

            // Butterfly Calculation
            DTYPE temp_R = X_R[i_lower] * c - X_I[i_lower] * s;
            DTYPE temp_I = X_I[i_lower] * c + X_R[i_lower] * s;

            Out_R[i_lower] = X_R[i] - temp_R;
            Out_I[i_lower] = X_I[i] - temp_I;
            Out_R[i]       = X_R[i] + temp_R;
            Out_I[i]       = X_I[i] + temp_I;
        }
    }
}

void fft(DTYPE X_R[SIZE], DTYPE X_I[SIZE])
{
    // PARTITIONING IS CRITICAL
    // Without this, you only have 1 or 2 memory ports, limiting throughput.
    //#pragma HLS ARRAY_PARTITION variable=X_R complete dim=1
    //#pragma HLS ARRAY_PARTITION variable=X_I complete dim=1
//#pragma HLS RESOURCE variable=X_R core=RAM_2P_BRAM
//#pragma HLS RESOURCE variable=X_I core=RAM_2P_BRAM
    // Intermediate storage for Dataflow
    // We need M stages of storage.
    // Partitioning creates registers instead of BRAM, allowing simultaneous access.
    DTYPE Stage_R[M+1][SIZE];
    DTYPE Stage_I[M+1][SIZE];
 //   #pragma HLS ARRAY_PARTITION variable=Stage_R dim=0 complete
 //   #pragma HLS ARRAY_PARTITION variable=Stage_I dim=0 complete

    // DATAFLOW ARCHITECTURE
    // This allows fft_stage(1) to run at the same time as fft_stage(2), etc.
    #pragma HLS DATAFLOW

    // 1. Bit Reverse (Input -> Stage 0)
    bit_reverse(X_R, X_I, Stage_R[0], Stage_I[0]);

    // 2. Unroll the stages
    // By unrolling, we create physical instances of fft_stage for every iteration.
    // The compiler sees: fft_stage(1, ...); fft_stage(2, ...); etc.
    stage_loop:
    for (int stage = 1; stage <= M; stage++) {
        //#pragma HLS UNROLL
            //#pragma HLS PIPELINE II=1
        // Input from previous stage buffer, Output to current stage buffer
        fft_stage(stage, 
                  Stage_R[stage-1], Stage_I[stage-1], 
                  Stage_R[stage],   Stage_I[stage]);
    }

    // 3. Copy final result back to X_R/X_I
    // (Or you can just use Stage_R[M] as your output if interface allows)
    write_out_loop:
    for (int i = 0; i < SIZE; i++) {
        //#pragma HLS UNROLL
		#pragma HLS PIPELINE II=1
        X_R[i] = Stage_R[M][i];
        X_I[i] = Stage_I[M][i];
    }
}
#endif