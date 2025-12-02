// --- fft.h (Conceptual) ---
#include "hls_stream.h" // ⭐ New inclusion
#include "fft.h"
// Define stream types for clarity
typedef hls::stream<DTYPE> DTYPE_STREAM;

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

// New function prototypes using streams:
void bit_reverse_stream(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE_STREAM& O_X_R, DTYPE_STREAM& O_X_I);
void fft_stages_stream(DTYPE_STREAM& X_R, DTYPE_STREAM& X_I, int STAGES, DTYPE_STREAM& OUT_R, DTYPE_STREAM& OUT_I);
void fft(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]); // Top level interface remains the same

void fft(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE])
{
#pragma HLS dataflow    

    // ⭐ Stage Interconnection using hls::stream
    // These now represent FIFOs connecting the modules.
    DTYPE_STREAM Stage0_R("S0R"), Stage0_I("S0I");   
    DTYPE_STREAM Stage1_R("S1R"), Stage1_I("S1I");
    DTYPE_STREAM Stage2_R("S2R"), Stage2_I("S2I");
    DTYPE_STREAM Stage3_R("S3R"), Stage3_I("S3I");
    DTYPE_STREAM Stage4_R("S4R"), Stage4_I("S4I");
    DTYPE_STREAM Stage5_R("S5R"), Stage5_I("S5I");
    DTYPE_STREAM Stage6_R("S6R"), Stage6_I("S6I");
    DTYPE_STREAM Stage7_R("S7R"), Stage7_I("S7I");
    DTYPE_STREAM Stage8_R("S8R"), Stage8_I("S8I");
    DTYPE_STREAM Stage9_R("S9R"), Stage9_I("S9I");
    DTYPE_STREAM Final_R("FR"), Final_I("FI"); // Buffer for the last stage output

    // 1. Bit Reverse (Array Input -> Stream Output)
    bit_reverse_stream(X_R, X_I, Stage0_R, Stage0_I);

    // 2. FFT Stages (Stream to Stream)
    fft_stages_stream(Stage0_R, Stage0_I, 1, Stage1_R, Stage1_I);
    fft_stages_stream(Stage1_R, Stage1_I, 2, Stage2_R, Stage2_I);
    fft_stages_stream(Stage2_R, Stage2_I, 3, Stage3_R, Stage3_I);
    fft_stages_stream(Stage3_R, Stage3_I, 4, Stage4_R, Stage4_I);
    fft_stages_stream(Stage4_R, Stage4_I, 5, Stage5_R, Stage5_I);
    fft_stages_stream(Stage5_R, Stage5_I, 6, Stage6_R, Stage6_I);
    fft_stages_stream(Stage6_R, Stage6_I, 7, Stage7_R, Stage7_I);
    fft_stages_stream(Stage7_R, Stage7_I, 8, Stage8_R, Stage8_I);
    fft_stages_stream(Stage8_R, Stage8_I, 9, Stage9_R, Stage9_I);
    fft_stages_stream(Stage9_R, Stage9_I, 10, Final_R, Final_I);

    // 3. Final Output (Stream to Array Output)
    for(int i = 0; i < SIZE; i++) {
        #pragma HLS PIPELINE II=1
        OUT_R[i] = Final_R.read();
        OUT_I[i] = Final_I.read();
    }
}

void bit_reverse_stream(
    DTYPE X_R_in[SIZE], DTYPE X_I_in[SIZE],
    DTYPE_STREAM& X_R_out, DTYPE_STREAM& X_I_out
) {
    // 1. Buffer the output of bit reversal
    DTYPE Temp_R[SIZE], Temp_I[SIZE];
    #pragma HLS ARRAY_PARTITION variable=Temp_R block factor=2 dim=1
    #pragma HLS ARRAY_PARTITION variable=Temp_I block factor=2 dim=1

    // Perform bit reversal (standard array operations)
    for (int i = 0; i < SIZE; i++) {
        #pragma HLS PIPELINE II=1
        int rev = 0;
        int temp = i;
        for (int j = 0; j < M; j++) {
            rev = (rev << 1) | (temp & 1);
            temp >>= 1;
        }
        Temp_R[i] = X_R_in[rev];
        Temp_I[i] = X_I_in[rev];
    }

    // 2. Transfer data sequentially to the output stream
    for(int i = 0; i < SIZE; i++) {
        #pragma HLS PIPELINE II=1
        X_R_out.write(Temp_R[i]);
        X_I_out.write(Temp_I[i]);
    }
}

void fft_stages_stream(
    DTYPE_STREAM& X_R_in, DTYPE_STREAM& X_I_in, 
    int stage, 
    DTYPE_STREAM& OUT_R_out, DTYPE_STREAM& OUT_I_out
) {
    // ⭐ 1. Local Arrays for Random Access
    DTYPE X_R[SIZE], X_I[SIZE];
    DTYPE OUT_R[SIZE], OUT_I[SIZE];
    
    // Maintain partitioning within the local arrays
    #pragma HLS ARRAY_PARTITION variable=X_R block factor=2 dim=1 // Increased factor for II=2
    #pragma HLS ARRAY_PARTITION variable=X_I block factor=2 dim=1
    #pragma HLS ARRAY_PARTITION variable=OUT_R block factor=2 dim=1
    #pragma HLS ARRAY_PARTITION variable=OUT_I block factor=2 dim=1

    // 2. Read all data sequentially from the input streams
    for(int i = 0; i < SIZE; i++) {
        #pragma HLS PIPELINE II=1
        X_R[i] = X_R_in.read();
        X_I[i] = X_I_in.read();
    }

    // 3. Perform FFT Stage (Original logic - operates on local arrays)
    int DFTpts = 1 << stage; 
    int numBF = DFTpts >> 1; 
    int step = SIZE >> stage;
    DTYPE e = e_LUT[stage-1];
    DTYPE a = 0.0;
    
    // Butterfly and DFT loops (No change to logic, just array names)
    butterfly_loop:
    for (int j = 0; j < numBF; j++) {
        //#pragma HLS PIPELINE II=2
        int tw_idx = j * step;
        DTYPE c = W_real[tw_idx];
        DTYPE s = W_imag[tw_idx];
        a = a + e;
        
        dft_loop:
        for (int i = j; i < SIZE; i += DFTpts) {
            #pragma HLS PIPELINE II=2 
            int i_lower = i + numBF; 
            DTYPE temp_R = X_R[i_lower] * c - X_I[i_lower] * s;
            DTYPE temp_I = X_I[i_lower] * c + X_R[i_lower] * s;
            OUT_R[i_lower] = X_R[i] - temp_R;
            OUT_I[i_lower] = X_I[i] - temp_I;
            OUT_R[i] = X_R[i] + temp_R;
            OUT_I[i] = X_I[i] + temp_I;
        }
    }
    
    // 4. Write all data sequentially to the output streams
    for(int i = 0; i < SIZE; i++) {
        #pragma HLS PIPELINE II=1
        OUT_R_out.write(OUT_R[i]);
        OUT_I_out.write(OUT_I[i]);
    }
}