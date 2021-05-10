/* Edge Impulse inferencing library
 * Copyright (c) 2021 EdgeImpulse Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _EIDSP_SPEECHPY_FEATURE_H_
#define _EIDSP_SPEECHPY_FEATURE_H_

#include <vector>
#include <stdint.h>
#include "functions.hpp"
#include "processing.hpp"
#include "../memory.hpp"

namespace ei {
namespace speechpy {
    // Actually only about half of this table is needed due to symmetry
//    float window[256] = {
//        0.01227154, 0.02454123, 0.03680722, 0.04906767, 0.06132074,
//       0.07356456, 0.08579731, 0.09801714, 0.11022221, 0.12241068,
//       0.13458071, 0.14673047, 0.15885814, 0.17096189, 0.18303989,
//       0.19509032, 0.20711138, 0.21910124, 0.23105811, 0.24298018,
//       0.25486566, 0.26671276, 0.27851969, 0.29028468, 0.30200595,
//       0.31368174, 0.32531029, 0.33688985, 0.34841868, 0.35989504,
//       0.37131719, 0.38268343, 0.39399204, 0.40524131, 0.41642956,
//       0.42755509, 0.43861624, 0.44961133, 0.46053871, 0.47139674,
//       0.48218377, 0.49289819, 0.50353838, 0.51410274, 0.52458968,
//       0.53499762, 0.54532499, 0.55557023, 0.56573181, 0.57580819,
//       0.58579786, 0.5956993 , 0.60551104, 0.61523159, 0.62485949,
//       0.63439328, 0.64383154, 0.65317284, 0.66241578, 0.67155895,
//       0.680601  , 0.68954054, 0.69837625, 0.70710678, 0.71573083,
//       0.72424708, 0.73265427, 0.74095113, 0.74913639, 0.75720885,
//       0.76516727, 0.77301045, 0.78073723, 0.78834643, 0.7958369 ,
//       0.80320753, 0.8104572 , 0.81758481, 0.8245893 , 0.83146961,
//       0.83822471, 0.84485357, 0.85135519, 0.85772861, 0.86397286,
//       0.87008699, 0.87607009, 0.88192126, 0.88763962, 0.8932243 ,
//       0.89867447, 0.90398929, 0.90916798, 0.91420976, 0.91911385,
//       0.92387953, 0.92850608, 0.9329928 , 0.93733901, 0.94154407,
//       0.94560733, 0.94952818, 0.95330604, 0.95694034, 0.96043052,
//       0.96377607, 0.96697647, 0.97003125, 0.97293995, 0.97570213,
//       0.97831737, 0.98078528, 0.98310549, 0.98527764, 0.98730142,
//       0.98917651, 0.99090264, 0.99247953, 0.99390697, 0.99518473,
//       0.99631261, 0.99729046, 0.99811811, 0.99879546, 0.99932238,
//       0.99969882, 0.9999247 , 1.        , 0.9999247 , 0.99969882,
//       0.99932238, 0.99879546, 0.99811811, 0.99729046, 0.99631261,
//       0.99518473, 0.99390697, 0.99247953, 0.99090264, 0.98917651,
//       0.98730142, 0.98527764, 0.98310549, 0.98078528, 0.97831737,
//       0.97570213, 0.97293995, 0.97003125, 0.96697647, 0.96377607,
//       0.96043052, 0.95694034, 0.95330604, 0.94952818, 0.94560733,
//       0.94154407, 0.93733901, 0.9329928 , 0.92850608, 0.92387953,
//       0.91911385, 0.91420976, 0.90916798, 0.90398929, 0.89867447,
//       0.8932243 , 0.88763962, 0.88192126, 0.87607009, 0.87008699,
//       0.86397286, 0.85772861, 0.85135519, 0.84485357, 0.83822471,
//       0.83146961, 0.8245893 , 0.81758481, 0.8104572 , 0.80320753,
//       0.7958369 , 0.78834643, 0.78073723, 0.77301045, 0.76516727,
//       0.75720885, 0.74913639, 0.74095113, 0.73265427, 0.72424708,
//       0.71573083, 0.70710678, 0.69837625, 0.68954054, 0.680601  ,
//       0.67155895, 0.66241578, 0.65317284, 0.64383154, 0.63439328,
//       0.62485949, 0.61523159, 0.60551104, 0.5956993 , 0.58579786,
//       0.57580819, 0.56573181, 0.55557023, 0.54532499, 0.53499762,
//       0.52458968, 0.51410274, 0.50353838, 0.49289819, 0.48218377,
//       0.47139674, 0.46053871, 0.44961133, 0.43861624, 0.42755509,
//       0.41642956, 0.40524131, 0.39399204, 0.38268343, 0.37131719,
//       0.35989504, 0.34841868, 0.33688985, 0.32531029, 0.31368174,
//       0.30200595, 0.29028468, 0.27851969, 0.26671276, 0.25486566,
//       0.24298018, 0.23105811, 0.21910124, 0.20711138, 0.19509032,
//       0.18303989, 0.17096189, 0.15885814, 0.14673047, 0.13458071,
//       0.12241068, 0.11022221, 0.09801714, 0.08579731, 0.07356456,
//       0.06132074, 0.04906767, 0.03680722, 0.02454123, 0.01227154,
//       0.
//    };
    float window[512] = {
      0.00613588, 0.01227154, 0.01840673, 0.02454123, 0.0306748 ,
       0.03680722, 0.04293826, 0.04906767, 0.05519524, 0.06132074,
       0.06744392, 0.07356456, 0.07968244, 0.08579731, 0.09190896,
       0.09801714, 0.10412163, 0.11022221, 0.11631863, 0.12241068,
       0.12849811, 0.13458071, 0.14065824, 0.14673047, 0.15279719,
       0.15885814, 0.16491312, 0.17096189, 0.17700422, 0.18303989,
       0.18906866, 0.19509032, 0.20110463, 0.20711138, 0.21311032,
       0.21910124, 0.22508391, 0.23105811, 0.23702361, 0.24298018,
       0.24892761, 0.25486566, 0.26079412, 0.26671276, 0.27262136,
       0.27851969, 0.28440754, 0.29028468, 0.29615089, 0.30200595,
       0.30784964, 0.31368174, 0.31950203, 0.32531029, 0.33110631,
       0.33688985, 0.34266072, 0.34841868, 0.35416353, 0.35989504,
       0.365613  , 0.37131719, 0.37700741, 0.38268343, 0.38834505,
       0.39399204, 0.3996242 , 0.40524131, 0.41084317, 0.41642956,
       0.42200027, 0.42755509, 0.43309382, 0.43861624, 0.44412214,
       0.44961133, 0.45508359, 0.46053871, 0.4659765 , 0.47139674,
       0.47679923, 0.48218377, 0.48755016, 0.49289819, 0.49822767,
       0.50353838, 0.50883014, 0.51410274, 0.51935599, 0.52458968,
       0.52980362, 0.53499762, 0.54017147, 0.54532499, 0.55045797,
       0.55557023, 0.56066158, 0.56573181, 0.57078075, 0.57580819,
       0.58081396, 0.58579786, 0.5907597 , 0.5956993 , 0.60061648,
       0.60551104, 0.61038281, 0.61523159, 0.62005721, 0.62485949,
       0.62963824, 0.63439328, 0.63912444, 0.64383154, 0.6485144 ,
       0.65317284, 0.65780669, 0.66241578, 0.66699992, 0.67155895,
       0.6760927 , 0.680601  , 0.68508367, 0.68954054, 0.69397146,
       0.69837625, 0.70275474, 0.70710678, 0.7114322 , 0.71573083,
       0.72000251, 0.72424708, 0.72846439, 0.73265427, 0.73681657,
       0.74095113, 0.74505779, 0.74913639, 0.7531868 , 0.75720885,
       0.76120239, 0.76516727, 0.76910334, 0.77301045, 0.77688847,
       0.78073723, 0.7845566 , 0.78834643, 0.79210658, 0.7958369 ,
       0.79953727, 0.80320753, 0.80684755, 0.8104572 , 0.81403633,
       0.81758481, 0.82110251, 0.8245893 , 0.82804505, 0.83146961,
       0.83486287, 0.83822471, 0.84155498, 0.84485357, 0.84812034,
       0.85135519, 0.85455799, 0.85772861, 0.86086694, 0.86397286,
       0.86704625, 0.87008699, 0.87309498, 0.87607009, 0.87901223,
       0.88192126, 0.8847971 , 0.88763962, 0.89044872, 0.8932243 ,
       0.89596625, 0.89867447, 0.90134885, 0.90398929, 0.9065957 ,
       0.90916798, 0.91170603, 0.91420976, 0.91667906, 0.91911385,
       0.92151404, 0.92387953, 0.92621024, 0.92850608, 0.93076696,
       0.9329928 , 0.93518351, 0.93733901, 0.93945922, 0.94154407,
       0.94359346, 0.94560733, 0.94758559, 0.94952818, 0.95143502,
       0.95330604, 0.95514117, 0.95694034, 0.95870347, 0.96043052,
       0.9621214 , 0.96377607, 0.96539444, 0.96697647, 0.96852209,
       0.97003125, 0.97150389, 0.97293995, 0.97433938, 0.97570213,
       0.97702814, 0.97831737, 0.97956977, 0.98078528, 0.98196387,
       0.98310549, 0.98421009, 0.98527764, 0.9863081 , 0.98730142,
       0.98825757, 0.98917651, 0.99005821, 0.99090264, 0.99170975,
       0.99247953, 0.99321195, 0.99390697, 0.99456457, 0.99518473,
       0.99576741, 0.99631261, 0.9968203 , 0.99729046, 0.99772307,
       0.99811811, 0.99847558, 0.99879546, 0.99907773, 0.99932238,
       0.99952942, 0.99969882, 0.99983058, 0.9999247 , 0.99998118,
       1.        , 0.99998118, 0.9999247 , 0.99983058, 0.99969882,
       0.99952942, 0.99932238, 0.99907773, 0.99879546, 0.99847558,
       0.99811811, 0.99772307, 0.99729046, 0.9968203 , 0.99631261,
       0.99576741, 0.99518473, 0.99456457, 0.99390697, 0.99321195,
       0.99247953, 0.99170975, 0.99090264, 0.99005821, 0.98917651,
       0.98825757, 0.98730142, 0.9863081 , 0.98527764, 0.98421009,
       0.98310549, 0.98196387, 0.98078528, 0.97956977, 0.97831737,
       0.97702814, 0.97570213, 0.97433938, 0.97293995, 0.97150389,
       0.97003125, 0.96852209, 0.96697647, 0.96539444, 0.96377607,
       0.9621214 , 0.96043052, 0.95870347, 0.95694034, 0.95514117,
       0.95330604, 0.95143502, 0.94952818, 0.94758559, 0.94560733,
       0.94359346, 0.94154407, 0.93945922, 0.93733901, 0.93518351,
       0.9329928 , 0.93076696, 0.92850608, 0.92621024, 0.92387953,
       0.92151404, 0.91911385, 0.91667906, 0.91420976, 0.91170603,
       0.90916798, 0.9065957 , 0.90398929, 0.90134885, 0.89867447,
       0.89596625, 0.8932243 , 0.89044872, 0.88763962, 0.8847971 ,
       0.88192126, 0.87901223, 0.87607009, 0.87309498, 0.87008699,
       0.86704625, 0.86397286, 0.86086694, 0.85772861, 0.85455799,
       0.85135519, 0.84812034, 0.84485357, 0.84155498, 0.83822471,
       0.83486287, 0.83146961, 0.82804505, 0.8245893 , 0.82110251,
       0.81758481, 0.81403633, 0.8104572 , 0.80684755, 0.80320753,
       0.79953727, 0.7958369 , 0.79210658, 0.78834643, 0.7845566 ,
       0.78073723, 0.77688847, 0.77301045, 0.76910334, 0.76516727,
       0.76120239, 0.75720885, 0.7531868 , 0.74913639, 0.74505779,
       0.74095113, 0.73681657, 0.73265427, 0.72846439, 0.72424708,
       0.72000251, 0.71573083, 0.7114322 , 0.70710678, 0.70275474,
       0.69837625, 0.69397146, 0.68954054, 0.68508367, 0.680601  ,
       0.6760927 , 0.67155895, 0.66699992, 0.66241578, 0.65780669,
       0.65317284, 0.6485144 , 0.64383154, 0.63912444, 0.63439328,
       0.62963824, 0.62485949, 0.62005721, 0.61523159, 0.61038281,
       0.60551104, 0.60061648, 0.5956993 , 0.5907597 , 0.58579786,
       0.58081396, 0.57580819, 0.57078075, 0.56573181, 0.56066158,
       0.55557023, 0.55045797, 0.54532499, 0.54017147, 0.53499762,
       0.52980362, 0.52458968, 0.51935599, 0.51410274, 0.50883014,
       0.50353838, 0.49822767, 0.49289819, 0.48755016, 0.48218377,
       0.47679923, 0.47139674, 0.4659765 , 0.46053871, 0.45508359,
       0.44961133, 0.44412214, 0.43861624, 0.43309382, 0.42755509,
       0.42200027, 0.41642956, 0.41084317, 0.40524131, 0.3996242 ,
       0.39399204, 0.38834505, 0.38268343, 0.37700741, 0.37131719,
       0.365613  , 0.35989504, 0.35416353, 0.34841868, 0.34266072,
       0.33688985, 0.33110631, 0.32531029, 0.31950203, 0.31368174,
       0.30784964, 0.30200595, 0.29615089, 0.29028468, 0.28440754,
       0.27851969, 0.27262136, 0.26671276, 0.26079412, 0.25486566,
       0.24892761, 0.24298018, 0.23702361, 0.23105811, 0.22508391,
       0.21910124, 0.21311032, 0.20711138, 0.20110463, 0.19509032,
       0.18906866, 0.18303989, 0.17700422, 0.17096189, 0.16491312,
       0.15885814, 0.15279719, 0.14673047, 0.14065824, 0.13458071,
       0.12849811, 0.12241068, 0.11631863, 0.11022221, 0.10412163,
       0.09801714, 0.09190896, 0.08579731, 0.07968244, 0.07356456,
       0.06744392, 0.06132074, 0.05519524, 0.04906767, 0.04293826,
       0.03680722, 0.0306748 , 0.02454123, 0.01840673, 0.01227154,
       0.00613588, 0.  
    };
class feature {
public:
    /**
     * Compute the Mel-filterbanks. Each filter will be stored in one rows.
     * The columns correspond to fft bins.
     *
     * @param filterbanks Matrix of size num_filter * coefficients
     * @param num_filter the number of filters in the filterbank
     * @param coefficients (fftpoints//2 + 1)
     * @param sampling_freq  the samplerate of the signal we are working
     *                       with. It affects mel spacing.
     * @param low_freq lowest band edge of mel filters, default 0 Hz
     * @param high_freq highest band edge of mel filters, default samplerate / 2
     * @param output_transposed If set to true this will transpose the matrix (memory efficient).
     *                          This is more efficient than calling this function and then transposing
     *                          as the latter requires the filterbank to be allocated twice (for a short while).
     * @returns EIDSP_OK if OK
     */
    static int filterbanks(
#if EIDSP_QUANTIZE_FILTERBANK
        quantized_matrix_t *filterbanks,
#else
        matrix_t *filterbanks,
#endif
        uint16_t num_filter, int coefficients, uint32_t sampling_freq,
        uint32_t low_freq, uint32_t high_freq,
        bool output_transposed = false
        )
    {
        const size_t mels_mem_size = (num_filter + 2) * sizeof(float);
        const size_t hertz_mem_size = (num_filter + 2) * sizeof(float);
        const size_t freq_index_mem_size = (num_filter + 2) * sizeof(int);

        float *mels = (float*)ei_dsp_malloc(mels_mem_size);
        if (!mels) {
            EIDSP_ERR(EIDSP_OUT_OF_MEM);
        }

        if (filterbanks->rows != num_filter || filterbanks->cols != static_cast<uint32_t>(coefficients)) {
            EIDSP_ERR(EIDSP_MATRIX_SIZE_MISMATCH);
        }

#if EIDSP_QUANTIZE_FILTERBANK
        memset(filterbanks->buffer, 0, filterbanks->rows * filterbanks->cols * sizeof(uint8_t));
#else
        memset(filterbanks->buffer, 0, filterbanks->rows * filterbanks->cols * sizeof(float));
#endif

        // Computing the Mel filterbank
        // converting the upper and lower frequencies to Mels.
        // num_filter + 2 is because for num_filter filterbanks we need
        // num_filter+2 point.
        numpy::linspace(
            functions::frequency_to_mel(static_cast<float>(low_freq)),
            functions::frequency_to_mel(static_cast<float>(high_freq)),
            num_filter + 2,
            mels);

        // we should convert Mels back to Hertz because the start and end-points
        // should be at the desired frequencies.
        float *hertz = (float*)ei_dsp_malloc(hertz_mem_size);
        if (!hertz) {
            ei_dsp_free(mels, mels_mem_size);
            EIDSP_ERR(EIDSP_OUT_OF_MEM);
        }
        for (uint16_t ix = 0; ix < num_filter + 2; ix++) {
            hertz[ix] = functions::mel_to_frequency(mels[ix]);
            if (hertz[ix] < low_freq) {
                hertz[ix] = low_freq;
            }
            if (hertz[ix] > high_freq) {
                hertz[ix] = high_freq;
            }

            // here is a really annoying bug in Speechpy which calculates the frequency index wrong for the last bucket
            // the last 'hertz' value is not 8,000 (with sampling rate 16,000) but 7,999.999999
            // thus calculating the bucket to 64, not 65.
            // we're adjusting this here a tiny bit to ensure we have the same result
            if (ix == num_filter + 2 - 1) {
                hertz[ix] -= 0.001;
            }
        }
        ei_dsp_free(mels, mels_mem_size);

        // The frequency resolution required to put filters at the
        // exact points calculated above should be extracted.
        //  So we should round those frequencies to the closest FFT bin.
        int *freq_index = (int*)ei_dsp_malloc(freq_index_mem_size);
        if (!freq_index) {
            ei_dsp_free(hertz, hertz_mem_size);
            EIDSP_ERR(EIDSP_OUT_OF_MEM);
        }
        for (uint16_t ix = 0; ix < num_filter + 2; ix++) {
            freq_index[ix] = static_cast<int>(floor((coefficients + 1) * hertz[ix] / sampling_freq));
        }
        ei_dsp_free(hertz, hertz_mem_size);

        for (size_t i = 0; i < num_filter; i++) {
            int left = freq_index[i];
            int middle = freq_index[i + 1];
            int right = freq_index[i + 2];

            EI_DSP_MATRIX(z, 1, (right - left + 1));
            if (!z.buffer) {
                ei_dsp_free(freq_index, freq_index_mem_size);
                EIDSP_ERR(EIDSP_OUT_OF_MEM);
            }
            numpy::linspace(left, right, (right - left + 1), z.buffer);
            functions::triangle(z.buffer, (right - left + 1), left, middle, right);

            // so... z now contains some values that we need to overwrite in the filterbank
            for (int zx = 0; zx < (right - left + 1); zx++) {
                size_t index = (i * filterbanks->cols) + (left + zx);

                if (output_transposed) {
                    index = ((left + zx) * filterbanks->rows) + i;
                }

#if EIDSP_QUANTIZE_FILTERBANK
                filterbanks->buffer[index] = numpy::quantize_zero_one(z.buffer[zx]);
#else
                filterbanks->buffer[index] = z.buffer[zx];
#endif
            }
        }

        if (output_transposed) {
            uint16_t r = filterbanks->rows;
            filterbanks->rows = filterbanks->cols;
            filterbanks->cols = r;
        }

        ei_dsp_free(freq_index, freq_index_mem_size);

        return EIDSP_OK;
    }

    /**
     * Compute Mel-filterbank energy features from an audio signal.
     * @param out_features Use `calculate_mfe_buffer_size` to allocate the right matrix.
     * @param out_energies A matrix in the form of Mx1 where M is the rows from `calculate_mfe_buffer_size`
     * @param signal: audio signal structure with functions to retrieve data from a signal
     * @param sampling_frequency (int): the sampling frequency of the signal
     *     we are working with.
     * @param frame_length (float): the length of each frame in seconds.
     *     Default is 0.020s
     * @param frame_stride (float): the step between successive frames in seconds.
     *     Default is 0.02s (means no overlap)
     * @param num_filters (int): the number of filters in the filterbank,
     *     default 40.
     * @param fft_length (int): number of FFT points. Default is 512.
     * @param low_frequency (int): lowest band edge of mel filters.
     *     In Hz, default is 0.
     * @param high_frequency (int): highest band edge of mel filters.
     *     In Hz, default is samplerate/2
     * @EIDSP_OK if OK
     */
    static int mfe(matrix_t *out_features, matrix_t *out_energies,
        signal_t *signal,
        uint32_t sampling_frequency,
        float frame_length, float frame_stride, uint16_t num_filters,
        uint16_t fft_length, uint32_t low_frequency, uint32_t high_frequency,
        uint16_t version
        )
    {
        int ret = 0;

        if (high_frequency == 0) {
            high_frequency = sampling_frequency / 2;
        }

        if (low_frequency == 0) {
            low_frequency = 300;
        }

        stack_frames_info_t stack_frame_info = { 0 };
        stack_frame_info.signal = signal;

        ret = processing::stack_frames(
            &stack_frame_info,
            sampling_frequency,
            frame_length,
            frame_stride,
            false,
            version
        );
        if (ret != 0) {
            EIDSP_ERR(ret);
        }

        if (stack_frame_info.frame_ixs->size() != out_features->rows) {
            EIDSP_ERR(EIDSP_MATRIX_SIZE_MISMATCH);
        }

        if (num_filters != out_features->cols) {
            EIDSP_ERR(EIDSP_MATRIX_SIZE_MISMATCH);
        }

        if (stack_frame_info.frame_ixs->size() != out_energies->rows || out_energies->cols != 1) {
            EIDSP_ERR(EIDSP_MATRIX_SIZE_MISMATCH);
        }

        for (uint32_t i = 0; i < out_features->rows * out_features->cols; i++) {
            *(out_features->buffer + i) = 0;
        }

        uint16_t coefficients = fft_length / 2 + 1;

        // calculate the filterbanks first... preferably I would want to do the matrix multiplications
        // whenever they happen, but OK...
#if EIDSP_QUANTIZE_FILTERBANK
        EI_DSP_QUANTIZED_MATRIX(filterbanks, num_filters, coefficients, &numpy::dequantize_zero_one);
#else
        EI_DSP_MATRIX(filterbanks, num_filters, coefficients);
#endif
        if (!filterbanks.buffer) {
            EIDSP_ERR(EIDSP_OUT_OF_MEM);
        }

        ret = feature::filterbanks(
            &filterbanks, num_filters, coefficients, sampling_frequency, low_frequency, high_frequency, true);
        if (ret != 0) {
            EIDSP_ERR(ret);
        }
        for (size_t ix = 0; ix < stack_frame_info.frame_ixs->size(); ix++) {
            size_t power_spectrum_frame_size = (fft_length / 2 + 1);

            EI_DSP_MATRIX(power_spectrum_frame, 1, power_spectrum_frame_size);
            if (!power_spectrum_frame.buffer) {
                EIDSP_ERR(EIDSP_OUT_OF_MEM);
            }

            // get signal data from the audio file
            EI_DSP_MATRIX(signal_frame, 1, stack_frame_info.frame_length);

            // don't read outside of the audio buffer... we'll automatically zero pad then
            size_t signal_offset = stack_frame_info.frame_ixs->at(ix);
            size_t signal_length = stack_frame_info.frame_length;
            if (signal_offset + signal_length > stack_frame_info.signal->total_length) {
                signal_length = signal_length -
                    (stack_frame_info.signal->total_length - (signal_offset + signal_length));
            }

            ret = stack_frame_info.signal->get_data(
                signal_offset,
                signal_length,
                signal_frame.buffer
            );
            if (ret != 0) {
                EIDSP_ERR(ret);
            }
            
            for (int i=0; i < signal_length; i++) {
                signal_frame.buffer[i] *= window[i];
            }

            ret = processing::power_spectrum(
                signal_frame.buffer,
                stack_frame_info.frame_length,
                power_spectrum_frame.buffer,
                power_spectrum_frame_size,
                fft_length
            );

            if (ret != 0) {
                EIDSP_ERR(ret);
            }

            float energy = numpy::sum(power_spectrum_frame.buffer, power_spectrum_frame_size);
            if (energy == 0) {
                energy = FLT_EPSILON;
            }

            out_energies->buffer[ix] = energy;

            // calculate the out_features directly here
            ret = numpy::dot_by_row(
                ix,
                power_spectrum_frame.buffer,
                power_spectrum_frame_size,
                &filterbanks,
                out_features
            );

            if (ret != 0) {
                EIDSP_ERR(ret);
            }
        }

        functions::zero_handling(out_features);
        
        // *TJG* Apply log to MFE
        ret = numpy::log(out_features);
        if (ret != EIDSP_OK) {
            EIDSP_ERR(ret);
        }        

        return EIDSP_OK;
    }

    /**
     * Compute spectrogram from a sensor signal.
     * @param out_features Use `calculate_mfe_buffer_size` to allocate the right matrix.
     * @param signal: audio signal structure with functions to retrieve data from a signal
     * @param sampling_frequency (int): the sampling frequency of the signal
     *     we are working with.
     * @param frame_length (float): the length of each frame in seconds.
     *     Default is 0.020s
     * @param frame_stride (float): the step between successive frames in seconds.
     *     Default is 0.02s (means no overlap)
     * @param fft_length (int): number of FFT points. Default is 512.
     * @EIDSP_OK if OK
     */
    static int spectrogram(matrix_t *out_features,
        signal_t *signal, uint32_t sampling_frequency,
        float frame_length, float frame_stride, uint16_t fft_length,
        uint16_t version
        )
    {
        int ret = 0;

        stack_frames_info_t stack_frame_info = { 0 };
        stack_frame_info.signal = signal;

        ret = processing::stack_frames(
            &stack_frame_info,
            sampling_frequency,
            frame_length,
            frame_stride,
            false,
            version
        );
        if (ret != 0) {
            EIDSP_ERR(ret);
        }

        if (stack_frame_info.frame_ixs->size() != out_features->rows) {
            EIDSP_ERR(EIDSP_MATRIX_SIZE_MISMATCH);
        }

        uint16_t coefficients = fft_length / 2 + 1;

        if (coefficients != out_features->cols) {
            EIDSP_ERR(EIDSP_MATRIX_SIZE_MISMATCH);
        }

        for (uint32_t i = 0; i < out_features->rows * out_features->cols; i++) {
            *(out_features->buffer + i) = 0;
        }

        for (size_t ix = 0; ix < stack_frame_info.frame_ixs->size(); ix++) {
            // get signal data from the audio file
            EI_DSP_MATRIX(signal_frame, 1, stack_frame_info.frame_length);

            // don't read outside of the audio buffer... we'll automatically zero pad then
            size_t signal_offset = stack_frame_info.frame_ixs->at(ix);
            size_t signal_length = stack_frame_info.frame_length;
            if (signal_offset + signal_length > stack_frame_info.signal->total_length) {
                signal_length = signal_length -
                    (stack_frame_info.signal->total_length - (signal_offset + signal_length));
            }

            ret = stack_frame_info.signal->get_data(
                signal_offset,
                signal_length,
                signal_frame.buffer
            );
            if (ret != 0) {
                EIDSP_ERR(ret);
            }

            ret = processing::power_spectrum(
                signal_frame.buffer,
                stack_frame_info.frame_length,
                out_features->buffer + (ix * coefficients),
                coefficients,
                fft_length
            );

            if (ret != 0) {
                EIDSP_ERR(ret);
            }
        }

        functions::zero_handling(out_features);

        return EIDSP_OK;
    }

    /**
     * Calculate the buffer size for MFE
     * @param signal_length: Length of the signal.
     * @param sampling_frequency (int): The sampling frequency of the signal.
     * @param frame_length (float): The length of the frame in second.
     * @param frame_stride (float): The stride between frames.
     * @param num_filters
     */
    static matrix_size_t calculate_mfe_buffer_size(
        size_t signal_length,
        uint32_t sampling_frequency,
        float frame_length, float frame_stride, uint16_t num_filters,
        uint16_t version)
    {
        uint16_t rows = processing::calculate_no_of_stack_frames(
            signal_length,
            sampling_frequency,
            frame_length,
            frame_stride,
            false,
            version);
        uint16_t cols = num_filters;

        matrix_size_t size_matrix;
        size_matrix.rows = rows;
        size_matrix.cols = cols;
        return size_matrix;
    }

    /**
     * Compute MFCC features from an audio signal.
     * @param out_features Use `calculate_mfcc_buffer_size` to allocate the right matrix.
     * @param signal: audio signal structure from which to compute features.
     *     has functions to retrieve data from a signal lazily.
     * @param sampling_frequency (int): the sampling frequency of the signal
     *     we are working with.
     * @param frame_length (float): the length of each frame in seconds.
     *     Default is 0.020s
     * @param frame_stride (float): the step between successive frames in seconds.
     *     Default is 0.01s (means no overlap)
     * @param num_cepstral (int): Number of cepstral coefficients.
     * @param num_filters (int): the number of filters in the filterbank,
     *     default 40.
     * @param fft_length (int): number of FFT points. Default is 512.
     * @param low_frequency (int): lowest band edge of mel filters.
     *     In Hz, default is 0.
     * @param high_frequency (int): highest band edge of mel filters.
     *     In Hz, default is samplerate/2
     * @param dc_elimination Whether the first dc component should
     *     be eliminated or not.
     * @returns 0 if OK
     */
    static int mfcc(matrix_t *out_features, signal_t *signal,
        uint32_t sampling_frequency, float frame_length, float frame_stride,
        uint8_t num_cepstral, uint16_t num_filters, uint16_t fft_length,
        uint32_t low_frequency, uint32_t high_frequency, bool dc_elimination,
        uint16_t version)
    {
        if (out_features->cols != num_cepstral) {
            EIDSP_ERR(EIDSP_MATRIX_SIZE_MISMATCH);
        }

        matrix_size_t mfe_matrix_size =
            calculate_mfe_buffer_size(
                signal->total_length,
                sampling_frequency,
                frame_length,
                frame_stride,
                num_filters,
                version);

        if (out_features->rows != mfe_matrix_size.rows) {
            EIDSP_ERR(EIDSP_MATRIX_SIZE_MISMATCH);
        }

        int ret = EIDSP_OK;

        // allocate some memory for the MFE result
        EI_DSP_MATRIX(features_matrix, mfe_matrix_size.rows, mfe_matrix_size.cols);
        if (!features_matrix.buffer) {
            EIDSP_ERR(EIDSP_OUT_OF_MEM);
        }

        EI_DSP_MATRIX(energy_matrix, mfe_matrix_size.rows, 1);
        if (!energy_matrix.buffer) {
            EIDSP_ERR(EIDSP_OUT_OF_MEM);
        }

        ret = mfe(&features_matrix, &energy_matrix, signal,
            sampling_frequency, frame_length, frame_stride, num_filters, fft_length,
            low_frequency, high_frequency, version);
        if (ret != EIDSP_OK) {
            EIDSP_ERR(ret);
        }

        // ok... now we need to calculate the MFCC from this...
        // first do log() over all features...
        ret = numpy::log(&features_matrix);
        if (ret != EIDSP_OK) {
            EIDSP_ERR(ret);
        }

        // now do DST type 2
        ret = numpy::dct2(&features_matrix, DCT_NORMALIZATION_ORTHO);
        if (ret != EIDSP_OK) {
            EIDSP_ERR(ret);
        }

        // replace first cepstral coefficient with log of frame energy for DC elimination
        if (dc_elimination) {
            for (size_t row = 0; row < features_matrix.rows; row++) {
                features_matrix.buffer[row * features_matrix.cols] = numpy::log(energy_matrix.buffer[row]);
            }
        }

        // copy to the output...
        for (size_t row = 0; row < features_matrix.rows; row++) {
            for(int i = 0; i < num_cepstral; i++) {
                *(out_features->buffer + (num_cepstral * row) + i) = *(features_matrix.buffer + (features_matrix.cols * row) + i);
            }
        }

        return EIDSP_OK;
    }

    /**
     * Calculate the buffer size for MFCC
     * @param signal_length: Length of the signal.
     * @param sampling_frequency (int): The sampling frequency of the signal.
     * @param frame_length (float): The length of the frame in second.
     * @param frame_stride (float): The stride between frames.
     * @param num_cepstral
     */
    static matrix_size_t calculate_mfcc_buffer_size(
        size_t signal_length,
        uint32_t sampling_frequency,
        float frame_length, float frame_stride, uint16_t num_cepstral,
        uint16_t version)
    {
        uint16_t rows = processing::calculate_no_of_stack_frames(
            signal_length,
            sampling_frequency,
            frame_length,
            frame_stride,
            false,
            version);
        uint16_t cols = num_cepstral;

        matrix_size_t size_matrix;
        size_matrix.rows = rows;
        size_matrix.cols = cols;
        return size_matrix;
    }
};

} // namespace speechpy
} // namespace ei

#endif // _EIDSP_SPEECHPY_FEATURE_H_
