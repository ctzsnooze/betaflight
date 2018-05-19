/*
 * This file is part of Cleanflight and Betaflight.
 *
 * Cleanflight and Betaflight are free software. You can redistribute
 * this software and/or modify this software under the terms of the
 * GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Cleanflight and Betaflight are distributed in the hope that they
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.
 *
 * If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "platform.h"

#include "common/filter.h"
#include "common/maths.h"
#include "common/utils.h"

#define M_LN2_FLOAT 0.69314718055994530942f
#define M_PI_FLOAT  3.14159265358979323846f
#define BIQUAD_BANDWIDTH 1.9f     /* bandwidth in octaves */
#define BIQUAD_Q 1.0f / sqrtf(2.0f)     /* quality factor - 2nd order butterworth*/

// NULL filter

FAST_CODE float nullFilterApply(filter_t *filter, float input)
{
    UNUSED(filter);
    return input;
}


// PT1 Low Pass filter

float pt1FilterGain(uint16_t f_cut, float dT)
{
    float RC = 1 / ( 2 * M_PI_FLOAT * f_cut);
    return dT / (RC + dT);
}

void pt1FilterInit(pt1Filter_t *filter, float k)
{
    filter->k = k;
}

FAST_CODE float pt1FilterApply(pt1Filter_t *filter, float input)
{
    filter->state = filter->state + filter->k * (input - filter->state);
    return filter->state;
}

// Slew filter with limit

void slewFilterInit(slewFilter_t *filter, float slewLimit, float threshold)
{
    filter->state = 0.0f;
    filter->slewLimit = slewLimit;
    filter->threshold = threshold;
}

FAST_CODE float slewFilterApply(slewFilter_t *filter, float input)
{
    if (filter->state >= filter->threshold) {
        if (input >= filter->state - filter->slewLimit) {
            filter->state = input;
        }
    } else if (filter->state <= -filter->threshold) {
        if (input <= filter->state + filter->slewLimit) {
            filter->state = input;
        }
    } else {
        filter->state = input;
    }
    return filter->state;
}

// get notch filter Q given center frequency (f0) and lower cutoff frequency (f1)
// Q = f0 / (f2 - f1) ; f2 = f0^2 / f1
float filterGetNotchQ(float centerFreq, float cutoffFreq) {
    return centerFreq * cutoffFreq / (centerFreq * centerFreq - cutoffFreq * cutoffFreq);
}

/* sets up a biquad Filter */
void biquadFilterInitLPF(biquadFilter_t *filter, float filterFreq, uint32_t refreshRate)
{
    biquadFilterInit(filter, filterFreq, refreshRate, BIQUAD_Q, FILTER_LPF);
}

// Q coefficients for Butterworth filter of given order. Implicit odd section coefficient of 0.5
// coefficients should be in ascending order
// generated by http://www.earlevel.com/main/2016/09/29/cascading-filters/
static const float butterworthLpfQ[] = {
#if BIQUAD_LPF_ORDER_MAX >= 2
    0.70710678,                         // 2nd
#endif
#if BIQUAD_LPF_ORDER_MAX >= 3
    1.0,                                // 3rd
#endif
#if BIQUAD_LPF_ORDER_MAX >= 4
    0.54119610, 1.3065630,              // 4th
#endif
#if BIQUAD_LPF_ORDER_MAX >= 5
    0.61803399, 1.6180340,              // 5th
#endif
#if BIQUAD_LPF_ORDER_MAX >= 6
    0.51763809, 0.70710678, 1.9318517   // 6th
#endif
};
#define BUTTERWORTH_QINDEX(o) (((o) - 1) * ((o) - 1) / 4)

// make sure that we have correct number of coefficients
STATIC_ASSERT(BUTTERWORTH_QINDEX(BIQUAD_LPF_ORDER_MAX + 1) == ARRAYLEN(butterworthLpfQ), butterworthLpfQ_mismatch);

// setup cascade of biquad sections for Butterworth LPF filter
// sections is pointer to array of biquad sections, it must be long enough
// return number of sections used
int biquadFilterLpfCascadeInit(biquadFilter_t *sections, int order, float filterFreq, uint32_t refreshRate)
{
    biquadFilter_t *section = sections;
    // single-pole section first
    if (order % 2 == 1) {
        biquadFilterInit(section, filterFreq, refreshRate, 0.5, FILTER_LPF1);
        section++;
    }
    const float *Qptr =  butterworthLpfQ + BUTTERWORTH_QINDEX(order);  // first coefficient for given order
    // 2 poles per section
    for (int i = order; i >= 2; i -= 2) {
        biquadFilterInit(section, filterFreq, refreshRate, *Qptr, FILTER_LPF);
        Qptr++; section++;
    }
    return section - sections;
}

void biquadFilterInit(biquadFilter_t *filter, float filterFreq, uint32_t refreshRate, float Q, biquadFilterType_e filterType)
{
    // setup variables
    const float omega = 2.0f * M_PI_FLOAT * filterFreq * refreshRate * 0.000001f;
    const float sn = sin_approx(omega);
    const float cs = cos_approx(omega);
    const float alpha = sn / (2.0f * Q);

    float b0 = 0, b1 = 0, b2 = 0, a0 = 0, a1 = 0, a2 = 0;

    switch (filterType) {
    case FILTER_LPF:
        // 2nd order Butterworth (with Q=1/sqrt(2)) / Butterworth biquad section with Q
        // described in http://www.ti.com/lit/an/slaa447/slaa447.pdf
        b0 = (1 - cs) * 0.5f;
        b1 = 1 - cs;
        b2 = (1 - cs) * 0.5f;
        a0 = 1 + alpha;
        a1 = -2 * cs;
        a2 = 1 - alpha;
        break;
    case FILTER_LPF1: {
        // 1st order Butterworth, H(s) = 1 / (1 + s),
        // transformed according to http://www.iowahills.com/A4IIRBilinearTransform.html
        const float T = 2.0f * sn / (cs + 1.0f);   // T = 2 * tan(omega / 2)
        b0 = T;
        b1 = T;
        b2 = 0;
        a0 = T + 2.0f;
        a1 = T - 2.0f;
        a2 = 0;
        break;
    }
    case FILTER_NOTCH:
        b0 =  1;
        b1 = -2 * cs;
        b2 =  1;
        a0 =  1 + alpha;
        a1 = -2 * cs;
        a2 =  1 - alpha;
        break;
    case FILTER_BPF:
        b0 = alpha;
        b1 = 0;
        b2 = -alpha;
        a0 = 1 + alpha;
        a1 = -2 * cs;
        a2 = 1 - alpha;
        break;
    }

    // precompute the coefficients
    filter->b0 = b0 / a0;
    filter->b1 = b1 / a0;
    filter->b2 = b2 / a0;
    filter->a1 = a1 / a0;
    filter->a2 = a2 / a0;

    // zero initial samples
    filter->x1 = filter->x2 = 0;
    filter->y1 = filter->y2 = 0;
}

FAST_CODE void biquadFilterUpdate(biquadFilter_t *filter, float filterFreq, uint32_t refreshRate, float Q, biquadFilterType_e filterType)
{
    // backup state
    float x1 = filter->x1;
    float x2 = filter->x2;
    float y1 = filter->y1;
    float y2 = filter->y2;

    biquadFilterInit(filter, filterFreq, refreshRate, Q, filterType);

    // restore state
    filter->x1 = x1;
    filter->x2 = x2;
    filter->y1 = y1;
    filter->y2 = y2;
}

/* Computes a biquadFilter_t filter on a sample (slightly less precise than df2 but works in dynamic mode) */
FAST_CODE float biquadFilterApplyDF1(biquadFilter_t *filter, float input)
{
    /* compute result */
    const float result = filter->b0 * input + filter->b1 * filter->x1 + filter->b2 * filter->x2 - filter->a1 * filter->y1 - filter->a2 * filter->y2;

    /* shift x1 to x2, input to x1 */
    filter->x2 = filter->x1;
    filter->x1 = input;

    /* shift y1 to y2, result to y1 */
    filter->y2 = filter->y1;
    filter->y1 = result;

    return result;
}

/* Computes a biquadFilter_t filter in direct form 2 on a sample (higher precision but can't handle changes in coefficients */
FAST_CODE float biquadFilterApply(biquadFilter_t *filter, float input)
{
    const float result = filter->b0 * input + filter->x1;
    filter->x1 = filter->b1 * input - filter->a1 * result + filter->x2;
    filter->x2 = filter->b2 * input - filter->a2 * result;
    return result;
}

FAST_CODE float biquadCascadeFilterApply(biquadFilterCascade_t *filter, float input)
{
    for (int i = 0; i < filter->sections; i++ ) {
        input = biquadFilterApply(&filter->biquad[i], input);
    }
    return input;
}

/*
 * FIR filter
 */
void firFilterInit2(firFilter_t *filter, float *buf, uint8_t bufLength, const float *coeffs, uint8_t coeffsLength)
{
    filter->buf = buf;
    filter->bufLength = bufLength;
    filter->coeffs = coeffs;
    filter->coeffsLength = coeffsLength;
    filter->movingSum = 0.0f;
    filter->index = 0;
    filter->count = 0;
    memset(filter->buf, 0, sizeof(float) * filter->bufLength);
}

/*
 * FIR filter initialisation
 * If the FIR filter is just to be used for averaging, then coeffs can be set to NULL
 */
void firFilterInit(firFilter_t *filter, float *buf, uint8_t bufLength, const float *coeffs)
{
    firFilterInit2(filter, buf, bufLength, coeffs, bufLength);
}

void firFilterUpdate(firFilter_t *filter, float input)
{
    filter->buf[filter->index++] = input; // index is at the first empty buffer positon
    if (filter->index >= filter->bufLength) {
        filter->index = 0;
    }
}

/*
 * Update FIR filter maintaining a moving sum for quick moving average computation
 */
void firFilterUpdateAverage(firFilter_t *filter, float input)
{
    filter->movingSum += input; // sum of the last <count> items, to allow quick moving average computation
    filter->movingSum -=  filter->buf[filter->index]; // subtract the value that "drops off" the end of the moving sum
    filter->buf[filter->index++] = input; // index is at the first empty buffer positon
    if (filter->index >= filter->bufLength) {
        filter->index = 0;
    }
    if (filter->count < filter->bufLength) {
        ++filter->count;
    }
}

FAST_CODE float firFilterApply(const firFilter_t *filter)
{
    float ret = 0.0f;
    int ii = 0;
    int index;
    for (index = filter->index - 1; index >= 0; ++ii, --index) {
        ret += filter->coeffs[ii] * filter->buf[index];
    }
    for (index = filter->bufLength - 1; ii < filter->coeffsLength; ++ii, --index) {
        ret += filter->coeffs[ii] * filter->buf[index];
    }
    return ret;
}

FAST_CODE float firFilterUpdateAndApply(firFilter_t *filter, float input)
{
    firFilterUpdate(filter, input);
    return firFilterApply(filter);
}

/*
 * Returns average of the last <count> items.
 */
float firFilterCalcPartialAverage(const firFilter_t *filter, uint8_t count)
{
    float ret = 0.0f;
    int index = filter->index;
    for (int ii = 0; ii < filter->coeffsLength; ++ii) {
        --index;
        if (index < 0) {
            index = filter->bufLength - 1;
        }
        ret += filter->buf[index];
    }
    return ret / count;
}

float firFilterCalcMovingAverage(const firFilter_t *filter)
{
    return filter->movingSum / filter->count;
}

float firFilterLastInput(const firFilter_t *filter)
{
    // filter->index points to next empty item in buffer
    const int index = filter->index == 0 ? filter->bufLength - 1 : filter->index - 1;
    return filter->buf[index];
}

#if defined(USE_FIR_FILTER_DENOISE)
void firFilterDenoiseInit(firFilterDenoise_t *filter, uint8_t gyroSoftLpfHz, uint16_t targetLooptime)
{
    memset(filter, 0, sizeof(firFilterDenoise_t));
    filter->targetCount = constrain(lrintf((1.0f / (0.000001f * (float)targetLooptime)) / gyroSoftLpfHz), 1, MAX_FIR_DENOISE_WINDOW_SIZE);
}

// prototype function for denoising of signal by dynamic moving average. Mainly for test purposes
float firFilterDenoiseUpdate(firFilterDenoise_t *filter, float input)
{
    filter->state[filter->index] = input;
    filter->movingSum += filter->state[filter->index++];
    if (filter->index == filter->targetCount) {
        filter->index = 0;
    }
    filter->movingSum -= filter->state[filter->index];

    if (filter->targetCount >= filter->filledCount) {
        return filter->movingSum / filter->targetCount;
    } else {
        return filter->movingSum / ++filter->filledCount + 1;
    }
}
#endif

// ledvinap's proposed RC+FIR2 Biquad-- 1st order IIR, RC filter k
void biquadRCFIR2FilterInit(biquadFilter_t *filter, float k)
{
    filter->b0 = k / 2;
    filter->b1 = k / 2;
    filter->b2 = 0;
    filter->a1 = -(1 - k);
    filter->a2 = 0;
}

// rs2k's fast "kalman" filter
void fastKalmanInit(fastKalman_t *filter, float k)
{
    filter->x     = 0.0f;  // set initial value, can be zero if unknown
    filter->lastX = 0.0f;  // set initial value, can be zero if unknown
    filter->k     = k;     // "kalman" gain - half of RC coefficient, calculated externally
}

FAST_CODE float fastKalmanUpdate(fastKalman_t *filter, float input)
{
    filter->x += (filter->x - filter->lastX);
    filter->lastX = filter->x;
    filter->x += filter->k * (input - filter->x);
    return filter->x;
}

void lmaSmoothingInit(laggedMovingAverage_t *filter, uint8_t windowSize, float weight)
{
    filter->movingWindowIndex = 0;
    filter->windowSize = windowSize;
    filter->weight = weight;
}

FAST_CODE float lmaSmoothingUpdate(laggedMovingAverage_t *filter, float input)
{

    filter->movingSum -= filter->buf[filter->movingWindowIndex];
    filter->buf[filter->movingWindowIndex] = input;
    filter->movingSum += input;

    if (++filter->movingWindowIndex == filter->windowSize) {
        filter->movingWindowIndex = 0;
    }

    return input + (((filter->movingSum  / filter->windowSize) - input) * filter->weight);
}

void smaSmoothingInit(simpleMovingAverage_t *filter, uint8_t windowSize)
{
    filter->movingWindowIndex = 0;
    filter->windowSize = windowSize;
}

FAST_CODE float smaSmoothingUpdate(simpleMovingAverage_t *filter, float input)
{

    filter->movingSum -= filter->buf[filter->movingWindowIndex];
    filter->buf[filter->movingWindowIndex] = input;
    filter->movingSum += input;

    if (++filter->movingWindowIndex == filter->windowSize) {
        filter->movingWindowIndex = 0;
    }

    return filter->movingSum / filter->windowSize;
}

