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

#include <math.h>
#include "platform.h"

#ifdef USE_FEEDFORWARD

#include "build/debug.h"

#include "common/maths.h"

#include "fc/rc.h"

#include "flight/pid.h"

#include "feedforward.h"

static float setpointDelta[XYZ_AXIS_COUNT];

typedef struct laggedMovingAverageCombined_s {
     laggedMovingAverage_t filter;
     float buf[4];
} laggedMovingAverageCombined_t;

laggedMovingAverageCombined_t  setpointDeltaAvg[XYZ_AXIS_COUNT];

static float prevSetpoint[XYZ_AXIS_COUNT];
static float prevSetpointSpeed[XYZ_AXIS_COUNT];
static float prevAcceleration[XYZ_AXIS_COUNT];
static float prevRcCommandDelta[XYZ_AXIS_COUNT];
static bool prevDuplicatePacket[XYZ_AXIS_COUNT];
static uint8_t averagingCount;
static float feedforwardMaxRateLimit[XYZ_AXIS_COUNT];
static float feedforwardMaxRate[XYZ_AXIS_COUNT];

void feedforwardInit(const pidProfile_t *pidProfile) {
    const float feedforwardMaxRateScale = pidProfile->feedforward_max_rate_limit * 0.01f;
    averagingCount = pidProfile->feedforward_averaging + 1;
    for (int i = 0; i < XYZ_AXIS_COUNT; i++) {
        feedforwardMaxRate[i] = applyCurve(i, 1.0f);
        feedforwardMaxRateLimit[i] = feedforwardMaxRate[i] * feedforwardMaxRateScale;
        laggedMovingAverageInit(&setpointDeltaAvg[i].filter, averagingCount, (float *)&setpointDeltaAvg[i].buf[0]);
    }
}

FAST_CODE_NOINLINE float feedforwardApply(int axis, bool newRcFrame, feedforwardAveraging_t feedforwardAveraging) {

    if (newRcFrame) {
        float rcCommandDelta = getRcCommandDelta(axis);
        float setpoint = getRawSetpoint(axis);
        const float rxInterval = getCurrentRxRefreshRate() * 1e-6f;
        const float rxRate = 1.0f / rxInterval;
        float setpointSpeed = (setpoint - prevSetpoint[axis]) * rxRate;
        float absPrevSetpointSpeed = fabsf(prevSetpointSpeed[axis]);
        float setpointAcceleration = 0.0f;
        const float feedforwardSmoothFactor = pidGetFeedforwardSmoothFactor();
        const float feedforwardJitterFactor = pidGetFeedforwardJitterFactor();
        float feedforward;

        // calculate an attenuator from average of two most recent rcCommand deltas vs jitter threshold
        if (axis == FD_ROLL) {
            DEBUG_SET(DEBUG_FEEDFORWARD, 3, lrintf(rcCommandDelta * 100.0f)); // rcCommand packet difference = steps of 50 mean 2000 RC steps
        }
        rcCommandDelta = fabsf(rcCommandDelta);
        float boostAttenuator = 1.0f;
        if (feedforwardJitterFactor) {
            if (rcCommandDelta < feedforwardJitterFactor) {
                boostAttenuator = MAX(1.0f - ((rcCommandDelta + prevRcCommandDelta[axis]) / 2.0f) / feedforwardJitterFactor, 0.0f);
                boostAttenuator = 1.0f - boostAttenuator * boostAttenuator;
            }
        }

        const float setpointPercent = fabsf(setpoint) / feedforwardMaxRate[axis];
        float absSetpointSpeed = fabsf(setpointSpeed); // unsmoothed for kick prevention

        // interpolate setpoint if necessary
        if (rcCommandDelta == 0.0f) {
            if (prevDuplicatePacket[axis] == false) {
                // first duplicate after movement
                // don't interpolate if sticks close to centre or max
                if (setpointPercent > 0.02f && setpointPercent < 0.95f) {
                    // simple setpoint interpolation
                    setpoint = prevSetpoint[axis] + prevSetpointSpeed[axis] * rxInterval;
                    // recalculate setpointSpeed and acceleration from this new setpoint value
                    setpointSpeed = (setpoint - prevSetpoint[axis]) * rxRate;
                    setpointAcceleration = prevAcceleration[axis];
                }
            } else {
                // force to zero
                boostAttenuator *= 0.0f;
            }
            prevDuplicatePacket[axis] = true;
        } else {
            // movement!
            // zero boost for the first move after a duplicate or flat spot, but not when coming back from max
            if (prevDuplicatePacket[axis] == true && setpointPercent < 0.95f) {
                boostAttenuator *= 0.0f;
            }
            prevDuplicatePacket[axis] = false;
        }

        prevSetpoint[axis] = setpoint;

        if (axis == FD_ROLL) {
            DEBUG_SET(DEBUG_FEEDFORWARD, 0, lrintf(setpoint)); // setpoint after interpolations
        }

        // first order type smoothing for derivative
        setpointSpeed = prevSetpointSpeed[axis] + feedforwardSmoothFactor * (setpointSpeed - prevSetpointSpeed[axis]);

        if (setpointAcceleration == 0.0f) {
            // use previous acceleration for duplicate to avoid sudden drop to zero
            setpointAcceleration = setpointSpeed - prevSetpointSpeed[axis];
        }
        // second order smoothing for for acceleration
        setpointAcceleration = prevAcceleration[axis] + feedforwardSmoothFactor * (setpointAcceleration - prevAcceleration[axis]);

        prevSetpointSpeed[axis] = setpointSpeed;
        prevAcceleration[axis] = setpointAcceleration;
        prevRcCommandDelta[axis] = rcCommandDelta;

        setpointAcceleration *= pidGetDT();
        feedforward = setpointSpeed * pidGetDT();

        // calculate boost and prevent kick-back spike at max deflection
        const float feedforwardBoostFactor = pidGetFeedforwardBoostFactor();
        float boostAmount = 0.0f;
        if (feedforwardBoostFactor) {
            // allows boost when returning from max, but not when hitting max on the way up
            if (setpointPercent < 0.95f || absSetpointSpeed > 3.0f * absPrevSetpointSpeed) {
                boostAmount = feedforwardBoostFactor * setpointAcceleration;
            }
        }

        if (axis == FD_ROLL) {
            DEBUG_SET(DEBUG_FEEDFORWARD, 1, lrintf(feedforward * 100.0f)); // delta after interpolating duplicates and smoothing
            DEBUG_SET(DEBUG_FEEDFORWARD, 2, lrintf(boostAmount * 100.0f)); // boost amount after jitter reduction and smoothing
            // debug 2 is interpolated setpoint, above
            // debug 3 is rcCommand delta, above
        }

        // add attenuated boost to base feedforward
        feedforward += boostAmount * boostAttenuator;

        // apply averaging, if enabled
        if (feedforwardAveraging) {
            setpointDelta[axis] = laggedMovingAverageUpdate(&setpointDeltaAvg[axis].filter, feedforward);
        } else {
            setpointDelta[axis] = feedforward;
        }
    }
    return setpointDelta[axis];
}

FAST_CODE_NOINLINE float applyFeedforwardLimit(int axis, float value, float Kp, float currentPidSetpoint) {
    switch (axis) {
    case FD_ROLL:
        DEBUG_SET(DEBUG_FEEDFORWARD_LIMIT, 0, value);

        break;
    case FD_PITCH:
        DEBUG_SET(DEBUG_FEEDFORWARD_LIMIT, 1, value);

        break;
    }

    if (value * currentPidSetpoint > 0.0f) {
        if (fabsf(currentPidSetpoint) <= feedforwardMaxRateLimit[axis]) {
            value = constrainf(value, (-feedforwardMaxRateLimit[axis] - currentPidSetpoint) * Kp, (feedforwardMaxRateLimit[axis] - currentPidSetpoint) * Kp);
        } else {
            value = 0;
        }
    }

    if (axis == FD_ROLL) {
        DEBUG_SET(DEBUG_FEEDFORWARD_LIMIT, 2, value);
    }

    return value;
}

bool shouldApplyFeedforwardLimits(int axis)
{
    return feedforwardMaxRateLimit[axis] != 0.0f && axis < FD_YAW;
}
#endif
