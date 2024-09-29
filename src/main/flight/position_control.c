/*
 * This file is part of Betaflight.
 *
 * Betaflight is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Betaflight is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Betaflight. If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "platform.h"

#include "build/debug.h"

#include "common/maths.h"
#include "flight/pid.h" // for pidCoefficient_t
#include "flight/imu.h" // for getCosTiltAngle
#include "rx/rx.h" // for PWM Ranges

#include "position.h"
#include "position_control.h"

#define ALTITUDE_P_SCALE  0.01f
#define ALTITUDE_I_SCALE  0.003f
#define ALTITUDE_D_SCALE  0.01f
#define ALTITUDE_F_SCALE  0.01f

static pidCoefficient_t altitudePidCoeffs;
static float altitudeI = 0.0f;
static float throttleOut = 0.0f;

void positionControlInit(void)
{
    altitudePidCoeffs.Kp = positionControlConfig()->altitude_P * ALTITUDE_P_SCALE;
    altitudePidCoeffs.Ki = positionControlConfig()->altitude_I * ALTITUDE_I_SCALE;
    altitudePidCoeffs.Kd = positionControlConfig()->altitude_D * ALTITUDE_D_SCALE;
    altitudePidCoeffs.Kf = positionControlConfig()->altitude_F * ALTITUDE_F_SCALE;
}

const pidCoefficient_t *getAltitudePidCoeffs(void)
{
    return &altitudePidCoeffs;
}

void resetAltitudeControl (void) {
    altitudeI = 0.0f;
}

void altitudeControl(float targetAltitudeCm, float taskIntervalS, float verticalVelocity, float targetAltitudeStep) {

    float altitudeErrorCm = targetAltitudeCm - getAltitudeCm();
    const float altitudeP = altitudeErrorCm * altitudePidCoeffs.Kp;

    // reduce the iTerm gain for errors greater than 2m, otherwise it winds up too much
    const float largeError = 200.0f; // 2m
    const float itermRelax = (fabsf(altitudeErrorCm) < largeError) ? 1.0f : 0.1f;
    altitudeI += altitudeErrorCm * altitudePidCoeffs.Ki * itermRelax * taskIntervalS;
    const float iTermLimit = 200.0f;
    altitudeI = constrainf(altitudeI, -iTermLimit, iTermLimit); 

    const float altitudeD = -verticalVelocity * altitudePidCoeffs.Kd;

    const float altitudeF = targetAltitudeStep * altitudePidCoeffs.Kf;

    const float throttleAdjustment = altitudeP + altitudeI + altitudeD + altitudeF;
    const float tiltMultiplier = 2.0f - fmaxf(getCosTiltAngle(), 0.5f);

    // 1 = flat, 1.24 at 40 degrees, max 1.5 around 60 degrees, the default limit of Angle Mode
    // 2 - cos(x) is between 1/cos(x) and 1/sqrt(cos(x)) in this range

    const float hoverThrottle = positionControlConfig()->hover_throttle - PWM_RANGE_MIN;

    float newThrottle = (hoverThrottle + throttleAdjustment) * tiltMultiplier + PWM_RANGE_MIN;
    newThrottle = constrainf(newThrottle, positionControlConfig()->alt_control_throttle_min, positionControlConfig()->alt_control_throttle_max);
    throttleOut = scaleRangef(throttleOut, MAX(rxConfig()->mincheck, PWM_RANGE_MIN), PWM_RANGE_MAX, 0.0f, 1.0f);
}

float positionControlThrottle(void)
{
    return throttleOut;
}
