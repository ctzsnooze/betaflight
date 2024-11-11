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
#include "common/axis.h"
#include "common/filter.h"
#include "common/maths.h"
#include "common/vector.h"
#include "fc/rc.h"
#include "fc/runtime_config.h"

#include "flight/imu.h"
#include "flight/position.h"
#include "flight/pos_hold.h"
#include "rx/rx.h"
#include "sensors/gyro.h"

#include "autopilot.h"

#define ALTITUDE_P_SCALE  0.01f
#define ALTITUDE_I_SCALE  0.003f
#define ALTITUDE_D_SCALE  0.01f
#define ALTITUDE_F_SCALE  0.01f
#define POSITION_P_SCALE  0.0012f
#define POSITION_I_SCALE  0.0001f
#define POSITION_D_SCALE  0.0015f
#define POSITION_A_SCALE  0.0008f
#define UPSAMPLING_CUTOFF_HZ 5.0f

static pidCoefficient_t altitudePidCoeffs;
static pidCoefficient_t positionPidCoeffs;

static float altitudeI = 0.0f;
static float throttleOut = 0.0f;


typedef enum {
    EW = 0,
    NS
} axisEF_t;

// apply PIDs on orthodonal axes
typedef struct {
    bool isStarting;
    float previousDistance;
    float previousVelocity;
    float integral;
    pt1Filter_t velocityLpf;
    pt1Filter_t accelerationLpf;
} pidAxis_t;

static struct posHoldState {
    float lpfCutoff;
    float maxAngle;
    gpsLocation_t targetLocation;           // current target
    float sanityCheckDistance;              // maximum allowed distance from target
//    float gpsDataInterval;                  // actual GPS interval, [s]
    bool sticksActive;
    vector2_t pidAngles;                    // pid output, rotated to body frame
    pt3Filter_t upsample[RP_AXIS_COUNT];    // upsampling filter
    pidAxis_t pidAxis[XY_AXIS_COUNT];
} posHold = {
    .lpfCutoff = 1.0f,
    .sanityCheckDistance = 1000.0f,
    .sticksActive = false,
};

float autopilotAngle[RP_AXIS_COUNT];

static void resetPidAxisFilters(pidAxis_t* axis, float filterGain) {
    pt1FilterInit(&axis->velocityLpf, filterGain);
    pt1FilterInit(&axis->accelerationLpf, filterGain);
}

static void resetPidAxisParams(pidAxis_t* axis, float filterGain)
{
    // at start only
    resetPidAxisFilters(axis, filterGain);
    axis->isStarting = true; // Enter starting 'phase'
    axis->integral = 0.0f;
}

void resetPositionControl(const gpsLocation_t* initialTargetLocation)
{
    // from pos_hold.c when initiating position hold at target location
    posHold.targetLocation = *initialTargetLocation;
    posHold.sticksActive = false;
    // set sanity check distance according to groundspeed at start
    posHold.sanityCheckDistance = fmaxf(gpsSol.groundSpeed * 2.0f, 1000.0f);
    for (unsigned i = 0; i < ARRAYLEN(posHold.pidAxis); i++) {
        // disable filter until first call
        resetPidAxisParams(&posHold.pidAxis[i], 1.0f);
    }
}

void autopilotInit(const autopilotConfig_t *config)
{
    posHold.sticksActive = false;
    posHold.maxAngle = autopilotConfig()->max_angle;
    // TODO - separate altitude PID to function
    altitudePidCoeffs.Kp = config->altitude_P * ALTITUDE_P_SCALE;
    altitudePidCoeffs.Ki = config->altitude_I * ALTITUDE_I_SCALE;
    altitudePidCoeffs.Kd = config->altitude_D * ALTITUDE_D_SCALE;
    altitudePidCoeffs.Kf = config->altitude_F * ALTITUDE_F_SCALE;
    positionPidCoeffs.Kp = config->position_P * POSITION_P_SCALE;
    positionPidCoeffs.Ki = config->position_I * POSITION_I_SCALE;
    positionPidCoeffs.Kd = config->position_D * POSITION_D_SCALE;
    positionPidCoeffs.Kf = config->position_A * POSITION_A_SCALE; // Kf used for acceleration
    // initialise filters with approximate filter gain
    const float upsampleCutoff = pt3FilterGain(UPSAMPLING_CUTOFF_HZ, HZ_TO_INTERVAL(POSHOLD_TASK_RATE_HZ)); // 5Hz
    for (unsigned i = 0; i < ARRAYLEN(posHold.upsample); i++) {
        pt3FilterInit(&posHold.upsample[i], upsampleCutoff);
    }

    posHold.lpfCutoff = config->position_cutoff * 0.01f;
    // Initialise PT1 filters for earth frame axes EW and NS
    const float lpfGain = pt1FilterGain(posHold.lpfCutoff,  0.1f); // assume 10Hz GPS connection at start; value is overwritten before first filter usek
    for (unsigned i = 0; i < ARRAYLEN(posHold.pidAxis); i++) {
        resetPidAxisFilters(&posHold.pidAxis[i], lpfGain);
    }
}

void resetAltitudeControl (void) {
    altitudeI = 0.0f;
}

void altitudeControl(float targetAltitudeCm, float taskInterval, float verticalVelocity, float targetAltitudeStep)
{
    const float altitudeErrorCm = targetAltitudeCm - getAltitudeCm();
    const float altitudeP = altitudeErrorCm * altitudePidCoeffs.Kp;

    // reduce the iTerm gain for errors greater than 200cm (2m), otherwise it winds up too much
    const float itermRelax = (fabsf(altitudeErrorCm) < 200.0f) ? 1.0f : 0.1f;
    altitudeI += altitudeErrorCm * altitudePidCoeffs.Ki * itermRelax * taskInterval;
    // limit iTerm to not more than 200 throttle units
    altitudeI = constrainf(altitudeI, -200.0f, 200.0f);

    const float altitudeD = verticalVelocity * altitudePidCoeffs.Kd;

    const float altitudeF = targetAltitudeStep * altitudePidCoeffs.Kf;

    const float hoverOffset = autopilotConfig()->hover_throttle - PWM_RANGE_MIN;
    float throttleOffset = altitudeP + altitudeI - altitudeD + altitudeF + hoverOffset;

    const float tiltMultiplier = 1.0f / fmaxf(getCosTiltAngle(), 0.5f);
    // 1 = flat, 1.3 at 40 degrees, 1.56 at 50 deg, max 2.0 at 60 degrees or higher
    // note: the default limit of Angle Mode is 60 degrees

    throttleOffset *= tiltMultiplier;

    float newThrottle = PWM_RANGE_MIN + throttleOffset;
    newThrottle = constrainf(newThrottle, autopilotConfig()->throttle_min, autopilotConfig()->throttle_max);
    DEBUG_SET(DEBUG_AUTOPILOT_ALTITUDE, 0, lrintf(newThrottle)); // normal range 1000-2000 but is before constraint

    newThrottle = scaleRangef(newThrottle, MAX(rxConfig()->mincheck, PWM_RANGE_MIN), PWM_RANGE_MAX, 0.0f, 1.0f);

    throttleOut = constrainf(newThrottle, 0.0f, 1.0f);

    DEBUG_SET(DEBUG_AUTOPILOT_ALTITUDE, 1, lrintf(tiltMultiplier * 100));
    DEBUG_SET(DEBUG_AUTOPILOT_ALTITUDE, 3, lrintf(targetAltitudeCm));
    DEBUG_SET(DEBUG_AUTOPILOT_ALTITUDE, 4, lrintf(altitudeP));
    DEBUG_SET(DEBUG_AUTOPILOT_ALTITUDE, 5, lrintf(altitudeI));
    DEBUG_SET(DEBUG_AUTOPILOT_ALTITUDE, 6, lrintf(-altitudeD));
    DEBUG_SET(DEBUG_AUTOPILOT_ALTITUDE, 7, lrintf(altitudeF));
}

void setSticksActiveStatus(bool areSticksActive)
{
    posHold.sticksActive = areSticksActive;
}


static void setTargetLocationAxis(const gpsLocation_t* newTargetLocation, axisEF_t efAxisIdx)
{
    if (efAxisIdx == EW) {
        posHold.targetLocation.lon = newTargetLocation->lon; // update East-West / / longitude position
    } else {
        posHold.targetLocation.lat = newTargetLocation->lat; // update North-South / latitude position
    }
    posHold.pidAxis[efAxisIdx].previousDistance = 0.0f; // and reset the previous distance to avoid D and A spikes
}

static void setTargetLocation(const gpsLocation_t* newTargetLocation)
{
    for (unsigned i = XY_AXIS_COUNT; i < 2; i++) {
        setTargetLocationAxis(newTargetLocation, i);
    }
    // function is intended for only small changes in position
    // for example, where the step distance change reflects an intended velocity, determined by a client function
    // if we had a 'target_ground_speed' value, like in gps_rescue, we can make a function that starts and stops smoothly and targets that velocity
}

bool positionControl(void)
{
    unsigned debugAxis = gyroConfig()->gyro_filter_debug_axis;

    static uint16_t previousGpsStamp = 0;
    if (getGpsStamp() != previousGpsStamp) {
        previousGpsStamp = getGpsStamp();
        float gpsDataInterval = getGpsDataIntervalSeconds(); // interval for current GPS data value 0.01s to 1.0s
        float gpsDataFreq = 1.0f / gpsDataInterval;

        // first get NS and EW distances from current location (gpsSol.llh) to target location
        vector2_t distance;
        GPS_distance2d(&gpsSol.llh, &posHold.targetLocation, &distance); // X is EW, Y is NS
        const float distanceNorm = vector2Norm(&distance);

        // ** Sanity check **
        // primarily to detect flyaway from no Mag or badly oriented Mag
        // must accept some overshoot at the start, especially if entering at high speed
        if (distanceNorm > posHold.sanityCheckDistance) {
            return false;
        }

        // update position and velocity filter coefs to actual GPS rate
        const float lpfGain = pt1FilterGain(posHold.lpfCutoff, gpsDataInterval);

        vector2_t pidSum = { 0 };       // P+I in loop, D+A added after limiting it
        vector2_t pidDA;                // D+A
        for (axisEF_t pidAxisIdx = 0; pidAxisIdx < XY_AXIS_COUNT; pidAxisIdx++) {
            pidAxis_t *pidAxis = &posHold.pidAxis[pidAxisIdx];
            // separate PID controllers for longitude (EastWest or EW, X) and latitude (NorthSouth or NS, Y)
            const float axisDistance = distance.v[pidAxisIdx];

            // ** P **
            const float pidP = axisDistance * positionPidCoeffs.Kp;
            pidSum.v[pidAxisIdx] += pidP;

            // ** I **
            // accumulate iTerm only while in hold phase
            pidAxis->integral += pidAxis->isStarting ? 0.0f : axisDistance * gpsDataInterval;
            const float pidI = pidAxis->integral * positionPidCoeffs.Ki;
            pidSum.v[pidAxisIdx] += pidI;

            // ** D ** //
            // Velocity derived from GPS position works better than module supplied GPS Speed and Heading information
            const float deltaPos = (axisDistance - pidAxis->previousDistance) * gpsDataFreq; // cm/s
            pidAxis->previousDistance = axisDistance;
            // filter delta to get velocity
            pt1FilterUpdateCutoff(&pidAxis->velocityLpf, lpfGain);
            const float velocity = pt1FilterApply(&pidAxis->velocityLpf, deltaPos);
            float pidD = velocity * positionPidCoeffs.Kd;

            // differentiate velocity another time to get acceleration
            float deltaVel = (velocity - pidAxis->previousVelocity) / gpsDataFreq;
            pidAxis->previousVelocity = velocity;
            // apply second filter to acceleration
            pt1FilterUpdateCutoff(&pidAxis->accelerationLpf, lpfGain);
            const float acceleration = pt1FilterApply(&pidAxis->accelerationLpf, deltaVel);
            const float pidA = acceleration * positionPidCoeffs.Kf;

            if (!posHold.sticksActive && pidAxis->isStarting) {
                // 'phase' after sticks stop, but before craft has stopped in this axis
                pidD *= 1.6f; // aribitrary D boost to stop more quickly when sticks are centered
                // detect deltaPos zero crossing (velocity is delayed by filter)
                if (deltaPos * velocity < 0.0f) {
                    // when an axis has nearly stopped moving, reset target position for it
                    setTargetLocationAxis(&gpsSol.llh, pidAxisIdx);
                    pidAxis->isStarting = false;
                }
            }
            pidDA.v[pidAxisIdx] = pidD + pidA;
            if (debugAxis == pidAxisIdx) {
                DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 0, lrintf(axisDistance));
                DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 4, lrintf(pidP * 10));
                DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 5, lrintf(pidI * 10));
                DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 6, lrintf(pidD * 10));
                DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 7, lrintf(pidA * 10));
            }
        } // pidAxis loop
        {
            // limit sum of D and A per axis based on total DA vector length, otherwise can be too aggressive when starting at speed
            // limit is 35 degrees from D and A alone, arbitrary value.  20 is a bit too low, allows a lot of overshoot
            // note: an angle of more than 35 degrees can still be achieved as P and I grow
            const float maxDAAngle = 35.0f; // D+A limit in degrees; arbitrary angle
            const float mag = vector2Norm(&pidDA);
            if (mag > maxDAAngle) {
                vector2Scale(&pidDA, &pidDA, maxDAAngle / mag);
            }
        }
        // add constrained DA to sum
        vector2Add(&pidSum, &pidSum, &pidDA);

        // calculate pid output - desired ROLL/PITCH angles, in body frame
        vector2_t anglesBF;
        if (posHold.sticksActive) {             // sticks active 'phase'
            // while sticks are moving, reset target on each cycle, to maintain a usable D value
            setTargetLocation(&gpsSol.llh);
            // update sanity check distance while sticks are out
            posHold.sanityCheckDistance = fmaxf(gpsSol.groundSpeed * 2.0f, 1000.0f);
            // leak iTerm while sticks are centered, 2.4s time constant
            const float leak = 1 - pt1FilterGainFromDelay(2.4, gpsDataInterval);
            for (unsigned i = 0; i < ARRAYLEN(posHold.pidAxis); i++) {
                // deceleration phase
                posHold.pidAxis[i].isStarting = true;
                // slowly leak iTerm away
                posHold.pidAxis[i].integral *= leak;
            }
            // if a Position Hold deadband is set, and sticks are outside deadband, allow pilot control in angle mode
            anglesBF = (vector2_t){{0,0}};             // set output from of PIDS to 0; upsampling filter will smooth this
        } else {
            // ** Rotate pid Sum to body frame, and convert it into pitch and roll **
            // attitude.values.yaw is clockwise from north
            // PID is running in ENU, adapt angle (to 0deg = EAST);
            //  rotation is from EF to BF, no change of sign from heading
            const float angle = DECIDEGREES_TO_RADIANS(attitude.values.yaw - 900);
            vector2_t pidBF;   // pid output in body frame; X is forward, Y is left
            vector2Rotate(&pidBF, &pidSum, angle);  // rotate by angle counterclockwise
            anglesBF.v[AI_ROLL] = -pidBF.y;         // need negative roll to fly left
            anglesBF.v[AI_PITCH] = pidBF.x;         // positive pitch for forward
             // limit angle vector to maxAngle
            const float mag = vector2Norm(&anglesBF);
            if (mag > posHold.maxAngle && mag > 0.0f) {
                vector2Scale(&anglesBF, &anglesBF, posHold.maxAngle / mag);
            }
        }
        posHold.pidAngles = anglesBF;    // this value will be upsampled
    } // new GPS position

    // Final output to pid.c Angle Mode at 100Hz with PT3 upsampling
    for (unsigned i = 0; i < RP_AXIS_COUNT; i++) {
        // note: upsampling should really be done in earth frame, to avoid 10Hz wobbles if pilot yaws and the controller is applying significant pitch or roll
        autopilotAngle[i] = pt3FilterApply(&posHold.upsample[i], posHold.pidAngles.v[i]);
    }

    if (debugAxis < 2) {
        // debug[1] is not set
        DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 2, lrintf(posHold.pidAngles.v[debugAxis] * 10)); // deg * 10 ; can be set in GPS change only
        DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 3, lrintf(autopilotAngle[debugAxis] * 10));   // deg * 10
    }
    return true;
}

bool isBelowLandingAltitude(void)
{
    return getAltitudeCm() < 100.0f * autopilotConfig()->landing_altitude_m;
}

float getAutopilotThrottle(void)
{
    return throttleOut;
}

bool isAutopilotActive(void)
{
    return !posHold.sticksActive;
}
