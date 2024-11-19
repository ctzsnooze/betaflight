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
#include "flight/alt_hold.h" // only to get ALTHOLD_TASK_RATE_HZ
#include "rx/rx.h"
#include "sensors/gyro.h"

#include "pg/autopilot.h"
#include "autopilot.h"

#define ALTITUDE_P_SCALE  0.005f
#define ALTITUDE_I_SCALE  0.00015f
#define ALTITUDE_D_SCALE  0.005f
#define ALTITUDE_F_SCALE  0.004f
#define POSITION_P_SCALE  0.001f
#define POSITION_I_SCALE  0.0001f
#define POSITION_D_SCALE  0.0015f
#define POSITION_A_SCALE  0.0008f
#define UPSAMPLING_CUTOFF_HZ 5.0f

static pidCoefficient_t altitudePidCoeffs;
static pidCoefficient_t positionPidCoeffs;

static float altitudeI = 0.0f;
static float throttleOut = 0.0f;
pt1Filter_t altitudePHpf;

typedef struct earthFrame_s {
    bool isStopping;
    float distance;
    float previousDistance;
    float previousVelocity;
    float integral;
    float pidSum;
    pt1Filter_t velocityLpf;
    pt1Filter_t accelerationLpf;
} earthFrame_t;

typedef enum {
    EW = 0,
    NS
} axisEF_t;

typedef struct autopilotState_s {
    float gpsDataIntervalS;
    float gpsDataFreqHz;
    float sanityCheckDistance;
    float upsampleCutoff;
    float pt1Cutoff;
    float pt1Gain;
    bool sticksActive;
    float maxAngle;
    float distanceCm;
    float pidSumRP[EF_AXIS_COUNT];
    pt3Filter_t upsample[EF_AXIS_COUNT];
    earthFrame_t efAxis[EF_AXIS_COUNT];
} autopilotState_t;

static autopilotState_t ap = {
    .gpsDataIntervalS = 0.1f,
    .gpsDataFreqHz = 10.0f,
    .sanityCheckDistance = 1000.0f,
    .pt1Cutoff = 1.0f,
    .pt1Gain = 1.0f,
    .sticksActive = false,
    .distanceCm = 0.0f,
    .pidSumRP = { 0.0f, 0.0f },
    .upsample = { {0}, {0} },
    .efAxis = { {0} }
};

static gpsLocation_t currentTargetLocation = {0, 0, 0};
float autopilotAngle[RP_AXIS_COUNT];

void initializeEfAxisFilters(earthFrame_t *efAxis, float filterGain) {
    pt1FilterInit(&efAxis->velocityLpf, filterGain);
    pt1FilterInit(&efAxis->accelerationLpf, filterGain);
}

void resetPt3UpsampleFilters(void)
{
    pt3FilterInit(&ap.upsample[AI_ROLL], ap.upsampleCutoff);
    pt3FilterInit(&ap.upsample[AI_PITCH], ap.upsampleCutoff);
}

void autopilotInit(void)
{
    ap.sticksActive = false;
    ap.maxAngle = autopilotConfig()->ap_max_angle;
    altitudePidCoeffs.Kp = autopilotConfig()->ap_altitude_P * ALTITUDE_P_SCALE;
    altitudePidCoeffs.Ki = autopilotConfig()->ap_altitude_I * ALTITUDE_I_SCALE;
    altitudePidCoeffs.Kd = autopilotConfig()->ap_altitude_D * ALTITUDE_D_SCALE;
    altitudePidCoeffs.Kf = autopilotConfig()->ap_altitude_F * ALTITUDE_F_SCALE;
    positionPidCoeffs.Kp = autopilotConfig()->ap_position_P * POSITION_P_SCALE;
    positionPidCoeffs.Ki = autopilotConfig()->ap_position_I * POSITION_I_SCALE;
    positionPidCoeffs.Kd = autopilotConfig()->ap_position_D * POSITION_D_SCALE;
    positionPidCoeffs.Kf = autopilotConfig()->ap_position_A * POSITION_A_SCALE; // Kf used for acceleration
    // initialise PT3 filters with approximate filter gain
    ap.upsampleCutoff = pt3FilterGain(UPSAMPLING_CUTOFF_HZ, 0.01f); // assuming 100Hz task rate
    resetPt3UpsampleFilters();
    // Initialise PT1 filters for earth frame axes NS and EW
    ap.pt1Cutoff = autopilotConfig()->ap_position_cutoff * 0.01f;
    ap.pt1Gain = pt1FilterGain(ap.pt1Cutoff, 0.1f); // assume 10Hz GPS connection at start
    initializeEfAxisFilters(&ap.efAxis[EW], ap.pt1Gain);
    initializeEfAxisFilters(&ap.efAxis[NS], ap.pt1Gain);

    const float altPHpfGain = pt1FilterGain(0.2f, HZ_TO_INTERVAL(ALTHOLD_TASK_RATE_HZ)); // Approx 1s time constant
    pt1FilterInit(&altitudePHpf, altPHpfGain);
}

void resetAltitudeControl(void) {
    altitudeI = 0.0f;
}

void altitudeControl(float targetAltitudeCm, float taskIntervalS, float targetAltitudeStep)
{
    const float altitudeErrorCm = targetAltitudeCm - getAltitudeCm();
    float altitudeP = altitudeErrorCm * altitudePidCoeffs.Kp;

    // add highpass of P to compensate for motor to altitude lag; D isn't quite the same and is needed as well
    const float highpassAltitudeP = altitudeP - pt1FilterApply(&altitudePHpf, altitudeP);
    altitudeP += highpassAltitudeP;

    // reduce the iTerm gain for errors greater than 200cm (2m), otherwise it winds up too much
    const float itermRelax = (fabsf(altitudeErrorCm) < 200.0f) ? 1.0f : 0.1f;
    altitudeI += altitudeErrorCm * altitudePidCoeffs.Ki * itermRelax * taskIntervalS;
    // limit iTerm to not more than 200 throttle units
    altitudeI = constrainf(altitudeI, -200.0f, 200.0f);

    // increase D when velocity is high, typically when initiating hold at high vertical speeds
    // 1.0 when less than 5 m/s, 2x at 10m/s, 2.5 at 20 m/s, 2.8 at 50 m/s, asymptotes towards max 3.0.
    float dBoost = 1.0f;
    const float startValue = 500.0f; // velocity in cm/s at which D should start to increase
    const float altDeriv = fabsf(getAltitudeDerivative());
    if (altDeriv > startValue) {
        const float ratio = altDeriv / startValue;
        dBoost = (3.0f * ratio - 2.0f) / ratio;
    }

    const float altitudeD = getAltitudeDerivative() * dBoost * altitudePidCoeffs.Kd;

    const float altitudeF = targetAltitudeStep * altitudePidCoeffs.Kf;

    const float hoverOffset = autopilotConfig()->ap_hover_throttle - PWM_RANGE_MIN;
    float throttleOffset = altitudeP + altitudeI - altitudeD + altitudeF + hoverOffset;

    const float tiltMultiplier = 1.0f / fmaxf(getCosTiltAngle(), 0.5f);
    // 1 = flat, 1.3 at 40 degrees, 1.56 at 50 deg, max 2.0 at 60 degrees or higher
    // note: the default limit of Angle Mode is 60 degrees

    throttleOffset *= tiltMultiplier;

    float newThrottle = PWM_RANGE_MIN + throttleOffset;
    newThrottle = constrainf(newThrottle, autopilotConfig()->ap_throttle_min, autopilotConfig()->ap_throttle_max);
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

void resetPositionControlEFParams(earthFrame_t *efAxis)
{
    // at start only
    pt1FilterInit(&efAxis->velocityLpf, ap.pt1Gain);
    pt1FilterInit(&efAxis->accelerationLpf, ap.pt1Gain);
    efAxis->isStopping = true; // Enter starting 'phase'
    efAxis->integral = 0.0f;
}

void resetPositionControl(gpsLocation_t *initialTargetLocation)
{
    // from pos_hold.c when initiating position hold at target location
    currentTargetLocation = *initialTargetLocation;
    ap.sticksActive = false;
    // set sanity check distance according to groundspeed at start
    ap.sanityCheckDistance = gpsSol.groundSpeed > 500 ? gpsSol.groundSpeed * 2.0f : 1000.0f;
    ap.velocityI = 0.0f;
    resetPositionControlEFParams(&ap.efAxis[EW]); // includes reset of PT1 filters
    resetPositionControlEFParams(&ap.efAxis[NS]);
    resetPt3UpsampleFilters(); // clear residual in filter from a previous run
}

void setSticksActiveStatus(bool areSticksActive)
{
    ap.sticksActive = areSticksActive;
}

void moveTargetLocation(int32_t latStep, int32_t lonStep)
{
    currentTargetLocation.lon += ap.efAxis[EW].isStopping ? 0 : lonStep;
    currentTargetLocation.lat += ap.efAxis[NS].isStopping ? 0 : latStep;
    // client runs this function, once position hold is initiated, to provide stepwise position change ie velocity control
    // in GPS Rescue, steps reflect requested velocity, which will be zero until rotation is complete
}

void updateLocation(earthFrame_t *efAxis, axisEF_t loopAxis)
{
    if (loopAxis == EW) {
        currentTargetLocation.lon = gpsSol.llh.lon; // update East-West / / longitude position
    } else {
        currentTargetLocation.lat = gpsSol.llh.lat; // update North-South / latitude position
    }
    efAxis->previousDistance = 0.0f; // and reset the previous distance
}

void (posControlOnNewGpsData) (void)
{
    ap.gpsDataIntervalS = getGpsDataIntervalSeconds(); // interval for current GPS data value 0.01s to 1.0s
    ap.gpsDataFreqHz = getGpsDataFrequencyHz();

    // get NS and EW distances from current location (gpsSol.llh) to target location in cm
    vector2_t gpsDistance;
    GPS_latLongVectors(&gpsSol.llh, &currentTargetLocation, &gpsDistance.x, &gpsDistance.y); // X is EW, Y is NS
    ap.efAxis[EW].distance = gpsDistance.x;
    ap.efAxis[NS].distance = gpsDistance.y;

    ap.distanceCm = vector2Norm(&gpsDistance);

    const float leak = 1.0f - 0.4f * ap.gpsDataIntervalS;
    // leak iTerm while sticks are centered, 2s time constant approximately
    const float lpfGain = pt1FilterGain(ap.pt1Cutoff, ap.gpsDataIntervalS);

    static float prevPidDASquared = 0.0f; // if we limit DA on true vector length
    const float maxDAAngle = 35.0f; // limit in degrees; arbitrary angle
    const float sqMaxDAAngle = sq(maxDAAngle);

    for (axisEF_t loopAxis = EW; loopAxis <= NS; loopAxis++) {
        earthFrame_t *efAxis = &ap.efAxis[loopAxis];
        // separate PID controllers for latitude (NorthSouth or NS) and longitude (EastWest or EW)

        // ** P **
        const float pidP = efAxis->distance * positionPidCoeffs.Kp;

        // ** I **
        // no accumulation while stopping,
        efAxis->integral += efAxis->isStopping ? 0.0f : efAxis->distance * ap.gpsDataIntervalS;

        // maybe limit I ? ... no more than 20 degrees of angle from iTerm alone
        // efAxis->integral += efAxis->isStopping || fabsf(efAxis->integral) > 200.0f ? 0.0f : efAxis->distance * ap.gpsDataIntervalS;

        // no iTerm in GPS Rescue mode for now but perhaps an effective limiter may work better
        const float pidI = FLIGHT_MODE(GPS_RESCUE_MODE) ? 0.0f : efAxis->integral * positionPidCoeffs.Ki;

        // ** D ** //
        // Velocity derived from GPS position works better than module supplied GPS Speed and Heading information
        float velocity = (efAxis->distance - efAxis->previousDistance) * ap.gpsDataFreqHz; // cm/s, minimum step 11.1 cm/s
        efAxis->previousDistance = efAxis->distance;
        pt1FilterUpdateCutoff(&efAxis->velocityLpf, lpfGain);
        const float velocityFiltered = pt1FilterApply(&efAxis->velocityLpf, velocity);
        float pidD = velocityFiltered * positionPidCoeffs.Kd;

        // ** A ** //
        float acceleration = (velocity - efAxis->previousVelocity) * ap.gpsDataFreqHz;
        efAxis->previousVelocity = velocity;
        pt1FilterUpdateCutoff(&efAxis->accelerationLpf, lpfGain);
        const float accelerationFiltered = pt1FilterApply(&efAxis->accelerationLpf, acceleration);
        const float pidA = accelerationFiltered * positionPidCoeffs.Kf;

        if (ap.sticksActive) {
            // sticks active 'phase'
            efAxis->isStopping = true;
            updateLocation(efAxis, loopAxis);
            // while sticks are moving, update the location on each axis, maintaining a usable D value
            // slowly leak iTerm away, approx 2s time constant
            efAxis->integral *= leak;
            // increase sanity check distance depending on speed, typically maximal when sticks stop
        } else if (efAxis->isStopping) {
            // 'phase' after sticks stop, but before craft has stopped
            pidD *= 1.6f; // aribitrary D boost to stop more quickly when sticks are centered
            if (velocity * velocityFiltered < 0.0f) {
                // when an axis has nearly stopped moving, reset it and end its 'stopping' phase
                updateLocation(efAxis, loopAxis);
                efAxis->isStopping = false;
            }
        }

       float pidDA = pidD + pidA;
 
        // limit sum of D and A per axis based on total DA vector length
        // limit is 35 degrees from D and A alone, arbitrary value.
        // note: an angle of more than 35 degrees can still be achieved as P and I grow
        const float pidDASquared = sq(pidD + pidA);
        float magSq = pidDASquared + prevPidDASquared;
        if (magSq > sqMaxDAAngle && magSq > 0.0f) {
            pidDA *= maxDAAngle / sqrtf(magSq);
        }
        prevPidDASquared = pidDASquared;

        // ** PID Sum **
        efAxis->pidSum = pidP + pidI + pidDA;

        if (gyroConfig()->gyro_filter_debug_axis == loopAxis) {
            DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 0, lrintf(ap.distanceCm));
            DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 4, lrintf(pidP * 10));
            DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 5, lrintf(pidI * 10));
            DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 6, lrintf(pidD * 10));
            DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 7, lrintf(pidA * 10));

            DEBUG_SET(DEBUG_GPS_RESCUE_VELOCITY, 4, lrintf(pidP * 10));
            DEBUG_SET(DEBUG_GPS_RESCUE_VELOCITY, 5, lrintf(pidD * 10));
            DEBUG_SET(DEBUG_GPS_RESCUE_VELOCITY, 6, lrintf(pidA * 10));
        }
    }

    DEBUG_SET(DEBUG_GPS_RESCUE_VELOCITY, 7, lrintf(ap.distanceCm));


    if (ap.sticksActive) {
        // keep updating sanity check
        ap.sanityCheckDistance = gpsSol.groundSpeed > 500 ? gpsSol.groundSpeed * 2.0f : 1000.0f;
        // allow pilot control, in angle mode
        ap.pidSumRP[AI_ROLL] = 0.0f;
        ap.pidSumRP[AI_PITCH] = 0.0f;
    } else {
        // ** Rotate pid Sum to quad frame of reference, into pitch and roll **
        const float headingRads = DECIDEGREES_TO_RADIANS(attitude.values.yaw);
        const float sinHeading = sin_approx(headingRads);
        const float cosHeading = cos_approx(headingRads);
        ap.pidSumRP[AI_ROLL] = -sinHeading * ap.efAxis[NS].pidSum + cosHeading * ap.efAxis[EW].pidSum;
        ap.pidSumRP[AI_PITCH] = cosHeading * ap.efAxis[NS].pidSum + sinHeading * ap.efAxis[EW].pidSum;

        // limit angle vector to maxAngle
        const float angleMagSq = sq(ap.pidSumRP[AI_ROLL]) + sq(ap.pidSumRP[AI_PITCH]);
        if (angleMagSq > sq(ap.maxAngle) && angleMagSq > 0.0f) {
            const float limiter = ap.maxAngle / sqrtf(angleMagSq);
            ap.pidSumRP[AI_ROLL] *= limiter;  // Scale the roll value
            ap.pidSumRP[AI_PITCH] *= limiter; // Scale the pitch value
        }
    }
}

void posControlOutput (void)
{
    // ** Final output to pid.c Angle Mode at 100Hz with primitive upsampling**
    autopilotAngle[AI_ROLL] = pt3FilterApply(&ap.upsample[AI_ROLL], ap.pidSumRP[AI_ROLL]);
    autopilotAngle[AI_PITCH] = pt3FilterApply(&ap.upsample[AI_PITCH], ap.pidSumRP[AI_PITCH]);
    // note: upsampling should really be done in earth frame, to avoid 10Hz wobbles if pilot yaws and the controller is applying significant pitch or roll

    if (gyroConfig()->gyro_filter_debug_axis == FD_ROLL) {
        DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 1, lrintf(ap.efAxis[EW].distance));    // cm
        DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 2, lrintf(ap.efAxis[EW].pidSum * 10)); // deg
        DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 3, lrintf(autopilotAngle[AI_ROLL] * 10));   // deg
    } else {
        DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 1, lrintf(ap.efAxis[NS].distance));
        DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 2, lrintf(ap.efAxis[NS].pidSum * 10));
        DEBUG_SET(DEBUG_AUTOPILOT_POSITION, 3, lrintf(autopilotAngle[AI_PITCH] * 10));
    }
}

bool positionControl(void) 
{
    static uint16_t gpsStamp = 0;
    if (gpsHasNewData(&gpsStamp)) {
        posControlOnNewGpsData();
    }
    posControlOutput();
    // ** Sanity check for Pos Hold failure **
    // primarily to detect flyaway from no Mag or badly oriented Mag in pos hold mode
    // GPS Rescue has its own sanity checks
    if (FLIGHT_MODE(POS_HOLD_MODE) && ap.distanceCm > ap.sanityCheckDistance) {
        return false;
    }
    return true;
}

void apVelocityControl(float targetVelocityCmS)
{
    // ultra basic DI velocity controller, intended to pitch forward at a set velocity, needs to be checked
    const float velocityError = targetVelocityCmS - gpsSol.groundSpeed;
    const float velocityD = velocityError * positionPidCoeffs.Kd;
    ap.velocityI += velocityError * positionPidCoeffs.Ki * ap.gpsDataIntervalS;
    const float velocityPidSum = velocityD + ap.velocityI;
    ap.pidSumRP[AI_ROLL] = 0.0f;
    ap.pidSumRP[AI_PITCH] = velocityPidSum;
}

bool isBelowLandingAltitude(void)
{
    return getAltitudeCm() < 100.0f * autopilotConfig()->ap_landing_altitude_m;
}

float getAutopilotThrottle(void)
{
    return throttleOut;
}

bool isAutopilotActive(void)
{
    return !ap.sticksActive;
}
