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
#include <math.h>

#include "platform.h"

#ifdef USE_GPS_RESCUE

#include "build/debug.h"

#include "common/axis.h"
#include "common/filter.h"
#include "common/maths.h"
#include "common/utils.h"

#include "config/config.h"
#include "drivers/time.h"

#include "fc/core.h"
#include "fc/rc_controls.h"
#include "fc/rc_modes.h"
#include "fc/runtime_config.h"

#include "flight/autopilot.h"
#include "flight/failsafe.h"
#include "flight/imu.h"
#include "flight/pid.h"
#include "flight/position.h"

#include "io/gps.h"
#include "rx/rx.h"
#include "pg/autopilot.h"
#include "sensors/acceleration.h"
#include "sensors/compass.h"

#include "gps_rescue.h"

typedef enum {
    RESCUE_IDLE,
    RESCUE_INITIALIZE,
    RESCUE_ATTAIN_ALT,
    RESCUE_ORIENT_IMU,
    RESCUE_ROTATE,
    RESCUE_FLY_HOME,
    RESCUE_DESCENT,
    RESCUE_LANDING,
    RESCUE_DO_NOTHING,
    RESCUE_ABORT
} rescuePhase_e;

typedef enum {
    RESCUE_HEALTHY,
    RESCUE_FLYAWAY,
    RESCUE_GPSLOST,
    RESCUE_LOWSATS,
    RESCUE_CRASHFLIP_DETECTED,
    RESCUE_STALLED,
    RESCUE_TOO_CLOSE,
    RESCUE_NO_HOME_POINT
} rescueFailureState_e;

typedef struct {
    float maxAltitudeCm;
    float returnAltitudeCm;
    float targetAltitudeCm;
    float targetAltitudeStepCm;
    float targetVelocityCmS;
    float descentDistanceCm;
    int8_t secondsFailing;
    float yawAttenuator;
    float velocityAttenuator;
    float proximityAttenuator;
    float disarmThreshold;
    vector2_t latLonSteps; // lat first to match gpsLocation_t which has latitude first
    float cmToEarthAngle;
    float initialClimbCm;
    bool forceDisableMag;
} rescueIntent_s;

typedef struct {
    uint16_t groundSpeedCmS;
    bool gpsHealthy;
    float errorAngleDeg;
    float velocityToHomeCmS;
    float imuYawCogGain;
} rescueSensorData_s;

typedef struct {
    rescuePhase_e phase;
    rescueFailureState_e failure;
    rescueSensorData_s sensor;
    rescueIntent_s intent;
    bool isAvailable;
} rescueState_s;

#define GPS_RESCUE_MAX_YAW_RATE          180    // deg/sec max yaw rate
#define GPS_RESCUE_ALLOWED_YAW_RANGE   30.0f   // yaw error must be less than this to enter fly home phase, and to pitch during descend()

static const float taskIntervalSeconds = HZ_TO_INTERVAL(TASK_GPS_RESCUE_RATE_HZ); // i.e. 0.01 s
static float rescueYaw;

rescueState_s rescueState;

void gpsRescueInit(void)
{
    rescueState.intent.cmToEarthAngle = 1.0f / EARTH_ANGLE_TO_CM; // approx 0.898 cm per unit lat at equator
    rescueState.intent.initialClimbCm = gpsRescueConfig()->initialClimbM * 100.0f;
    rescueState.intent.disarmThreshold = gpsRescueConfig()->disarmThreshold * 0.1f;
    rescueState.intent.descentDistanceCm = gpsRescueConfig()->descentDistanceM * 100.0f;
}

static void rescueStart(void)
{
    rescueState.phase = RESCUE_INITIALIZE;
}

static void rescueStop(void)
{
    rescueState.phase = RESCUE_IDLE;
}

// While idle, update the altitude targets 
static void setReturnAltitude(void)
{
    // While disarmed, force maxAltitude to zero , unless set_home_point_once is true, when maxAlt is only reset on a power cycle
    if (!ARMING_FLAG(ARMED) && !gpsConfig()->gps_set_home_point_once) {
        rescueState.intent.maxAltitudeCm = 0.0f;
        return;
    }
    // Otherwise store maxAltitude for this armed period
    rescueState.intent.maxAltitudeCm = fmaxf(getAltitudeCm(), rescueState.intent.maxAltitudeCm);
}

static void setLatLonSteps(void)
{
    // precalculate the latitude and longitude step per cm according to direction to home
    float directionToHomeDegrees = DECIDEGREES_TO_DEGREES(GPS_directionToHome);
    if (directionToHomeDegrees > 180) {
        directionToHomeDegrees -= 360;
    } else if (directionToHomeDegrees < -180) {
        directionToHomeDegrees += 360;
    }
    const float directionToHomeRadians = DEGREES_TO_RADIANS(directionToHomeDegrees);
    rescueState.intent.latLonSteps.v[0] = cos_approx(directionToHomeRadians) * rescueState.intent.cmToEarthAngle;  // Latitude (North)
    rescueState.intent.latLonSteps.v[1] = sin_approx(directionToHomeRadians) * rescueState.intent.cmToEarthAngle / getGpsCosLat(); // Longitude (East)
}

static bool isHeadingOK(void)
{
    return (
#ifdef USE_MAG
    compassEnabledAndCalibrated() ||
#endif
    canUseGPSHeading);
}

static void rescueAttainPosition(bool newGpsData)
{
    // runs at 100hz, but only updates RPYT settings when new GPS Data arrives and when not in idle phase.
    switch (rescueState.phase) {
    case RESCUE_IDLE:
        // do nothing
        return;

    case RESCUE_DO_NOTHING:
        // 20s of hover at current altitude, for switch induced sanity failures, to allow time to recover
        // don't hold position in case IMU disorientation caused a sanity failure
        return;

    case RESCUE_INITIALIZE:
        rescueState.intent.targetAltitudeCm = getAltitudeCm(); // initally current altitude
        switch (gpsRescueConfig()->altitudeMode) {
            case GPS_RESCUE_ALT_MODE_FIXED:
                rescueState.intent.returnAltitudeCm = gpsRescueConfig()->returnAltitudeM * 100.0f;
                break;
            case GPS_RESCUE_ALT_MODE_CURRENT:
                // climb above current altitude, but always return at least initial height above takeoff point, in case current altitude was negative
                rescueState.intent.returnAltitudeCm = fmaxf(rescueState.intent.initialClimbCm, getAltitudeCm() + rescueState.intent.initialClimbCm);
                break;
            case GPS_RESCUE_ALT_MODE_MAX:
            default:
                rescueState.intent.returnAltitudeCm = rescueState.intent.maxAltitudeCm + rescueState.intent.initialClimbCm;
                break;
        }
        // initialise the required autopilot functions
        resetAltitudeControl();
        resetPositionControl(&gpsSol.llh, TASK_GPS_RESCUE_RATE_HZ); // enables position hold at current location with hard stop
        return;

     default:
        break;
    }

    /**
        Altitude (throttle) controller
    */
    altitudeControl(rescueState.intent.targetAltitudeCm, taskIntervalSeconds, rescueState.intent.targetAltitudeStepCm);

    if (rescueState.phase == RESCUE_ORIENT_IMU) {
        // only handle altitude if IMU is not oriented
        return;
    }

    /**
        Heading / yaw controller
    */
    float yawRate = rescueState.sensor.errorAngleDeg * rescueState.intent.yawAttenuator * gpsRescueConfig()->yawP * 0.1f;
    yawRate = constrainf(yawRate, -GPS_RESCUE_MAX_YAW_RATE, GPS_RESCUE_MAX_YAW_RATE);
    rescueYaw = yawRate * GET_DIRECTION(rcControlsConfig()->yaw_control_reversed);
    // (extern) rescueYaw is the yaw rate in deg/s to correct the heading error
    DEBUG_SET(DEBUG_GPS_RESCUE_HEADING, 7, rescueYaw);

    /*
        Pitch / velocity controller
    */
    if (newGpsData) {
        if (rescueState.intent.targetVelocityCmS > 0.0f) {
            // only possible in fly home or descend modes
            // move target location along a path, step by step
            const float distanceToMove = rescueState.intent.targetVelocityCmS * getGpsDataIntervalSeconds();
            setLatLonSteps(); // update latitude and longitude step from current location to home at current target velocity
            vector2Scale(&rescueState.intent.latLonSteps, &rescueState.intent.latLonSteps, distanceToMove);
            // send steps to update the target location in autopilot.c 
            moveTargetLocation(rescueState.intent.latLonSteps);
            // run the autopilot function that calculates earth frame PID sums and converts to pitch and roll values
            // must have an accurate aircraft heading estimate from the IMU
        }
        posControlOnNewGpsData(); // hold position if zero set velocity, otherwise move at target velocity
    }
    // upsample the pitch and roll setpoints for pid.c, at gps_rescue task rate
    posControlOutput();
    DEBUG_SET(DEBUG_GPS_RESCUE_TRACKING, 1, lrintf(rescueState.intent.targetVelocityCmS));
}

bool oneSecondPassed(timeUs_t currentTimeUs, timeUs_t *lastTimeUs) {
    timeDelta_t deltaTime = cmpTimeUs(currentTimeUs, *lastTimeUs);
    if (deltaTime >= 1000000) {
        *lastTimeUs = currentTimeUs;
        return true;
    }
    return false;
}

static void performSanityChecks(void)
{
    static float prevAltitudeCm = 0.0f;            // to calculate ascent or descent change
    static float prevTargetAltitudeCm = 0.0f;      // to calculate ascent or descent target change
    static float previousDistanceToHomeCm = 0.0f;  // to check that we are returning
    static int8_t secondsLowSats = 0;              // Minimum sat detection
    static int8_t secondsDoingNothing;             // Limit on doing nothing
    const timeUs_t currentTimeUs = micros();

    if (rescueState.phase == RESCUE_IDLE) {
        rescueState.failure = RESCUE_HEALTHY;
        return;
    } else if (rescueState.phase == RESCUE_INITIALIZE) {
        // Initialize these variables each time a GPS Rescue is started
        prevAltitudeCm = getAltitudeCm();
        prevTargetAltitudeCm = rescueState.intent.targetAltitudeCm;
        previousDistanceToHomeCm = GPS_distanceToHomeCm;
        secondsLowSats = 0;
        secondsDoingNothing = 0;
    }

    // Handle events that set a failure mode to other than healthy.
    // Disarm via Abort when sanity on, or for hard Rx loss in FS_ONLY mode
    // Otherwise allow 20s of semi-controlled descent with impact disarm detection
    const bool hardFailsafe = !isRxReceivingSignal();

    if (rescueState.failure != RESCUE_HEALTHY) {
        // Default to 20s semi-controlled descent with impact detection, then abort
        rescueState.phase = RESCUE_DO_NOTHING;

        switch(gpsRescueConfig()->sanityChecks) {
        case RESCUE_SANITY_ON:
            rescueState.phase = RESCUE_ABORT;
            break;
        case RESCUE_SANITY_FS_ONLY:
            if (hardFailsafe) {
                rescueState.phase = RESCUE_ABORT;
            }
            break;
        default:
            // even with sanity checks off,
            // override when Allow Arming without Fix is enabled without GPS_FIX_HOME and no Control link available.
            if (gpsRescueConfig()->allowArmingWithoutFix && !STATE(GPS_FIX_HOME) && hardFailsafe) {
                rescueState.phase = RESCUE_ABORT;
            }
        }
    }

    // Crash detection is enabled in all rescues.  If triggered, immediately disarm.
    if (crashRecoveryModeActive()) {
        setArmingDisabled(ARMING_DISABLED_ARM_SWITCH);
        disarm(DISARM_REASON_CRASH_PROTECTION);
        rescueStop();
    }

    // Check if GPS comms are healthy
    if (!rescueState.sensor.gpsHealthy) {
        rescueState.failure = RESCUE_GPSLOST;
    }

    //  Things that should run at a low refresh rate (such as flyaway detection, etc) will be checked at 1Hz
    static timeUs_t lastSanityCheck = 0;
    if (!oneSecondPassed(currentTimeUs, &lastSanityCheck)) {
        return;
    }

    // checks that we are getting closer to home.
    // if the quad is stuck, or if GPS data packets stop, there will be no change in distance to home
    // we can't use rescueState.sensor.currentVelocity because it will be held at the last good value if GPS data updates stop
    if (rescueState.phase == RESCUE_FLY_HOME) {
        const float velocityToHomeCmS = previousDistanceToHomeCm - GPS_distanceToHomeCm; // cm/s
        previousDistanceToHomeCm = GPS_distanceToHomeCm;
        rescueState.intent.secondsFailing += (velocityToHomeCmS < 0.1f * rescueState.intent.targetVelocityCmS) ? 1 : -1;
        rescueState.intent.secondsFailing = constrain(rescueState.intent.secondsFailing, 0, 30);
        if (rescueState.intent.secondsFailing >= 30) {
#ifdef USE_MAG
            //If there is an active mag, try again without it
            if (compassEnabledAndCalibrated() && !rescueState.intent.forceDisableMag) {
                //Try again with mag disabled
                rescueState.intent.forceDisableMag = true;
                rescueState.intent.secondsFailing = 0;
            } else
#endif
            {
                rescueState.failure = RESCUE_FLYAWAY;
            }
        }
    }

    secondsLowSats += (!STATE(GPS_FIX) || (gpsSol.numSat < GPS_MIN_SAT_COUNT)) ? 1 : -1;
    secondsLowSats = constrain(secondsLowSats, 0, 10);

    if (secondsLowSats == 10) {
        rescueState.failure = RESCUE_LOWSATS;
    }

    // These conditions ignore sanity mode settings, and apply in all rescues, to handle getting stuck in a climb or descend

    const float actualAltitudeChange = getAltitudeCm() - prevAltitudeCm;
    // ** possibly could use getAltitudeDerivative() for  for actual altitude change, though it is smoothed
    const float targetAltitudeChange = rescueState.intent.targetAltitudeCm - prevTargetAltitudeCm;
    const float ratio = actualAltitudeChange / targetAltitudeChange;
    prevAltitudeCm = getAltitudeCm();
    prevTargetAltitudeCm = rescueState.intent.targetAltitudeCm;

    switch (rescueState.phase) {
    case RESCUE_LANDING:
        rescueState.intent.secondsFailing += ratio > 0.5f ? -1 : 1;
        rescueState.intent.secondsFailing = constrain(rescueState.intent.secondsFailing, 0, 10);
        if (rescueState.intent.secondsFailing >= 10) {
            rescueState.phase = RESCUE_ABORT;
            // Landing mode shouldn't take more than 10s
        }
        break;
    case RESCUE_ATTAIN_ALT:
        rescueState.intent.secondsFailing += ratio > 0.5f ? -1 : 1;
        rescueState.intent.secondsFailing = constrain(rescueState.intent.secondsFailing, 0, 10);
        if (rescueState.intent.secondsFailing >= 10) {
            rescueState.phase = RESCUE_DESCENT;
            rescueState.intent.secondsFailing = 0;
            // if can't climb, go to descend mode
        }
        break;
    case RESCUE_DESCENT:
        rescueState.intent.secondsFailing += ratio > 0.5f ? -1 : 1;
        rescueState.intent.secondsFailing = constrain(rescueState.intent.secondsFailing, 0, 10);
        if (rescueState.intent.secondsFailing >= 20) {
            rescueState.phase = RESCUE_LANDING;
            rescueState.intent.secondsFailing = 0;
            // if can't descend, enable impact detection and go to landing mode
        }
        break;
    case RESCUE_DO_NOTHING:
        secondsDoingNothing = MIN(secondsDoingNothing + 1, 20);
        if (secondsDoingNothing >= 20) {
            rescueState.phase = RESCUE_ABORT;
            // time-limited semi-controlled fall with impact detection
        }
        break;
    default:
        // do nothing
        break;
    }

    DEBUG_SET(DEBUG_RTH, 2, (rescueState.failure * 10 + rescueState.phase));
    DEBUG_SET(DEBUG_RTH, 3, (rescueState.intent.secondsFailing * 100 + secondsLowSats));
}

static void sensorUpdate(bool newGpsData)
{
    rescueState.sensor.gpsHealthy = gpsIsHealthy();

    static float prevDistanceToHomeCm = 0.0f;
    if (newGpsData) {
        rescueState.sensor.velocityToHomeCmS = (prevDistanceToHomeCm - GPS_distanceToHomeCm) * getGpsDataFrequencyHz();
        prevDistanceToHomeCm = GPS_distanceToHomeCm;
        // positive = towards home.  First value is useless since prevDistanceToHomeCm was zero.
    }

    // ** heading values **
    const float bearingToHomeDeg = DECIDEGREES_TO_DEGREES(GPS_directionToHome); // 0 to 360
    const float aircraftHeadingDeg = DECIDEGREES_TO_DEGREES(attitude.values.yaw); // 0 to 360

    rescueState.sensor.errorAngleDeg = aircraftHeadingDeg - bearingToHomeDeg;
    // normalise to -180 ... + 180 for plus and minus yaw corrections
    if (rescueState.sensor.errorAngleDeg <= -180) {
        rescueState.sensor.errorAngleDeg += 360;
    } else if (rescueState.sensor.errorAngleDeg > 180) {
        rescueState.sensor.errorAngleDeg -= 360;
    }

    DEBUG_SET(DEBUG_ATTITUDE, 0, lrintf(aircraftHeadingDeg));
    DEBUG_SET(DEBUG_ATTITUDE, 2, lrintf(rescueState.sensor.velocityToHomeCmS));

    DEBUG_SET(DEBUG_GPS_RESCUE_VELOCITY, 0, lrintf(rescueState.intent.targetVelocityCmS)); // target velocity to home
    DEBUG_SET(DEBUG_GPS_RESCUE_VELOCITY, 1, lrintf(rescueState.sensor.velocityToHomeCmS)); // target velocity to home
    DEBUG_SET(DEBUG_GPS_RESCUE_VELOCITY, 2, lrintf(rescueState.intent.latLonSteps.v[0] * 100.0f)); // latitude step
    DEBUG_SET(DEBUG_GPS_RESCUE_VELOCITY, 3, lrintf(rescueState.intent.latLonSteps.v[1] * 100.0f)); // longitude step

    DEBUG_SET(DEBUG_GPS_RESCUE_HEADING, 0, lrintf(rescueState.sensor.velocityToHomeCmS));
    DEBUG_SET(DEBUG_GPS_RESCUE_HEADING, 1, gpsSol.groundCourse);            // deg * 10
    DEBUG_SET(DEBUG_GPS_RESCUE_HEADING, 2, attitude.values.yaw);            // deg * 10
    DEBUG_SET(DEBUG_GPS_RESCUE_HEADING, 3, GPS_directionToHome);            // deg * 10

    const float currentAltitudeCm = getAltitudeCm();
    DEBUG_SET(DEBUG_GPS_RESCUE_TRACKING, 0, lrintf(rescueState.sensor.velocityToHomeCmS));
    DEBUG_SET(DEBUG_GPS_RESCUE_TRACKING, 2, lrintf(currentAltitudeCm));
    DEBUG_SET(DEBUG_GPS_RESCUE_TRACKING, 4, lrintf(aircraftHeadingDeg));    // estimated heading of the quad (direction nose is pointing in)
    DEBUG_SET(DEBUG_GPS_RESCUE_TRACKING, 5, lrintf(bearingToHomeDeg));      // angle to home derived from GPS location and home position
}

// This function runs every iteration and flashes "RESCUE N/A" in the OSD if any of:
// 1. sensor not healthy - GPS data is not being received.
// 2. GPS has no Home or 3D fix.
// 3. We are relying on GPS for heading, and the IMU is not oriented by prior forward flight
// 4. Sat count < minimum configured or no 3D fix for at least 2s..
// Note 1: cannot arm without the required number of sats, unless bypassed in CLI
// Note 2: the sanity checks are separate, this just provides the OSD warning - maybe can be merged?

static bool checkGPSRescueIsAvailable(void)
{
    if (!gpsIsHealthy() || !STATE(GPS_FIX_HOME) || !isHeadingOK()) {
        return false;
    }

    static int8_t secondsFailing = 0;
    const timeUs_t currentTimeUs = micros();

    //  Check low sats / loss of 3D fix, once per second
    static timeUs_t lastGpsRescueCheck = 0;
    if (!oneSecondPassed(currentTimeUs, &lastGpsRescueCheck)) {
        return secondsFailing < 3;
    }

    const bool failure = gpsSol.numSat < GPS_MIN_SAT_COUNT || !STATE(GPS_FIX);
    secondsFailing = constrain(secondsFailing + (failure ? 1 : -1), 0, 3);
    return secondsFailing < 3;
}

void disarmOnImpact(void)
{
    if (acc.accMagnitude > rescueState.intent.disarmThreshold) {
        setArmingDisabled(ARMING_DISABLED_ARM_SWITCH);
        disarm(DISARM_REASON_GPS_RESCUE);
        rescueStop();
    }
}

void initDescent(void)
{
    rescueState.intent.yawAttenuator = 0.0f; // block yaw while descending
    rescueState.intent.proximityAttenuator = GPS_distanceToHomeCm / rescueState.intent.descentDistanceCm;
}

void descend(bool newGpsData)
{
    if (rescueState.intent.velocityAttenuator < 1.0f) { // acquire velocity over one second if not yet there
        rescueState.intent.velocityAttenuator += taskIntervalSeconds;
    }
    if (newGpsData) {
        rescueState.intent.proximityAttenuator = GPS_distanceToHomeCm / rescueState.intent.descentDistanceCm;
    }
    rescueState.intent.targetVelocityCmS = gpsRescueConfig()->groundSpeedCmS * rescueState.intent.proximityAttenuator * rescueState.intent.velocityAttenuator;

    float altitudeStepCm = taskIntervalSeconds * gpsRescueConfig()->descendRate;

    // at or below 10m: descend at 0.6x set value; above 10m, descend faster, to max 3.0x at 50m
    altitudeStepCm *= scaleRangef(constrainf(rescueState.intent.targetAltitudeCm, 1000, 5000), 1000, 5000, 0.6f, 3.0f);

    rescueState.intent.targetAltitudeStepCm = -altitudeStepCm;
    rescueState.intent.targetAltitudeCm -= altitudeStepCm;
}

void initialiseRescueValues (void)
{
    if (GPS_distanceToHomeCm < gpsRescueConfig()->minStartDistM * 100.0f) {
        rescueState.intent.returnAltitudeCm = fmaxf(500.0f, getAltitudeCm() + (gpsRescueConfig()->initialClimbM * 100.0f));
        // climb above current height by buffer height, to at least 5m altitude
        // set the descent distance to current distance, noting this could be zero
    }

    rescueState.sensor.imuYawCogGain = 1.0f;
    rescueState.intent.secondsFailing = 0; // reset the sanity check timer
    rescueState.intent.yawAttenuator = 0.0f; // no yaw in the climb
    rescueState.intent.velocityAttenuator = 0.0f; // control velocity acquisition
    rescueState.intent.targetVelocityCmS = 0.0f; // stop the quad immediately
    rescueState.intent.targetAltitudeStepCm = 0.0f;
    rescueState.sensor.velocityToHomeCmS = 0.0f;
    vector2Zero(&rescueState.intent.latLonSteps);
    rescueState.intent.forceDisableMag = false; // re-enable Mag on next rescue start even if it failed on a previous rescue
}

void gpsRescueUpdate(void)
// runs at gpsRescueTaskIntervalSeconds, and runs whether or not rescue is active
{
    static uint16_t gpsStamp = 0;
    bool newGpsData = gpsHasNewData(&gpsStamp);

    if (!FLIGHT_MODE(GPS_RESCUE_MODE)) {
        rescueStop(); // sets phase to RESCUE_IDLE; does nothing else.  RESCUE_IDLE tasks still run.
    } else if (FLIGHT_MODE(GPS_RESCUE_MODE) && rescueState.phase == RESCUE_IDLE) {
        rescueStart(); // sets phase to rescue_initialise if we enter GPS Rescue mode while idle
        rescueAttainPosition(false); // Initialise basic parameters when a Rescue starts (can't initialise sensor data reliably)
        performSanityChecks(); // Initialises sanity check values when a Rescue starts
    }
    // Will now be in RESCUE_INITIALIZE mode, if just entered Rescue while IDLE, otherwise stays IDLE

    sensorUpdate(newGpsData); // always get latest GPS and Altitude data, update ascend and descend rates

    static bool returnAltitudeLow = true;
    rescueState.isAvailable = checkGPSRescueIsAvailable();

    switch (rescueState.phase) {
    case RESCUE_IDLE:
        // in Idle phase = NOT in GPS Rescue
        // update the return altitude and descent distance values, to have valid settings immediately they are needed
        setReturnAltitude();
        break;
        // sanity checks are bypassed in IDLE mode; instead, failure state is always initialised to HEALTHY
        // target altitude is always set to current altitude.

    case RESCUE_INITIALIZE:
        // Things that should be done at the start of a Rescue
        if (!STATE(GPS_FIX_HOME)) {
            // we didn't get a home point on arming
            rescueState.failure = RESCUE_NO_HOME_POINT;
            // will result in a disarm via the sanity check system, or float around if switch induced
            break;
        }

        if (GPS_distanceToHome < 5 && isBelowLandingAltitude()) {
            // attempted initiation within 5m of home, and 'on the ground' -> instant disarm, for safety reasons
            rescueState.phase = RESCUE_ABORT;
            break;
        }

        initialiseRescueValues(); // initialise the control related values
        returnAltitudeLow = getAltitudeCm() < rescueState.intent.returnAltitudeCm;
        rescueState.phase = RESCUE_ATTAIN_ALT;
        rescueState.failure = RESCUE_HEALTHY;
        break;

    case RESCUE_ATTAIN_ALT:
        // gradually increment the target altitude until the craft reaches target altitude
        // note that this can mean the target altitude may increase above returnAltitude if the craft lags target
        // sanity check will abort if altitude gain is blocked for a cumulative period
        if (returnAltitudeLow == (getAltitudeCm() < rescueState.intent.returnAltitudeCm)) {
            // we started low, and still are low; also true if we started high, and still are too high
            rescueState.intent.targetAltitudeStepCm = (returnAltitudeLow ? gpsRescueConfig()->ascendRate : -1.0f * gpsRescueConfig()->descendRate) * taskIntervalSeconds;
            rescueState.intent.targetAltitudeCm += rescueState.intent.targetAltitudeStepCm;
        } else {
            // target altitude achieved - move on to ROTATE phase, returning at target altitude
            rescueState.intent.targetAltitudeCm = rescueState.intent.returnAltitudeCm;
            rescueState.intent.targetAltitudeStepCm = 0.0f;
            if (isHeadingOK()) {
                rescueState.phase = RESCUE_ROTATE;
            } else {
                rescueState.phase = RESCUE_ORIENT_IMU; // only for emergency imu orientation
            }
        }
        break;

    case RESCUE_ORIENT_IMU:
        if (isHeadingOK()) {
            rescueState.phase = RESCUE_ROTATE;
        } else {
            // fixed pitch forward for emergency IMU orientation to GPS course over ground
            autopilotAngle[AI_ROLL] = 0.0f;
            autopilotAngle[AI_PITCH] = 35.0f;
            rescueYaw = 0.0f;
        }
        break;

    case RESCUE_ROTATE:
        if (rescueState.intent.yawAttenuator < 1.0f) { // acquire yaw authority over one second
            rescueState.intent.yawAttenuator += taskIntervalSeconds;
        }
        if (fabsf(rescueState.sensor.errorAngleDeg) < GPS_RESCUE_ALLOWED_YAW_RANGE) {
            // yaw angle error is small enough to allow us to enter fly home or descend phase
            // but first check that, if there is no mag, the imu is oriented
            if (GPS_distanceToHomeCm < rescueState.intent.descentDistanceCm) {
                initDescent();
                rescueState.phase = RESCUE_DESCENT; // enter descend phase
            } else {
                rescueState.phase = RESCUE_FLY_HOME; // enter fly home phase
            }
            rescueState.intent.secondsFailing = 0; // reset sanity timer for flight home
        }
        break;

    case RESCUE_FLY_HOME:
        if (rescueState.intent.yawAttenuator < 1.0f) { // be sure to accumulate full yaw authority
            rescueState.intent.yawAttenuator += taskIntervalSeconds;
        }
        if (rescueState.intent.velocityAttenuator < 1.0f) { // acquire velocity over one second
            rescueState.intent.velocityAttenuator += taskIntervalSeconds;
        }
        rescueState.intent.targetVelocityCmS = gpsRescueConfig()->groundSpeedCmS * rescueState.intent.velocityAttenuator;
        if (newGpsData) {
            if (GPS_distanceToHomeCm < rescueState.intent.descentDistanceCm) {
                initDescent();
                rescueState.phase = RESCUE_DESCENT; // enter descend phase
                rescueState.intent.secondsFailing = 0; // reset sanity timer for descent
            }
        }
        break;

    case RESCUE_DESCENT:
        if (isBelowLandingAltitude()) {
            // enter landing mode once below landing altitude
            rescueState.phase = RESCUE_LANDING;
            rescueState.intent.secondsFailing = 0; // reset sanity timer for landing
        }
        descend(newGpsData);
        break;

    case RESCUE_LANDING:
        // Reduce altitude target steadily until impact, then disarm.
        // control yaw angle and throttle and pitch, attenuate velocity, roll and pitch iTerm
        // increase velocity smoothing cutoff as we get closer to ground
        descend(newGpsData);
        disarmOnImpact();
        break;

    case RESCUE_ABORT:
        setArmingDisabled(ARMING_DISABLED_ARM_SWITCH);
        disarm(DISARM_REASON_FAILSAFE);
        rescueState.intent.secondsFailing = 0; // reset sanity timers so we can re-arm
        rescueStop();
        break;

    case RESCUE_DO_NOTHING:
        disarmOnImpact();
        break;

    default:
        break;
    }

    DEBUG_SET(DEBUG_GPS_RESCUE_VELOCITY, 1, rescueState.phase);
    DEBUG_SET(DEBUG_GPS_RESCUE_TRACKING, 3, lrintf(rescueState.intent.targetAltitudeCm));
    DEBUG_SET(DEBUG_RTH, 0, lrintf(rescueState.intent.maxAltitudeCm / 10.0f));

    performSanityChecks();
    rescueAttainPosition(newGpsData);
}

float gpsRescueGetYawRate(void)
{
    return rescueYaw; // the control yaw value for rc.c to be used while flightMode gps_rescue is active.
}

bool gpsRescueIsConfigured(void)
{
    return failsafeConfig()->failsafe_procedure == FAILSAFE_PROCEDURE_GPS_RESCUE || isModeActivationConditionPresent(BOXGPSRESCUE);
}

bool gpsRescueIsAvailable(void)
{
    return rescueState.isAvailable; // flashes the warning when not available (low sats, )
}

bool gpsRescueIsDisabled(void)
// used for OSD warning
{
    return (!STATE(GPS_FIX_HOME));
}

#ifdef USE_MAG
bool gpsRescueDisableMag(void)
{
    // Enable mag on user request, but don't use it during fly home or if force disabled
    // Note that while flying home the course over ground from GPS provides a heading that is less affected by wind
    return rescueState.intent.forceDisableMag;
}
#endif
#endif
