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
#include "sensors/acceleration.h"
#include "sensors/compass.h"

#include "gps_rescue.h"

typedef enum {
    RESCUE_IDLE,
    RESCUE_INITIALIZE,
    RESCUE_ATTAIN_ALT,
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
    float descentDistanceM;
    int8_t secondsFailing;
    float verticalVelocityMultiplier;
    float yawAttenuator;
    float disarmThreshold;
    float latFactor;
    float lonFactor;
} rescueIntent_s;

typedef struct {
    uint16_t groundSpeedCmS;
    bool healthy;
    float errorAngleDeg;
    float gpsDataIntervalSeconds;
    float velocityToHomeCmS;
    float imuYawCogGain;
    gpsLocation_t previousLocation;
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
bool magForceDisable = false;

rescueState_s rescueState;

void gpsRescueInit(void)
{
    // currently no filters to initialise
}

// check for new GPS Data
static bool newGPSData = false;
static uint16_t previousGpsStamp = ~0;
static void gpsRescueNewGpsData(void)
{
    if (getGpsStamp() != previousGpsStamp) {
        previousGpsStamp = getGpsStamp();
        newGPSData = true;
    }
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
    // Hold maxAltitude at zero while disarmed, but if set_home_point_once is true, hold maxAlt until power cycled
    if (!ARMING_FLAG(ARMED) && !gpsConfig()->gps_set_home_point_once) {
        rescueState.intent.maxAltitudeCm = 0.0f;
        return;
    }

    // While armed...
    rescueState.intent.maxAltitudeCm = fmaxf(getAltitudeCm(), rescueState.intent.maxAltitudeCm);

    if (newGPSData) {
        // set the target altitude to the current altitude.
        rescueState.intent.targetAltitudeCm = getAltitudeCm();

        // Set descend distance to the minimum of the user's preference or half the distance to home
        rescueState.intent.descentDistanceM = fminf(0.5f * GPS_distanceToHome, gpsRescueConfig()->descentDistanceM);

        const float initialClimbCm = gpsRescueConfig()->initialClimbM * 100.0f;
        switch (gpsRescueConfig()->altitudeMode) {
            case GPS_RESCUE_ALT_MODE_FIXED:
                rescueState.intent.returnAltitudeCm = gpsRescueConfig()->returnAltitudeM * 100.0f;
                break;
            case GPS_RESCUE_ALT_MODE_CURRENT:
                // climb above current altitude, but always return at least initial height above takeoff point, in case current altitude was negative
                rescueState.intent.returnAltitudeCm = fmaxf(initialClimbCm, getAltitudeCm() + initialClimbCm);
                break;
            case GPS_RESCUE_ALT_MODE_MAX:
            default:
                rescueState.intent.returnAltitudeCm = rescueState.intent.maxAltitudeCm + initialClimbCm;
                break;
        }
    }
}

static void rescueAttainPosition(void)
{
    // runs at 100hz, but only updates RPYT settings when new GPS Data arrives and when not in idle phase.
    switch (rescueState.phase) {
    case RESCUE_IDLE:
        // do nothing
        return;
    case RESCUE_INITIALIZE:
        // initialise the positioning related settings
        rescueState.intent.disarmThreshold = gpsRescueConfig()->disarmThreshold * 0.1f;
        // initialise the required autopilot functions
        resetAltitudeControl();
        resetPositionControl(gpsSol.llh);
        // pre-calculate the latitude and longitude adjustment factors
        rescueState.intent.latFactor = cos_approx(DEGREES_TO_RADIANS(DECIDEGREES_TO_DEGREES(GPS_directionToHome))) / EARTH_ANGLE_TO_CM ;
        rescueState.intent.lonFactor = sin_approx(DEGREES_TO_RADIANS(DECIDEGREES_TO_DEGREES(GPS_directionToHome))) / (EARTH_ANGLE_TO_CM * getGpsCosLat());
        rescueState.sensor.imuYawCogGain = 1.0f;
        rescueState.sensor.previousLocation = gpsSol.llh;
        return;
    case RESCUE_DO_NOTHING:
        // 20s of hover at current altitude, for switch induced sanity failures, to allow time to recover
        // do nothing
        setTargetLocation(gpsSol.llh);
        return;
     default:
        break;
    }

    /**
        Altitude (throttle) controller
    */
    const float verticalVelocity = getAltitudeDerivative() * rescueState.intent.verticalVelocityMultiplier;
    altitudeControl(rescueState.intent.targetAltitudeCm, taskIntervalSeconds, verticalVelocity, rescueState.intent.targetAltitudeStepCm);

    /**
        Heading / yaw controller
    */
    // simple yaw P controller with roll mixed in.
    // attitude.values.yaw is set by imuCalculateEstimatedAttitude() and is updated from GPS while groundspeed exceeds 2 m/s
    // below 2m/s groundspeed, the IMU uses gyro to estimate yaw attitude change from previous values
    // above 2m/s, GPS course over ground us ysed to 'correct' the IMU heading
    // if the course over ground, due to wind or pre-exiting movement, is different from the attitude of the quad, the GPS correction will be less accurate
    // the craft should not return much less than 5m/s during the rescue or the GPS corrections may be inaccurate.
    // the faster the return speed, the more accurate the IMU will be, but the consequences of IMU error at the start are greater
    // A compass (magnetometer) is vital for accurate GPS rescue at slow speeds, but must be calibrated and validated
    // WARNING:  Some GPS units give false Home values!  Always check the arrow points to home on leaving home.
    rescueYaw = rescueState.sensor.errorAngleDeg * gpsRescueConfig()->yawP * rescueState.intent.yawAttenuator / 10.0f;
    rescueYaw = constrainf(rescueYaw, -GPS_RESCUE_MAX_YAW_RATE, GPS_RESCUE_MAX_YAW_RATE);
    // rescueYaw is the yaw rate in deg/s to correct the heading error

    rescueYaw *= GET_DIRECTION(rcControlsConfig()->yaw_control_reversed);
    // rescueYaw is the yaw rate in deg/s to correct the heading error

    DEBUG_SET(DEBUG_GPS_RESCUE_HEADING, 7, rescueYaw);                  // the yaw rate in deg/s to correct a yaw error

    /*
        Pitch / velocity controller
    */
    if (newGPSData) {
        // TODO / queries
        // limit to pitch only responses in autopilot when IMU is disoriented?
        // check modifying D at start and on descent
        // check no iTerm at start

        float distanceToMove = rescueState.intent.targetVelocityCmS * rescueState.sensor.gpsDataIntervalSeconds;
        static bool arrivedHome = false;
        gpsLocation_t newLocation;

        if (!arrivedHome) {
            if (GPS_distanceToHomeCm < 10) {
                // within 10cm of home
                rescueState.intent.targetVelocityCmS = 0.0f;
                arrivedHome = true;
                // lock autopilot target at home location
                setTargetLocation(GPS_home_llh);
            }
        } else if (distanceToMove > 0.0f) {
            // Calculate the required change in latitude and longitude based on the bearing
            // if no requested distance, target location doesn't change
            // Update the new location so that we move along the intended flightpath
            newLocation.lat = rescueState.sensor.previousLocation.lat + (int32_t)(distanceToMove * rescueState.intent.latFactor);
            newLocation.lon = rescueState.sensor.previousLocation.lon + (int32_t)(distanceToMove * rescueState.intent.lonFactor);
            rescueState.sensor.previousLocation = newLocation;
            setTargetLocation(newLocation); // update the autopilot target location to the new location

       }
        DEBUG_SET(DEBUG_GPS_RESCUE_VELOCITY, 1, lrintf(distanceToMove));
        DEBUG_SET(DEBUG_GPS_RESCUE_VELOCITY, 3, lrintf(rescueState.intent.targetVelocityCmS)); // target velocity to home
        DEBUG_SET(DEBUG_GPS_RESCUE_VELOCITY, 4, lrintf(distanceToMove * rescueState.intent.lonFactor)); // all unused at present
        DEBUG_SET(DEBUG_GPS_RESCUE_VELOCITY, 2, lrintf(distanceToMove * rescueState.intent.latFactor)); // all unused at present

        posControlOnNewGPSData();
    }
    posControlOutput(); // upsample the setpoints for pid.c
    DEBUG_SET(DEBUG_GPS_RESCUE_VELOCITY, 0, GPS_distanceToHomeCm);
    DEBUG_SET(DEBUG_GPS_RESCUE_TRACKING, 1, lrintf(rescueState.intent.targetVelocityCmS)); // target velocity to home
}

static void performSanityChecks(void)
{
    static timeUs_t previousTimeUs = 0; // Last time Stalled/LowSat was checked
    static float prevAltitudeCm = 0.0f; // to calculate ascent or descent change
    static float prevTargetAltitudeCm = 0.0f; // to calculate ascent or descent target change
    static float previousDistanceToHomeCm = 0.0f; // to check that we are returning
    static int8_t secondsLowSats = 0; // Minimum sat detection
    static int8_t secondsDoingNothing; // Limit on doing nothing
    const timeUs_t currentTimeUs = micros();

    if (rescueState.phase == RESCUE_IDLE) {
        rescueState.failure = RESCUE_HEALTHY;
        return;
    } else if (rescueState.phase == RESCUE_INITIALIZE) {
        // Initialize these variables each time a GPS Rescue is started
        previousTimeUs = currentTimeUs;
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
    // ToDo - check if we have an altitude reading; if we have Baro, we can use Landing mode for controlled descent without GPS
    if (!rescueState.sensor.healthy) {
        rescueState.failure = RESCUE_GPSLOST;
    }

    //  Things that should run at a low refresh rate (such as flyaway detection, etc) will be checked at 1Hz
    const timeDelta_t dTime = cmpTimeUs(currentTimeUs, previousTimeUs);
    if (dTime < 1000000) { //1hz
        return;
    }
    previousTimeUs = currentTimeUs;

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
            //If there is a mag and has not been disabled, try again without the mag
if (sensors(SENSOR_MAG) && gpsRescueConfig()->useMag && !magForceDisable) {
                //Try again with mag disabled
                magForceDisable = true;
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
    case RESCUE_DESCENT:
        rescueState.intent.secondsFailing += ratio > 0.5f ? -1 : 1;
        rescueState.intent.secondsFailing = constrain(rescueState.intent.secondsFailing, 0, 10);
        if (rescueState.intent.secondsFailing >= 10) {
            rescueState.phase = RESCUE_LANDING;
            rescueState.intent.secondsFailing = 0;
            // if can't climb, or slow descending, enable impact detection and time out in 10s
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

static void sensorUpdate(void)
{
    rescueState.sensor.healthy = gpsIsHealthy();
    gpsRescueNewGpsData(); // set newGPSData if there is new GPS data
    const float currentAltitudeCm = getAltitudeCm();

    const float bearingToHomeDeg = DECIDEGREES_TO_DEGREES(GPS_directionToHome); // 0 to 360
    const float aircraftHeadingDeg = DECIDEGREES_TO_DEGREES(attitude.values.yaw); // 0 to 360
    const float groundCourseDeg = DECIDEGREES_TO_DEGREES(gpsSol.groundCourse);    // 0 to 360

    // assess relationship between aircraft heading and GPS ground course
    // to determine which heading to use in case of IMU disorientation
    float headingErrorDeg = aircraftHeadingDeg - groundCourseDeg;
    // normalise to -180 ... + 180
    if (headingErrorDeg > 180) {
        headingErrorDeg -= 360;
    } else if (headingErrorDeg < -180) {
        headingErrorDeg += 360;
    }

    float headingVsCogFactor = fabsf(cos_approx(DEGREES_TO_RADIANS(headingErrorDeg)) - 1.0f);
    // 0 = aligned (0 degrees error)
    // 1 = 90 degrees off
    // 2 = 180 degrees off

    float headingSpeedFactor = gpsRescueConfig()->groundSpeedCmS ? (gpsSol.groundSpeed / gpsRescueConfig()->groundSpeedCmS) : 0.0f;
    headingSpeedFactor = fmaxf(headingSpeedFactor - 1.0f, 0.0f);
    // 0.0 until target groundspeed, 1.0 at double target groundspeed, 2.0 at three times
    // above zero strongly implies IMU error or strong drift

    headingVsCogFactor *= headingSpeedFactor; // zero is well aligned at inside normal speed
    // above 1 implies IMU error, requiring an angle error and a speed overshoot;
    // 4 at double target speed and 180 error
    // note : this factor may be useful for cogYawGain and to suppress Roll in IMU error situations

    const float headingVsCogFader = fminf(headingVsCogFactor, 1.0f);

    // when closer to 0, use GPS groundcourse for heading, since there must be drift or the IMU is incorrect or groundspeed is low
    float headingToUse = aircraftHeadingDeg * (1.0f - headingVsCogFader) + groundCourseDeg * headingVsCogFader;
    // prefer the IMU heading when aligned with groundcourse, or moving slowly, otherwise use groundcourse


// don't think this is needed

//     if (headingToUse < 0) {
//         headingToUse += 360;
//     }
//     if (headingToUse > 360) {
//         headingToUse -= 360;
//     }



    // for now, do nothing with headingToUse
    //    rescueState.sensor.errorAngleDeg = (rescueState.phase == RESCUE_FLY_HOME) ? (headingToUse - bearingToHomeDeg) : aircraftHeadingDeg - bearingToHomeDeg;

    rescueState.sensor.errorAngleDeg = aircraftHeadingDeg - bearingToHomeDeg;
    // normalise to -180 ... + 180
    if (rescueState.sensor.errorAngleDeg <= -180) {
        rescueState.sensor.errorAngleDeg += 360;
    } else if (rescueState.sensor.errorAngleDeg > 180) {
        rescueState.sensor.errorAngleDeg -= 360;
    }

    static float prevDistanceToHomeCm = 0.0f;
    if (newGPSData) {
        rescueState.sensor.gpsDataIntervalSeconds = getGpsDataIntervalSeconds();
        // Range from 10ms (100hz) to 1000ms (1Hz). Intended to cover common GPS data rates and exclude unusual values.

        // *** all the following is now only needed to calculate imuYawCogGain for GPS-only rescues
        rescueState.sensor.velocityToHomeCmS = ((prevDistanceToHomeCm - GPS_distanceToHomeCm) / rescueState.sensor.gpsDataIntervalSeconds);
        prevDistanceToHomeCm = GPS_distanceToHomeCm;
        // positive = towards home.  First value is useless since prevDistanceToHomeCm was zero.

        // Handle IMU disorientation
        // If disoriented, thinks it is flying in the direction of home, pitches forward, but distance to home never closes
        // iTerm accumulates and PIDs pitch progressively further forward so velocity increases
        // we need to detect this and encourage the IMU to correct.
        if (gpsRescueConfig()->groundSpeedCmS) {
            const float sensitivity = gpsRescueConfig()->groundSpeedCmS;
            const float groundspeedErrorRatio = fabsf(gpsSol.groundSpeed - rescueState.sensor.velocityToHomeCmS) / sensitivity;
            // 0 if groundspeed = velocity to home, or both are zero, meaning OK
            // 1 if forward velocity is zero but sideways speed equal to return speed
            // 2 if moving backwards at expected return speed, 4 if moving backwards at 2* expected return speed, etc

            const float pitchForwardFactor = (attitude.values.pitch > 0.0f) ? fminf(attitude.values.pitch / 30.0f, 2.0f) : 0.0f;
            // 0 when flat or pitched back, 1.0 at 30 degrees, 2.0 at 60 degrees

            if (rescueState.phase == RESCUE_FLY_HOME) {
                rescueState.sensor.imuYawCogGain = pitchForwardFactor + fminf(groundspeedErrorRatio, 3.5f);
                // imuYawCogGain will be more positive at higher pitch angles and higher groundspeed errors
                // imuYawCogGain will be lowest (close to zero) at lower pitch angles and when flying straight towards home
            } else {
                rescueState.sensor.imuYawCogGain = 0.0f;
                // accept whatever IMU state we had, until fly home phase
            }
            DEBUG_SET(DEBUG_ATTITUDE, 5, lrintf(groundspeedErrorRatio * 100));
            DEBUG_SET(DEBUG_ATTITUDE, 6, lrintf(rescueState.sensor.imuYawCogGain * 10));
        }
    }

    DEBUG_SET(DEBUG_ATTITUDE, 4, rescueState.sensor.velocityToHomeCmS); // velocity to home

    DEBUG_SET(DEBUG_GPS_RESCUE_HEADING, 0, lrintf(headingErrorDeg)); 
    DEBUG_SET(DEBUG_GPS_RESCUE_HEADING, 1, lrintf(headingVsCogFactor * 100.0f)); 
    DEBUG_SET(DEBUG_GPS_RESCUE_HEADING, 2, lrintf(aircraftHeadingDeg));        // degrees * 10
    DEBUG_SET(DEBUG_GPS_RESCUE_HEADING, 3, lrintf(bearingToHomeDeg)); // computed from current GPS position in relation to home
    DEBUG_SET(DEBUG_GPS_RESCUE_HEADING, 4, lrintf(headingToUse)); // computed from current GPS position in relation to home
    DEBUG_SET(DEBUG_GPS_RESCUE_HEADING, 5, gpsSol.groundCourse / 10);  // heading value used for yaw angle
    DEBUG_SET(DEBUG_GPS_RESCUE_HEADING, 6, lrintf(rescueState.sensor.imuYawCogGain * 10.0f));  // heading value used for yaw angle

    DEBUG_SET(DEBUG_GPS_RESCUE_TRACKING, 2, lrintf(currentAltitudeCm));
    DEBUG_SET(DEBUG_GPS_RESCUE_TRACKING, 4, lrintf(aircraftHeadingDeg));                 // estimated heading of the quad (direction nose is pointing in)
    DEBUG_SET(DEBUG_GPS_RESCUE_TRACKING, 5, lrintf(bearingToHomeDeg));  // angle to home derived from GPS location and home position
    DEBUG_SET(DEBUG_GPS_RESCUE_TRACKING, 0, lrintf(rescueState.sensor.velocityToHomeCmS));
}

// This function runs every iteration and flashes "RESCUE N/A" in the OSD if any of:
// 1. sensor not healthy - GPS data is being received.
// 2. GPS has no 3D fix.
// 3. GPS number of satellites is greater than or equal to the minimum configured satellite count.
// Note 1: cannot arm without the required number of sats
// hence this flashing indicates that after having enough sats, we now have below the minimum and the rescue likely would fail
// Note 2: this function does not take into account the distance from home
// The sanity checks are independent, this just provides the OSD warning
static bool checkGPSRescueIsAvailable(void)
{
    static timeUs_t previousTimeUs = 0; // Last time LowSat was checked
    const timeUs_t currentTimeUs = micros();
    static int8_t secondsLowSats = 0; // Minimum sat detection
    static bool lowsats = false;
    static bool noGPSfix = false;
    bool result = true;

    if (!gpsIsHealthy() || !STATE(GPS_FIX_HOME)) {
        return false;
    }

    //  Things that should run at a low refresh rate >> ~1hz
    const timeDelta_t dTime = cmpTimeUs(currentTimeUs, previousTimeUs);
    if (dTime < 1000000) { //1hz
        if (noGPSfix || lowsats) {
            result = false;
        }
        return result;
    }

    previousTimeUs = currentTimeUs;

    if (!STATE(GPS_FIX)) {
        result = false;
        noGPSfix = true;
    } else {
        noGPSfix = false;
    }

    secondsLowSats = constrain(secondsLowSats + ((gpsSol.numSat < GPS_MIN_SAT_COUNT) ? 1 : -1), 0, 2);
    if (secondsLowSats == 2) {
        lowsats = true;
        result = false;
    } else {
        lowsats = false;
    }

    return result;
}

void disarmOnImpact(void)
{
    if (acc.accMagnitude > rescueState.intent.disarmThreshold) {
        setArmingDisabled(ARMING_DISABLED_ARM_SWITCH);
        disarm(DISARM_REASON_GPS_RESCUE);
        rescueStop();
    }
}

void descend(void)
{
    if (newGPSData) {
        const float velocityAttenuator = fminf(GPS_distanceToHome / rescueState.intent.descentDistanceM, 1.0f);
        rescueState.intent.targetVelocityCmS = gpsRescueConfig()->groundSpeedCmS * velocityAttenuator;
    }

    // ensure we have full yaw authority in case we entered descent mode without enough time in fly home to acquire it gracefully
    rescueState.intent.yawAttenuator = 1.0f;

    // set the altitude step, considering the interval between altitude readings and the descent rate
    float altitudeStepCm = taskIntervalSeconds * gpsRescueConfig()->descendRate;

    // descend more slowly if the return home altitude was less than 20m
    const float descentRateAttenuator = constrainf (rescueState.intent.returnAltitudeCm / 2000.0f, 0.25f, 1.0f);
    altitudeStepCm *= descentRateAttenuator;
    // slowest descent rate will be 1/4 of normal when we start descending at or below 5m above take-off point.
    // otherwise a rescue initiated very low and close may not get all the way home

    // descend faster while the quad is at higher altitudes
    const float descentRateMultiplier = constrainf(rescueState.intent.targetAltitudeCm / 5000.0f, 0.0f, 1.0f);
    altitudeStepCm *= 1.0f + (2.0f * descentRateMultiplier);
    // maximum descent rate increase is 3x default above 50m, 2x above 25m, 1.2x at 5m, default by ground level

    // also increase throttle D up to 2x in the descent phase when altitude descent rate is faster, for better control
    rescueState.intent.verticalVelocityMultiplier = 1.0f + descentRateMultiplier;

    rescueState.intent.targetAltitudeStepCm = -altitudeStepCm;
    rescueState.intent.targetAltitudeCm += rescueState.intent.targetAltitudeStepCm;
}

void initialiseRescueValues (void)
{
    rescueState.intent.secondsFailing = 0; // reset the sanity check timer
    rescueState.intent.yawAttenuator = 0.0f; // no yaw in the climb
    rescueState.intent.targetVelocityCmS = 0.0f; // might as well stop the quad immediately
    rescueState.intent.verticalVelocityMultiplier = 1.0f;
    rescueState.intent.targetAltitudeStepCm = 0.0f;
}

void gpsRescueUpdate(void)
// runs at gpsRescueTaskIntervalSeconds, and runs whether or not rescue is active
{
    if (!FLIGHT_MODE(GPS_RESCUE_MODE)) {
        rescueStop(); // sets phase to RESCUE_IDLE; does nothing else.  RESCUE_IDLE tasks still run.
    } else if (FLIGHT_MODE(GPS_RESCUE_MODE) && rescueState.phase == RESCUE_IDLE) {
        rescueStart(); // sets phase to rescue_initialise if we enter GPS Rescue mode while idle
        rescueAttainPosition(); // Initialise basic parameters when a Rescue starts (can't initialise sensor data reliably)
        performSanityChecks(); // Initialises sanity check values when a Rescue starts
    }

    // Will now be in RESCUE_INITIALIZE mode, if just entered Rescue while IDLE, otherwise stays IDLE

    sensorUpdate(); // always get latest GPS and Altitude data, update ascend and descend rates

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
        } else {
            if (GPS_distanceToHome < 5.0f && isBelowLandingAltitude()) {
                // attempted initiation within 5m of home, and 'on the ground' -> instant disarm, for safety reasons
                rescueState.phase = RESCUE_ABORT;
            } else {
                if (GPS_distanceToHome < gpsRescueConfig()->minStartDistM) {
                    rescueState.intent.returnAltitudeCm = fmaxf(500.0f, getAltitudeCm() + (gpsRescueConfig()->initialClimbM * 100.0f));
                    // climb above current height by buffer height, to at least 5m altitude
                    rescueState.intent.descentDistanceM = GPS_distanceToHome;
                    // set the descent distance to current distance, whatever that is
                }
                // otherwise behave as for a normal rescue
                initialiseRescueValues(); // initialise the control related values
                returnAltitudeLow = getAltitudeCm() < rescueState.intent.returnAltitudeCm;
                rescueState.phase = RESCUE_ATTAIN_ALT;
                rescueState.failure = RESCUE_HEALTHY;
            }
        }
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


            // **** TO DO:  only applies if no Mag, or no GPS lock
            // if initiated too close, do not rotate or do anything else until sufficiently far away that a 'normal' rescue can happen
            if (GPS_distanceToHome > gpsRescueConfig()->minStartDistM) {
                rescueState.phase = RESCUE_ROTATE;
            }



        }

        // stop the quad
        rescueState.intent.targetVelocityCmS = 0.0f;
        break;

    case RESCUE_ROTATE:
        if (rescueState.intent.yawAttenuator < 1.0f) { // acquire yaw authority over one second
            rescueState.intent.yawAttenuator += taskIntervalSeconds;
        }
        if (fabsf(rescueState.sensor.errorAngleDeg) < GPS_RESCUE_ALLOWED_YAW_RANGE) {
            // enter fly home phase, and enable pitch, when the yaw angle error is small enough
            rescueState.phase = RESCUE_FLY_HOME; // enter fly home phase
            rescueState.intent.secondsFailing = 0; // reset sanity timer for flight home
        }
        rescueState.intent.targetVelocityCmS = 0.0f;
        break;

    case RESCUE_FLY_HOME:
        if (rescueState.intent.yawAttenuator < 1.0f) { // be sure to accumulate full yaw authority
            rescueState.intent.yawAttenuator += taskIntervalSeconds;
        }
        if (newGPSData) {
            rescueState.intent.targetVelocityCmS = gpsRescueConfig()->groundSpeedCmS;
            if (GPS_distanceToHome <= rescueState.intent.descentDistanceM) {
                rescueState.phase = RESCUE_DESCENT;
                rescueState.intent.secondsFailing = 0; // reset sanity timer for descent
            }
        }
        break;

    case RESCUE_DESCENT:
        // attenuate velocity and altitude targets while updating the heading to home
        if (isBelowLandingAltitude()) {
            // enter landing mode once below landing altitude
            rescueState.phase = RESCUE_LANDING;
            rescueState.intent.secondsFailing = 0; // reset sanity timer for landing
        }
        descend();
        break;

    case RESCUE_LANDING:
        // Reduce altitude target steadily until impact, then disarm.
        // control yaw angle and throttle and pitch, attenuate velocity, roll and pitch iTerm
        // increase velocity smoothing cutoff as we get closer to ground
        descend();
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

    DEBUG_SET(DEBUG_GPS_RESCUE_TRACKING, 3, lrintf(rescueState.intent.targetAltitudeCm));
    DEBUG_SET(DEBUG_RTH, 0, lrintf(rescueState.intent.maxAltitudeCm / 10.0f));

    performSanityChecks();
    rescueAttainPosition();

    newGPSData = false;
}

float gpsRescueGetYawRate(void)
{
    return rescueYaw; // the control yaw value for rc.c to be used while flightMode gps_rescue is active.
}

float gpsRescueGetImuYawCogGain(void)
{
    return rescueState.sensor.imuYawCogGain; // to speed up the IMU orientation to COG when needed
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
    return !(gpsRescueConfig()->useMag && rescueState.phase != RESCUE_FLY_HOME && !magForceDisable);
}
#endif
#endif
