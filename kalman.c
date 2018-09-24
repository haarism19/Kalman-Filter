 /*
  * The library is not extensively tested and only
  * meant as a simple explanation and for inspiration.
  * NO WARRANTY of ANY KIND is provided.
  */

// Original code in C++. Ported to C, to be used by common MCU


#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include "kalman.h"


void kalman_init(kalman * p_kalman) {

    /* We will set the variables like so, these can also be tuned by the user */
    
    p_kalman->Q_angle = 0.001f;
    p_kalman->Q_bias = 0.003f;
    p_kalman->R_measure = 0.03f;

    p_kalman->angle = 0.0f; // Reset the angle
    p_kalman->bias = 0.0f; // Reset bias

    p_kalman->P[0][0] = 0.0f; // Since we assume that the bias is 0 and we know the starting angle (use setAngle), the error covariance matrix is set like so - see: http://en.wikipedia.org/wiki/Kalman_filter#Example_application.2C_technical
    p_kalman->P[0][1] = 0.0f;
    p_kalman->P[1][0] = 0.0f;
    p_kalman->P[1][1] = 0.0f;
	
}


float kalman_get_angle(kalman * p_kalman, float newAngle, float newRate, float dt){

    p_kalman->rate = newRate - p_kalman->bias;
    p_kalman->angle += dt * p_kalman->rate;

    // Update estimation error covariance - Project the error covariance ahead
    
    p_kalman->P[0][0] += dt * (dt*p_kalman->P[1][1] - p_kalman->P[0][1] - p_kalman->P[1][0] + p_kalman->Q_angle);
    p_kalman->P[0][1] -= dt * p_kalman->P[1][1];
    p_kalman->P[1][0] -= dt * p_kalman->P[1][1];
    p_kalman->P[1][1] += p_kalman->Q_bias * dt;

    // Discrete Kalman filter measurement update equations - Measurement Update ("Correct")
    // Calculate Kalman gain - Compute the Kalman gain
    
    float S = p_kalman->P[0][0] + p_kalman->R_measure; // Estimate error
    
    float K[2]; // Kalman gain - This is a 2x1 vector
    K[0] = p_kalman->P[0][0] / S;
    K[1] = p_kalman->P[1][0] / S;

    // Calculate angle and bias - Update estimate with measurement zk (newAngle)
    
    float y = newAngle - p_kalman->angle; // Angle difference
    
    p_kalman->angle += K[0] * y;
    p_kalman->bias += K[1] * y;

    // Calculate estimation error covariance - Update the error covariance
    
    float P00_temp = p_kalman->P[0][0];
    float P01_temp = p_kalman->P[0][1];

    p_kalman->P[0][0] -= K[0] * P00_temp;
    p_kalman->P[0][1] -= K[0] * P01_temp;
    p_kalman->P[1][0] -= K[1] * P00_temp;
    p_kalman->P[1][1] -= K[1] * P01_temp;

    return p_kalman->angle;
}



void setAngle(kalman * p_kalman , float angle) {
p_kalman->angle=angle;
}