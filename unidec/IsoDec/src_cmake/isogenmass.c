//
// Created by Michael Marty on 10/17/2024.
//
#include "isogenmass.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "isogenmass_model_64.h"


void mass_to_vector(const float mass, float *vector) {
    // const float x0 = mass / 10000000.0f;
    const float x1 = mass / 1000000.0f;
    const float x2 = mass / 100000.0f;
    const float x3 = mass / 10000.0f;
    const float x4 = mass / 1000.0f;
    const float x5 = mass / 100.0f;
    vector[0] = fmaxf(0.0f, fminf(1.0f, x1));
    vector[1] = fmaxf(0.0f, fminf(1.0f, x2));
    vector[2] = fmaxf(0.0f, fminf(1.0f, x3));
    vector[3] = fmaxf(0.0f, fminf(1.0f, x4));
    vector[4] = fmaxf(0.0f, fminf(1.0f, x5));
    // vector[5] = fmaxf(0.0f, fminf(1.0f, x5));
    // printf("Vector: %f %f %f %f %f\n", vector[0], vector[1], vector[2], vector[3], vector[4]);
}

// Set up the weights
struct IsoGenWeights SetupWeights(const int veclen, const int outlen) {
    struct IsoGenWeights weights;
    weights.vl1 = veclen;
    weights.vl2 = outlen;
    weights.vl3 = outlen;
    weights.vl4 = outlen;
    weights.ml1 = weights.vl1 * weights.vl2;
    weights.ml2 = weights.vl2 * weights.vl3;
    weights.ml3 = weights.vl3 * weights.vl4;
    weights.tot = weights.ml1 + weights.vl2 + weights.ml2 + weights.vl3 + weights.ml3 + weights.vl4;
    // Set up the weights and biases
    weights.w1 = (float *) calloc(weights.ml1, sizeof(float));
    weights.b1 = (float *) calloc(weights.vl2, sizeof(float));
    weights.w2 = (float *) calloc(weights.ml2, sizeof(float));
    weights.b2 = (float *) calloc(weights.vl3, sizeof(float));
    weights.w3 = (float *) calloc(weights.ml3, sizeof(float));
    weights.b3 = (float *) calloc(weights.vl4, sizeof(float));

    // Check for null
    if (weights.w1 == NULL || weights.b1 == NULL || weights.w2 == NULL || weights.b2 == NULL || weights.w3 == NULL || weights.b3 == NULL) {
        printf("Error allocating memory for weights\n");
        exit(100);
    }
    return weights;
}

// Free Weights
void FreeWeights(const struct IsoGenWeights weights) {
    free(weights.w1);
    free(weights.b1);
    free(weights.w2);
    free(weights.b2);
    free(weights.w3);
    free(weights.b3);
}

// Load the weights
struct IsoGenWeights LoadWeights(const struct IsoGenWeights weights, const unsigned char *model_weights) {
    // Allocate memory for the default weights
    float *default_weights = calloc(weights.tot, sizeof(float));
    if (default_weights == NULL) {
        printf("Error allocating memory for default weights\n");
        exit(100);
    }
    // Copy the default weights
    memcpy(default_weights, model_weights, weights.tot * sizeof(float));

    // Copy the weights into the object
    memcpy(weights.w1, default_weights, weights.ml1 * sizeof(float));
    memcpy(weights.b1, default_weights + weights.ml1, weights.vl2 * sizeof(float));
    memcpy(weights.w2, default_weights + weights.ml1 + weights.vl2, weights.ml2 * sizeof(float));
    memcpy(weights.b2, default_weights + weights.ml1 + weights.vl2 + weights.ml2, weights.vl3 * sizeof(float));
    memcpy(weights.w3, default_weights + weights.ml1 + weights.vl2 + weights.ml2 + weights.vl3, weights.ml3 * sizeof(float));
    memcpy(weights.b3, default_weights + weights.ml1 + weights.vl2 + weights.ml2 + weights.vl3 + weights.ml3, weights.vl4 * sizeof(float));
    // Free Memory
    free(default_weights);
    return weights;
}

// Softsign Function
float softsign(const float x) {return x / (1.0f + fabsf(x));}

// Sigmoid Function
float sigmoid(const float x) {return 1.0f / (1.0f + expf(-x));}

void neural_net(const float *vector, float *isodist, const struct IsoGenWeights weights) {
    //Allocate the memory for the intermediate values
    float *l1 = (float *) calloc(weights.vl2, sizeof(float));
    float *l2 = (float *) calloc(weights.vl3, sizeof(float));

    // Calculate the First Layer
    for (int i = 0; i < weights.vl2; i++) {
        l1[i] = weights.b1[i];
        for (int j = 0; j < weights.vl1; j++) {
            l1[i] += vector[j] * weights.w1[i * weights.vl1 + j];
        }
        l1[i] = softsign(l1[i]);
    }

    // Calculate the second layer
    for (int i = 0; i < weights.vl3; i++) {
        l2[i] = weights.b2[i];
        for (int j = 0; j < weights.vl2; j++) {
            l2[i] += l1[j] * weights.w2[i * weights.vl2 + j];
        }
        l2[i] = softsign(l2[i]);
    }

    // Calculate the third layer
    for (int i = 0; i < weights.vl4; i++) {
        isodist[i] = weights.b3[i];
        for (int j = 0; j < weights.vl3; j++) {
            isodist[i] += l2[j] * weights.w3[i * weights.vl3 + j];
        }
        isodist[i] = sigmoid(isodist[i]);
    }

    // Free Memory
    free(l1);
    free(l2);
}

void isogenmass(const float mass, float *isodist) {
    // Create Encoding Vector
    float vector[5];
    mass_to_vector(mass, vector);
    // Setup and Load Weights
    struct IsoGenWeights weights = SetupWeights(5,64);
    weights = LoadWeights(weights, isogenmass_model_64_bin);
    // Run the Neural Net
    neural_net(vector, isodist, weights);
    // Free Memory
    FreeWeights(weights);
}

float isogenmass_fancy(const float mass, float *isodist, const int isolen, const int offset) {
    // Run Isogen mass
    isogenmass(mass, isodist);
    // Return the max value
    float max = 0.0f;
    for (int i = isolen-offset-1; i >= 0 ; i--) {
        isodist[i+offset] = isodist[i];
        if (isodist[i] > max) {
            max = isodist[i];
        }
        if (i<=offset){isodist[i]=0.0f;}
    }
    return max;
}

// Isotope Parameters
float isoparams[10] = {
    1.00840852e+00f, 1.25318718e-03f, 2.37226341e+00f, 8.19178000e-04f, -4.37741951e-01f,
    6.64992972e-04f, 9.94230511e-01f, 4.64975237e-01f, 1.00529041e-02f, 5.81240792e-01f
};

// Function used in isotope distribution calculation
float isotopemid(const float mass, const float *isoparams) {
    const float a = isoparams[4];
    const float b = isoparams[5];
    const float c = isoparams[6];
    return a + b * powf(mass, c);
}

// Function used in isotope distribution calculation
float isotopesig(const float mass, const float *isoparams) {
    const float a = isoparams[7];
    const float b = isoparams[8];
    const float c = isoparams[9];
    return a + b * powf(mass, c);
}

// Function used in isotope distribution calculation
float isotopealpha(const float mass, const float *isoparams) {
    const float a = isoparams[0];
    const float b = isoparams[1];
    return a * expf(-mass * b);
}

// Function used in isotope distribution calculation
float isotopebeta(const float mass, const float *isoparams) {
    const float a = isoparams[2];
    const float b = isoparams[3];
    return a * expf(-mass * b);
}

float isomike(const float mass, float * isodist, const int isolen, const int offset) {
    const float mid = isotopemid(mass, isoparams);
    const float sig = isotopesig(mass, isoparams);
    if (sig == 0) {
        printf("Error: Sigma Isotope Parameter is 0");
        exit(102);
    }
    const float alpha = isotopealpha(mass, isoparams);
    const float amp = 1.0f - alpha;
    const float beta = isotopebeta(mass, isoparams);

    float max = 0.0f;
    for (int k = 0; k < isolen-offset; k++) {
        const float kfloat = (float) k;
        const float e = alpha * expf(-kfloat * beta);
        const float g = amp / (sig * 2.50662827f) * expf(-powf(kfloat - mid, 2) / (2 * powf(sig, 2)));
        const float val = e + g;
        if(val> max) {max = val;}
        isodist[k+offset] = val;
    }
    return max;
}