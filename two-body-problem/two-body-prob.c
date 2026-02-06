/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2026, Abhijnan Saraswat Gogoi
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
usage: ./two-body-prob.o [--help] [--force] --m1 --m2 --r --v --vt --k --n --name [--step_size] [--steps] [--read_steps]

Options:
 --help                         Show this help message and exit
 --force        Optional        Force overwrite of existing simulation data
 --m1           Required        Mass of object 1
 --m2           Required        Mass of object 2
 --r            Required        Initial distance between the objects
 --v            Required        Magnitude of initial relative velocity of object 2 with respect to object 1
 --vt           Required        Angle of orientation of the initial relative velocity of object 2 with respect to object 1 (in radians)
 --k            Required        Force law proportionality constant
 --n            Required        Force law power (an integer)
 --name         Required        Name of the simulation
 --step_size    Optional        Step size (in position space) for each iteration
 --steps        Optional        Total number of steps to be taken
 --read_steps   Optional        Number of iterations after which a step is sampled
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include <time.h>

#define PI 3.1415926535

double U(double r); // Potential energy function.
void iterate_next(double *t, double *r, double *th, double *pr, double L, double m, double step);
void argparse(int argc, char **argv, short *force, double *m1, double *m2, double *r, double *v, double *vt, char **name, double *step_size, long *steps, long *read_steps);    // To parse CLI arguments.

int n;
double k;

int main(int argc, char **argv) {
    time_t time_ini, time_fin;

    FILE *fp;
    struct stat statbuf;

    char *name = (char *)malloc(100*sizeof(char));
    char file_str[210];
    short force = 0;

    double step_size;
    long steps;
    long read_steps;

    double m1, m2;                  // Masses.
    double m;                       // Reduced mass.
    double r, th, pr, L, E, t;      // Generalised coordinates in phase space (..., in radians, ..., ...) and total energy and time.
    double v, v_ang;                // Initial velocity parameters (..., in radians).

    th = 0; // Initial theta is 0.
    t = 0;  // Initial time is 0.

    // CLI argument parsing
    argparse(argc, argv, &force, &m1, &m2, &r, &v, &v_ang, &name, &step_size, &steps, &read_steps);

    // Check/create output directories:
    if(stat(name, &statbuf)!=0) {
        if(mkdir(name, 0755)!=0) {
            printf("Error: Could not create %s directory. Please create it manually.\n", name);
            printf("Exiting code\n");
            exit(EXIT_FAILURE);
        }
    }
    else if(!force) {
        printf("Simulation with the same name already exists.\nUse --force to overwrite existing data.\n");
        printf("Exiting code\n");
        exit(EXIT_FAILURE);
    }

    // Obtained derived parameters:
    m = m1*m2/(m1+m2);
    pr = m*v*sin(v_ang);
    L = r*m*v*cos(v_ang);
    E = U(r) + 0.5*m*pow(v,2);

    // Display on stdout:
    printf("2-body central force problem simulation\n");
    printf("Simulation name: %s\n\n", name);

    printf("Potential function:\tU(r) = k*r^n\n");
    printf("k\t= %g\nn\t= %d\n\n", k, n);

    printf("Input parameters:\n");
    printf("m1\t= %g\nm2\t= %g\nr\t= %g\nv\t= %g\nv_ang\t= %g\n\n", m1, m2, r, v, v_ang);

    printf("Derived parameters:\n");
    printf("m\t= %g\t(reduced mass)\npr\t= %g\t(canonical momentum for r)\nL\t= %g\t(canonical momentum for theta)(conserved)\nE\t= %g\t(total energy)(conserved)\n\n", m, pr, L, E);

    printf("Number of steps: %ld\nStep size: %g\nSampling every %ld steps\n", steps, step_size, read_steps);

    // Print parameter data into params file:
    sprintf(file_str, "%s/%s_params.txt", name, name);
    fp = fopen(file_str, "w");
    
    fprintf(fp, "# 2-body central force problem simulation\n");
    fprintf(fp, "# Simulation name: %s\n\n", name);
    fprintf(fp, "# Potential function:\tU(r) = k*r^n\n");
    fprintf(fp, "# k\t= %g\n", k);
    fprintf(fp, "# n\t= %d\n\n", n);

    fprintf(fp, "# Input parameters:\n");
    fprintf(fp, "# m1\t= %g\n", m1);
    fprintf(fp, "# m2\t= %g\n", m2);
    fprintf(fp, "# r\t\t= %g\n", r);
    fprintf(fp, "# v\t\t= %g\n", v);
    fprintf(fp, "# v_ang\t= %g\n\n", v_ang);

    fprintf(fp, "# Derived parameters:\n");
    fprintf(fp, "# m\t\t= %g\t(reduced mass)\n", m);
    fprintf(fp, "# pr\t= %g\t(canonical momentum for r)\n", pr);
    fprintf(fp, "# L\t\t= %g\t(canonical momentum for theta)(conserved)\n", L);
    fprintf(fp, "# E\t\t= %g\t(total energy)(conserved)\n\n", E);

    fprintf(fp, "# Number of steps: %ld\n", steps);
    fprintf(fp, "# Step size: %g\n", step_size);
    fprintf(fp, "# Sampling every %ld steps\n", read_steps);

    fclose(fp);

    // Simulation related:
    sprintf(file_str, "%s/%s_data.txt", name, name);
    fp = fopen(file_str, "w");

    fprintf(fp, "# 2-body central force problem simulation\n");
    fprintf(fp, "# Simulation name: %s\n", name);
    fprintf(fp, "# Number of steps: %ld\n", steps);
    fprintf(fp, "# Step size: %g\n", step_size);
    fprintf(fp, "# Reading every %ld steps\n", read_steps);
    fprintf(fp, "# \t\tstep\t\tt\t\t\tr\t\ttheta\t\tpr\n");
    fprintf(fp, "\t%10d\t%g\t%g\t%g\t%g\n", 0, t, r, th, pr);

    time_ini = time(NULL);
    printf("Starting simulation:\n");
    printf("[%s]\tsteps taken=%10d\t\tsteps recorded=%10d\tt=%g\t\t\tr=%g\tth=%g\tpr=%g\n", argv[0], 0, 0, t, r, th, pr);
    for(long int i=1; i<=steps; i++) {
        iterate_next(&t, &r, &th, &pr, L, m, step_size);
        if(!(i%read_steps)) {
            fprintf(fp, "\t%ld\t%g\t%g\t%g\t%g\n", i, t, r, th, pr);
        }
        if(!(i%1000000)) {
            printf("[%s]\tsteps taken=%10ld\t\tsteps read=%10d\tt=%g\tr=%g\tth=%g\tpr=%g\n", argv[0], i, (int)floor(i/100), t, r, th, pr);
        }
    }
    fclose(fp);
    time_fin = time(NULL);
    printf("Complete after %g seconds!\n", difftime(time_fin, time_ini));
    return 0;
}

void iterate_next(double *t, double *r, double *th, double *pr, double L, double m, double step) {
    double t0 = *t;
    double r0 = *r;
    double th0 = *th;
    double pr0 = *pr;

    double dUdr = k*n*pow(r0,n-1);
    double A, B;    // Some expressions for easy coding

    double dt, dr, dth, dpr;        // Steps
    double Dt, Dr, Dth, Dpr;        // First order derivatives
    double D2t, D2r, D2th, D2pr;    // Second order derivatives

    if(abs(pr0)<1e-10 && abs(L)<1e-10) {    // if pr0==0 && L==0
        dth = 0;
        
        if(-1*dUdr>0) dr = 1;
        else if(-1*dUdr<0) dr = -1;
        else dr = 0;

        if(dUdr!=0) dt = sqrt((2*m*step)/abs(dUdr));
        else dt = 1;
        dpr = -1*dUdr*dt;
    }
    else {
        A = pow(r0*pr0,2)+pow(L,2);
        B = (pow(L,2)/(m*pow(r0,3)))-dUdr;

        Dt = (m*r0)/sqrt(A);
        Dr = (r0*pr0)/sqrt(A);
        Dth = L/(r0*sqrt(A));
        Dpr = ((m*r0)/sqrt(A))*B;

        D2t = m*r0*((pr0/A) - ((pow(r0,2)*pow(pr0,3))/pow(A,2)) - ((pr0*pow(L,2))/pow(A,2)) + ((m*pow(r0,3)*pr0)/pow(A,2))*dUdr);
        D2r = r0*pr0*((pr0/A) + ((m*r0)/(pr0*A))*B - ((m*pow(r0,3)*pr0)/pow(A,2))*B - ((pow(r0,2)*(pow(pr0,3)))/pow(A,2)));
        D2th = ((-1*L)/(r0))*((pr0/A) + ((pow(r0,2)*pow(pr0,3))/pow(A,2)) + ((m*pow(r0,3)*pr0)/pow(A,2))*B);
        D2pr = m*r0*B*((pr0/A) - ((pow(r0,2)*pow(pr0,3))/pow(A,2)) - ((m*pow(r0,3)*pr0)/pow(A,2))*B - ((3*pow(L,2)*pr0)/(m*pow(r0,3)*A)) - (n-1)*(pr0/A)*dUdr);

        dt = step*Dt + 0.5*pow(step,2)*D2t;
        dr = step*Dr + 0.5*pow(step,2)*D2r;
        dth = step*Dth + 0.5*pow(step,2)*D2th;
        dpr = step*Dpr + 0.5*pow(step,2)*D2pr;
    }

    *t = t0 + dt;
    *r = r0 + dr;
    *th = th0 + dth;
    *pr = pr0 + dpr;
}

double U(double r) {
    return(k * pow(r,n));
}

void argparse(int argc, char **argv, short *force, double *m1, double *m2, double *r, double *v, double *v_ang, char **name, double *step_size, long *steps, long *read_steps) {
    char *end_ptr;
    char *token, arg_opt[100];
    char *help = (char *) malloc(1000*sizeof(char));
    sprintf(help, "usage: %s [--help] [--force] --m1 --m2 --r --v --vt --k --n --name [--step_size] [--steps] [--read_steps]\n\nOptions:\n --help\t\t\t\tShow this help message and exit\n --force\tOptional\tForce overwrite of existing simulation data\n --m1\t\tRequired\tMass of object 1\n --m2\t\tRequired\tMass of object 2\n --r\t\tRequired\tInitial distance between the objects\n --v\t\tRequired\tMagnitude of initial relative velocity of object 2 with respect to object 1\n --vt\t\tRequired\tAngle of orientation of the initial relative velocity of object 2 with respect to object 1 (in radians)\n --k\t\tRequired\tForce law proportionality constant\n --n\t\tRequired\tForce law power (an integer)\n --name\t\tRequired\tName of the simulation\n --step_size\tOptional\tStep size (in position space) for each iteration\n --steps\tOptional\tTotal number of steps to be taken\n --read_steps\tOptional\tNumber of iterations after which a step is read", argv[0]);

    struct cli_checklist {
        short m1;
        short m2;
        short r;
        short v;
        short vt;
        short k;
        short n;
        short name;
        short step_size;
        short steps;
        short read_steps;
    } args_check;
    args_check.m1 = 0;
    args_check.m2 = 0;
    args_check.r = 0;
    args_check.v = 0;
    args_check.vt = 0;
    args_check.k = 0;
    args_check.n = 0;
    args_check.name = 0;
    args_check.step_size = 0;
    args_check.steps = 0;
    args_check.read_steps = 0;

    *force = 0; // Set default first.

    // Obtain input parameters from CLI args:
    for(int i=1; i<argc; i++) {
        strcpy(arg_opt, argv[i]);
        token = strtok(arg_opt, "=");
        token = strtok(NULL, "=");

        if(!strcmp(arg_opt, "--help")) {
            printf("Ignoring value %s for %s\n", token, arg_opt);
            printf("%s\n", help);
            exit(EXIT_SUCCESS);
        }
        else if(!strcmp(arg_opt, "--force")) {
            printf("Ignoring value %s for %s\n", token, arg_opt);
            *force = 1;
        }
        else if(!strcmp(arg_opt, "--m1")) {
            *m1 = strtod(token, &end_ptr);
            args_check.m1 = 1;
        }
        else if(!strcmp(arg_opt, "--m2")) {
            *m2 = strtod(token, &end_ptr);
            args_check.m2 = 1;
        }
        else if(!strcmp(arg_opt, "--r")) {
            *r = strtod(token, &end_ptr);
            args_check.r = 1;
        }
        else if(!strcmp(arg_opt, "--v")) {
            *v = strtod(token, &end_ptr);
            args_check.v = 1;
        }
        else if(!strcmp(arg_opt, "--vt")) {
            *v_ang = strtod(token, &end_ptr);
            args_check.vt = 1;
        }
        else if(!strcmp(arg_opt, "--k")) {
            *(&k) = strtod(token, &end_ptr);
            args_check.k = 1;
        }
        else if(!strcmp(arg_opt, "--n")) {
            *(&n) = atoi(token);
            args_check.n = 1;
        }
        else if(!strcmp(arg_opt, "--name")) {
            strcpy(*name, token);
            args_check.name = 1;
        }
        else if(!strcmp(arg_opt, "--step_size")) {
            *step_size = strtod(token, &end_ptr);
            args_check.step_size = 1;
        }
        else if(!strcmp(arg_opt, "--steps")) {
            *steps = atol(token);
            args_check.steps = 1;
        }
        else if(!strcmp(arg_opt, "--read_steps")) {
            *read_steps = atol(token);
            args_check.read_steps = 1;
        }
        else {
            printf("Error: invalid option %s\n", arg_opt);
            printf("%s\n", help);
            printf("Exiting code\n");
            exit(EXIT_FAILURE);
        }
    }

    // Check for required values:
    if(!(args_check.m1 && args_check.m2 && args_check.r && args_check.v && args_check.vt && args_check.k && args_check.n && args_check.name)) {
        printf("Error: some required options are missing.\n");
        printf("%s\n", help);
        printf("Exiting code\n");
        exit(EXIT_FAILURE);
    }

    // Set default values for optional parameters if unfilled:
    if(!args_check.step_size) {
        *step_size = 1e-6;
    }
    if(!args_check.steps) {
        *steps = 10000000;
    }
    if(!args_check.read_steps) {
        *read_steps = 40000;
    }
    return;
}
