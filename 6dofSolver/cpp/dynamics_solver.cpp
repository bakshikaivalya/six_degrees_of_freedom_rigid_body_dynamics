// Author: Kaivalya Bakshi (kaivalya.bakshi@zeusnumerix.com)

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>

using namespace std;

#include "dynamics_solver.h"
#include "define_variables.cpp"

// _Function to calculate the transpose of a 3 by 3 matrix_
void DynamicsSolver::transposeMatrix(double a[][3], double b[][3]) {
	
	for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            b[i][j] = a[j][i];
        }
    }
	
}

// _Function to calculate the product of two 3 by 3 matrices_
void DynamicsSolver::matrixMultiply3_by_3(double c[][3], double d[][3], double e[][3]) {
	
	int i, j, k;
    for(i=0; i<3; i++) {
        for(j=0; j<3; j++) {
            e[i][j] = 0.0;
            for(k=0; k<3; k++) {
                e[i][j] =   e[i][j] + c[i][k] * d[k][j];
            }
        }
    }

}

// _Function to calculate the inverse of a 3 by 3 matrix_
void DynamicsSolver::matrixInverse (double A[][3] ,double B[][3]) {
	
	double d1,d2,d3,det,invdet;

    d1=  A[0][0]*(A[1][1]*A[2][2]-A[2][1]*A[1][2]);
    d2 = A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]);
    d3 = A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);

    det = d1-d2+d3;

    invdet = 1/det;

    B[0][0] =  (A[1][1]*A[2][2]-A[2][1]*A[1][2])*invdet;
    B[0][1] = -(A[0][1]*A[2][2]-A[0][2]*A[2][1])*invdet;
    B[0][2] =  (A[0][1]*A[1][2]-A[0][2]*A[1][1])*invdet;
    B[1][0] = -(A[1][0]*A[2][2]-A[1][2]*A[2][0])*invdet;
    B[1][1] =  (A[0][0]*A[2][2]-A[0][2]*A[2][0])*invdet;
    B[1][2] = -(A[0][0]*A[1][2]-A[1][0]*A[0][2])*invdet;
    B[2][0] =  (A[1][0]*A[2][1]-A[2][0]*A[1][1])*invdet;
    B[2][1] = -(A[0][0]*A[2][1]-A[2][0]*A[0][1])*invdet;
    B[2][2] =  (A[0][0]*A[1][1]-A[1][0]*A[0][1])*invdet;
	
}

// _The 6 DoF solver which is used to solve the Dynamics Equations for x, y and theta    
//  which are three degrees of freedom in all, for an aerospace vehicle with 
//  no roll, yaw and forces normal to its plane of motion_
double DynamicsSolver::DynamicsSolver_ThreeDoF () {

    // _Loading the Aerodynamic, Thrust and Structural property tables in *.dat format
    //  as arrays which can be utilized later_
    double EAT_Output_Table[1000][7], Thrust_Table[1000][2], Structural_Output[2][5];
    int lenEAT_Output_Table, lenThrust_Table;
    double tbr;
    double Range, Modulus_v_E, Altitude;
    
    // Loading EAT_Output_Table.dat as an array of the data
    ifstream file1 ("EAT_Output_Table.dat");
    
    int row1 = 0; 
    
    while (!file1.eof()) {
		file1 >> EAT_Output_Table[row1][0] >> EAT_Output_Table[row1][1] >> EAT_Output_Table[row1][2] >>
	     EAT_Output_Table[row1][3] >> EAT_Output_Table[row1][4] >> EAT_Output_Table[row1][5] >>
	     EAT_Output_Table[row1][6];
	    row1 ++; 	
	}
	
	lenEAT_Output_Table = row1 - 1;

	// Loading Thrust_Table.dat as an array of the data
	ifstream file2 ("Thrust_Table.dat");
	
	int row2 = 0;
	
	while (!file2.eof()) {
		file2 >> Thrust_Table[row2][0] >> Thrust_Table[row2][1];
	    row2 ++;
	}
	lenThrust_Table = row2 - 1;
	tbr = Thrust_Table[lenThrust_Table - 1][0];
	
	// Loading Structural_Output.dat as an array of the data
	ifstream file3 ("Structural_Output.dat");
	
	int row3 = 0;
	
	while (!file3.eof()) {
		file3 >> Structural_Output[row3][0] >> Structural_Output[row3][1] >> Structural_Output[row3][2] >>
	     Structural_Output[row3][3] >> Structural_Output[row3][4];
	    row3 ++; 	
	}
		
	// Initialization of state and other variables
	int IntervalsofAlpha;
	double t, x, y, z, u, v, w, phi, theta, psi, p, q, r;
	double LaunchHeadingAngle, LaunchFlightPathAngle;
	double h;
	
	x = 0;
	y = 0;
	z = 0;
	u = - V;
	v = 0;
	w = 0;
	phi = 0;
	LaunchFlightPathAngle = LaunchAngle * pi/180; theta = LaunchFlightPathAngle;
	LaunchHeadingAngle = 0; psi = pi + LaunchHeadingAngle;
	p = 0;
	q = 0;
	r = 0;
	
	IntervalsofAlpha = 1.0/(IntervalLengthAlpha-0.00000001) * UpperBoundAlpha;
	
	// Vehicle geometry parameters
	double A_cs;
	double massi, Ixxi, Ixyi, Ixzi, Iyxi, Iyyi, Iyzi, Izxi, Izyi, Izzi, xcmi;
	double massf, Ixxf, Ixyf, Ixzf, Iyxf, Iyyf, Iyzf, Izxf, Izyf, Izzf, xcmf;
	
	A_cs = (3.14159 * d * d)/4 ;
	massi = Structural_Output[0][0];  massf = Structural_Output[1][0]; 
    Ixxi =  Structural_Output[0][2]/1000000;  Ixxf = Structural_Output[1][2]/1000000;
    Ixyi =  0;  Ixyf = 0;
    Ixzi =  0;  Ixzf = 0;
    Iyxi =  0;  Iyxf = 0;
    Iyyi =  Structural_Output[0][3]/1000000;  Iyyf = Structural_Output[1][3]/1000000;
    Iyzi =  0;  Iyzf = 0;
    Izxi =  0;  Izxf = 0;
    Izyi =  0;  Izyf = 0;
    Izzi =  Structural_Output[0][4]/1000000;  Izzf = Structural_Output[1][4]/1000000;
    xcmi = Structural_Output[0][1]/1000;  xcmf = Structural_Output[1][1]/1000;
    
    // Initialization of time variable
    t = 0; // In seconds
     
    // Time step for Runge Kutta Method used for solution of ODEs
    h = TimeIncrement;
    
    // Initialization of Kutta Matrix
    double k [5][12];
    
    for (int i=0; i<12; i++) {
        k[0][i] = 0;	
	}
    
    // _Solution of ODEs by Runge Kutta 4 Method_
    double xr, yr, zr, ur, vr, wr, phir, thetar, psir, pr, qr, rr;
    double timestep;
    
    // Creating an output data file Dynamics_Trajectory_Output.dat in which the generated
    //  Dynamics Trajectory data can be recorded
    ofstream fwrite;
    fwrite.open ( "Dynamics_Trajectory_Output.dat", ios :: trunc );
            
    while (z >= 0) {	
		
		for (int j=1; j<5; j++) {
			
			if (j==4) {
            timestep=h;
            }
            else {
            timestep=h/2.0;
            }
        
        ur = u + timestep * k[j-1][0];
        vr = v + timestep * k[j-1][1];
        wr = w + timestep * k[j-1][2];
        pr = p + timestep * k[j-1][3];
        qr = q + timestep * k[j-1][4];
        rr = r + timestep * k[j-1][5];
        phir = phi + timestep * k[j-1][6];
        thetar = theta + timestep * k[j-1][7];
        psir = psi + timestep * k[j-1][8];
        xr = x + timestep * k[j-1][9];
        yr = y + timestep * k[j-1][10];
        zr = z + timestep * k[j-1][11];
        
        // Calculation of Matrix of Inertia Tensor by Linear Interpolation
        //  of comprising elements with time and Calculation of Thrust
        //  by Linear Interpolation over time variable
        double Inertia_Matrix[3][3];
        double mass, Ixx, Ixy, Ixz, Iyx, Iyy, Iyz, Izx, Izy, Izz, xcm, ycm, zcm;
        double Thrust;
        
        if (t<=tbr) {
            mass = massi - (massi - massf) * (t/tbr);
            Ixx = Ixxi - (Ixxi-Ixxf) * (t/tbr);
            Ixy = Ixyi - (Ixyi-Ixyf) * (t/tbr);
            Ixz = Ixzi - (Ixzi-Ixzf) * (t/tbr);
            Iyx = Iyxi - (Iyxi-Iyxf) * (t/tbr);
            Iyy = Iyyi - (Iyyi-Iyyf) * (t/tbr);
            Iyz = Iyzi - (Iyzi-Iyzf) * (t/tbr);
            Izx = Izxi - (Izxi-Izxf) * (t/tbr);
            Izy = Izyi - (Izyi-Izyf) * (t/tbr);
            Izz = Izzi - (Izzi-Izzf) * (t/tbr);
            xcm = xcmi - (xcmi-xcmf) * (t/tbr);
            ycm = 0;
            zcm = 0;
		
            int thrustvar;
            
            for (int c=0; c<lenThrust_Table; c++) {
                double C;
                C = Thrust_Table[c][0] - t;
                if (C>0) {
                    thrustvar = c - 1;
                    break;
                }
			}
			
            Thrust = Thrust_Table[thrustvar][1] + (Thrust_Table[thrustvar+1][1] -
             Thrust_Table[thrustvar][1])/(Thrust_Table[thrustvar+1][0] -
             Thrust_Table[thrustvar][0]) * (t - Thrust_Table[thrustvar][0]);
		}
		else {
		    mass = massf;
            Ixx = Ixxf;
            Ixy = Ixyf;
            Ixz = Ixzf;
            Iyx = Iyxf;
            Iyy = Iyyf;
            Iyz = Iyzf;
            Izx = Izxf;
            Izy = Izyf;
            Izz = Izzf;
            xcm = xcmf;
            ycm = 0;
            zcm = 0;
            
            Thrust = 0;
		}
		
		Thrust = 9.81 * Thrust;
		
		Inertia_Matrix[0][0] = Ixx;
        Inertia_Matrix[0][1] = Ixy;
        Inertia_Matrix[0][2] = Ixz;
        Inertia_Matrix[1][0] = Iyx;
        Inertia_Matrix[1][1] = Iyy;
        Inertia_Matrix[1][2] = Iyz;
        Inertia_Matrix[2][0] = Izx;
        Inertia_Matrix[2][1] = Izy;
        Inertia_Matrix[2][2] = Izz;
        
        // Calculation of Coordinates Transformation Matrix
		double T_BE[3][3], T_EB[3][3];
		
		T_BE[0][0] = cos(psir) * cos(thetar);
        T_BE[0][1] = (sin(psir)) * cos(thetar);
        T_BE[0][2] = - sin(thetar);
        T_BE[1][0] = (cos(psir)) * (sin(thetar)) * sin(phir) - (sin(psir)) * cos(phir);
        T_BE[1][1] = (sin(psir)) * (sin(thetar)) * sin(phir) + (cos(psir)) * cos(phir);
        T_BE[1][2] = (cos(thetar)) * sin(phir);
        T_BE[2][0] = (cos(psir)) * (sin(thetar)) * cos(phir) + (sin(psir)) * sin(phir);
        T_BE[2][1] = (sin(psir)) * (sin(thetar)) * cos(phir) - (cos(psir)) * sin(phir);
        T_BE[2][2] = (cos(thetar)) * cos(phir);
        transposeMatrix(T_BE, T_EB);
        
        // Calculation of rho, Ma, alpha, qbar
        double rho, qbar, vs, Ma, alpha_radians, alpha;
        
        if (zr<11000) {
            rho = 1.225 * pow((1.0 - 0.000022557695 * zr), (4.25587));
		}
        else {
            rho = 0.36392 * exp(- 1.576584e-04 * (zr - 11000));
		}
            
        vs = 340.249 * sqrt((288.15 - 0.0065 * zr)/288.15);
        Ma = sqrt (ur * ur + vr * vr + wr * wr)/vs;
        
        alpha_radians = atan (wr/ur);
        alpha = alpha_radians * 180/pi;
        
        qbar = 0.5 * rho * (ur * ur + vr * vr + wr * wr);
        
        // Calculation of Cx, Cz, Cm
        int n, m;
        double A, B;
        
        for (int a=0; a<lenEAT_Output_Table; a++) {
            A = EAT_Output_Table[a][1] - fabs(alpha); 
            if (A>0) {
                n = a - 1;
                break;
			}
		}
        
        for (int b=0; b<lenEAT_Output_Table; b++) {
            B = EAT_Output_Table[b][0] - Ma;
            if (B>0) {
                m = b - (IntervalsofAlpha + 1);
                break;
			}
		}
		
		if (Ma<IntervalLengthAlpha) {
			m = m + (IntervalsofAlpha + 1);
		}

		// Interpolation of Cx, Cz, Cm0, Cmq with y as independent variable at x = xi, x= xi+1
		double Cx, Cz, Cm0, Cnalpha, xcp, Cmq;
		double CX1, CX2, CZ1, CZ2, CM01, CM02, CNALPHA1, CNALPHA2, XCP1, XCP2;
		
		CX1 = EAT_Output_Table[m+n][6] + (EAT_Output_Table[m+n+1][6] - 
				EAT_Output_Table[m+n][6])/(EAT_Output_Table[n+1][1] - 
				EAT_Output_Table[n][1]) * (fabs(alpha) - EAT_Output_Table[n][1]);
				
		
        CX2 = EAT_Output_Table[m+(IntervalsofAlpha+1)+n][6] + 
				(EAT_Output_Table[m+(IntervalsofAlpha+1)+n+1][6] - 
				EAT_Output_Table[m+(IntervalsofAlpha+1)+n][6])/
				(EAT_Output_Table[n+1][1] - EAT_Output_Table[n][1]) * 
				(fabs(alpha) - EAT_Output_Table[n][1]);
		
        CZ1 = EAT_Output_Table[m+n][3] + (EAT_Output_Table[m+n+1][3] -
				EAT_Output_Table[m+n][3])/(EAT_Output_Table[n+1][1] -
				EAT_Output_Table[n][1]) * (fabs(alpha) - EAT_Output_Table[n][1]);

        CZ2 = EAT_Output_Table[m+(IntervalsofAlpha+1)+n][3] +
				(EAT_Output_Table[m+(IntervalsofAlpha+1)+n+1][3] -
				EAT_Output_Table[m+(IntervalsofAlpha+1)+n][3])/
				(EAT_Output_Table[n+1][1] - EAT_Output_Table[n][1]) *
				(fabs(alpha) - EAT_Output_Table[n][1]);

        CM01 = EAT_Output_Table[m+n][5] + (EAT_Output_Table[m+n+1][5] - 
				EAT_Output_Table[m+n][5])/(EAT_Output_Table[n+1][1] - 
				EAT_Output_Table[n][1]) * (fabs(alpha) - EAT_Output_Table[n][1]);
    
        CM02 = EAT_Output_Table[m+(IntervalsofAlpha+1)+n][5] +
				(EAT_Output_Table[m+(IntervalsofAlpha+1)+n+1][5] -
				EAT_Output_Table[m+(IntervalsofAlpha+1)+n][5])/
				(EAT_Output_Table[n+1][1] - EAT_Output_Table[n][1]) *
				(fabs(alpha) - EAT_Output_Table[n][1]);
        
        CNALPHA1 = EAT_Output_Table[m+n][2] + (EAT_Output_Table[m+n+1][2] -
				EAT_Output_Table[m+n][2])/(EAT_Output_Table[n+1][1] -
				EAT_Output_Table[n][1]) * (fabs(alpha) - EAT_Output_Table[n][1]);
        
        CNALPHA2 = EAT_Output_Table[m+(IntervalsofAlpha+1)+n][2] +
				(EAT_Output_Table[m+(IntervalsofAlpha+1)+n+1][2] -
				EAT_Output_Table[m+(IntervalsofAlpha+1)+n][2])/
				(EAT_Output_Table[n+1][1] - EAT_Output_Table[n][1]) *
				(fabs(alpha) - EAT_Output_Table[n][1]);
        
        XCP1 = EAT_Output_Table[m+n][4] + (EAT_Output_Table[m+n+1][4] -
				EAT_Output_Table[m+n][4])/(EAT_Output_Table[n+1][1] -
				EAT_Output_Table[n][1]) * (fabs(alpha) - EAT_Output_Table[n][1]);
        
        XCP2 = EAT_Output_Table[m+(IntervalsofAlpha+1)+n][4] +
				(EAT_Output_Table[m+(IntervalsofAlpha+1)+n+1][4] -
				 EAT_Output_Table[m+(IntervalsofAlpha+1)+n][4])/
				 (EAT_Output_Table[n+1][1] - EAT_Output_Table[n][1]) *
				 (fabs(alpha) - EAT_Output_Table[n][1]);
        
        // Interpolation of Cx, Cz, Cm0 with x as independent variable at y = yi, y= yi+1
        Cx = CX1 + (CX2 - CX1)/(EAT_Output_Table[m+(IntervalsofAlpha+1)][0] - EAT_Output_Table[m][0]) *
				(Ma - EAT_Output_Table[m][0]);
		
        Cz = CZ1 + (CZ2 - CZ1)/(EAT_Output_Table[m+(IntervalsofAlpha+1)][0] - EAT_Output_Table[m][0]) *
				(Ma - EAT_Output_Table[m][0]);
		
        Cm0 = CM01 + (CM02 - CM01)/(EAT_Output_Table[m+(IntervalsofAlpha+1)][0] - EAT_Output_Table[m][0]) *
				(Ma - EAT_Output_Table[m][0]);
        
        Cnalpha = CNALPHA1 + (CNALPHA2 - CNALPHA1)/(EAT_Output_Table[m+(IntervalsofAlpha+1)][0] - EAT_Output_Table[m][0]) *
				(Ma - EAT_Output_Table[m][0]);
        
        xcp = (XCP1 + (XCP2 - XCP1)/(EAT_Output_Table[m+(IntervalsofAlpha+1)][0] - EAT_Output_Table[m][0]) *
				(Ma - EAT_Output_Table[m][0]))/1000;
        
        Cmq =   (- Cnalpha) * ((xcp - xcm) * (xcp - xcm))/(d * d);
        
        // Calculation of Aerodynamic Forces and Moments
        double Fx_aero, Fy_aero, Fz_aero, Mx_aero, My_aero, Mz_aero;
        
        Fx_aero = Cx * qbar * A_cs;
		
        Fy_aero = 0;

        if (alpha>0) {
            Fz_aero = Cz * qbar * A_cs;
		}
        else {
            Fz_aero = - Cz * qbar * A_cs;
		}
		
		Mx_aero = 0;
		
		if (alpha>0) {
            My_aero = (- Cm0 +  Cmq * ((qr * d)/(2 * sqrt (ur * ur + vr * vr + wr * wr)))) * qbar * A_cs * d;
		}
        else {
            My_aero = (Cm0 +  Cmq * ((qr * d)/(2 * sqrt (ur * ur + vr * vr + wr * wr)))) * qbar * A_cs * d;
		}
		
		Mz_aero = 0;
		
		// Calculation of Moment of Force due to Displacement of Aero Data
        //  Reference point from the CM at instant of time t
        double Mx_aeroref_cm, My_aeroref_cm, Mz_aeroref_cm;
        
        Mx_aeroref_cm = ycm * Fz_aero - zcm * Fy_aero;
    
        My_aeroref_cm = zcm * Fx_aero -  (- xcm) * Fz_aero;
    
        Mz_aeroref_cm = xcm * Fy_aero - ycm * Fx_aero;
        
        // Calculation of Forces and Moments due to Propulsion.
        //  Note that components of Moment due to Thrust Force are Zero since
        //  the Thrust Force is aligned with X axis of Body Coordinate System
        double Fx_thrust, Fy_thrust, Fz_thrust;
          
        Fx_thrust = - Thrust;
    
        Fy_thrust = 0;
    
        Fz_thrust = 0;
        
        double Fx_total, Fy_total, Fz_total;
        
        Fx_total = Fx_aero + Fx_thrust + - sin(thetar) * (mass * (- 9.81));
        Fy_total = Fy_aero + Fy_thrust + (cos(thetar)) * sin(phir) * (mass * (- 9.81));
        Fz_total = Fz_aero + Fz_thrust + (cos(thetar)) * cos(phir) * (mass * (- 9.81));
         
        double Mx_total, My_total, Mz_total;
        
        Mx_total = Mx_aero + Mx_aeroref_cm;
        My_total = My_aero + My_aeroref_cm;
        Mz_total = Mz_aero + Mz_aeroref_cm;
        
        // Calculation of k1, k2, k3, k4 for u, v, w
        double omega_BE_B_repeatedcolumnvectors[3][3], Omega_BE_B[3][3],
         vel_E_B_repeatedcolumnvectors[3][3], temp_RotVel[3][3];
        
        omega_BE_B_repeatedcolumnvectors[0][0] = pr;
        omega_BE_B_repeatedcolumnvectors[0][1] = pr;
        omega_BE_B_repeatedcolumnvectors[0][2] = pr;
        omega_BE_B_repeatedcolumnvectors[1][0] = qr;
        omega_BE_B_repeatedcolumnvectors[1][1] = qr;
        omega_BE_B_repeatedcolumnvectors[1][2] = qr;
        omega_BE_B_repeatedcolumnvectors[2][0] = rr;
        omega_BE_B_repeatedcolumnvectors[2][1] = rr;
        omega_BE_B_repeatedcolumnvectors[2][2] = rr;
        
        Omega_BE_B[0][0] = 0;
        Omega_BE_B[0][1] = - rr;
        Omega_BE_B[0][2] = qr;
        Omega_BE_B[1][0] = rr;
        Omega_BE_B[1][1] = 0;
        Omega_BE_B[1][2] = - pr;
        Omega_BE_B[2][0] = - qr;
        Omega_BE_B[2][1] = pr;
        Omega_BE_B[2][2] = 0;
        
        vel_E_B_repeatedcolumnvectors[0][0] = ur;
        vel_E_B_repeatedcolumnvectors[0][1] = ur;
        vel_E_B_repeatedcolumnvectors[0][2] = ur;
        vel_E_B_repeatedcolumnvectors[1][0] = vr;
        vel_E_B_repeatedcolumnvectors[1][1] = vr;
        vel_E_B_repeatedcolumnvectors[1][2] = vr;
        vel_E_B_repeatedcolumnvectors[2][0] = wr;
        vel_E_B_repeatedcolumnvectors[2][1] = wr;
        vel_E_B_repeatedcolumnvectors[2][2] = wr;
        
        
        
        matrixMultiply3_by_3 (Omega_BE_B, vel_E_B_repeatedcolumnvectors, temp_RotVel);
        
        k[j][0] = (Fx_total/mass) - temp_RotVel[0][0];
        k[j][1] = (Fy_total/mass) - temp_RotVel[1][0];
        k[j][2] = (Fz_total/mass) - temp_RotVel[2][0];
        
        // Calculation of k1, k2, k3, k4 for p, q, r         
        double inv_Inertia_Matrix[3][3], temp_Matrix1[3][3], temp_Matrix2[3][3],
         temp_Matrix_repeatedcolumnvectors[3][3], omegadot_E_B[3][3];
        
        matrixInverse (Inertia_Matrix, inv_Inertia_Matrix);
        
        matrixMultiply3_by_3 (Omega_BE_B, Inertia_Matrix, temp_Matrix1);
        
        matrixMultiply3_by_3 (temp_Matrix1, omega_BE_B_repeatedcolumnvectors, temp_Matrix2);
        
        temp_Matrix_repeatedcolumnvectors[0][0] =  Mx_total - temp_Matrix2[0][0];
        temp_Matrix_repeatedcolumnvectors[0][1] =  Mx_total - temp_Matrix2[0][0];
        temp_Matrix_repeatedcolumnvectors[0][2] =  Mx_total - temp_Matrix2[0][0];
        temp_Matrix_repeatedcolumnvectors[1][0] =  My_total - temp_Matrix2[1][0];
        temp_Matrix_repeatedcolumnvectors[1][1] =  My_total - temp_Matrix2[1][0];
        temp_Matrix_repeatedcolumnvectors[1][2] =  My_total - temp_Matrix2[1][0];
        temp_Matrix_repeatedcolumnvectors[2][0] =  Mz_total - temp_Matrix2[2][0];
        temp_Matrix_repeatedcolumnvectors[2][1] =  Mz_total - temp_Matrix2[2][0];
        temp_Matrix_repeatedcolumnvectors[2][2] =  Mz_total - temp_Matrix2[2][0];
        
        matrixMultiply3_by_3 (inv_Inertia_Matrix, temp_Matrix_repeatedcolumnvectors, omegadot_E_B);
        
        k[j][3] = omegadot_E_B[0][0];
        k[j][4] = omegadot_E_B[1][0];
        k[j][5] = omegadot_E_B[2][0];
        
        // Calculation of k1, k2, k3, k4 for psi, theta, phi
        k[j][6] = pr + qr * sin(phir) * tan(thetar) +
                        rr * cos(phir) * tan(thetar);
        k[j][7] = qr * cos(phir) - rr * sin(phir);
        k[j][8] = qr * sin(phir)/cos(thetar) +
                        rr * cos(phir)/cos(thetar);
                        
        // Calculation of k1, k2, k3, k4 for x, y, z
        k[j][9] = T_EB[0][0] * ur + T_EB[0][1] * vr + T_EB[0][2] * wr;
        k[j][10] = T_EB[1][0] * ur + T_EB[1][1] * vr + T_EB[1][2] * wr;
        k[j][11] = T_EB[2][0] * ur + T_EB[2][1] * vr + T_EB[2][2] * wr;
               
	}
	
	// Calculation of delta's of RK4 method for the 12 State Variables
	double RK[12];
	
	for (int incrementingvariable=0; incrementingvariable<12; incrementingvariable++) {
        RK[incrementingvariable] = (h/6.0) * (k[1][incrementingvariable] +
         2 * k[2][incrementingvariable] + 2 * k[3][incrementingvariable] +
          k[4][incrementingvariable]);
		}
	
	// Calculating Range and Velocity as functions of time
	
	Range = sqrt (x * x + y * y);
    Modulus_v_E = sqrt (u * u + v * v + w * w);
    Altitude = z;
	
	// Writing the generated Dynamics Trajectory data in the Dynamics_Trajectory_Output.dat output file
	fwrite << t << setw(15) << x << setw(15) << z << setw(15) 
		<< Modulus_v_E << endl;
	
    // Updating values of variables to obtain their values at time t
    u = u + RK[0];
	v = v + RK[1];
	w = w + RK[2];
	p = p + RK[3];
	q = q + RK[4];
	r = r + RK[5];
	phi = phi + RK[6];
	theta = theta + RK[7];
	psi = psi + RK[8];
	x = x + RK[9];
	y = y + RK[10];
    z = z + RK[11];	  
	
    // Incrementing time variable
    t = t + h;

}

return(Range);

}

// _Function which runs either the Optimizer OR the Dynamics Solver_
void DynamicsSolver::runSolverOROptimizer (){
	
    // Creating an output data file Dynamics_Range_Optimization_Output.dat in which the
    //  Dynamics Range Optimization results can be recorded
    ofstream fwrite;
    fwrite.open ( "Dynamics_Range_Optimization_Output.dat", ios :: trunc );
    
    if(optimizerON == 1){
		
	double range1 = 0, range2 = 1;
	LaunchAngle = initialLaunchAngle;
	
	while (range2 - range1 > 0) {
		range1 = DynamicsSolver_ThreeDoF();
		LaunchAngle += angleIncrement;
		range2 = DynamicsSolver_ThreeDoF();
		
		cout << "Range1 = " << range1 << "   " << "Range2 = " << range2 << endl;
		
	}
	
    LaunchAngle -= angleIncrement;
    range1 = DynamicsSolver_ThreeDoF();
    
    cout << "Optimal Launch Angle = " << LaunchAngle << endl;
	cout << "Maximum Range = " << range1 << endl;
	
	// Writing the Dynamics Range Optimization results data
    //  in the Dynamics_Range_Optimization_Output.dat output file
	fwrite << LaunchAngle << setw(15) << range1 << endl;
    
    }
    
    else {
		
    cout << "Range = " << DynamicsSolver_ThreeDoF () << endl;
   
    }
    
}

// Main function
int main() {

    DynamicsSolver obj;

    obj.takeInputs ();
    
    obj.runSolverOROptimizer ();

    cout << "THE DYNAMICS SOLVER AND RANGE OPTIMIZER ARE UP AND RUNNING."<< endl;

    return 0;
    
}
