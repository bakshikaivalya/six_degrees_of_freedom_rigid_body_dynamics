function [Range] = DynamicsSolver (LaunchAngle)

% _User Inputs_
V = 50;
LaunchHeadingAngle = 0;
LaunchFlightPathAngle = LaunchAngle * pi/180;
TimeIncrement = 0.001;

% _Input files and User Inputs related to Creation of Input Files_
load EAT_Output_Table.dat;
load Thrust_Table.dat;
IntervalLengthAlpha = 0.2;
UpperBoundAlpha = 3;
IntervalsofAlpha = 1/IntervalLengthAlpha * UpperBoundAlpha;

% _Initialization of vairables_

% Initialization of s_BE_E
x = 0;
y = 0;
z = 0;
s_BE_E = [x; y; z];

% Initialization of v_E_B
u =  - V;
v = 0;
w = 0;
v_E_B = [u; v; w];
     
% Initialization of phi, theta, psi
phi = 0;
theta = LaunchFlightPathAngle;
psi = pi + LaunchHeadingAngle;
     
% Initialization of omega vector and Omega tensor
p = 0;
q = 0;
r = 0;

% Vehicle Geometry Parameters
d = 0.214; 
A_cs = (3.14159*d*d)/4 ;
massi = 291.736597336487;  massf = 154.30324858168; 
Ixxi =  1.78147455320142;  Ixxf =  0.994737348254526;
Ixyi =  0;  Ixyf = 0;
Ixzi =  0;  Ixzf = 0;
Iyxi =  0;  Iyxf = 0;
Iyyi =  483.245811997219;  Iyyf = 270.097539984996;
Iyzi =  0;  Iyzf = 0;
Izxi =  0;  Izxf = 0;
Izyi =  0;  Izyf = 0;
Izzi =  483.245811997219;  Izzf = 270.097539984996;
xcmi =  2.48499781117085;  xcmf = 1.88691063312857;

% Initialization of gravitational acceleration vector
g = [0; 0; -9.81];

% Initialization of time variable and specification of burn time
t = 0;
tbr = 4.9;

% Time Step 
h = TimeIncrement;

% Initialization of Kutta Matrix
for (i=1:12)
    k(1,i) = 0;
end;

% _Solution of ODEs by Runge Kutta 4 Method_
while (z >= 0)
    
    for j=2:1:5
        
        if (j==5);
            timestep=h;
        else
            timestep=h/2.0;
        end
            
        ur = u + timestep * k(j-1,1);
        vr = v + timestep * k(j-1,2);
        wr = w + timestep * k(j-1,3);
        pr = p + timestep * k(j-1,4);
        qr = q + timestep * k(j-1,5);
        rr = r + timestep * k(j-1,6);
        phir = phi + timestep * k(j-1,7);
        thetar = theta + timestep * k(j-1,8);
        psir = psi + timestep * k(j-1,9);
        xr = x + timestep * k(j-1,10);
        yr = y + timestep * k(j-1,11);
        zr = z + timestep * k(j-1,12);
    
        % Calculation of Matrix of Inertia Tensor by Linear Interpolation ...
        % of comprising elements with time and Calculation of Thrust ...
        % by Linear Interpolation over time variable
        if (t<=tbr);
            mass = massi - (massi - massf)*(t/tbr);
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
   
            lenThrust_Table = length (Thrust_Table);

            for (c=1:1:lenThrust_Table)
                C = Thrust_Table(c,1) - t;
                if (C>0)
                    thrustvar = c - 1;
                    break;
                end;
            end;

            Thrust = Thrust_Table(thrustvar,2) + (Thrust_Table(thrustvar+1,2) - Thrust_Table(thrustvar,2))/(Thrust_Table(thrustvar+1,1) - Thrust_Table(thrustvar,1)) * (t - Thrust_Table(thrustvar,1));
            
        else
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
        end;      
        Thrust = 9.81 * Thrust;        
        
        Inertia_Matrix = [Ixx Ixy Ixz; Iyx Iyy Iyz; Izx Izy Izz]; 
 
        % Calculation of Coordinates Transformation Matrix  
        t11 = (cos(psir)) * cos(thetar);
        t12 = (sin(psir)) * cos(thetar);
        t13 = - sin(thetar);
        t21 = (cos(psir)) * (sin(thetar)) * sin(phir) - (sin(psir)) * cos(phir);
        t22 = (sin(psir)) * (sin(thetar)) * sin(phir) + (cos(psir)) * cos(phir);
        t23 = (cos(thetar)) * sin(phir);
        t31 = (cos(psir)) * (sin(thetar)) * cos(phir) + (sin(psir)) * sin(phir);
        t32 = (sin(psir)) * (sin(thetar)) * cos(phir) - (cos(psir)) * sin(phir);
        t33 = (cos(thetar)) * cos(phir);
        T_BE = [t11 t12 t13; t21 t22 t23; t31 t32 t33];
        T_EB = transpose(T_BE);
%__________________________________________________________________________    
        % _Calculation of rho, Ma, alpha, qbar_
        if (zr<11000)
            rho = 1.225 * (1.0 - 0.000022557695 * zr) ^ (4.25587);
        else
            rho = 0.36392 * exp(- 1.576584e-04 * (zr - 11000));
        end;
    
        vs = 340.249 * sqrt((288.15 - 0.0065 * zr)/288.15);
        Ma = sqrt (ur^2 + vr^2 + wr^2)/vs;
    
        %if (t<tbr)
            alpha_radians = atan (wr/ur);
            alpha = alpha_radians * 180/pi;%0.3 * sin (2000*t);
         %else
            %alpha = 0;
        %end;
    
        qbar = 0.5 * rho * (ur^2 + vr^2 + wr^2);
    
        % _Calculation of Cx, Cy, Cz, Cm_
        lenEAT_Output_Table = length (EAT_Output_Table);

        for (a=1:1:lenEAT_Output_Table)
            A = EAT_Output_Table(a,2) - abs(alpha);
            if (A>0)
                n = a - 1;
                break;
            end;
        end;
        
        for (b=1:1:lenEAT_Output_Table)
            B =  EAT_Output_Table(b,1) - Ma;
            if (B>0)
                m = b - (IntervalsofAlpha+1);
                break;
            end;
        end;

        % Interpolation of Cx, Cz, Cm0 with y as independent variable at x = xi, x= xi+1
        CX1 = EAT_Output_Table(m+n-1,7) + (EAT_Output_Table(m+n-1+1,7) - EAT_Output_Table(m+n-1,7))/(EAT_Output_Table(n+1,2) - EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));

        CX2 = EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,7) + (EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1+1,7) - EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,7))/(EAT_Output_Table(n+1,2) - EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));

        CZ1 = EAT_Output_Table(m+n-1,4) + (EAT_Output_Table(m+n-1+1,4) - EAT_Output_Table(m+n-1,4))/(EAT_Output_Table(n+1,2) - EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));

        CZ2 = EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,4) + (EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1+1,4) - EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,4))/(EAT_Output_Table(n+1,2) - EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));

        CM01 = EAT_Output_Table(m+n-1,6) + (EAT_Output_Table(m+n-1+1,6) - EAT_Output_Table(m+n-1,6))/(EAT_Output_Table(n+1,2)-EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));
    
        CM02 = EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,6) + (EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1+1,6) - EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,6))/(EAT_Output_Table(n+1,2) - EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));
        
        CNALPHA1 = EAT_Output_Table(m+n-1,3) + (EAT_Output_Table(m+n-1+1,3) - EAT_Output_Table(m+n-1,3))/(EAT_Output_Table(n+1,2) - EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));
        
        CNALPHA2 = EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,3) + (EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1+1,3) - EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,3))/(EAT_Output_Table(n+1,2) - EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));
        
        XCP1 = EAT_Output_Table(m+n-1,5) + (EAT_Output_Table(m+n-1+1,5) - EAT_Output_Table(m+n-1,5))/(EAT_Output_Table(n+1,2) - EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));
        
        XCP2 = EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,5) + (EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1+1,5) - EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,5))/(EAT_Output_Table(n+1,2) - EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));

        % Interpolation of Cx, Cz, Cm0 with x as independent variable at y = yi, y= yi+1
        Cx = CX1 + (CX2 - CX1)/(EAT_Output_Table(m+(IntervalsofAlpha+1),1) - EAT_Output_Table(m,1)) * (Ma - EAT_Output_Table(m,1));

        Cz = CZ1 + (CZ2 - CZ1)/(EAT_Output_Table(m+(IntervalsofAlpha+1),1) - EAT_Output_Table(m,1)) * (Ma - EAT_Output_Table(m,1));

        Cm0 = CM01 + (CM02 - CM01)/(EAT_Output_Table(m+(IntervalsofAlpha+1),1) - EAT_Output_Table(m,1)) * (Ma - EAT_Output_Table(m,1));
        
        Cnalpha = CNALPHA1 + (CNALPHA2 - CNALPHA1)/(EAT_Output_Table(m+(IntervalsofAlpha+1),1) - EAT_Output_Table(m,1)) * (Ma - EAT_Output_Table(m,1));
        
        xcp = (XCP1 + (XCP2 - XCP1)/(EAT_Output_Table(m+(IntervalsofAlpha+1),1) - EAT_Output_Table(m,1)) * (Ma - EAT_Output_Table(m,1)))/1000;
        
        Cmq =   - Cnalpha * (xcp - xcm)^2/d^2;
%__________________________________________________________________________    
        % Calculation of Aerodynamic Forces and Moments
    
        Fx_aero = Cx * qbar * A_cs;

        Fy_aero = 0;

        if (alpha>0)
            Fz_aero = Cz * qbar * A_cs;
        else
            Fz_aero = - Cz * qbar * A_cs;
        end;

       %Cm = Fz_aero * (xcp - xcm)/(qbar * A_cs * d) + Cmq * qr * d/(2 * sqrt (ur^2 + vr^2 + wr^2));
        
        Mx_aero = 0;
        
        if (alpha>0)
            %My_aero =  - ((Fz_aero * (xcp - xcm)/(qbar * A_cs * d)) + Cmq * ((q * d)/(2 * sqrt (ur^2 + vr^2 + wr^2))))* qbar * A_cs * d;
            My_aero = (- Cm0 +  Cmq * ((qr * d)/(2 * sqrt (ur^2 + vr^2 + wr^2)))) * qbar * A_cs * d;
        else 
            My_aero = (Cm0 +  Cmq * ((qr * d)/(2 * sqrt (ur^2 + vr^2 + wr^2)))) * qbar * A_cs * d;
        end;
        
        Mz_aero = 0;
    
        % Calculation of Moment of Force due to Displacement of Aero Data ...
        % Reference point from the CM at instant of time t
        Mx_aeroref_cm = ycm * Fz_aero - zcm * Fy_aero;
    
        My_aeroref_cm = zcm * Fx_aero -  (- xcm) * Fz_aero;
    
        Mz_aeroref_cm = xcm * Fy_aero - ycm * Fx_aero;
    
        % Calculation of Forces and Moments due to Propulsion ...
        % Note that components of Moment due to Thrust Force are Zero since ...
        %  the Thrust Force is aligned with X axis of Body Coordinate System
        Fx_thrust = - Thrust;
    
        Fy_thrust = 0;
    
        Fz_thrust = 0;
    
        % Calculation of Forces due to Gravitational Force
        % Note that local level reference frame is used so that ...
        % the Earth is not rotating and is flat
        F_g_matrix = T_BE * (mass * g);
    
        Fx_total = Fx_aero + Fx_thrust + F_g_matrix (1);
        Fy_total = Fy_aero + Fy_thrust + F_g_matrix (2);
        Fz_total = Fz_aero + Fz_thrust + F_g_matrix (3);
  
        Mx_total = Mx_aero + Mx_aeroref_cm;
        My_total = My_aero + My_aeroref_cm;
        Mz_total = Mz_aero + Mz_aeroref_cm;
    
        % Calculation of k2, k3, k4, k5 for u, v, w
        omega_BE_B = [pr; qr; rr];
        Omega_BE_B = [0 -rr qr; rr 0 -pr; -qr pr 0];
        veldot_E_B = [(Fx_total/mass); (Fy_total/mass); (Fz_total/mass)] - Omega_BE_B * [ur; vr; wr];
         
        k(j,1) = veldot_E_B (1);
        k(j,2) = veldot_E_B (2);
        k(j,3) = veldot_E_B (3);
     
        % Calculation of k2, k3, k4, k5 for p, q, r
        omegadot_E_B = inv(Inertia_Matrix) * ([Mx_total; My_total; Mz_total] - Omega_BE_B * Inertia_Matrix * omega_BE_B);
    
        k(j,4) = omegadot_E_B (1);
        k(j,5) = omegadot_E_B (2);
        k(j,6) = omegadot_E_B (3);
    
        % Calculation of k2, k3, k4, k5 for psi, theta, phi
        psithetaphidot = [1 (sin(phir) * tan(thetar)) (cos(phir) * tan(thetar)); 0 cos(phir) -sin(phir); 0 (sin(phir)/cos(thetar)) (cos(phir)/cos(thetar))] * [pr; qr; rr];
    
        k(j,7) = psithetaphidot (1);
        k(j,8) = psithetaphidot (2);    
        k(j,9) = psithetaphidot (3);
    
        % Calculation of k2, k3, k4, k5 for x, y, z
        xyzdot = T_EB * [ur; vr; wr];
    
        k(j,10) = xyzdot (1);
        k(j,11) = xyzdot (2);
        k(j,12) = xyzdot (3);
            
        end;
    
    
    % Calculation of delta's of RK4 method for the 12 State Variables
    for incrementingvariable=1:12;
          RK(incrementingvariable,1) = (h/6.0) * (k(2,incrementingvariable) + 2*k(3,incrementingvariable) + 2*k(4,incrementingvariable) + k(5,incrementingvariable));
    end;
    
    % Calculating Range and Velocity as functions of time
        
    Range = sqrt (x^2 + y^2);
    Modulus_v_E = sqrt (u^2 + v^2 + w^2);
    Altitude = z;
    
    % Updating values of variables to obtain their values at time t    
    u = u + RK(1,1);
    v = v + RK(2,1);
    w = w + RK(3,1);
    p = p + RK(4,1);
    q = q + RK(5,1);
    r = r + RK(6,1);
    phi = phi + RK(7,1);
    theta = theta + RK(8,1);
    psi = psi + RK(9,1);
    x = x + RK(10,1);
    y = y + RK(11,1);
    z = z + RK(12,1);
           
    % Incrementing time variable
    t = t + h;
        
    plot (Range, Altitude);
    %plot (t, Modulus_v_E);
    %plot (t, Altitude);
    %plot (t, theta);
    hold on;
end;