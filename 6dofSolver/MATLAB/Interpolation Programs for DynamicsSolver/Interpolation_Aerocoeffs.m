clc;
clear all;
close all;

load EAT_Output_Table.dat;
load Thrust_Table.dat;
IntervalLengthAlpha = 0.2;
UpperBoundAlpha = 3;
IntervalsofAlpha = 1/IntervalLengthAlpha * UpperBoundAlpha;

Ma = 0;
alpha = 0.4;

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
        
        %Cmq =   - Cnalpha * (xcp - xcm)^2/d^2;