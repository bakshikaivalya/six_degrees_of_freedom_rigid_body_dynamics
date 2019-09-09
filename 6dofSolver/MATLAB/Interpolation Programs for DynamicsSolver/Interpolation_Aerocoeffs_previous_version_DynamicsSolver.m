clc
clear

Ma = 0;
alpha = 2.68;

load EAT_Output_Table.dat
IntervalLengthAlpha = 0.25;
UpperBoundAlpha = 3;
IntervalsofAlpha = 1/IntervalLengthAlpha * UpperBoundAlpha;

lenEAT_Output_Table = length (EAT_Output_Table);

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

% Interpolation of Cx, Cz, Cm with y as independent variable at x = xi, x= xi+1
CX1 = EAT_Output_Table(m+n-1,6) + (EAT_Output_Table(m+n-1+1,6) - EAT_Output_Table(m+n-1,6))/(EAT_Output_Table(n+1,2) - EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));

CX2 = EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,6) + (EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1+1,6) - EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,6))/(EAT_Output_Table(n+1,2) - EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));

CZ1 = EAT_Output_Table(m+n-1,3) + (EAT_Output_Table(m+n-1+1,3) - EAT_Output_Table(m+n-1,3))/(EAT_Output_Table(n+1,2) - EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));

CZ2 = EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,3) + (EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1+1,3) - EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,3))/(EAT_Output_Table(n+1,2) - EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));

CM1 = EAT_Output_Table(m+n-1,5) + (EAT_Output_Table(m+n-1+1,5) - EAT_Output_Table(m+n-1,5))/(EAT_Output_Table(n+1,2)-EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));
   
CM2 = EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,5) + (EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1+1,5) - EAT_Output_Table(m+(IntervalsofAlpha+1)+n-1,5))/(EAT_Output_Table(n+1,2) - EAT_Output_Table(n,2)) * (abs(alpha) - EAT_Output_Table(n,2));

% Interpolation of Cx, Cz, Cm with x as independent variable at y = yi, y= yi+1
Cx = CX1 + (CX2 - CX1)/(EAT_Output_Table(m+(IntervalsofAlpha+1),1) - EAT_Output_Table(m,1)) * (Ma - EAT_Output_Table(m,1));

Cz = CZ1 + (CZ2 - CZ1)/(EAT_Output_Table(m+(IntervalsofAlpha+1),1) - EAT_Output_Table(m,1)) * (Ma - EAT_Output_Table(m,1));

Cm = CM1 + (CM2 - CM1)/(EAT_Output_Table(m+(IntervalsofAlpha+1),1) - EAT_Output_Table(m,1)) * (Ma - EAT_Output_Table(m,1));