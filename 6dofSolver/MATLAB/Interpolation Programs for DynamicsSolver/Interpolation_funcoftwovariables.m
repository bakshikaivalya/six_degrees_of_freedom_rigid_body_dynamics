clc
Ma = 0.375;
alpha = 1.125;

load EAT_Output_Table.dat

len = length (EAT_Output_Table);

for (s=1:1:len)
    B = EAT_Output_Table(s,2) - alpha;
    if (B>=0)
        n = s;
       break;
    end;
end;
        
for (r=1:1:len)
    B = EAT_Output_Table(r,1) - Ma;
    if (B>=0)
        m = r;
        break;
    end;
end;

% Interpolation of Cx, Cz, Cm with y as independent variable at x = xi, x= xi+1
CX1 = EAT_Output_Table(m+n-1,6) + (EAT_Output_Table(m+n-1+1,6) - EAT_Output_Table(m+n-1,6))/(EAT_Output_Table(n+1,2)-EAT_Output_Table(n,2));

CX2 = EAT_Output_Table(m+13+n-1,6) + (EAT_Output_Table(m+13+n-1+1,6) - EAT_Output_Table(m+13+n-1,6))/(EAT_Output_Table(n+1,2)-EAT_Output_Table(n,2));

CZ1 = EAT_Output_Table(m+n-1,3) + (EAT_Output_Table(m+n-1+1,3) - EAT_Output_Table(m+n-1,3))/(EAT_Output_Table(n+1,2)-EAT_Output_Table(n,2));

CZ2 = EAT_Output_Table(m+13+n-1,3) + (EAT_Output_Table(m+13+n-1+1,3) - EAT_Output_Table(m+13+n-1,3))/(EAT_Output_Table(n+1,2)-EAT_Output_Table(n,2));

CM1 = EAT_Output_Table(m+n-1,5) + (EAT_Output_Table(m+n-1+1,5) - EAT_Output_Table(m+n-1,5))/(EAT_Output_Table(n+1,2)-EAT_Output_Table(n,2));

CM2 = EAT_Output_Table(m+13+n-1,5) + (EAT_Output_Table(m+13+n-1+1,5) - EAT_Output_Table(m+13+n-1,5))/(EAT_Output_Table(n+1,2)-EAT_Output_Table(n,2));

% Interpolation of Cx, Cz, Cm with x as independent variable at y = yi, y= yi+1
Cx = CX1 + (CX2 - CX1)/(EAT_Output_Table(m+13,1) - EAT_Output_Table(m,1))

Cz = CZ1 + (CZ2 - CZ1)/(EAT_Output_Table(m+13,1) - EAT_Output_Table(m,1))

Cm = CM1 + (CM2 - CM1)/(EAT_Output_Table(m+13,1) - EAT_Output_Table(m,1))






















