clc

t = 0.000;

load Thrust_Table.dat;

lenThrust_Table = length (Thrust_Table);

for (c=1:1:lenThrust_Table)
    C = Thrust_Table(c,1) - t;
    if (C>0)
        thrustvar = c - 1;
       break;
    end;
end;

Thrust = Thrust_Table(thrustvar,2) + (Thrust_Table(thrustvar+1,2) - Thrust_Table(thrustvar,2))/(Thrust_Table(thrustvar+1,1) - Thrust_Table(thrustvar,1)) * (t - Thrust_Table(thrustvar,1));