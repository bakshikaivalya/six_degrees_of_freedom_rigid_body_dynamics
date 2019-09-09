clc
clear all
close all

LaunchAngle = 25;

approxgradientfactor = 1;
variable = 1;

while (approxgradientfactor>0)
    
    variable
    
    y1 = DynamicsSolver (LaunchAngle);
    y2 = DynamicsSolver (LaunchAngle + 1);
    
    approxgradientfactor = y2 - y1;
    
    LaunchAngle = LaunchAngle + 1;

    variable = variable + 1;
        
end;
    
Optimal_LaunchAngle = LaunchAngle
Maximum_Range = y1