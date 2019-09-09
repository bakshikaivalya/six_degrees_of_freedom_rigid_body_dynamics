clc
clear all
close all

LaunchAngle = 25;

gradRange = 1;

variable = 1;

while (abs(gradRange)>0.5)
    
    variable
    
    y1 = DynamicsSolver (LaunchAngle);
    y2 = DynamicsSolver (LaunchAngle + 0.01);    
    
    gradRange = (y2 - y1)/(0.01);
    
    LaunchAngle
    ImprovedRange = y1
    LaunchAngle = LaunchAngle + gradRange * 0.5;
    
    variable = variable + 1;
        
end;
    
    Optimal_LaunchAngle = LaunchAngle
    Maximum_Range = y1