// Author: Kaivalya Bakshi (kaivalya.bakshi@zeusnumerix.com)

#include <stdio.h>
#include <iostream>
#include <math.h>

using namespace std;

#include "define_variables.h"

#define pi 3.141592653589793

/// Initializes the global variables used in the code.
void DefineInputs::takeInputs(){

    // _Specification of values of user input variables_
    // User inputs for DynamicsSolverforRangeOptimization function
    V = 50; // In m/s
	LaunchAngle = 59; // In degrees
    TimeIncrement = 0.001; // In seconds. Some constraint such as
    //  0.01<TimeIncrement should be imposed on user input
    IntervalLengthAlpha = 0.2; // In degrees.
    UpperBoundAlpha = 3; // In degrees. Some constraint such as
    //  3<UpperBoundAlpha should be imposed on user input
    
    // User inputs for Optimizer function
    optimizerON = 1; // optimizerON = 1 => Optimizer function is run, 
    //  optimizerON = Any other integer => DynamicsSolverforRangeOptimization is run
    //  to calculate Range and trajectory data 
    initialLaunchAngle = 58.5; // In degrees. This is the initial guess
    //  given to the optimizer to obtain the optimal LaunchAngle. Note that the 
    //  user must provide an initial guess which is lower in value than the optimal angle
    // for this optimization algorithm to work
    angleIncrement = 0.1; // In degrees. This is the angle increment used to estimate the gradient 
    //  of the Range as a function of the LaunchAngle at any LaunchAngle
    
    d = 0.214; // In metres. The same input should be provided to the EAT
    
}
