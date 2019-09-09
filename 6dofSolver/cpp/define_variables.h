// Copyright 2011 Zeus Numerix Pvt. Ltd.
// All rights reserved

// Author: Kaivalya Bakshi (kaivalya.bakshi@zeusnumerix.com)

#ifndef DEFINE_VARIABLES_H
#define DEFINE_VARIABLES_H

class DefineInputs {
        
    public :

        virtual ~DefineInputs() { }

        void takeInputs();

        // _Declaration of all the user input variables_
        // Variables for the DynamicsSolverforRangeOptimization function
        double V, initialLaunchAngle, TimeIncrement, IntervalLengthAlpha, 
			UpperBoundAlpha, LaunchAngle; 
			
		// Structural variable given as input to the structural module
		double d; 
			
		// Variables concerning optimization function of code
		int optimizerON; // Variable used to choose between running the
		//  optimizer function or the DynamicsSolverforRangeOptimization function
		double angleIncrement; // Angle increment used to estimate the gradient 
		//  of the Range as a function of the LaunchAngle at any LaunchAngle

};
#endif
