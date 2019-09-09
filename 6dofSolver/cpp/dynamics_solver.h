// Copyright 2011 Zeus Numerix Pvt. Ltd.
// All rights reserved

// Author: Kaivalya Bakshi (kaivalya.bakshi@zeusnumerix.com)


#ifndef DYNAMICS_SOLVER_H
#define DYNAMICS_SOLVER_H

#include "define_variables.h"

class DynamicsSolver : public DefineInputs{

    public:
    
        virtual ~DynamicsSolver() { }

        // Declaration of all the C++ functions present in dynamics_solver.cpp
        void transposeMatrix (double input[][3], double output[][3]),
		 matrixMultiply3_by_3 (double input1[][3], double input2[][3], double output[][3]),
		 matrixInverse (double input[][3], double output[][3]);
		 
        double DynamicsSolver_ThreeDoF ();        
        void runSolverOROptimizer ();
		
    protected:    
    
};
#endif
