// Copyright 2011 Zeus Numerix Pvt. Ltd.
// All rights reserved

#ifndef DEFINE_VARIABLES_H
#define DEFINE_VARIABLES_H

/// Class for defining the variables 
///
/// \author Anant Diwakar <anant.diwakar@zeusnumerix.com>


class DefineInputs {
        
    public :

        virtual ~DefineInputs() { }

        void TakeInputs();

        int Nfin, noseType, conditionPropellant; 
        double lengthNose, diameter, lengthBody1, thickness, tNoseExt;  
        double Xt, Cr, Ct, s, nFin, GammaC, tFin;
        double AR, cant_angle;
        double K, C;

        double yieldStress, fos, internalPress, density, massWarhead, loadfactorPropellant, tBody;

};
#endif
