// Copyright 2011 Zeus Numerix Pvt. Ltd.
// All rights reserved

#ifndef CALCULATE_STRUCTURAL_PROP_H
#define CALCULATE_STRUCTURAL_PROP_H

/// Class for calculating structural properties and center of 
/// gravity for a given rocket assembly
///
/// \author Anant Diwakar <anant.diwakar@zeusnumerix.com>

#include"define_variables.h"

class CalculateStructuralProp : public DefineInputs{

    public:
    
        virtual ~CalculateStructuralProp() { }

        double CalcNoseVolume( double, double );
        double NoseShapeFunc( double, double, double );
        double CalculateNozzleThickness();
        double CalculateNoseProp();
        double CalculateNoseExtnProp();
        double CalculateBodyProp();
        double CalculateNozzleProp();
        double CalculateFinProp();
        double CalculateAssemblyProp();

        double intDia, lengthNose, centroidNose, massNose, intVolumeNose, centroidNoseExt, massNoseExt, massratioNose, massratioNoseExt, massWarheadAssembly, lengthWarhead, izNose, ixNose, izNoseExt, ixNoseExt, lengthBody, massBodyAssembly, centroidBody, izBody, ixBody, lengthNozzle, massNozzle,  centroidNozzle, izNozzle, ixNozzle, massFin, centroidFin, izFin, ixFin, lengthAssembly, massAssembly, centroidCompleteAssembly, CM, totalMass, Ixx, Iyy, Izz;

    protected:    
};
#endif
