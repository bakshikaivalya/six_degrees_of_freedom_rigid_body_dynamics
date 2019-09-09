// Copyright 2011 Zeus Numerix Pvt. Ltd.
// All rights reserved

/// Class for calculating structural properties and center of 
/// gravity for a given rocket assembly
///
/// \author Sudhir Muthyala <sudhir.m@zeusnumerix.com>

/// Converted into c++ by Anant Diwakar <anant.diwakar@zeusnumerix.com>


#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<math.h>

using namespace std;

#define pi 3.142857143

#include"calculate_structural_prop.h"
#include"define_variables.cpp"

double CalculateStructuralProp::NoseShapeFunc( double rd, double ln, double val ){

    double r, y, x, h, rn, K, C, theta;

    x = val; r = rd; h = ln; 

    switch (noseType) {

        case 1:  // Conical nose cone

            y =  ( x * r / h );

            break;

        case 2: // Elliptic nose cone

            y = r * pow (( 1 - ( x * x / ( h * h ))), 0.5 );

            break;

        case 3: // Power series 1/2 parabola

            y = r * pow ( ( x / h ), 0.5 );

            break;

        case 4: // Parabolic

            y = r * ( ( 2 * ( x / h ) - K * ( x / h ) * 
                        ( x / h ) ) / ( 2 - K ) );

            break;

        case 5: // Haack series 

            theta = acos( 1 - 2 * x / h );

            y = r * pow ( ( theta - 0.5 * sin ( 2 * theta ) + 
                        C * pow ( sin ( theta ), 3 ) ), 0.5 ) / sqrt ( pi );

            break;

        case 6: // Ogive

            rn = ( r * r + h * h ) / ( 2 * r );

            y = ( sqrt ( rn * rn - ( x - h ) * ( x - h ) ) 
                    + ( r - rn ) );

            break;

        default:

            cout << " nose type not defined " << endl;

            break;
    }

    return ( y );
}

double CalculateStructuralProp::CalcNoseVolume( double dia, double len ){

    double dh, a, b, s1, s2, x, integralValue;

    double rad = 0.5 * dia; 
    
    int i, nInterval = 10;

    a = 0; b = len;

    dh = ( b - a ) / ( 2 * nInterval );

    s1 = 0.0; s2 = 0.0; integralValue = 0.0;
    
    for ( i = 1; i <= nInterval; i++ ){
    
        x = a + dh * ( 2 * i - 1 );

        s1 += pi * pow (NoseShapeFunc( rad, len, x ), 2);
    }

    for ( i = 1; i < nInterval; i++ ){
        
        x = a + dh * 2 * i;
        
        s2 += pi * pow (NoseShapeFunc( rad, len, x ), 2);
    }

    integralValue  = dh * ( pi * pow (NoseShapeFunc( rad, len, a ), 2) + 
            pi * pow (NoseShapeFunc( rad, len, b ), 2) + 4 * s1 + 2 * s2 ) / 3;

    cout << " Volume = " << integralValue << endl;

    return (integralValue);
}

double CalculateStructuralProp::CalculateNozzleThickness(){

    double tol = 1e-5, err = 1;

    double tNozzleN, tNozzle; 

    double maxStress = yieldStress / fos;

    double meanDia = diameter;

    while (err > tol) {

        tNozzleN = (internalPress * meanDia) / (2 * maxStress); 

        intDia = diameter - 2 * tNozzleN;

        meanDia = ( intDia + diameter )/2;

        tNozzle = (internalPress * meanDia) / (2 * maxStress);

        err = abs(tNozzle - tNozzleN);
    }

    cout << "nozzle thickness : " << tNozzle << endl;
}

double CalculateStructuralProp::CalculateNoseProp(){

    double extVolume, tNose, intDiaNose;
    double rad_n, rad_n1, Mom_nose, sa_nose;

    if (diameter <= 22.5) {

        tNose = 0.3; 

    } else {

        tNose = 0.4;
    }

    
    lengthNose = 2.8*diameter;
    intDiaNose = diameter - 2 * tNose;

    switch ( noseType ){

        case 1:  // Conical

            extVolume = (1.0/3) * pi * pow (( 0.5 * diameter ), 2) * lengthNose;

            intVolumeNose = (1.0/3) * pi * pow (( 0.5 * intDiaNose ), 2) * 
                (lengthNose - tNose * 5); 

            centroidNose = (2.0/3) * lengthNose;

            izNose = (3.0/10) * (extVolume - intVolumeNose ) * density * 
                diameter * diameter;

            ixNose = (3.0/80) * ( extVolume - intVolumeNose ) * density * 
                ( 4 * ( diameter / 2) * ( diameter / 2) + 
                  lengthNose * lengthNose );
            
            break;

        case 2:  // Elliptical

            extVolume = (2.0/3) * pi * pow (( 0.5 * diameter ), 2) * lengthNose;

            intVolumeNose = (2.0/3) * pi * pow (( 0.5 * intDiaNose ), 2) * 
                (lengthNose - tNose * 3);
          
            centroidNose =  lengthNose - ((4/(3*pi)) * lengthNose);

            izNose = (extVolume - intVolumeNose ) * density * 
                pow (( 0.5 * diameter ), 2) / 5;
    
            ixNose = (extVolume - intVolumeNose ) * density * (( diameter / 2) * 
                    ( diameter / 2) + lengthNose * lengthNose )/ 10;

            break;

        case 3:   // Power series - 1/2 Parabola

            extVolume = (1.0/2) * pi * pow (( 0.5 * diameter ), 2) * lengthNose;

            intVolumeNose = (1.0/2) * pi * pow (( 0.5 * intDiaNose ), 2) * 
                (lengthNose - tNose * 3);
          
            centroidNose = (5.0/8) * lengthNose;
            
            izNose = (3.0/10) * (extVolume - intVolumeNose ) * density * 
                pow (( 0.5 * diameter ), 2); 
    
            ixNose = (3.0/80) * (extVolume - intVolumeNose ) * density * 
                ( 4 * ( diameter / 2) * ( diameter / 2) + 
                  lengthNose * lengthNose );

            break;

        case 4:    // Parabolic 

            extVolume = (8.0/15) * pi * pow (( 0.5 * diameter ), 2) * lengthNose;

            intVolumeNose = (8.0/15) * pi * pow (( 0.5 * intDiaNose ), 2) * 
                (lengthNose - tNose * 4);
            
            centroidNose = (1.0/4) * lengthNose;
        
            izNose = (3.0/10) * (extVolume - intVolumeNose ) * density * 
                pow (( 0.5 * diameter ), 2); 
    
            ixNose = (3.0/80) * (extVolume - intVolumeNose ) * density * 
                ( 4 * ( diameter / 2) * ( diameter / 2) + 
                  lengthNose * lengthNose );

            break;

        case 5:   // Haack series

            extVolume = CalcNoseVolume( diameter, lengthNose );
                      
            intVolumeNose = CalcNoseVolume( intDiaNose, lengthNose - tNose * 4 );
            
            centroidNose = ((5.0/8) * lengthNose + (1/4) * lengthNose )/2 ; 
        
            izNose = (3.0/10) * (extVolume - intVolumeNose ) * density * 
                pow (( 0.5 * diameter ), 2); 
    
            ixNose = (3.0/80) * (extVolume - intVolumeNose ) * density * 
                ( 4 * ( diameter / 2) * ( diameter / 2) + 
                  lengthNose * lengthNose );

            break;

        case 6:  // Ogive

            rad_n = ( (diameter/2) * (diameter/2) + lengthNose * lengthNose )
                / diameter;
        
            extVolume = pi * (lengthNose * rad_n * rad_n - 
                    pow (lengthNose, 3) / 3 - (rad_n - diameter/2) * 
                    rad_n * rad_n * asin (lengthNose/rad_n));
        
            rad_n1 = ( pow((intDiaNose/2), 2) + 
                    pow((lengthNose-tNose * 4), 2) ) / intDiaNose;
        
            intVolumeNose = pi * ((lengthNose-tNose * 4) * rad_n1 * rad_n1 - 
                    pow((lengthNose-tNose * 4), 3)/3 - (rad_n1 - intDiaNose/2) * 
                    rad_n1 * rad_n1 * asin((lengthNose-tNose * 4)/rad_n1));
        
            Mom_nose = 2 * pi * rad_n * ( 0.5 * lengthNose * lengthNose + 
                    (diameter/2 - rad_n) * ( rad_n - 
                        sqrt( rad_n * rad_n - lengthNose * lengthNose) ));
        
            sa_nose = 2 * pi * rad_n * ( ( (diameter/2 - rad_n) * 
                        asin(lengthNose/rad_n) ) + lengthNose );
                
            centroidNose =  lengthNose - (Mom_nose/sa_nose);
                
            izNose = (3.0/10) * (extVolume - intVolumeNose ) * density * 
                diameter * diameter; 

            ixNose = (3.0/80) * (extVolume - intVolumeNose ) * density * 
                ( 4 * ( diameter / 2) * ( diameter / 2) + 
                  lengthNose * lengthNose );

            break;

        default:

            cout << " This nose type not defined " << endl;

            break;
    }

    massNose = (extVolume - intVolumeNose) * density;

    cout << "Nose mass: " << massNose << endl;
}

double CalculateStructuralProp::CalculateNoseExtnProp(){

    double extVolume, intVolumeNoseExt, lengthNoseExt, intDiaNoseExt;

    intDiaNoseExt = diameter - 2 * tNoseExt;
    
    lengthNoseExt = 5.2*diameter;
    
    extVolume = pi * pow( (diameter / 2), 2) * lengthNoseExt;

    intVolumeNoseExt = pi * pow ((intDiaNoseExt / 2), 2) * lengthNoseExt;

    massNoseExt = (extVolume - intVolumeNoseExt) * density;

    centroidNoseExt = (1.0/2) * lengthNoseExt;

    massratioNose = (intVolumeNose / (intVolumeNoseExt + intVolumeNose) ) * (massWarhead - massNose - massNoseExt);

    massratioNoseExt = (intVolumeNoseExt / (intVolumeNoseExt + intVolumeNose) ) * (massWarhead - massNose - massNoseExt) ;


    massWarheadAssembly = massWarhead;
    
    lengthWarhead = lengthNose + lengthNoseExt;

    izNoseExt = (massNoseExt + massratioNoseExt) * pow( (diameter / 2), 2) / 2;

    ixNoseExt = (1.0/12) * (massNoseExt + massratioNoseExt) * (3 * pow ((diameter / 2), 2) + pow (lengthNoseExt, 2));
    }

double CalculateStructuralProp::CalculateBodyProp() {

    double extVolume, intVolume, intDiaBody, massPropellant, massBodyExt, densityPropellant;

    lengthBody = 13.5 * diameter;

    intDiaBody = diameter - 2 * tBody;

    extVolume = pi * pow ( (diameter/2), 2) * lengthBody;

    intVolume = pi * pow ((intDiaBody/2), 2) * lengthBody;

    switch (conditionPropellant) {

        case 1:  // Initial  

            densityPropellant = 1750e-06;

            break;

        case 2: // Burnout

            densityPropellant = 0;

            break;

            default:

        cout << " This propellant condition is not defined " << endl;
    }

    massPropellant = densityPropellant * intVolume * loadfactorPropellant;

    massBodyExt = (extVolume - intVolume )* density;

    massBodyAssembly = massPropellant + massBodyExt;

    centroidBody = (1.0/2) * lengthBody;

    izBody = (massBodyAssembly) * pow ((diameter/2), 2) / 2;

    ixBody = (1.0/12) * (massBodyAssembly) * (3 * pow ((diameter/2), 2) + pow (lengthBody, 2) );
}

double CalculateStructuralProp::CalculateNozzleProp() {

    double lengthNozzle, extVolume, intVolume;

    lengthNozzle = 1.8 * diameter;

    extVolume = pi * pow ((diameter/2), 2) * lengthNozzle;

    intVolume = pi * pow ((intDia/2), 2) * lengthNozzle;

    massNozzle = (extVolume - intVolume )* density;

    centroidNozzle = (1.0/2) * lengthNozzle;

    izNozzle = massNozzle * pow ((diameter/2), 2) / 2;

    ixNozzle = (1.0/12) * (massNozzle) * (3 * pow ((diameter/2), 2) + pow (lengthNozzle, 2));
    }

double CalculateStructuralProp::CalculateFinProp() {

    double lengthFin, heightFin, Volume;

    lengthFin = 1.4 * diameter;

    heightFin = 0.7 * diameter;

    Volume = lengthFin * heightFin * tFin * nFin * 1.1;

    massFin = Volume * density;

    centroidFin = (1.0/2) * lengthFin;

    izFin = (1.0/12) * massFin * (pow (tFin, 2) + pow (heightFin, 2)) + massFin * pow ((diameter/2 + heightFin/2), 2);

    if (nFin == 0){

    nFin=1;
    }

    else {}

    ixFin = (1.0/12) * (massFin/nFin) * 2 * (2 * pow (lengthFin, 2) + pow (tFin, 2) + pow (heightFin, 2)) + 2 * (massFin/nFin) * pow ((diameter/2 + heightFin/2), 2);
}

double CalculateStructuralProp::CalculateAssemblyProp() {

    double centroidNoseExtO, centroidBodyO, centroidNozzleO, centroidFinO; 

    lengthAssembly = lengthWarhead + lengthBody + lengthNozzle;

    massAssembly = massWarheadAssembly + massBodyAssembly + massNozzle + massFin;

    totalMass = massAssembly;

    centroidNoseExtO = centroidNoseExt + lengthNose;

    centroidBodyO = centroidBody + lengthWarhead;

    centroidNozzleO = centroidNozzle + lengthWarhead + lengthBody;

    centroidFinO = centroidFin + lengthWarhead + lengthBody;

    centroidCompleteAssembly = ( (centroidNose * (massNose + massratioNose)) + (centroidNoseExtO * (massNoseExt + massratioNoseExt)) + (centroidBodyO * massBodyAssembly) + (centroidNozzleO * massNozzle) + (centroidFinO * massFin) ) / massAssembly;

    CM = centroidCompleteAssembly;

    Ixx = ixNose + (massNose + massratioNose ) * pow ((centroidCompleteAssembly - centroidNose), 2) + ixNoseExt + (massNoseExt + massratioNoseExt) * pow ((centroidCompleteAssembly - centroidNoseExtO), 2) + ixBody + (massBodyAssembly) * pow ((centroidCompleteAssembly - centroidBodyO), 2) + ixNozzle + (massNozzle) * pow ((centroidCompleteAssembly - centroidNozzleO), 2) + ixFin + (massFin) * pow ((centroidCompleteAssembly - centroidFinO), 2);

    Iyy = Ixx;

    Izz = izNose + izNoseExt + izBody + izNozzle + izFin;

    cout << "Total mass of vehicle is " << totalMass << endl;

    cout << "Position of Centroid of vehicle relative to body coordinate system with origin at nose tip is " << CM << endl;

    cout << "Ixx of vehicle is " << Ixx << endl;

    cout << "Iyy of vehicle is " << Iyy << endl;

    cout << "Izz of vehicle is " << Izz << endl;
}

int main (){

    CalculateStructuralProp obj;

    obj.TakeInputs();
    obj.CalculateNozzleThickness();
    obj.CalculateNoseProp();
    obj.CalculateNoseExtnProp();
    obj.CalculateBodyProp();
    obj.CalculateNozzleProp();
    obj.CalculateFinProp();
    obj.CalculateAssemblyProp();

    cout << "Hello " << endl;

    return 0;
}
