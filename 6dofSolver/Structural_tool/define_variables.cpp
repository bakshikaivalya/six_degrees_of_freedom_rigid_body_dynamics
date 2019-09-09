

#include<stdio.h>
#include<iostream>
#include<math.h>

using namespace std;

#include"define_variables.h"

#define pi 3.142857143

void DefineInputs::TakeInputs(){

    //-------Nose and Body dimensions ( All dimensions in cms )-----------

    noseType = 3;  // 1 - Conical, 2 - Elliptic, 3 - Power series - 1/2 Parabolic

    K = 1;  // Full parabolic nose cone
    C = 0;  // LD haack - von karman

    
    diameter = 21.4;
    lengthBody1 = 438.1;
    tNoseExt = 0.5;
    tBody = 0.3;

    //-------------Fin parameters-----------------

    nFin = 4;
    Cr = 27.299;
    Ct = 20.295;
    s = 14.38;
    tFin = 0.3;                 // fin thickness
    cant_angle = 0.83;       // degrees   
    GammaC = 15;

    GammaC = ( pi / 180) * GammaC;

    cant_angle = ( pi / 180) * cant_angle;

    //---------------Structural properties---------------

    yieldStress = 850;
    fos = 1.4;
    internalPress = 14;
    density = 7800e-06;
    massWarhead = 100;

    //--------------Propellant properties---------------
    loadfactorPropellant = 0.8;
    conditionPropellant = 1;     //1 - Initial, 2 - Burnout
}
