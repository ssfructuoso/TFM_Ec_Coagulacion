/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: salvador
 *
 */

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <array>
#include <cstdio>
#include <math.h>  
#include <fstream>
#include <string>
#include "header.h"
#include <vector>
using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    /*
    cout<< "Seccional: (MidPoint Rule y RK3)" << endl;
    Seccional secRK3 = Seccional("./simulaciones/seccional/seccional_nx100_midpoint_1");
    
    //insertarGrid(double v0, double R, int m, int numParticionesIntegrales, bool dominioEquiespaciado, bool pasoAdaptativo)
    secRK3.insertarGrid(1e-4,500,100,400,false, true);
    //insertarTiempo(double t0, double tFinal, double incTiempo)
    secRK3.insertarTiempo(0.001, 2, 1e-2);
    secRK3.calcular();
    */
    
//////////////////////////////////////////    
//////////////////////////////////////////
    
    
    cout<< "Volumenes finitos: (MidPoint Rule y RK3)" << endl;
    VolFinitos volRK3=VolFinitos("./simulaciones/vol_finitos/volFinitos_nx100_midpoint_1");
    
    //insertarGrid(double xInicio, double R, double Nx, int numParticionesIntegrales, bool dominioEquiespaciado, bool pasoAdaptativo)
    volRK3.insertarGrid(1e-4,500,100,400,false, true);
    
    //insertarTiempo(double t0, double tFinal, double incTiempo)
    volRK3.insertarTiempo(0.001,2,1e-2);
    volRK3.calcular();
    
     
    
    return 0;
}
