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
    
    cout<< "Seccional:---" << endl;
    Seccional secMidPoint = Seccional("./simulaciones/seccional/2s-seccional_RK4_SUMA-mid_point_200");
    //insertarGrid(double v0, double R, int m, int numParticionesIntegrales, bool dominioEquiespaciado, bool pasoAdaptativo)
    secMidPoint.insertarGrid(1e-3,500,200,400,false,true);
    //insertarTiempo(double t0, double tFinal, double incTiempo)
    secMidPoint.insertarTiempo(0, 2, 1e-2);
    secMidPoint.calcular();
    
    
    /*
    cout<< "Volumenes finitos:---" << endl;
    VolFinitos volMidPoint=VolFinitos("./simulaciones/vol_finitos/2s-vol_finitos_RK4_ID_mid_point_50");
    //insertarGrid(double xInicio, double R, double Nx, int numParticionesIntegrales, bool dominioEquiespaciado, bool pasoAdaptativo)
    volMidPoint.insertarGrid(1e-10,50,50,400,false,true);
    volMidPoint.insertarTiempo(0,2,1e-2);
    volMidPoint.calcular();
    */
    return 0;
}
