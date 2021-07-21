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
    Seccional secRK3 = Seccional("./simulaciones/seccional/SECC-PROD-finish_100");
    //insertarGrid(double v0, double R, int m, int numParticionesIntegrales, bool dominioEquiespaciado,double logPuntosGrid)
    secRK3.insertarGrid(1e-4,500,100,400,false,0.01);
    //insertarTiempo(double t0, double tFinal, double incTiempo)
    secRK3.insertarTiempo(0.001, 2, 1e-2);
    secRK3.calcular();
    
    /*
    cout<< "Volumenes finitos:---" << endl;
    VolFinitos vol1=VolFinitos("./simulaciones/vol_finitos/2s-vol_finitos_RK3_PROD_trapecios_compuesto_50");
    //insertarGrid(double xInicio, double R, double Nx, int numParticionesIntegrales, bool dominioEquiespaciado, double kdominio)
    vol1.insertarGrid(1e-4,500,400,400,false,1.1);
    vol1.insertarTiempo(0.001,2,1e-2);
    vol1.calcular();
     * */
    
    
    return 0;
}