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
        
    //Testeo de los dos metodos para Oceanografia
    /*    
    cout<< "Seccional:---" << endl;
    Seccional secRK3 = Seccional("./simulaciones/SECC_test-100_e-13mol_Z3_6m_EXPONENCIAL_1m");
    //insertarGrid(double v0, double R, int m, int numParticionesIntegrales, bool dominioEquiespaciado, bool pasoAdaptativo)
    secRK3.insertarGrid(1.9268e-15,5.13925e-6,50,400,false,false);
    //insertarTiempo(double t0, double tFinal, double incTiempo)
    secRK3.insertarTiempo(0,2600000,20); //1 mes
    secRK3.calcular();
     */
    
    
    
    
    cout<< "Volumenes finitos:---" << endl;
    VolFinitos volRK3=VolFinitos("./simulaciones/VOL_TESTs_e-13mol_Z3_6m_EXPONENCIAL_1m");
    //insertarGrid(double xInicio, double R, double Nx, int numParticionesIntegrales, bool dominioEquiespaciado,bool pasoAdaptativo)
    volRK3.insertarGrid(1.9268e-15,5.13925e-6,49,400,false,false);//SI (2.45044e-15, 4.23436e-6)
    volRK3.insertarTiempo(0,2600000,20); //1 mes
    volRK3.calcular();
   
    
    
    
    //////Fin testeo
    
    
    //Simulacion para el caso i) del capitulo 3 de la memoria del TFM.
    /*
    cout<< "Volumenes finitos:---" << endl;
    VolFinitos volRK3=VolFinitos("./simulaciones/VOL-200_e-13mol_Z3_6m_EXPONENCIAL");
    //insertarGrid(double xInicio, double R, double Nx, int numParticionesIntegrales, bool dominioEquiespaciado, bool pasoAdaptativo)
    volRK3.insertarGrid(2.32026e-15,4.46e-6,200,400,false,1.1);//SI (2.45044e-15, 4.23436e-6)
    volRK3.insertarTiempo(0,15600000,20); //6 meses
    volRK3.calcular();
     * */
    
    
    
    return 0;
}
