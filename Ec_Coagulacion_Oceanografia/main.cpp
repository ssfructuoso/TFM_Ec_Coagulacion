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
    cout<< "Seccional:---" << endl;
    Seccional secRK3 = Seccional("./simulaciones/seccional/2s-seccional_RK3_PROD_trapecios_compuesto_50_OK");
    //insertarGrid(double v0, double R, int m, int numParticionesIntegrales, bool dominioEquiespaciado,double logPuntosGrid)
    secRK3.insertarGrid(1e-4,1000,50,400,false,0.01);
    //insertarTiempo(double t0, double tFinal, double incTiempo)
    secRK3.insertarTiempo(0.0001, 2, 1e-2);
    secRK3.calcular();
   
    cout<< "Volumenes finitos:---" << endl;
    VolFinitos vol1=VolFinitos("./simulaciones/vol_finitos/2s-vol_finitos_RK3_PROD_trapecios_compuesto_50");
    //insertarGrid(double xInicio, double R, double Nx, int numParticionesIntegrales, bool dominioEquiespaciado, double kdominio)
    vol1.insertarGrid(1e-4,500,400,400,false,1.1);
    vol1.insertarTiempo(0.001,2,1e-2);
    vol1.calcular();
     * */
        
    //Testeo de los dos metodos para Oceanografia
    
    
    cout<< "Seccional:---" << endl;
    Seccional secRK3 = Seccional("./simulaciones/v5/SECC_test-100_e-13mol_Z3_6m_EXPONENCIAL_1m");
    //insertarGrid(double v0, double R, int m, int numParticionesIntegrales, bool dominioEquiespaciado,double logPuntosGrid)
    secRK3.insertarGrid(1.9268e-15,5.13925e-6,50,400,false,1.1);
    //insertarTiempo(double t0, double tFinal, double incTiempo)
    secRK3.insertarTiempo(0,2600000,20);
    secRK3.calcular();
    
    
    /*
    cout<< "Volumenes finitos:---" << endl;
    VolFinitos vol1=VolFinitos("./simulaciones/v5/VOL_TESTs_e-13mol_Z3_6m_EXPONENCIAL_1m");
    //insertarGrid(double xInicio, double R, double Nx, int numParticionesIntegrales, bool dominioEquiespaciado, double kdominio)
    vol1.insertarGrid(1.9268e-15,5.13925e-6,49,400,false,1.1);//SI (2.45044e-15, 4.23436e-6)
    //vol1.insertarGrid(0.936e-16,247.17,99,400,false,1.1); //Test de exp-x con el seccional
    vol1.insertarTiempo(0,2600000,20); //2 semanas
    vol1.calcular();
     * */
    
    
    
    //FIn testeo
    
    /*
    cout<< "Volumenes finitos:---" << endl;
    VolFinitos vol1=VolFinitos("./simulaciones/v5/VOL-100_e-13mol_Z3_6m_EXPONENCIAL");
    //insertarGrid(double xInicio, double R, double Nx, int numParticionesIntegrales, bool dominioEquiespaciado, double kdominio)
    //vol1.insertarGrid(1.084e-16,222,250,400,false,1.1); // 1nm a 1.2mm (en centigramos)
    //vol1.insertarGrid(5.8e-8,111.5,200,400,false,1.1); e-5 g
    vol1.insertarGrid(2.32026e-15,4.46e-6,200,400,false,1.1);//SI (2.45044e-15, 4.23436e-6)
    //vol1.insertarGrid(0.936e-16,247.17,99,400,false,1.1); //Test de exp-x con el seccional
    vol1.insertarTiempo(0,15600000,20); //2 meses
    vol1.calcular();
     * */
    
    
    
    return 0;
}
