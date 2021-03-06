/*
 * Clase para la simulacion de soluciones de la ecuacion de Smoluchowski
 * utilizando un metodo seccional.
 * En este código primeramente se definen distintas funciones (calculo de coeficientes seccionales,
 * definicion del nucleo, etc) para luego poder plasmar y calcular el esquema numerico en la 
 * funcion "calcular()".
 * 
 * Los resultados obtenidos se guardan en un fichero de texto con un formato que se puede consultar en la 
 * funcion "calcular()".
 * 
 * En las formulas de cuadratura se aplica la formula de cuadratura de Simpson compuesto y un metodo de Runge-Kutta 
 * de orden 3 para la aproximacion temporal.
 * 
 */

/* 
 * File:   Seccional.cpp
 * Author: Salvador Fructuoso
 *
 */

#include <cstdlib>
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
#include <ctime> 
#include "header.h"
#include <vector> 


using namespace std;

/*
 * 
 * 
 */

const long double PI = std::atan(1.0)*4;
/*
long double Seccional::densidadN0(long double v, long double t) {
    auto bessel = [&] (long double s) {
        long double limInf = 0;
        long double limSup = PI;

        int nParticiones = minimoParticiones + (int) (numParticionesIntegrales * (limSup - limInf));
        if (nParticiones % 3 == 1) {
            nParticiones = nParticiones + 2;
        }
        if (nParticiones % 3 == 2) {
            nParticiones = nParticiones + 1;
        }

        long double precision = (limSup - limInf) / nParticiones;
        //int numParticiones = (int) ((limSup - limInf) / precision + 1);
        long double xRange [nParticiones + 1];
        for (int i = 0; i < nParticiones + 1; i++) {
            xRange[i] = limInf + i*precision;
        }
        long double suma = 0;
        long double sum1 = 0;
        long double sum2 = 0;
        int k = nParticiones / 3;
        for (int i = 1; i <= k; i++) {
            //sum1 = sum1 + gIni(xRange[3 * i - 2], t) + gIni(xRange[3 * i - 1], t);
            sum1 = sum1 + exp(s * cos(xRange[3 * i - 2])) * cos(xRange[3 * i - 2]) + exp(s * cos(xRange[3 * i - 1])) * cos(xRange[3 * i - 1]);
            if (i < k) {
                sum2 = sum2 + exp(s * cos(xRange[3 * i])) * cos(xRange[3 * i]);
            }
        }
        suma = (3 * precision / 8)*(exp(s * cos(xRange[0])) * cos(xRange[0]) + 3 * sum1 + 2 * sum2 + exp(s * cos(xRange[nParticiones])) * cos(xRange[nParticiones]));
        suma = suma / PI;
        return suma;
    };

    long double n;
    if (t <= 1) {
        n = exp(-(1 + t) * v) * bessel(2 * v * pow(t, 0.5)) / (v * v * pow(t, 0.5));
    } else {
        n = exp(-2 * pow(t, 0.5) * v) * bessel(2 * v * pow(t, 0.5)) / (v * v * pow(t, 0.5));
    }
    return n;
}
*/

long double Seccional::densidadN0(long double v, long double t) {
    double M00=1;
    long double n = pow(((2 * M00) / (2 + M00 * t)), 2) * exp(-2 * M00 * v / (2 + M00 * t));
    //long double n = pow(2 * PI, -0.5) * exp(-t) * pow(v, -1.5) * exp(-v * 0.5 * exp(-2 * t));
    //n = 100 * exp(-v);
    return n;
}



long double Seccional::densidad_media0(int i, long double t) {

    long double limInf = v[i];
    long double limSup = v[i + 1];

    int nParticiones = numParticionesIntegrales;
    if (nParticiones % 3 == 1) {
        nParticiones = nParticiones + 2;
    }
    if (nParticiones % 3 == 2) {
        nParticiones = nParticiones + 1;
    }

    long double precision = (limSup - limInf) / nParticiones;
    //int numParticiones = (int) ((limSup - limInf) / precision + 1);
    long double xRange [nParticiones + 1];
    for (int i = 0; i < nParticiones + 1; i++) {
        xRange[i] = limInf + i*precision;
    }
    
    //Simpson compuesto
    auto integrando = [&](long double x, long double t) {
        return densidadN0(x,t);
    };
    
    long double suma; 
    long double sum_aux;
    suma=integrando(xRange[0],t);
    sum_aux=0;
    for(int j=1;j<=(nParticiones/2-1);j++){
        sum_aux+=integrando(xRange[2*j],t);
    }
    sum_aux*=2;
    suma+=sum_aux;
    sum_aux=0;
    for(int j=1;j<=(nParticiones/2);j++){
        sum_aux+=integrando(xRange[2*j-1],t);
    }
    sum_aux*=4;
    suma+=sum_aux;
    suma+=integrando(xRange[nParticiones],t);
    suma*=precision/3;
    
    suma = suma / (limSup - limInf);
    
    return suma;    
}

/**
long double Seccional::ker(long double x, long double y) {
    long double exponente = 1. / D;
    long double ker_brow;
    long double ker_sh;
    long double ker_ds;

    ker_brow = 2 * kBolztmann * T * pow(pow(x, exponente) + pow(y, exponente), 2) / (3 * din_viscosity * pow(x*y, exponente));
    ker_sh = (4. / 3) * sh_rate * pow(a_0, 3) / pow(vi[0], 3 * exponente) * pow(pow(x, exponente) + pow(y, exponente), 3);

    long double wx = (g / (6 * PI * (x / p_0)*(a_0 / vi[0])))*(1. / p_w - 1. / p_0) * pow(x, 1 - exponente);
    long double wy = (g / (6 * PI * (y / p_0)*(a_0 / vi[0])))*(1. / p_w - 1. / p_0) * pow(y, 1 - exponente);
    ker_ds = PI * pow(a_0, 2) / pow(vi[0], 2 * exponente) * pow(pow(x, exponente) + pow(y, exponente), 2) * abs(wx - wy);

    return ker_brow + ker_sh + ker_ds;
}
 * */



long double Seccional::ker(long double x, long double y) {
    return 1.;
}


long double Seccional::fInverse_point(long double x) {
    //long double v = x; //En este caso se toma la funcion identidad
    //long double v = 4 * PI * pow(x / 2, 3.0) / 3.0;
    return x;
}

long double Seccional::f_point(long double v) {
    //long double exponente = 1.0 / 3;
    //long double x = pow(6 * v / PI, exponente); //x es el diametro de un volumen esferico
    return v;
}

/*
long double Seccional::df(long double v) {
    //long double exponente=1./3;
    //long double df=(1./3)*pow(6/PI,exponente)*pow(1/pow(0.5*(v[l+1]+v[l]),2),exponente);
    return 1;
}
 */
long double Seccional::thetaInfSup(long double limInf, long double limSup, long double valor) {
    if (valor < limSup && valor > limInf) {
        return 1;
    } else {
        return 0;
    }
}

long double Seccional::thetaInf(long double limInf, long double valor) {
    if (valor > limInf) {
        return 1;
    } else {
        return 0;
    }
}

long double Seccional::thetaSup(long double limSup, long double valor) {
    if (valor < limSup) {
        return 1;
    } else {
        return 0;
    }
}



/*
 *Coeficiente seccional C_1
 */
long double Seccional::coefC1(int i, int j, int l) {
    long double vInf = v[l];
    long double vSup = v[l + 1];
    long double limInfX = x[i];
    long double limSupX = x[i + 1];
    long double limInfY = x[j];
    long double limSupY = x[j + 1];

    //9JUN
    //int nParticionesV=minimoParticiones+numParticionesIntegrales*((int)((vSup-vInf)*m/R));
    int nParticionesX = minimoParticiones + (int) (numParticionesIntegrales * (limSupX - limInfX));
    int nParticionesY = minimoParticiones + (int) (numParticionesIntegrales * (limSupY - limInfY));

    //

    //int nParticiones = numParticionesIntegrales;

    long double precisionX = (limSupX - limInfX) / nParticionesX;
    long double precisionY = (limSupY - limInfY) / nParticionesY;
    long double xRange [nParticionesX + 1];
    for (int i = 0; i < nParticionesX + 1; i++) {
        xRange[i] = limInfX + i*precisionX;
    }
    long double yRange [nParticionesY + 1];
    for (int i = 0; i < nParticionesY + 1; i++) {
        yRange[i] = limInfY + i*precisionY;
    }
    
    auto integrando = [&](long double x, long double y) {
        return thetaInfSup(vInf, vSup, x+y)*ker(x,y);
    };
    
    long double sum_aux=0;
    //###### Bloque 1
    long double suma_bloque1=0;
    suma_bloque1 += integrando(xRange[0],yRange[0]);
    
    for(int i=1; i<= (nParticionesX/2-1);i++){
        sum_aux += integrando(xRange[2*i],yRange[0]);
    }
    suma_bloque1 += 2*sum_aux;
    sum_aux=0;
    for(int i=1; i<=(nParticionesX/2);i++){
        sum_aux += integrando(xRange[2*i-1],yRange[0]);
    }
    suma_bloque1+=4*sum_aux;
    suma_bloque1+=integrando(xRange[nParticionesX],yRange[0]);
    
    //#####Bloque 2
    long double suma_bloque2=0;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        sum_aux+=integrando(xRange[0],yRange[2*j]);
    }
    suma_bloque2+=sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        for(int i=1; i<= (nParticionesX/2-1);i++){
            sum_aux+=integrando(xRange[2*i],yRange[2*j]);
        }
    }
    suma_bloque2+=2*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        for(int i=1; i<= (nParticionesX/2);i++){
            sum_aux+=integrando(xRange[2*i-1],yRange[2*j]);
        }
    }
    suma_bloque2+=4*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        sum_aux+=integrando(xRange[nParticionesX],yRange[2*j]);
    }
    suma_bloque2+=sum_aux;
    suma_bloque2 *= 2;
    
    //## Bloque 3:
    long double suma_bloque3=0;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        sum_aux+=integrando(xRange[0],yRange[2*j-1]);
    }
    suma_bloque3+=sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        for(int i=1; i<= (nParticionesX/2-1);i++){
            sum_aux+=integrando(xRange[2*i],yRange[2*j-1]);
        }
    }
    suma_bloque3+=2*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        for(int i=1; i<= (nParticionesX/2);i++){
            sum_aux+=integrando(xRange[2*i-1],yRange[2*j-1]);
        }
    }
    suma_bloque3+=4*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        sum_aux+=integrando(xRange[nParticionesX],yRange[2*j-1]);
    }
    suma_bloque3+=sum_aux;
    suma_bloque3*=4;
    
    //## Bloque 4:
    long double suma_bloque4=0;
    suma_bloque4=integrando(xRange[0],yRange[nParticionesY]);
    sum_aux=0;
    for(int i=1; i<=(nParticionesX/2-1);i++){
        sum_aux+=integrando(xRange[2*i],yRange[nParticionesY]);
    }
    suma_bloque4+=2*sum_aux;
    sum_aux=0;
    for(int i=1; i<=(nParticionesX/2);i++){
        sum_aux+=integrando(xRange[2*i-1],yRange[nParticionesY]);
    }
    suma_bloque4+=4*sum_aux;
    sum_aux=0;
    suma_bloque4+=integrando(xRange[nParticionesX],yRange[nParticionesY]);
    
    ///////
    
    long double suma;
    
    suma=precisionX*precisionY/9*(suma_bloque1+suma_bloque2+suma_bloque3+suma_bloque4);
    
    suma = suma / ((limSupX - limInfX)*(limSupY - limInfY));

    return suma;
    
}


/*
 *Coeficiente seccional C_2
 */
long double Seccional::coefC2(int i, int l) {
    long double vInf = v[l + 1];
    long double vSup = v[l + 1];
    long double limInfX = x[i];
    long double limSupX = x[i + 1];
    long double limInfY = x[l];
    long double limSupY = x[l + 1];

    //9JUN
    //int nParticionesV=minimoParticiones+numParticionesIntegrales*((int)((vSup-vInf)*m/R));
    int nParticionesX = minimoParticiones + (int) (numParticionesIntegrales * (limSupX - limInfX));
    int nParticionesY = minimoParticiones + (int) (numParticionesIntegrales * (limSupY - limInfY));

    //

    //int nParticiones = numParticionesIntegrales;

    long double precisionX = (limSupX - limInfX) / nParticionesX;
    long double precisionY = (limSupY - limInfY) / nParticionesY;
    long double xRange [nParticionesX + 1];
    for (int i = 0; i < nParticionesX + 1; i++) {
        xRange[i] = limInfX + i*precisionX;
    }
    long double yRange [nParticionesY + 1];
    for (int i = 0; i < nParticionesY + 1; i++) {
        yRange[i] = limInfY + i*precisionY;
    }

    auto integrando = [&](long double x, long double y) {
        return thetaInf(vInf, x+y)*ker(x,y);
    };

    
    long double sum_aux=0;
    //###### Bloque 1
    long double suma_bloque1=0;
    suma_bloque1 += integrando(xRange[0],yRange[0]);
    
    for(int i=1; i<= (nParticionesX/2-1);i++){
        sum_aux += integrando(xRange[2*i],yRange[0]);
    }
    suma_bloque1 += 2*sum_aux;
    sum_aux=0;
    for(int i=1; i<=(nParticionesX/2);i++){
        sum_aux += integrando(xRange[2*i-1],yRange[0]);
    }
    suma_bloque1+=4*sum_aux;
    suma_bloque1+=integrando(xRange[nParticionesX],yRange[0]);
    
    //#####Bloque 2
    long double suma_bloque2=0;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        sum_aux+=integrando(xRange[0],yRange[2*j]);
    }
    suma_bloque2+=sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        for(int i=1; i<= (nParticionesX/2-1);i++){
            sum_aux+=integrando(xRange[2*i],yRange[2*j]);
        }
    }
    suma_bloque2+=2*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        for(int i=1; i<= (nParticionesX/2);i++){
            sum_aux+=integrando(xRange[2*i-1],yRange[2*j]);
        }
    }
    suma_bloque2+=4*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        sum_aux+=integrando(xRange[nParticionesX],yRange[2*j]);
    }
    suma_bloque2+=sum_aux;
    suma_bloque2 *= 2;
    
    //## Bloque 3:
    long double suma_bloque3=0;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        sum_aux+=integrando(xRange[0],yRange[2*j-1]);
    }
    suma_bloque3+=sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        for(int i=1; i<= (nParticionesX/2-1);i++){
            sum_aux+=integrando(xRange[2*i],yRange[2*j-1]);
        }
    }
    suma_bloque3+=2*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        for(int i=1; i<= (nParticionesX/2);i++){
            sum_aux+=integrando(xRange[2*i-1],yRange[2*j-1]);
        }
    }
    suma_bloque3+=4*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        sum_aux+=integrando(xRange[nParticionesX],yRange[2*j-1]);
    }
    suma_bloque3+=sum_aux;
    suma_bloque3*=4;
    
    //## Bloque 4:
    long double suma_bloque4=0;
    suma_bloque4=integrando(xRange[0],yRange[nParticionesY]);
    sum_aux=0;
    for(int i=1; i<=(nParticionesX/2-1);i++){
        sum_aux+=integrando(xRange[2*i],yRange[nParticionesY]);
    }
    suma_bloque4+=2*sum_aux;
    sum_aux=0;
    for(int i=1; i<=(nParticionesX/2);i++){
        sum_aux+=integrando(xRange[2*i-1],yRange[nParticionesY]);
    }
    suma_bloque4+=4*sum_aux;
    sum_aux=0;
    suma_bloque4+=integrando(xRange[nParticionesX],yRange[nParticionesY]);
    
    ///////
    
    long double suma;
    
    suma=precisionX*precisionY/9*(suma_bloque1+suma_bloque2+suma_bloque3+suma_bloque4);
    
    suma = suma / ((limSupX - limInfX)*(limSupY - limInfY));

    return suma;
    
}


/*
 *Coeficiente seccional C_3
 */
long double Seccional::coefC3(int l) {
    long double vInf = v[l + 1];
    long double vSup = v[l + 1];
    long double limInfX = x[l];
    long double limSupX = x[l + 1];
    long double limInfY = x[l];
    long double limSupY = x[l + 1];

    //9JUN
    //int nParticionesV=minimoParticiones+numParticionesIntegrales*((int)((vSup-vInf)*m/R));
    int nParticionesX = minimoParticiones + (int) (numParticionesIntegrales * (limSupX - limInfX));
    int nParticionesY = minimoParticiones + (int) (numParticionesIntegrales * (limSupY - limInfY));

    //

    //int nParticiones = numParticionesIntegrales;

    long double precisionX = (limSupX - limInfX) / nParticionesX;
    long double precisionY = (limSupY - limInfY) / nParticionesY;
    long double xRange [nParticionesX + 1];
    for (int i = 0; i < nParticionesX + 1; i++) {
        xRange[i] = limInfX + i*precisionX;
    }
    long double yRange [nParticionesY + 1];
    for (int i = 0; i < nParticionesY + 1; i++) {
        yRange[i] = limInfY + i*precisionY;
    }

    auto integrando = [&](long double x, long double y) {
        return (thetaInf(vInf, x+y)*2 + thetaSup(vSup, x+y))*ker(x,y);
    };
  
    long double sum_aux=0;
    //###### Bloque 1
    long double suma_bloque1=0;
    suma_bloque1 += integrando(xRange[0],yRange[0]);
    
    for(int i=1; i<= (nParticionesX/2-1);i++){
        sum_aux += integrando(xRange[2*i],yRange[0]);
    }
    suma_bloque1 += 2*sum_aux;
    sum_aux=0;
    for(int i=1; i<=(nParticionesX/2);i++){
        sum_aux += integrando(xRange[2*i-1],yRange[0]);
    }
    suma_bloque1+=4*sum_aux;
    suma_bloque1+=integrando(xRange[nParticionesX],yRange[0]);
    
    //#####Bloque 2
    long double suma_bloque2=0;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        sum_aux+=integrando(xRange[0],yRange[2*j]);
    }
    suma_bloque2+=sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        for(int i=1; i<= (nParticionesX/2-1);i++){
            sum_aux+=integrando(xRange[2*i],yRange[2*j]);
        }
    }
    suma_bloque2+=2*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        for(int i=1; i<= (nParticionesX/2);i++){
            sum_aux+=integrando(xRange[2*i-1],yRange[2*j]);
        }
    }
    suma_bloque2+=4*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        sum_aux+=integrando(xRange[nParticionesX],yRange[2*j]);
    }
    suma_bloque2+=sum_aux;
    suma_bloque2 *= 2;
    
    //## Bloque 3:
    long double suma_bloque3=0;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        sum_aux+=integrando(xRange[0],yRange[2*j-1]);
    }
    suma_bloque3+=sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        for(int i=1; i<= (nParticionesX/2-1);i++){
            sum_aux+=integrando(xRange[2*i],yRange[2*j-1]);
        }
    }
    suma_bloque3+=2*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        for(int i=1; i<= (nParticionesX/2);i++){
            sum_aux+=integrando(xRange[2*i-1],yRange[2*j-1]);
        }
    }
    suma_bloque3+=4*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        sum_aux+=integrando(xRange[nParticionesX],yRange[2*j-1]);
    }
    suma_bloque3+=sum_aux;
    suma_bloque3*=4;
    
    //## Bloque 4:
    long double suma_bloque4=0;
    suma_bloque4=integrando(xRange[0],yRange[nParticionesY]);
    sum_aux=0;
    for(int i=1; i<=(nParticionesX/2-1);i++){
        sum_aux+=integrando(xRange[2*i],yRange[nParticionesY]);
    }
    suma_bloque4+=2*sum_aux;
    sum_aux=0;
    for(int i=1; i<=(nParticionesX/2);i++){
        sum_aux+=integrando(xRange[2*i-1],yRange[nParticionesY]);
    }
    suma_bloque4+=4*sum_aux;
    sum_aux=0;
    suma_bloque4+=integrando(xRange[nParticionesX],yRange[nParticionesY]);
    
    ///////
    
    long double suma;
    
    suma=precisionX*precisionY/9*(suma_bloque1+suma_bloque2+suma_bloque3+suma_bloque4);
    
    suma = suma / ((limSupX - limInfX)*(limSupY - limInfY));

    return suma;
}

/*
 *Coeficiente seccional C_4
 */
long double Seccional::coefC4(int i, int l) {
    long double limInfX = x[i];
    long double limSupX = x[i + 1];
    long double limInfY = x[l];
    long double limSupY = x[l + 1];

    //9JUN
    //int nParticionesV=minimoParticiones+numParticionesIntegrales*((int)((vSup-vInf)*m/R));
    int nParticionesX = minimoParticiones + (int) (numParticionesIntegrales * (limSupX - limInfX));
    int nParticionesY = minimoParticiones + (int) (numParticionesIntegrales * (limSupY - limInfY));

    //

    //int nParticiones = numParticionesIntegrales;

    long double precisionX = (limSupX - limInfX) / nParticionesX;
    long double precisionY = (limSupY - limInfY) / nParticionesY;
    long double xRange [nParticionesX + 1];
    for (int i = 0; i < nParticionesX + 1; i++) {
        xRange[i] = limInfX + i*precisionX;
    }
    long double yRange [nParticionesY + 1];
    for (int i = 0; i < nParticionesY + 1; i++) {
        yRange[i] = limInfY + i*precisionY;
    }

    auto integrando = [&](long double x, long double y) {
        return ker(x,y);
    };
    
    long double sum_aux=0;
    //###### Bloque 1
    long double suma_bloque1=0;
    suma_bloque1 += integrando(xRange[0],yRange[0]);
    
    for(int i=1; i<= (int)(nParticionesX/2-1);i++){
        sum_aux += integrando(xRange[2*i],yRange[0]);
    }
    suma_bloque1 += 2*sum_aux;
    sum_aux=0;
    for(int i=1; i<=(int)(nParticionesX/2);i++){
        sum_aux += integrando(xRange[2*i-1],yRange[0]);
    }
    suma_bloque1+=4*sum_aux;
    suma_bloque1+=integrando(xRange[nParticionesX],yRange[0]);
    
    //#####Bloque 2
    long double suma_bloque2=0;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        sum_aux+=integrando(xRange[0],yRange[2*j]);
    }
    suma_bloque2+=sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        for(int i=1; i<= (nParticionesX/2-1);i++){
            sum_aux+=integrando(xRange[2*i],yRange[2*j]);
        }
    }
    suma_bloque2+=2*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        for(int i=1; i<= (nParticionesX/2);i++){
            sum_aux+=integrando(xRange[2*i-1],yRange[2*j]);
        }
    }
    suma_bloque2+=4*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2-1);j++){
        sum_aux+=integrando(xRange[nParticionesX],yRange[2*j]);
    }
    suma_bloque2+=sum_aux;
    suma_bloque2 *= 2;
    
    //## Bloque 3:
    long double suma_bloque3=0;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        sum_aux+=integrando(xRange[0],yRange[2*j-1]);
    }
    suma_bloque3+=sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        for(int i=1; i<= (nParticionesX/2-1);i++){
            sum_aux+=integrando(xRange[2*i],yRange[2*j-1]);
        }
    }
    suma_bloque3+=2*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        for(int i=1; i<= (nParticionesX/2);i++){
            sum_aux+=integrando(xRange[2*i-1],yRange[2*j-1]);
        }
    }
    suma_bloque3+=4*sum_aux;
    sum_aux=0;
    for(int j=1; j<=(nParticionesY/2);j++){
        sum_aux+=integrando(xRange[nParticionesX],yRange[2*j-1]);
    }
    suma_bloque3+=sum_aux;
    suma_bloque3*=4;
    
    //## Bloque 4:
    long double suma_bloque4=0;
    suma_bloque4=integrando(xRange[0],yRange[nParticionesY]);
    sum_aux=0;
    for(int i=1; i<=(nParticionesX/2-1);i++){
        sum_aux+=integrando(xRange[2*i],yRange[nParticionesY]);
    }
    suma_bloque4+=2*sum_aux;
    sum_aux=0;
    for(int i=1; i<=(nParticionesX/2);i++){
        sum_aux+=integrando(xRange[2*i-1],yRange[nParticionesY]);
    }
    suma_bloque4+=4*sum_aux;
    sum_aux=0;
    suma_bloque4+=integrando(xRange[nParticionesX],yRange[nParticionesY]);
    
    ///////
    
    long double suma;
    
    suma=precisionX*precisionY/9*(suma_bloque1+suma_bloque2+suma_bloque3+suma_bloque4);
    
    suma = suma / ((limSupX - limInfX)*(limSupY - limInfY));

    return suma;
}

/*
 *Constructor de la clase Seccional
 */
Seccional::Seccional(const char * nombreTxtResultado) {
    this->nombreTxtResultado = nombreTxtResultado;
}

void Seccional::insertarGrid(long double v0, long double R, int m, int numParticionesIntegrales, bool dominioEquiespaciado, bool pasoAdaptativo) {
    this->m = m;
    this->v0 = v0;
    this->R = R;
    this->numParticionesIntegrales = numParticionesIntegrales;
    this->dominioEquiespaciado = dominioEquiespaciado;
    this->pasoAdaptativo=pasoAdaptativo;
}

void Seccional::insertarTiempo(long double t0, long double tFinal, long double incTiempo) {
    this->t0 = t0;
    this->tFinal = tFinal;
    this->incTiempo = incTiempo;
    this->numIntervalosTiempo = (int) ((tFinal - t0) / incTiempo) + 1;
}


/*
 * dQ/dt: se invoca en el metodo de Runge-Kutta posteriormente
 */
long double Seccional::dQ(int l, long double add) {
    long double QAux [m];
    for (int i = 0; i < m; i++) {
        QAux[i] = Q[0][i] + add;
    }

    long double suma1 = 0;
    long double suma2 = 0;
    long double suma3 = 0;
    long double suma4 = 0;
    long double k = 0;
    if (l > 1) {
        for (int i = 1; i < l; i++) {
            for (int j = 1; j < l; j++) {
                //17May : cambio l-1 por l-2
                suma1 = suma1 + C1[i - 1][j - 1][l - 2] * QAux[j - 1] * QAux[i - 1];
            }
        }
        suma1 = 0.5 * suma1;

        for (int i = 1; i < l; i++) {
            //!!
            suma2 = suma2 + C2[i - 1][l - 2] * QAux[i - 1];

        }
        suma2 = suma2 * QAux[l - 1];
    }

    suma3 = 0.5 * C3[l - 1] * pow(QAux[l - 1], 2);

    if (l < m) {
        //cambio el i=2 por l+1
        for (int i = l + 1; i < m + 1; i++) {
            //17May : cambio l-1 por l-2
            suma4 = suma4 + C4[i - 2][l - 1] * QAux[i - 1];

        }
        suma4 = suma4 * QAux[l - 1];
    }
    k = suma1 - suma2 - suma3 - suma4;

    return k;

}



/*
 * Mediante esta funcion se invocan a las funciones definidas anteriormente y se 
 * ejecutan las simulaciones
 */
void Seccional::calcular() {
    time_t start, end;
    time(&start);//Inicio del calculo del tiempo de ejecucion

    alfa = 1.0;
    gamma = 0;

    cout << numIntervalosTiempo << endl;

    v = new long double [m + 1];
    v[m] = v0 + R;

    long double vv [m + 1];
    vv[m] = v[m];
    if (dominioEquiespaciado == false) {
        long double k_esc = pow((v0 + R) / v0, 1. / (m));
        cout << "K escala=";
        cout << k_esc << endl;
        for (int i = 0; i < m; i++) {
            v[m - i - 1] = v[m - i] / k_esc;
        }
    } else {
        long double h = R / m;
        for (int i = 0; i < m + 1; i++) {
            v[i] = v0 + h*i;
        }
    }

    //////////////////

    vi = new long double [m];
    for (int i = 0; i < m; i++) {
        vi[i] = 0.5 * (v[i] + v[i + 1]);
    }

    ////////////////////////////

    //cout << "X------" << endl;
    x = new long double [m + 1];
    for (int i = 0; i < m + 1; i++) {
        x[i] = f_point(v[i]);
    }

    Q = new long double * [2];
    for (int i = 0; i < 2; i++) {
        Q[i] = new long double [m];
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < m; j++) {
            Q[i][j] = 0;
        }
    }

    factorQ = new long double [m];
    cout << "Q:----" << endl;
    for (int l = 0; l < m; l++) {
        factorQ[l] = alfa * pow(vi[l], gamma) * (f_point(v[l + 1]) - f_point(v[l]));
        Q[0][l] = factorQ[l]*densidad_media0(l, t0);
        //cout << Q[0][l] << endl;
    }

    
    ////////////////////////////
    ///Generacion de un fichero de texto en el cual guardar los resultados 
    ///de la simulacion
    ofstream fichero;
    char buf[120];
    snprintf(buf, sizeof (buf), "%s.txt", this->nombreTxtResultado);
    fichero.open(buf);
    fichero.precision(20);
    //Escritura de variables en el fichero de texto
    fichero << m;
    fichero << ";";
    fichero << v0;
    fichero << ";";
    fichero << R;
    fichero << ";";
    fichero << t0;
    fichero << ";";
    fichero << tFinal;
    fichero << ";";
    fichero << numIntervalosTiempo;
    fichero << "\n";

    for (int i = 0; i < m; i++) {
        fichero << v[i];
        fichero << ";";
    }
    fichero << to_string(v[m]) + "\n";

    for (int j = 0; j < m; j++) {
        fichero << Q[0][j]/factorQ[j];//Escritura de la distribucion inicial del problema
        if (j == m - 1) {
            fichero << "\n";
        } else {
            fichero << ";";
        }
    }

    /////////
    ////////////////Coeficientes seccionales
    C1 = new long double **[m - 1];
    C2 = new long double *[m - 1];
    C3 = new long double [m];
    C3[m - 1] = 0;
    C4 = new long double *[m - 1];
    for (int i = 0; i < m - 1; i++) {
        C3[i] = 0;
        C2[i] = new long double [m - 1];
        C4[i] = new long double [m - 1];
        C1[i] = new long double *[m - 1];
        for (int j = 0; j < m - 1; j++) {
            C2[i][j] = 0;
            C4[i][j] = 0;
            C1[i][j] = new long double [m - 1];
            for (int k = 0; k < m - 1; k++) {
                C1[i][j][k] = 0;
            }
        }
    }


    //C1///////////////////////////////////
    long double porcentaje;
    long double porcentaje2 = 0;
    cout << "C1:----------" << endl;

    for (int l = 2; l < m + 1; l++) {
        for (int i = 1; i < l; i++) {
            for (int j = 1; j < i + 1; j++) {
                //
                C1[i - 1][j - 1][l - 2] = coefC1(i - 1, j - 1, l - 1);

                if (j != i) {
                    C1[j - 1][i - 1][l - 2] = C1[i - 1][j - 1][l - 2];
                }
                if (isnan(C1[j - 1][i - 1][l - 2]) == true) {
                    cout << "Error en C1" << endl;
                    std::exit(EXIT_FAILURE);
                }
            }
        }
        porcentaje = (int) (100 * l / m);
        if (porcentaje > porcentaje2) {
            porcentaje2 = porcentaje;
            cout << porcentaje2;
            cout << " %" << endl;
        }

    }



    ///C2///////////////////////////////////////////
    porcentaje2 = 0;
    cout << "\n";
    cout << "C2:----------" << endl;

    for (int i = 1; i < m; i++) {
        for (int l = i + 1; l < m + 1; l++) {
            C2[i - 1][l - 2] = coefC2(i - 1, l - 1);
            if (isnan(C2[i - 1][l - 2]) == true) {
                cout << "Error en C2" << endl;
                std::exit(EXIT_FAILURE);
            }
        }
        porcentaje = (int) (100 * i / (m - 1));
        if (porcentaje > porcentaje2) {
            porcentaje2 = porcentaje;
            cout << porcentaje2;
            cout << " %" << endl;
        }
    }


    //C3/////////////////////////////////////////////////////////
    porcentaje2 = 0;
    cout << "\n";
    cout << "C3:----------" << endl;

    for (int l = 1; l < m + 1; l++) {
        C3[l - 1] = coefC3(l - 1);
        if (isnan(C3[l - 1]) == true) {
            cout << "Error en C3" << endl;
            std::exit(EXIT_FAILURE);
        }
        porcentaje = (int) (100 * l / m);
        if (porcentaje > porcentaje2) {
            porcentaje2 = porcentaje;
            cout << porcentaje2;
            cout << " %" << endl;
        }
    }

    //C4/////////////////////////////////////////////
    porcentaje2 = 0;
    cout << "\n";
    cout << "C4:----------" << endl;

    for (int i = 2; i < m + 1; i++) {
        for (int l = 1; l < i; l++) {
            C4[i - 2][l - 1] = coefC4(i - 1, l - 1);
            if (isnan(C4[i - 2][l - 1]) == true) {
                cout << "Error en C4" << endl;
                std::exit(EXIT_FAILURE);
            }

            porcentaje = (int) (100 * i / m);
            if (porcentaje > porcentaje2) {
                porcentaje2 = porcentaje;
                cout << porcentaje2;
                cout << " %" << endl;
            }
        }
    }
    //}

    ///////////Inicio de la implementación del esquema numerico
    ///////////////////////

    int z = 1;
       //Pendientes del metodo de RK3
    long double k1 = 0;
    long double k2 = 0;
    long double k3 = 0;


    long double incTime = 0;//Variable que define \Delta t
    bool tiempoFino = false;//Variable booleana para comprobar si el 
                      //salto temporal cumple unas condiciones impuestas
                      //Se usa para el caso de un salto temporal variable

    long double vecDomTiempoAux = 0;
    vector<long double> domTiempo;//Vector donde guardar cada uno de los tiempos usados en las iteraciones
    domTiempo.push_back(t0);

    long double error = 0; //Error de aproximacion en cada salto temporal
                           //respecto a una solución exacta conocida.
    vector<long double> errores;//Aray donde guardar los errores de cada salto
                                //temporal
    errores.push_back(0);

    int contador = 0;
    int ll;
    long double value;

    vector<long double> domTiempo_printed;//Vector donde guardar los instantes de
                                          //tiempo de los puntos valores progresados
                                          //al fichero de texto.
    int indiceErrores = 1;
    long double tiempoAcumulado = 0;
    int porcentajeRK = 0;//Variable para mostrar en pantalla el porcentaje 
    //de la simulación completado.

    cout << ">>> RK3 >>>" << endl;

    while (domTiempo[z - 1] < tFinal) {

        incTime = incTiempo;//Tiempo fijo \Delta t definido en la clase main()
        value = 0;//Variable auxiar para evaluar el valor la diferencia entre 
                  //n(x,t^n) y n(x,t^{n+1})
        tiempoFino = false;//Boleano para ver si el mallado temporal variable es apto.
        
        
        long double vFinal;
        int contadorVFinalNegativo = 0;
        int binarioParidad = 2;
        int numValoresNegativos = 0;
        
        
        ///////////////////////////////////
        ///Bucle para evaluar el salto temporal en cada instante
        //Reduce \Delta t si no cumple las condiciones del bucle:
        //      1. |n(x,t^{n+1})-n(x,t^n})|<=0.0005 en este caso
        //      2. La aproximación tenga más de umbral de puntos definido que
        //         sean negativos. En este caso hay un límite de subdivisión de 
        //         \Delta t; si se alcaza ese límite el programa se para y muestra
        //         un Warning por pantalla.
        //
        //En caso de querer un dominio temporal fijo, basta con imponer una condición
        //falsa el bucle para que no se ejecute.
        while (tiempoFino == false && pasoAdaptativo== true) {
            contadorVFinalNegativo = 0;
            ll = 1;
            double incLL;
            if (binarioParidad % 2 == 0) {
                incLL = 1;
            } else {
                incLL = 4;
            }
            while (ll < m + 1) {
                k1 = dQ(ll, 0);
                k2 = dQ(ll, incTime * k1 / 3.0);
                k3 = dQ(ll, 2 * incTime * k2 / 3.0);

                value = incTime * (k1 + 3 * k3) / 4.0;

                vFinal = Q[0][ll - 1] + value;

                if (abs(value) > 0.0005) {
                    incTime = 0.5 * incTime;
                    break;
                } else if (vFinal < 0 && incTime > 1e-8) {
                    incTime = 0.5 * incTime;
                    break;
                } else if (vFinal < 0 && incTime < 1e-8) {
                    numValoresNegativos++;
                    int umbralMaximoValoresNegativos = 5;
                    if (vFinal > 1e-10 && numValoresNegativos < umbralMaximoValoresNegativos) {
                        cout << "[WARNING]: Valor negativo de la aproximación (;";
                        cout << vFinal;
                        cout << ") en el instante de tiempo t=";
                        cout << domTiempo[z - 1] << endl;
                    } else if (vFinal > 1e-10 || numValoresNegativos >= umbralMaximoValoresNegativos) {
                        cout << "[ERROR] Aproximación cancelada: Error en aproximación. Valores negativos." << endl;
                        std::exit(EXIT_FAILURE);
                    }
                }
                if (ll + incLL >= m + 1) {
                    tiempoFino = true;
                }
                ll = ll + incLL;
            }
        }
        binarioParidad++;
        //
        //Fin bucle de adaptación del dominio temporal.
        //////////////////////
        

        vecDomTiempoAux = domTiempo[z - 1] + incTime;
        domTiempo.resize(z + 1, vecDomTiempoAux);

        //Aplicacion del metodo de Runge-Kutta.
        for (int l = 1; l < m + 1; l++) {
            k1 = dQ(l, 0);
            k2 = dQ(l, incTime * k1 / 3.0);
            k3 = dQ(l, 2 * incTime * k2 / 3.0);

            Q[1][l - 1] = Q[0][l - 1] + incTime * (k1 + 3 * k3) / 4.0;

            if (isnan(Q[1][l - 1]) == true) {
                cout << "[ERROR] Aproximación cancelada: Nan." << endl;
                std::exit(EXIT_FAILURE);
            }
        }

        for (int i = 0; i < m; i++) {
            Q[0][i] = Q[1][i];
        }
        tiempoAcumulado = tiempoAcumulado + incTime;

        if (tiempoAcumulado >= 0.005) {//De esta forma se guarda en el fichero de texto
            //puntos con una distancia temporal mínima, para asi evitar que el fichero
            //generado sea de mucho tamaño.
            error = 0;
            domTiempo_printed.resize(indiceErrores, vecDomTiempoAux);
            indiceErrores++;
            tiempoAcumulado = 0;
            long double valorEscribir = 0;
            for (int j = 0; j < m; j++) {
                valorEscribir = Q[1][j]/factorQ[j];//Escritura de soluciones en el fichero de texto

                fichero << valorEscribir;
                if (j < (m - 1)) {
                    fichero << ";";
                }
                if (j == m - 1) {

                    fichero << "\n";
                }
                
                //Introduccion de los errores de calculo    
                //error = error + (v[j + 1] - v[j]) * abs(valorEscribir - 0.5 * (v[j + 1] + v[j]) * densidadN0(0.5 * (v[j + 1] + v[j]), domTiempo[z]));
                error = error + (v[j + 1] - v[j]) * vi[j] * abs(valorEscribir - densidadN0(vi[j], domTiempo[z]));
            }
            errores.resize(indiceErrores, error);
        }
        
        //
        //Mostrado por pantalla del porcentaje completado de la simulación
        if (((int) (100 * (domTiempo[z - 1] / tFinal))) > porcentajeRK) {
            porcentajeRK = ((int) (100 * (domTiempo[z - 1] / tFinal)));
            cout << to_string((int) (100 * (domTiempo[z - 1] / tFinal))) + " % ; t=";
            cout << domTiempo[z - 1] << endl;
        }
        z = z + 1;
    }

    //Guardado en el fichero de texto el array de errores respecto a soluciones conocidas.No se guardan todos, sino los que 
    //cumplan un mínimo \Delta t en función de la variable "tiempoAcumulado" para así evitar un 
    //exceso de tamaño del fichero.
    for (int i = 1; i < indiceErrores; i++) {
        fichero << errores[i];
        if (i + 1 < indiceErrores) {
            fichero << ";";
        }
    }
    ////////////
    fichero << "\n";
    
    //Guardado en el fichero de texto los instantes temporales. No se guardan todos, sino los que 
    //cumplan un mínimo \Delta t en función de la variable "tiempoAcumulado" para así evitar un 
    //exceso de tamaño del fichero.
    for (int i = 1; i < indiceErrores - 1; i++) {
        fichero << domTiempo_printed[i];
        fichero << ";";
    }
    //////////
    time(&end);//Fin del cálculo del tiempo de ejecución del programa

    double tiempoEjecucion = double(end - start);
    fichero << "\n";
    fichero << tiempoEjecucion;
    fichero.close();
    cout << "[INFO]: Simulación finalizada correctamente." << endl;
}