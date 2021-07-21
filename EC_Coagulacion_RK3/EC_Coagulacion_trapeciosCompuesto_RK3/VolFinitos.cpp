/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   volFinitos.cpp
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

#include <iomanip>  

#include "header.h"
#include <vector>
using namespace std;

const long double PI = std::atan(1.0)*4;


long double VolFinitos::fIni(long double x, long double t) {
    auto bessel = [&] (long double s) {
        long double limInf = 0;
        long double limSup = PI;

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
        //Trapecios compuesto
        auto integrando = [&](long double theta, long double ss) {
            return exp(ss* cos(theta))*cos(theta);
        };

        long double suma;
        long double sum_aux = 0;
        suma = integrando(xRange[0], s);
        for (int j = 1; j <= nParticiones - 1; j++) {
            sum_aux += integrando(xRange[j], s);
        }
        suma += 2 * sum_aux;
        suma += integrando(xRange[nParticiones], s);
        suma *= precision / 2;
        suma = suma / PI;
        return suma;
    };

    long double n;
    if (t <= 1) {
        n = exp(-(1 + t) * x) * bessel(2 * x * pow(t, 0.5)) / (x * x * pow(t, 0.5));
    } else {
        n = exp(-2 * pow(t, 0.5) * x) * bessel(2 * x * pow(t, 0.5)) / (x * x * pow(t, 0.5));
    }
    return n;
}


/*
 * 
 */


long double VolFinitos::f_media0(int i, long double t) {
    long double dom_Lambda[2];
    dom_Lambda [0] = x12[i];
    dom_Lambda [1] = x12[i + 1];

    long double limInf = dom_Lambda[0];
    long double limSup = dom_Lambda[1];

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
    
    //Simpson compuesto
    auto integrando = [&](long double x, long double t) {
        return fIni(x,t);
    };
    
    long double suma; 
    long double sum_aux=0;
    suma=integrando(xRange[0],t);
    for(int j=1;j<=nParticiones-1;j++){
        sum_aux+=integrando(xRange[j],t);
    }
    suma+=2*sum_aux;
    suma+=integrando(xRange[nParticiones],t);
    suma*=precision/2;
    suma = suma / (limSup-limInf);
    
    return suma;

}

long double VolFinitos::gIni(long double x, long double t) {
    //long double g = x * fIni(x, t);
    long double g = x * fIni(x, t);
    return g;
}

long double VolFinitos::g_media0(int i, long double t) {
    long double dom_Lambda[2];
    dom_Lambda [0] = x12[i];
    dom_Lambda [1] = x12[i + 1];

    long double limInf = dom_Lambda[0];
    long double limSup = dom_Lambda[1];

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
    //Simpson compuesto
    auto integrando = [&](long double x, long double t) {
        return gIni(x,t);
    };
    
    
    long double suma; 
    long double sum_aux=0;
    suma=integrando(xRange[0],t);
    for(int j=1;j<=nParticiones-1;j++){
        sum_aux+=integrando(xRange[j],t);
    }
    suma+=2*sum_aux;
    suma+=integrando(xRange[nParticiones],t);
    suma*=precision/2;
    suma = suma / (limSup-limInf);
    
    return suma;

}
/*
long double VolFinitos::ker(long double x, long double y) {
    long double exponente = 1. / D;
    long double ker_brow;
    long double ker_sh;
    long double ker_ds;

    ker_brow = 2 * kBolztmann * T * pow(pow(x, exponente) + pow(y, exponente), 2) / (3 * din_viscosity * pow(x*y, exponente));
    ker_sh = (4. / 3) * sh_rate * pow(a_0, 3) / pow(xi[0], 3 * exponente) * pow(pow(x, exponente) + pow(y, exponente), 3);

    long double wx = (g / (6 * PI * (x / p_0)*(a_0 / xi[0])))*(1. / p_w - 1. / p_0) * pow(x, 1 - exponente);
    long double wy = (g / (6 * PI * (y / p_0)*(a_0 / xi[0])))*(1. / p_w - 1. / p_0) * pow(y, 1 - exponente);
    ker_ds = PI * pow(a_0, 2) / pow(xi[0], 2 * exponente) * pow(pow(x, exponente) + pow(y, exponente), 2) * abs(wx - wy);

    return ker_brow + ker_sh + ker_ds;
}
 * */
long double VolFinitos::ker(long double x, long double y) {
    return x*y;
}

long double VolFinitos::A(int j, long double x_k) {
    long double dom_Lambda[2];
    //if (j + 1 == Nx + 1) {
    //  j = j - 1;
    //}
    dom_Lambda [0] = x12[j];
    dom_Lambda [1] = x12[j + 1];
    long double suma = 0;

    long double limInf = dom_Lambda[0];
    long double limSup = dom_Lambda[1];

    int nParticiones = minimoParticiones + (int) (numParticionesIntegrales * (limSup - limInf));
    if (nParticiones % 3 == 1) {
        nParticiones = nParticiones + 2;
    }
    if (nParticiones % 3 == 2) {
        nParticiones = nParticiones + 1;
    }

    long double precision = (limSup - limInf) / nParticiones;
    //int numParticiones = (int) ((limSup - limInf) / precision + 1);
    long double xRange [ nParticiones + 1];
    for (int i = 0; i < nParticiones + 1; i++) {
        xRange[i] = limInf + i*precision;
    }

    //Trapecios compuesto
    auto integrando = [&](long double x, long double x_k) {
        return ker(x, x_k) / x;
    };
    suma=integrando(xRange[0],x_k);
    long double sum_aux=0;
    for(int j=1;j<=nParticiones-1;j++){
        sum_aux+=integrando(xRange[j],x_k);
    }
    suma+=2*sum_aux;
    suma+=integrando(xRange[nParticiones],x_k);
    suma*=precision/2;
 
    return suma;
}

long double VolFinitos::B(int i, long double x_k, int alfa) {
    long double limInf = x12[i + 1] - x_k;
    long double limSup = x12[alfa];
    long double suma = 0;

    if (limInf < 0) {
        return suma;
    }
    if (limSup <= limInf) {
        return suma;
    }

    int nParticiones = minimoParticiones + (int) (numParticionesIntegrales * (limSup - limInf)); //minimoParticiones + numParticionesIntegrales * ((int) (0.1 + (limSup - limInf) * Nx / R));
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
  //Trapecios compuesto
    auto integrando = [&](long double x, long double x_k) {
        return ker(x, x_k) / x;
    };
    suma=integrando(xRange[0],x_k);
    long double sum_aux=0;
    for(int j=1;j<=nParticiones-1;j++){
        sum_aux+=integrando(xRange[j],x_k);
    }
    suma+=2*sum_aux;
    suma+=integrando(xRange[nParticiones],x_k);
    suma*=precision/2;

    return suma;
}

int VolFinitos::getAlfa(int i, int k) {
    bool encontrado = false;
    long double v12 = x12[i + 1];
    long double vi = xi[k];
    long double valor = x12[i + 1] - xi[k];
    int alfa = 1;
    if (valor < x12[0]) {
        return alfa;
    }
    int j = 0;
    while (encontrado == false) {
        if (valor < x12[j + 1] && valor >= x12[j]) {
            alfa = j + 1;
            encontrado = true;
        } else {
            j = j + 1;
        }
    }
    return alfa;
}

VolFinitos::VolFinitos(const char * nombreTxtResultado) {
    this->nombreTxtResultado = nombreTxtResultado;
}

void VolFinitos::insertarGrid(long double xInicio, long double R, long double Nx, int numParticionesIntegrales, bool dominioEquiespaciado, long double kdominio) {
    this->xInicio = xInicio;
    this->R = R;
    this->Nx = Nx;
    this->numParticionesIntegrales = numParticionesIntegrales;
    this->dominioEquiespaciado = dominioEquiespaciado;
    this->kdominio = kdominio;
}

void VolFinitos::insertarTiempo(long double t0, long double tFinal, long double incTiempo) {
    this->t0 = t0;
    this->tFinal = tFinal;
    this->incTiempo = incTiempo;
    this->numIntervalosTiempo = (int) ((tFinal - t0) / incTiempo) + 1;
}

long double VolFinitos::du(int i, long double add, int kS) {
    long double suma = 0;
    long double sumaK = 0;
    long double sumaJ = 0;
    long double x_k;
    int alfaValue;
    long double kReturn;
    long double uAux [Nx + 1];

    for (int i = 0; i < Nx + 1; i++) {
        uAux[i] = u[0][i] + add;
    }
    suma = 0;
    for (int k = 0; k <= i; k++) {//!!
        sumaK = 0;
        alfaValue = alfa[i][k];
        sumaJ = 0;
        x_k = xi[k];

        for (int j = alfaValue; j <= Nx; j++) {
            sumaJ = sumaJ + Aa[j][k] * uAux[j];
        }
        sumaK = sumaJ;
        sumaK = sumaK + Bb[i][k] * uAux[alfaValue - 1];

        sumaK = sumaK * (x12[k + 1] - x12[k]) * uAux[k];
        suma = suma + sumaK;
    }

    J[kS - 1][i + 1] = suma;
    kReturn = -(J[kS - 1][i + 1] - J[kS - 1][i]) / (x12[i + 1] - x12[i]);// - (w[i] / (Z)) * uAux[i] + I[i] * xi[i];

    return kReturn;
}

/*
 * 
 */
void VolFinitos::calcular() {
    time_t start, end;
    time(&start);

    x12 = new long double [Nx + 2];
    x12[0] = xInicio;
    x12[Nx + 1] = xInicio + R;
    xi = new long double [Nx + 1];


    if (dominioEquiespaciado == false) {
        long double k_esc = pow((xInicio + R) / xInicio, 1. / (Nx + 1));
        cout << "K escala=";
        cout << k_esc << endl;
        for (int i = 0; i < Nx + 1; i++) {
            x12[Nx - i] = x12[Nx + 1 - i] / k_esc;
        }
    } else {
        long double h = R / (Nx + 1);
        for (int i = 0; i < Nx + 2; i++) {
            x12[i] = xInicio + i*h;
        }
    }

    for (int i = 0; i < Nx + 1; i++) {
        xi[i] = (x12[i] + x12[i + 1]) / 2.0;
    }


    /*
    g=7.321e10;// m/d² 
    kBolztmann=1.380649e-23; // J/K
    T=298; //Kelvin
    din_viscosity=1.0325e-8; //N*d/m²
    sh_rate=86400; // N*s/m2
    p_0=2.25e12;// ug/m³
    p_w=1.027e12;// ug/m³
    Z = 30; //m
    a_0=pow(3*xi[0]/(4*PI*p_0),1./3); //Simetria esferica para las primeras particulas
    D = 2.6; //Dimension fractal
    w = new long double [Nx+1];
    I = new long double [Nx + 1];
    for (int i = 0; i < Nx + 1; i++) {
        I[i] = 0;
        w[i]=(g/(6*PI*(xi[i]/p_0)*(a_0/xi[0])))*(1./p_w-1./p_0)*pow(xi[i],1-1./D);
    }
    I[0] = 1000;//ug/m³*d^(-1)
    //
     * */

    g = 9.8; // m/d² 
    kBolztmann = 1.380649e-23; // J/K
    T = 298; //Kelvin
    din_viscosity = 0.0008921; //N*d/m²
    sh_rate = 1; // N*s/m2
    p_0 = 2.25e12; // ug/m³
    p_w = 1.027e12; // ug/m³
    Z = 30; //m
    a_0 = pow(3 * xi[0] / (4 * PI * p_0), 1. / 3); //Simetria esferica para las primeras particulas
    D = 2.6; //Dimension fractal
    w = new long double [Nx + 1];
    I = new long double [Nx + 1];
    for (int i = 0; i < Nx + 1; i++) {
        I[i] = 0;
        w[i] = (g / (6 * PI * (xi[i] / p_0)*(a_0 / xi[0])))*(1. / p_w - 1. / p_0) * pow(xi[i], 1 - 1. / D);
    }
    I[0] = 10e5/86400;//ug/m³*d^(-1 Z=30
    //I[0] = 10e4 / 86400;
    //

    /////////////
    //int numIntervalosT = (int) ((tFinal - t0) / estTiempo) + 1;
    cout << "--------------" << endl;
    cout << "Intervalos de tiempo:";
    cout << numIntervalosTiempo << endl;

    /////////////////////////////////////////////////

    ofstream fichero;
    char buf[120];
    snprintf(buf, sizeof (buf), "%s.txt", this->nombreTxtResultado);
    fichero.open(buf);

    fichero.precision(20);
    fichero << Nx;
    fichero << ";";
    fichero << xInicio;
    fichero << ";";
    fichero << R;
    fichero << ";";
    fichero << t0;
    fichero << ";";
    fichero << tFinal;
    fichero << ";";
    fichero << numIntervalosTiempo;
    fichero << "\n";
    for (int i = 0; i < Nx + 1; i++) {
        fichero << x12[i];
        fichero << ";";
    }
    fichero << x12[Nx + 1];
    fichero << "\n";


    Aa = new long double *[Nx + 1];
    Bb = new long double *[Nx + 1];
    for (int i = 0; i < Nx + 1; i++) {
        Aa[i] = new long double [Nx + 1];
        Bb[i] = new long double [Nx + 1];
    }


    int alfaValue;
    int t = 0;
    long double x_k = 0;

    long double porcentaje = 0;
    long double porcentaje2 = 0;
    long double auxB;
    ///////////

    alfa = new int *[Nx + 1]; //!!!
    for (int i = 0; i < Nx + 1; i++) {//!!
        alfa[i] = new int [Nx + 1]; //!!!
    }

    auto valores = [&]() {
        for (int i = 0; i < Nx + 1; i++) {//!!
            for (int j = 0; j < Nx + 1; j++) {//!!
                alfa[i][j] = getAlfa(i, j);

                if (isnan(alfa[i][j]) == 1) {
                    cout << "--NAN-- en alfa" << endl;
                    std::exit(EXIT_FAILURE);
                }
                /*
                if(j+1==Nx){
                    cout<< alfa[i][j]<< endl;
                }
                else{
                    cout << alfa[i][j];
                    cout<<";";
                }
                 * */
            }
        }
        for (int i = 0; i < Nx + 1; i++) {//!!
            for (int k = 0; k < Nx + 1; k++) {//!!
                alfaValue = alfa[i][k];
                x_k = xi[k];

                Aa[i][k] = A(i, x_k);

                Bb[i][k] = B(i, x_k, alfaValue);

                if (isnan(Aa[i][k]) == 1 || isnan(Bb[i][k]) == 1) {
                    cout << "--NAN-- en valores de la matriz A o de la matriz B" << endl;
                    std::exit(EXIT_FAILURE);
                }
            }

        }
    };
    valores();
    ///////////////////////////////////////////////////////////

    u = new long double *[2];
    for (int i = 0; i < 2; i++) {
        u[i] = new long double [Nx + 1];
    }
    for (int i = 0; i < Nx + 1; i++) {// Cambio Nx+1 por Nx
        u[0][i] = g_media0(i, t0);
    }

    J = new long double *[3];

    for (int i = 0; i < 3; i++) {
        J[i] = new long double [Nx + 2];
        J[i][0] = 0;
    }
    //////////
    int z = 1;

    long double suma = 0;
    long double sumaK = 0;
    long double sumaJ = 0;

    long double k1 = 0;
    long double k2 = 0;
    long double k3 = 0;

    long double incTime = 0;
    bool tiempoFino = false;

    long double vecDomTiempoAux = 0;
    vector<long double > domTiempo;
    domTiempo.push_back(t0);

    vector<long double > domTiempo_printed;

    long double error = 0;
    int indiceErrores = 1;
    vector<long double > errores;
    errores.push_back(0);

    long double tiempoAcumulado = 0;
    int porcentajeRK = 0;

    while (domTiempo[z - 1] < tFinal) {

        incTime = incTiempo;
        long double value = 0;
        tiempoFino = false;
        int ll;

        long double vFinal;
        int binarioParidad = 2;
        int numValoresNegativos = 0;
        while (tiempoFino == false) {
            ll = 0;
            double incLL;
            if (binarioParidad % 2 == 0) {
                incLL = 1;
            } else {
                incLL = 4;
            }
            while (ll <= Nx) {
                k1 = du(ll, 0, 1);
                k2 = du(ll, incTime * k1 / 3.0, 2);
                k3 = du(ll, 2 * incTime * k2 / 3.0, 3);

                value = incTime * (k1 + 3 * k3) / 4.0;

                vFinal = u[0][ll] + value;

                if (abs(value) > 0.0001) {
                    incTime = 0.5 * incTime;
                    break;
                }
                else if (vFinal < 0 && incTime > 1e-1) {
                    incTime = 0.5 * incTime;
                    break;
                } else if (vFinal < 0 && incTime < 1e-1) {
                    numValoresNegativos++;
                    int umbralMaximoValoresNegativos = 5;
                    if (vFinal > -1e-10 && numValoresNegativos < umbralMaximoValoresNegativos && 1<0) {
                        cout << "[WARNING]: Valor negativo de la aproximación (sección ";
                        cout << ll; cout << "; ";
                        cout << vFinal;
                        cout << ") en el instante de tiempo t=";
                        cout << domTiempo[z - 1] << endl;
                    } else if (vFinal < -1e-10 || numValoresNegativos >= umbralMaximoValoresNegativos) {
                        cout << "[ERROR] Aproximación cancelada: Error en aproximación. Valores negativos." << endl;
                        std::exit(EXIT_FAILURE);
                    }
                }
                

                if (ll + incLL > Nx) {
                    tiempoFino = true;
                }
                ll = ll + incLL;
            }
        }
        binarioParidad++;

        vecDomTiempoAux = domTiempo[z - 1] + incTime;
        domTiempo.resize(z + 1, vecDomTiempoAux);

        for (int i = 0; i <= Nx; i++) {
            k1 = du(i, 0, 1);
            k2 = du(i, incTime * k1 / 3.0, 2);
            k3 = du(i, 2 * incTime * k2 / 3.0, 3);

            u[1][i] = u[0][i] + incTime * (k1 + 3 * k3) / 4.0;

            if (isnan(u[1][i]) == true) {
                cout << "[ERROR] Aproximación cancelada: Nan." << endl;
                std::exit(EXIT_FAILURE);
            }
        }


        for (int i = 0; i < Nx + 1; i++) {
            u[0][i] = u[1][i];
        }

        tiempoAcumulado = tiempoAcumulado + incTime;

        if (tiempoAcumulado >= 0.001) {
            error = 0;
            domTiempo_printed.resize(indiceErrores, vecDomTiempoAux);
            indiceErrores++;
            tiempoAcumulado = 0;

            for (int j = 0; j < Nx + 1; j++) {
                fichero << u[1][j] / xi[j];

                if (j == Nx) {
                    fichero << "\n";
                } else {
                    fichero << ";";
                }

                error = error + (x12[j + 1] - x12[j]) * abs(u[1][j] - xi[j] * fIni(xi[j], domTiempo[z]));
            }
            errores.resize(indiceErrores, error);
        }

        if (((int) (100 * (domTiempo[z - 1] / tFinal))) > porcentajeRK) {
            porcentajeRK = ((int) (100 * (domTiempo[z - 1] / tFinal)));
            cout << std::to_string((int) (100 * (domTiempo[z - 1] / tFinal))) + " % ; t=";
            cout << domTiempo[z - 1] << endl;
        }
        cout << domTiempo[z - 1] << endl;
        z = z + 1;
    }

    for (int i = 1; i < indiceErrores; i++) {
        fichero << errores[i];
        if (i + 1 < indiceErrores) {
            fichero << ";";
        }
    }
    ////
    fichero << "\n";

    for (int i = 1; i < indiceErrores - 1; i++) {
        fichero << domTiempo_printed[i];
        fichero << ";";
    }

    time(&end);

    double tiempoEjecucion = double(end - start);
    fichero << "\n";
    fichero << tiempoEjecucion;

    fichero.close();
    cout << "[INFO]: Simulación finalizada correctamente." << endl;
}
