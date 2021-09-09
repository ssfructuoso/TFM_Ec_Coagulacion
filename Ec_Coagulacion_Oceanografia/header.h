/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   header.h
 * Author: salvador
 *
 * Created on 19 de abril de 2020, 0:55
 */

#ifndef HEADER_H
#define HEADER_H

#include <vector>

class Seccional {
    int m;
    long double v0;
    long double R;
    long double alfa;
    long double gamma;
    int numParticionesIntegrales;
    long double * v;
    long double * vi;
    long double * factorQ;
    long double * x;
    long double t0;
    long double tFinal;
    long double incTiempo;
    int numIntervalosTiempo;
    const char * nombreTxtResultado;
    const char * nombreTxtCoeficientesSeccionales;
    bool dominioEquiespaciado;
    bool pasoAdaptativo;
    
    int minimoParticiones = 25;

    long double ** Q;

    long double * I;
    long double * w;

    long double a_0;
    long double D;
    long double kBolztmann;
    long double T;
    long double din_viscosity;
    long double sh_rate;
    long double p_0;
    long double p_w;
    long double g;
    long double Z;

    long double *** C1;
    long double ** C2;
    long double * C3;
    long double ** C4;





public:
    Seccional(const char * nombreTxtResultado);

    void insertarGrid(long double v0, long double R, int m, int numParticionesIntegrales, bool dominioEquiespaciado, bool pasoAdaptativo);
    void insertarTiempo(long double t0, long double tFinal, long double incTiempo);

    void calcular();

    long double densidadN0(long double v, long double t);
    long double densidad_media0(int i, long double t);
    long double ker(long double x, long double y);
    long double f_point(long double v);
    long double df(long double v);
    long double fInverse_point(long double x);
    long double thetaInfSup(long double limInf, long double limSup, long double valor);
    long double thetaInf(long double limInf, long double valor);
    long double thetaSup(long double limSup, long double valor);
    //void vVector(long double v [], int lenV, long double v0,long  long double R,long  long double corte1,long  long double corte2);
    //void vKerId(long double v [], int lenV, long double v0, long double R);
    long double coefC1(int i, int j, int l);
    long double coefC2(int i, int l);
    long double coefC3(int l);
    long double coefC4(int i, int l);

    long double dQ(int l, long double add);
};
///////////////////////////////////////////////////
class VolFinitos {
    long double R;
    long double xInicio;
    int Nx;
    bool dominioAsimetrico;
    long double h;
    int numParticionesIntegrales;
    long double t0;
    long double tFinal;
    long double incTiempo;
    int numIntervalosTiempo;
    const char * nombreTxtResultado;
    long double kdominio;
    bool dominioEquiespaciado;
    bool pasoAdaptativo;
    int minimoParticiones = 25;

    //////////////////////

    int NxIni;

    long double * x12Ini;
    long double * xiIni;
    long double ** g;


    long double * I;
    long double * w;

    long double a_0;
    long double D;
    long double kBolztmann;
    long double T;
    long double din_viscosity;
    long double sh_rate;
    long double p_0;
    long double p_w;
    long double g_gravedad;
    long double Z;

    ////////////

    int contador;
    long double * x12;
    //std::vector <double> x12;
    long double * xi;


    //long double ** g;
    long double ** J;
    long double ** Aa;
    long double ** Bb;
    int ** alfa;

public:
    VolFinitos(const char * nombreTxtResultado);
    void insertarGrid(long double xInicio, long double R, long double Nx, int numParticionesIntegrales, bool dominioEquiespaciado, bool pasoAdaptativo);
    void insertarTiempo(long double t0, long double tFinal, long double incTiempo);

    void calcular();

    long double fIni(long double x, long double t);
    long double f_media0(int i, long double t);
    long double gIni(long double x, long double t);
    long double g_media0(int i, long double t);
    long double ker(long double x, long double y);
    long double A(int j, long double x_k);
    long double B(int i, long double x_k, int alfa);
    int getAlfa(int i, int k);
    
    long double dg(int i, long double add, int kS);

};


#endif /* HEADER_H */
