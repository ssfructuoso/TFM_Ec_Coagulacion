/*
 * Clase para la simulacion de soluciones de la ecuacion de Smoluchowski
 * utilizando un metodo de volumenes finitos.
 * En este codigo primeramente se definen distintas funciones (calculo de coeficientes seccionales,
 * definicion del nucleo, etc) para luego poder plasmar y calcular el esquema numerico en la 
 * funcion "calcular()".
 * 
 * Los resultados obtenidos se guardan en un fichero de texto con un formato que se puede consultar en la 
 * funcion "calcular()".
 * 
 * En las formulas de cuadratura se aplica el metodos de trapecios compuesto y un metodo de Runge-Kutta 
 * de orden 3 para la aproximacion temporal.
 * 
 */

/* 
 * File:   VolFinitos.cpp
 * Author: Salvador Fructuoso
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


/*
 * Distribucion inicial de puntos: en este caso es la solucion en t_0 para el nucleo producto.
 */
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

/*
 * Promedio de la distribuci??n inicial de puntos.
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
    
    //Trapecios compuesto
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

/*
 * Primera integral contenida en el flujo J_{i+1/2}
 * en el metodo de volumenes finitos
 */
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

/*
 * Segunda integral contenida en el flujo J_{i+1/2}
 * en el metodo de volumenes finitos
 */
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

/*
 Valor \alpha_{i,k} para los subindices de los distintas puntos en el dominio
 */
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


/*
 Constructor de la clase
 */
VolFinitos::VolFinitos(const char * nombreTxtResultado) {
    this->nombreTxtResultado = nombreTxtResultado;
}

void VolFinitos::insertarGrid(long double xInicio, long double R, long double Nx, int numParticionesIntegrales, bool dominioEquiespaciado, bool pasoAdaptativo) {
    this->xInicio = xInicio;
    this->R = R;
    this->Nx = Nx;
    this->numParticionesIntegrales = numParticionesIntegrales;
    this->dominioEquiespaciado = dominioEquiespaciado;
    this->pasoAdaptativo=pasoAdaptativo;
}

void VolFinitos::insertarTiempo(long double t0, long double tFinal, long double incTiempo) {
    this->t0 = t0;
    this->tFinal = tFinal;
    this->incTiempo = incTiempo;
    this->numIntervalosTiempo = (int) ((tFinal - t0) / incTiempo) + 1;
}


/*
 * dg/dt: se invoca en el metodo de Runge-Kutta posteriormente
 */
long double VolFinitos::dg(int i, long double add, int kS) {
    long double suma = 0;
    long double sumaK = 0;
    long double sumaJ = 0;
    long double x_k;
    int alfaValue;
    long double kReturn;
    long double uAux [Nx + 1];

    for (int i = 0; i < Nx + 1; i++) {
        uAux[i] = g[0][i] + add;
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
    kReturn = -(J[kS - 1][i + 1] - J[kS - 1][i]) / (x12[i + 1] - x12[i]);

    return kReturn;
}


/*
 * Mediante esta funcion se invocan a las funciones definidas anteriormente y se 
 * ejecutan las simulaciones
 */
void VolFinitos::calcular() {
    time_t start, end;
    time(&start);//Inicio del calculo del tiempo de ejecucion

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

    
    cout << "--------------" << endl;
    cout << "Intervalos de tiempo:";
    cout << numIntervalosTiempo << endl;

    /////////////////////////////////////////////////
    ////////////////////////////
    ///Generacion de un fichero de texto en el cual guardar los resultados 
    ///de la simulacion
    ofstream fichero;
    char buf[120];
    snprintf(buf, sizeof (buf), "%s.txt", this->nombreTxtResultado);
    fichero.open(buf);

    fichero.precision(20);
    //Escritura de variables en el fichero de texto
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
    valores(); //Invocacion de valores() para generar los valores de las distintas integrales presentes
               // en el flujo J_{i+1/2}
    
    
    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    ///////////Inicio de la implementaci??n del esquema numerico
    
    
    g = new long double *[2];
    for (int i = 0; i < 2; i++) {
        g[i] = new long double [Nx + 1];
    }
    for (int i = 0; i < Nx + 1; i++) {// Cambio Nx+1 por Nx
        g[0][i] = g_media0(i, t0);
    }

    J = new long double *[3];

    for (int i = 0; i < 3; i++) {
        J[i] = new long double [Nx + 2];
        J[i][0] = 0;
    }
    //////////
    int z = 1;

    long double k1 = 0;
    long double k2 = 0;
    long double k3 = 0;

    long double incTime = 0;//Variable que define \Delta t
    bool tiempoFino = false;//Variable booleana para comprobar si el 
                      //salto temporal cumple unas condiciones impuestas
                      //Se usa para el caso de un salto temporal variable.

    long double vecDomTiempoAux = 0;
    vector<long double > domTiempo;
    domTiempo.push_back(t0);

    vector<long double > domTiempo_printed;//Vector donde guardar cada uno de los tiempos usados en las iteraciones

    long double error = 0; //Error de aproximacion en cada salto temporal
                           //respecto a una soluci??n exacta conocida.
    int indiceErrores = 1;
    vector<long double > errores;//Aray donde guardar los errores de cada salto
                                //temporal
    errores.push_back(0);

    long double tiempoAcumulado = 0;
    int porcentajeRK = 0;//Variable para mostrar en pantalla el porcentaje 
    //de la simulaci??n completado.

    cout << ">>> RK3 >>>" << endl; 
    
    while (domTiempo[z - 1] < tFinal) {

        incTime = incTiempo;//Tiempo fijo \Delta t definido en la clase main()
        long double value = 0;//Variable auxiar para evaluar el valor la diferencia entre 
                              //n(x,t^n) y n(x,t^{n+1})
        tiempoFino = false;//Boleano para ver si el mallado temporal variable es apto.
        int ll;

        long double vFinal;
        int binarioParidad = 2;
        int numValoresNegativos = 0;
        
        ///////////////////////////////////
        ///Bucle para evaluar el salto temporal en cada instante
        //Reduce \Delta t si no cumple las condiciones del bucle:
        //      1. |n(x,t^{n+1})-n(x,t^n})|<=0.0005 en este caso
        //      2. La aproximaci??n tenga m??s de umbral de puntos definido que
        //         sean negativos. En este caso hay un l??mite de subdivisi??n de 
        //         \Delta t; si se alcaza ese l??mite el programa se para y muestra
        //         un Warning por pantalla.
        //
        //En caso de querer un dominio temporal fijo, basta con imponer una condici??n
        //falsa el bucle para que no se ejecute.
        while (tiempoFino == false && pasoAdaptativo==true) {
            ll = 0;
            double incLL;
            //Mediante la variable "binarioParidad" se testea puntos salteados mediante el siguiente
            //bloque IF, de modo que no se revisen todos los puntos del mallado computacional, sino
            //que estos cambien en cada iteracion y as?? la simulaci??n no se alargue mucho. 
            if (binarioParidad % 2 == 0) {
                incLL = 1;
            } else {
                incLL = 4;
            }
            while (ll <= Nx) {
                k1 = dg(ll, 0, 1);
                k2 = dg(ll, incTime * k1 / 3.0, 2);
                k3 = dg(ll, 2 * incTime * k2 / 3.0, 3);

                value = incTime * (k1 + 3 * k3) / 4.0;

                vFinal = g[0][ll] + value;

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
                        cout << "[WARNING]: Valor negativo de la aproximaci??n (secci??n ";
                        cout << ll; cout << "; ";
                        cout << vFinal;
                        cout << ") en el instante de tiempo t=";
                        cout << domTiempo[z - 1] << endl;
                    } else if (vFinal < -1e-10 || numValoresNegativos >= umbralMaximoValoresNegativos) {
                        cout << "[ERROR] Aproximaci??n cancelada: Error en aproximaci??n. Valores negativos." << endl;
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
        //
        //Fin bucle de adaptaci??n del dominio temporal.
        //////////////////////


        vecDomTiempoAux = domTiempo[z - 1] + incTime;
        domTiempo.resize(z + 1, vecDomTiempoAux);

        //Aplicacion del metodo de Runge-Kutta.
        for (int i = 0; i <= Nx; i++) {
            k1 = dg(i, 0, 1);
            k2 = dg(i, incTime * k1 / 3.0, 2);
            k3 = dg(i, 2 * incTime * k2 / 3.0, 3);

            g[1][i] = g[0][i] + incTime * (k1 + 3 * k3) / 4.0;

            if (isnan(g[1][i]) == true) {
                cout << "[ERROR] Aproximaci??n cancelada: Nan." << endl;
                std::exit(EXIT_FAILURE);
            }
        }


        for (int i = 0; i < Nx + 1; i++) {
            g[0][i] = g[1][i];
        }

        tiempoAcumulado = tiempoAcumulado + incTime;

        if (tiempoAcumulado >= 0.005) {//De esta forma se guarda en el fichero de texto
            //puntos con una distancia temporal m??nima, para asi evitar que el fichero
            //generado sea de mucho tama??o.
            error = 0;
            domTiempo_printed.resize(indiceErrores, vecDomTiempoAux);
            indiceErrores++;
            tiempoAcumulado = 0;

            for (int j = 0; j < Nx + 1; j++) {
                fichero << g[1][j] / xi[j];//Escritura de soluciones en el fichero de texto

                if (j == Nx) {
                    fichero << "\n";
                } else {
                    fichero << ";";
                }
                //Introduccion de los errores de calculo  
                error = error + (x12[j + 1] - x12[j]) * abs(g[1][j] - xi[j] * fIni(xi[j], domTiempo[z]));
            }
            errores.resize(indiceErrores, error);
        }
        //
        //Mostrado por pantalla del porcentaje completado de la simulaci??n 
        if (((int) (100 * (domTiempo[z - 1] / tFinal))) > porcentajeRK) {
            porcentajeRK = ((int) (100 * (domTiempo[z - 1] / tFinal)));
            cout << std::to_string((int) (100 * (domTiempo[z - 1] / tFinal))) + " % ; t=";
            cout << domTiempo[z - 1] << endl;
        }
        //cout << domTiempo[z - 1] << endl;
        z = z + 1;
    }
    
    
    //Guardado en el fichero de texto el array de errores respecto a soluciones conocidas.No se guardan todos, sino los que 
    //cumplan un m??nimo \Delta t en funci??n de la variable "tiempoAcumulado" para as?? evitar un 
    //exceso de tama??o del fichero.
    for (int i = 1; i < indiceErrores; i++) {
        fichero << errores[i];
        if (i + 1 < indiceErrores) {
            fichero << ";";
        }
    }
    ////
    fichero << "\n";
    //Guardado en el fichero de texto los instantes temporales. No se guardan todos, sino los que 
    //cumplan un m??nimo \Delta t en funci??n de la variable "tiempoAcumulado" para as?? evitar un 
    //exceso de tama??o del fichero.
    for (int i = 1; i < indiceErrores - 1; i++) {
        fichero << domTiempo_printed[i];
        fichero << ";";
    }

    //////////
    time(&end);//Fin del c??lculo del tiempo de ejecuci??n del programa

    double tiempoEjecucion = double(end - start);
    fichero << "\n";
    fichero << tiempoEjecucion;

    fichero.close();
    cout << "[INFO]: Simulaci??n finalizada correctamente." << endl;
}
