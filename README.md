# Métodos numéricos para aproximar soluciones de la ecuación de coagulación de Smoluchowski

En este repositorio están ubicados los códigos desarrollados para implementar los dos métodos numéricos expuestos en la memoria del Trabajo de Fin de Máster "Soluciones numéricas de la ecuación de coagulación de Smoluchowski: comparativa de dos métodos utilizados y una aplicación en Oceanografı́a.". El programado ha sido elaborado por mi, @ssfructuoso. La lógica de estos métodos fue publicada anteriormente por otros autores, la cual utilizo para el programado de los métodos; las fuentes consultadas son:

* El método de volúmenes finitos:

Filbet F, Laurençot P (2004) Numerical simulation of the Smoluchowski coagulation equation. SIAM J Sci Comput 25(6):2004–2028, doi:10.1137/S1064827503429132

* El método seccional:

Gelbard F, Tambour Y, Seinfeld JH (1980) Sectional representations for simulating aerosol dynamics. J Colloid Interface Sci 76(2):541-556. doi:10.1016/0021-9797(80)90394-X

------------------------------

Los métodos numéricos para aproximar soluciones de la ecuación de Smoluchowski en su forma continua están desarrollados en C++, cada uno en una clase distinta: Seccional.cpp y VolFinitos.cpp. La ejecución de estas simulaciones es registrada en un fichero de texto, en el directorio que se configure en la clase main() del script. Posteriormente, se procede con el graficado de los resultados obtenidos importando este fichero de texto como una matriz en unos scripts desarrollados en Python. Se ha optado por esta representación gráfica en Python por tener mayor conocimiento en el operado de ficheros y graficado con distintas configuraciones apoyado de paquetes de representación gráfica.

En este repositorio los directorios EC_Coagulacion_RK* contienen los métodos seccional y de volúmenes finitos en función del método de Runge-Kutta que se ha utilizado en ellos. Estos directorios a su vez se subdividen en función del tipo de fórmula de cuadratura utilizada para aproximar los valores de las distintas integrales definidas en los métodos. Se prueba con dos tipos de métodos de Runge-Kutta (el método de Heunn de orden 3 y Runge-Kutta clásico de orden 4) y con tres tipos de fórmulas de cuadratura (Regla del Punto Medio, Simpson y Trapecios). En la memoria del trabajo se concluye que para los dos núcleos testeados se obtienen mejores resultados utilizando el método de Runge-Kutta de orden 3 y la formula de cuadratura de Trapecios. 

A partir de esta version (RK3 y Trapecios) se clonan los scripts en el directorio [Ec_Coagulacion_Oceanografia](https://github.com/ssfructuoso/TFM_Ec_Coagulacion/tree/main/Ec_Coagulacion_Oceanografia) para adaptarlo al caso de Oceanografía expuesto en el capítulo tres de la memoria del TFM.

Respecto a los scripts de graficado de los resultados, se puede consultar el directorio [plot_graficas](https://github.com/ssfructuoso/TFM_Ec_Coagulacion/tree/main/plot_graficas).



A continuación se muestra un ejemplo de cómo utilizar los programas. En este caso se va a utilizar el método de seccional (con cuadratura de trapecios y RK3, es decir, la versión del directorio [EC_Coagulacion_trapeciosCompuesto_RK3](https://github.com/ssfructuoso/TFM_Ec_Coagulacion/tree/main/Ec_Coagulacion_RK3/EC_Coagulacion_trapeciosCompuesto_RK3) utilizando el núcleo suma Ker(x,y)=x+y; para el método de volúmenes finitos el procedimiento es análogo. Los puntos a seguir son los siguientes:

1. Definir la distribución inicial n(x_{0},t_{0}).
   En el caso del núcleo suma es conocida la existencia de solución. Por tanto, se opta por definir esta condición inicial en un instante t_{0} para así poder obtener    una comparativa de los errores de aproximación del método numérico bajo este núcleo. La distribución inicial de partículas n(x_{0},t_{0}) se define dentro de la    clase [Seccional.cpp](https://github.com/ssfructuoso/TFM_Ec_Coagulacion/blob/main/Ec_Coagulacion_RK3/EC_Coagulacion_trapeciosCompuesto_RK3/Seccional.cpp),          mediante la función "densidadN0()"; definimos esta función así:
   ```c++
      long double Seccional::densidadN0(long double v, long double t) {
          long double n = pow(2 * PI, -0.5) * exp(-t) * pow(v, -1.5) * exp(-v * 0.5 * exp(-2 * t));
          return n;
      }
   ```

2. Configurar el núcleo  Ker(x,y) de la ecuación de coagulación. En este caso Ker(x,y)=x+y se indica en la misma clase que antes (Seccional.cpp) en la función   "ker":
    ```c++
      long double Seccional::ker(long double x, long double y) {
         return x+y;
      }
    ```
    
3. Definir en la clase [main()](https://github.com/ssfructuoso/TFM_Ec_Coagulacion/blob/main/Ec_Coagulacion_RK3/EC_Coagulacion_trapeciosCompuesto_RK3/main.cpp) el resto de valores a tener en cuenta para ejecutar la simulación. Por ejemplo, pueden ser los siguientes:
    ```c++
      ...
      int main(int argc, char** argv) {

          //Metodo seccional
          cout<< "Seccional: (Trapecios y RK3)" << endl;
          //Constructor de la clase "Seccional", donde se especifica la ruta donde se genera el fichero de texto resultante de la simulacion
          Seccional secRK3 = Seccional("./simulaciones/seccional/seccional_nx100_trapeciosCompuesto_KerSUMA");
          
          /*####
          * insertarGrid(double v0, double R, int m, int numParticionesIntegrales, bool dominioEquiespaciado, bool pasoAdaptativo)
          * ####
          *
          * 1. v0: extremo inferior del dominio espacial.
          * 2. R: amplitud del dominio espacial. Por tanto el extremo derecho v_{final}=v0+R.
          * 3. numParticionesIntegrales: numero de particiones por unidad sobre el intervalo en la variable v definido en cada una de las secciones. 
          *  El particionado en las integrales es equiespaciado y variable.
          * 4. dominioEquiespaciado: determina si el dominio de la variable espacial es equiespaciado o logaritmico. Si dominioEquiespaciado=false entonces
          * el dominio es logaritmico.
          * 5. m: numero de puntos definidos en el mallado de la variable v
          * 6. pasoAdaptativo: determina si se aplica un paso adaptativo en la aproximacion temporal o no; en este caso si se ha aplicado.
          */
          secRK3.insertarGrid(1e-4,200,50,400,false,true);
            
          /*###
          * insertarTiempo(double t0, double tFinal, double incTiempo)
          * ###
          * En el caso de aplicarse un paso adaptativo, la variable incTiempo establece el máximo del salto temporal \Delta t aplicable en el metodo; 
          * en caso de que el dominio temporal sea fijo incTiempo es directamente \Delta t.
          */
          secRK3.insertarTiempo(0.001, 2, 1e-2);
          secRK3.calcular(); //Se inicia la ejecucion del metodo seccional
          
          return 0;
      }
      
    ```
    
 4. Ejecutamos el programa. Mientras esté corriendo va a ir mostrando por pantalla el porcentaje completado. Cuando finalice se habrá generado el fichero [seccional_nx50_trapeciosCompuesto_KerSUMA](https://github.com/ssfructuoso/TFM_Ec_Coagulacion/blob/main/Ec_Coagulacion_RK3/EC_Coagulacion_trapeciosCompuesto_RK3/simulaciones/seccional/seccional_nx50_trapeciosCompuesto_KerSUMA.txt).

 5. Enlazamos la ruta del fichero (o lo copiamos a una ruta que interese) en cualquiera de los scripts de representación gráfica de la carpeta [plot_graficas](https://github.com/ssfructuoso/TFM_Ec_Coagulacion/tree/main/plot_graficas). 
    Por ejemplo, si se busca graficar las soluciones obtenidas en las simulaciones se utiliza el script [plot_general.py](https://github.com/ssfructuoso/TFM_Ec_Coagulacion/blob/main/plot_graficas/plot_general.py). Hay que indicar la ruta donde esté localizado el fichero de texto resultante de la simulación
    ```python
       datosTxt='Ec_Coagulacion_RK3/EC_Coagulacion_trapeciosCompuesto_RK3/simulaciones/seccional/seccional_nx50_trapeciosCompuesto_KerSUMA.txt'
    ```
    y los instantes de tiempo de las soluciones que se quiere mostrar:
    ```python
       tiempo1=0.25
       tiempo2=0.8
       tiempo3=1
       tiempo4=1.75
    ```
    Es posible configurar más puntos relativos al graficado modificando valores de los métodos usados del paquete Matplotlib.
    
    Se obtiene la siguiente gráfica:
    
    ![Seccional_Kernel_Suma](https://github.com/ssfructuoso/TFM_Ec_Coagulacion/blob/main/plot_graficas/screenshots/Ejemplo_grafica_ker_suma_seccional.png)
     
