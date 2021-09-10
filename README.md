# Métodos numéricos para aproximar la ecuación de coagulación de Smoluchowski

En este repositorio están ubicados los códigos desarrollados para implementar los dos métodos numéricos expuestos en la memoria del Trabajo de Fin de Máster.
El programado ha sido llevado a cabo por mi, @ssfructuoso. La lógica de estos métodos fue publicada anteriormente por otros autores, la cual utilizo para el programado de los métodos; las fuentes consultadas son:

* El método de volúmenes finitos:
Filbet F, Laurençot P (2004) Numerical simulation of the Smoluchowski coagulation equation. SIAM J Sci Comput 25(6):2004–2028, doi:10.1137/S1064827503429132

* El método seccional:
Gelbard F, Tambour Y, Seinfeld JH (1980) Sectional representations for simulating aerosol dynamics. J Colloid Interface Sci 76(2):541-556. doi:10.1016/0021-9797(80)90394-X

------------------------------

Los métodos numéricos para aproximar soluciones de la ecuación de Smoluchowski están desarrollados en C++. La ejecución de estas simulaciones es registrada en un fichero de texto, en el directorio que se configure en la clase main() del script. Posteriormente, se procede con el graficado de los resultados obtenidos importando este fichero de texto como una matriz en unos scripts desarrollados en Python. Se ha optado por esta representación gráfica en Python por tener mayor conocimiento en el operado de ficheros y graficado con distintas configuraciones apoyado de paquetes de representación gráfica ampliamente conocidos.

En este repositorio los directorios EC_Coagulacion_RK* contienen los métodos seccional y de volúmenes finitos en función del método de Runge-Kutta que se ha utilizado en ellos. Estos directorios a su vez se subdividen en función del tipo de fórmula de cuadratura utilizada para aproximar los valores de las distintas integrales definidas en los métodos. Se prueba con dos tipos de métodos de Runge-Kutta (el método de Heunn de orden 3 y Runge-Kutta clásico de orden 4) y con tres tipos de fórmulas de cuadratura (Regla del Punto Medio, Simpson y Trapecios). En la memoria del trabajo se concluye que para los dos núcleos testeados se obtienen mejores resultados utilizando el método de Runge-Kutta de orden 3 y la formula de cuadratura de trapecios. 

A partir de esta version (RK3 y Trapecios) se clonan los scripts en el directorio "Ec_Coagulacion_Oceanografia" para adaptarlo al caso de Oceanografía expuesto en el el capítulo tres de la memoria del TFM.

Respecto a los scripts de graficado de los resultados, se puede consultar el directorio plot_graficas.





