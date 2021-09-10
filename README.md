# Métodos numéricos para aproximar la ecuación de coagulación de Smoluchowski

En este repositorio están ubicados los códigos desarrollados para implementar los dos métodos numéricos expuestos en la memoria del Trabajo de Fin de Máster.
El programado ha sido llevado a cabo por mi, @ssfructuoso. La lógica de estos métodos fue publicada anteriormente por otros autores, la cual utilizo para el programado de los métodos; las fuentes consultadas son:

* El método de volúmenes finitos:
Filbet F, Laurençot P (2004) Numerical simulation of the Smoluchowski coagulation equation. SIAM J Sci Comput 25(6):2004–2028, doi:10.1137/S1064827503429132

* El método seccional:
Gelbard F, Tambour Y, Seinfeld JH (1980) Sectional representations for simulating aerosol dynamics. J Colloid Interface Sci 76(2):541-556. doi:10.1016/0021-9797(80)90394-X

------------------------------

Los métodos numéricos para aproximar soluciones de la ecuación de Smoluchowski están desarrollados en C++. La ejecución de estas simulaciones es registrada en un fichero de texto, en el directorio que se configure en el script. Posteriormente, se procede con el graficado de los resultados obtenidos importando este fichero de texto como una matriz en unos scripts desarrollados en Python. Se ha optado por esta representación gráfica en Python por tener mayor conocimiento en el manejado del fichero y graficado con distintas configuraciones apoyado de paquetes de graficado ampliamente conocidos.



