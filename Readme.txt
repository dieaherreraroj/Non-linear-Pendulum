PROYECTO FINAL CURSO HERRAMIENTAS COMPUTACIONALES 2018 - I

FECHA: 20 de Mayo de 2018.
REALIZADOR: Diego Alejandro Herrera Rojas
	    Andrés Francisco Bohada.

DESCRIPCIÓN: El objetivo de este proyecto es estudiar. la relacion entre el rui-
	     do y el caos en un sistema de oscilador no lineal, forzado y con
	     fricción. Se utiliza la DFT para realizar un análisis del movimien-
	     to en el espacio de las frecuencias.

OBJETIVO PRINCIPAL:

	 Desarrollar una comprensión general de la relación entre los conceptos
	 de caos y ruido en un sistema mecánico simple.

OBJETIVOS ESPECÍFICOS:

	 Resolver numéricamente la ecuación de movimiento de un péndulo simple
	 no lineal, forzado y con fricción viscosa proporcional a la velocidad.

	 Realizar un análisis de Fourier de las señal generada por el desplaza-
	 miento angular del péndulo no lineal.

	 Estudiar el espectro de potencias del péndulo no lineal como función
	 del forzamiento externo. Haciéndo énfasis en el fenómeno de mezcla de
	 frecuencias.

	 Descomponer las señales generadas en componentes de ondas simples de a-
	 cuerdo con el modelo de Lorentz. Comparando los regímenes caóticos, no
	 caóticos, intermitentes y "period-doubled".

TAREAS: Las tareas principales se encuentran en el archivo pdf:

	    	   	     06-espectroOscFFT.pdf

        Así como una breve discución del problema a resolver. Sin embargo, se
	realizan unas tareas preliminares, y otras complementarias, con miras
	a producir las señales a descomponer, un informe en latex, y una expo-
	sición sobre el problema. La descripción de estas se encuentra en la
	bitácora del proyecto.

CONTENIDOS: A continuación se presenta un listado de los archivos principales
	    del proyecto:

	    ARCHIVOS DE SIMULACION:

	    En el presente proyecto se utilizaron programas en lenguaje C para
	    simular el sistema fisico considerado. Los archivos que contienen
	    estas rutinas son:

	    	  	         - rk4_integration.c
			         - routines.h

	    El primero es un archivo que al compilarse con la instruccion:
	    gcc rk4_integrarion -lfftw3, produce un ejecutable a.out que a
	    su vez crea archivos .txt con los datos necesarios para realizar
	    analisis grafico. El segundo es un fichero con la declaracion e
	    implementacion de las funciones necesarias para generar las simula-
	    ciones numericas.

	    ARCHIVOS DE DATOS:

	    El objetivo del proyecto es analizar el caos presente en un pendulo
	    simple con forzamiento, al variar la intensidad de la misma. Para
	    ello, se utilizan los parametros fijos:

	    	     	     q = 0.5	     wd = 2.0/3.0

	    Se varia la intensidad de la fuerza en el rango de (0.5:1.5). Para
	    cada una de las intensidades estudiadas se genera un archivo con los
	    datos de la dinamica del sistema, el espectro de potencias, y una a-
	    mnimacion para facilitar el analisis de los resultados. Para cada
	    intensidad se producen archivos con los siguientes prototipos de
	    nombres:

			         - nlp_Fd_*.gif
		                 - fft_Fd_*.txt
				 - ps_Fd_*.txt

	    El primer archivo es una animacion del movimiento que dura 3 ciclos
	    aproximadamente. El segundo contiene el espectro de potencias encon-
	    trado con la libreria fftw para C. El tercero, contiene todos los
	    datos necesarios para analizar la dinamica del sistema, incluyendo
	    la energia mecanica del sistema, y la potencia.
