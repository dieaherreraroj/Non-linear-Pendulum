PROYECTO FINAL CURSO HERRAMIENTAS COMPUTACIONALES 2018 - I

FECHA: 20 de Mayo de 2018.
REALIZADOR: Diego Alejandro Herrera Rojas
	    Andr�s Francisco Bohada.

DESCRIPCI�N: El objetivo de este proyecto es estudiar. la relacion entre el rui-
	     do y el caos en un sistema de oscilador no lineal, forzado y con
	     fricci�n. Se utiliza la DFT para realizar un an�lisis del movimien-
	     to en el espacio de las frecuencias.

OBJETIVO PRINCIPAL:

	 Desarrollar una comprensi�n general de la relaci�n entre los conceptos
	 de caos y ruido en un sistema mec�nico simple.

OBJETIVOS ESPEC�FICOS:

	 Resolver num�ricamente la ecuaci�n de movimiento de un p�ndulo simple
	 no lineal, forzado y con fricci�n viscosa proporcional a la velocidad.

	 Realizar un an�lisis de Fourier de las se�al generada por el desplaza-
	 miento angular del p�ndulo no lineal.

	 Estudiar el espectro de potencias del p�ndulo no lineal como funci�n
	 del forzamiento externo. Haci�ndo �nfasis en el fen�meno de mezcla de
	 frecuencias.

	 Descomponer las se�ales generadas en componentes de ondas simples de a-
	 cuerdo con el modelo de Lorentz. Comparando los reg�menes ca�ticos, no
	 ca�ticos, intermitentes y "period-doubled".

TAREAS: Las tareas principales se encuentran en el archivo pdf:

	    	   	     06-espectroOscFFT.pdf

        As� como una breve discuci�n del problema a resolver. Sin embargo, se
	realizan unas tareas preliminares, y otras complementarias, con miras
	a producir las se�ales a descomponer, un informe en latex, y una expo-
	sici�n sobre el problema. La descripci�n de estas se encuentra en la
	bit�cora del proyecto.

CONTENIDOS: A continuaci�n se presenta un listado de los archivos principales
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
