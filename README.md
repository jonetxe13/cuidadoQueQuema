# Proyecto MPI: Simulación de Difusión Térmica en Placas de Circuitos Impresos

## Introducción

Este repositorio contiene el proyecto desarrollado para la asignatura de Sistemas de Cómputo Paralelo (SCP) del Grado en Ingeniería Informática. El proyecto se centra en la paralelización de un algoritmo de simulación de difusión de calor en placas de circuitos impresos (PCB).

## Descripción del Problema

La empresa TXIPSA, dedicada a la fabricación de PCBs, enfrenta el desafío de la disipación de calor generado por los chips integrados. Una colocación inadecuada de estos componentes puede llevar a un sobrecalentamiento, afectando el rendimiento y la vida útil de la placa. Para mitigar este riesgo, se emplea un algoritmo de simulación de difusión del calor que evalúa diversas configuraciones de chips, con el fin de identificar aquella que minimice la temperatura media global de la tarjeta.

El método de simulación se basa en la discretización de la tarjeta en una rejilla bidimensional. La temperatura de cada punto de la rejilla se calcula de forma iterativa, considerando la temperatura de sus vecinos.

El proceso iterativo contempla la inyección de calor en las ubicaciones de los chips y la disipación en zonas designadas como ventiladas. La simulación converge cuando la variación de la temperatura media entre iteraciones sucesivas (evaluada cada 10 iteraciones) es inferior a un umbral (referido como `t_delta` en el proyecto) o se alcanza un número máximo de iteraciones preestablecido (referido como `max_iter`).

El programa original, `heat_s.c`, implementa esta simulación de forma secuencial. El objetivo fundamental del proyecto es desarrollar e implementar versiones paralelas de este algoritmo utilizando la interfaz de paso de mensajes (MPI) para su ejecución en un clúster de memoria distribuida.

## Fases del Desarrollo

La paralelización del algoritmo se aborda en dos fases principales:

### Fase 1: Paralelización Homogénea de Configuraciones

En esta primera fase, la simulación de cada configuración de chips es abordada concurrentemente por todos los procesos MPI. La rejilla bidimensional que representa la tarjeta se descompone en franjas horizontales, asignando la responsabilidad de la simulación de cada franja a un proceso MPI distinto. Esta descomposición requiere la comunicación de los valores de temperatura de las fronteras entre procesos adyacentes. Se han desarrollado dos variantes:
* **Versión con Comunicación Síncrona**: Empleando primitivas de MPI para la comunicación síncrona de datos (`MPI_Ssend` y `MPI_Recv`).
* **Versión con Comunicación Asíncrona (Inmediata)**: Optimizando la versión anterior mediante el uso de primitivas de MPI para la comunicación no bloqueante (`MPI_Isend` y `MPI_Irecv`).

### Fase 2: Paralelización con Planificación Dinámica (Manager-Worker)

La segunda fase introduce un paradigma de planificación dinámica basado en el modelo Manager-Worker. Un proceso MPI, designado como *manager*, se encarga de distribuir las tareas de simulación (configuraciones de chips) entre varios grupos de procesos *worker*. Cada grupo de *workers*, cuyo tamaño `P` se determina a partir de los resultados de la Fase 1, procesa una configuración en paralelo aplicando la estrategia de descomposición por franjas horizontales implementada en la Fase 1. Una vez completada la simulación, el grupo de *workers* devuelve el resultado al *manager*, quien procede a asignar una nueva configuración si quedan tareas pendientes. Dentro de cada grupo de `P` procesos, uno de ellos asume el rol de coordinador para la comunicación con el *manager*.

## Estructura del Repositorio

El proyecto se organiza en la siguiente estructura de directorios y ficheros:

* `Serie/`: Contiene la implementación secuencial de referencia.
    * `heat_s.c`: Programa principal.
    * `diffusion_s.c`: Funciones para la simulación de la difusión térmica.
    * `faux_s.c`: Funciones auxiliares (lectura de datos, gestión de resultados).
    * `defines.h`: Definiciones de constantes y estructuras de datos.
    * `script.sh`: Script para la compilación y ejecución.
* `Par_Sincrono/`: Implementación de la Fase 1 con comunicación MPI síncrona.
    * `heat_pSinc.c`: Programa principal paralelo.
    * `diffusion_pSinc.c`: Funciones de simulación adaptadas para MPI síncrono.
    * `faux_pSinc.c`: Funciones auxiliares adaptadas.
    * `defines.h`: Definiciones específicas.
    * `script.sh`: Script de compilación y ejecución.
* `Par_Inmediata/`: Implementación de la Fase 1 con comunicación MPI asíncrona (inmediata).
    * `heat_pAsinc.c`: Programa principal paralelo.
    * `diffusion_pAsinc.c`: Funciones de simulación adaptadas para MPI asíncrono.
    * `faux_pAsinc.c`: Funciones auxiliares adaptadas.
    * `defines.h`: Definiciones específicas.
    * `script.sh`: Script de compilación y ejecución.
* `Par_master-slave/`: Implementación de la Fase 2 (modelo Manager-Worker).
    * `heat_pMS.c`: Programa principal paralelo.
    * `diffusion_pMS.c`: Funciones de simulación adaptadas para el modelo manager-worker.
    * `faux_pMS.c`: Funciones auxiliares adaptadas.
    * `defines.h`: Definiciones específicas.
    * `script.sh`: Script de compilación y ejecución.
* `Resultados/`: Directorio destinado al almacenamiento de los ficheros de resultados generados por las distintas versiones del programa.
    * `resultados-serie.txt`: Salida de la versión secuencial.
    * `resultados-totales-sinc.txt`: Salida de la versión paralela síncrona.
    * `resultados-totales-asinc.txt`: Salida de la versión paralela asíncrona.
    * `resultados-totales-ms.txt`: Salida de la versión manager-worker.
* `Cards/`: Contiene los ficheros de entrada que definen las características de la tarjeta y las configuraciones de los chips a simular.

Cada directorio de implementación incluye, además de los ficheros fuente, un fichero `defines.h` con las estructuras de datos y constantes pertinentes, y un `script.sh` para facilitar la compilación y ejecución de las simulaciones.
