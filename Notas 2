=======================
#ARCHIVO DE OBJETIVOS DE foamToKIVA4Mesh
gedit /home/carlos/OpenFOAM/carlos-4.0/run/foamToKIVA4Mesh/Notas
pcmanfm -n /home/carlos/OpenFOAM/carlos-4.0/run/foamToKIVA4Mesh/ /home/carlos/OpenFOAM/carlos-4.0/platforms/linux64GccDPInt32Opt/bin
pcmanfm -n /home/carlos/OpenFOAM/carlos-4.0/platforms/linux64GccDPInt32Opt/bin /home/carlos/Descargas/TMP "/home/carlos/archivos/Colciencias/Jovenes Investigadores 2015/Mallas"


gedit --new-window /home/carlos/OpenFOAM/carlos-4.0/platforms/linux64GccDPInt32Opt/bin/kiva4grid





pcmanfm -n
gedit --new-window
Se compila con:    clear && clear && wmake 


carpetas: 
 source:       /home/carlos/OpenFOAM/carlos-4.0/run/foamToKIVA4Mesh/
 binarios:    /home/carlos/OpenFOAM/carlos-4.0/platforms/linux64GccDPInt32Opt/bin
 testeo cubo:  /home/carlos/archivos/Colciencias/Jovenes Investigadores 2015/Mallas /Malla Kiva4-Testing
 mallas varias: /home/carlos/archivos/Colciencias/Jovenes Investigadores 2015/Mallas /Malla Salome

convertir mallas:
     1. copiar carperas constant y system en una ruta sin espacios
     1.5. copiar la malla en formato unv:  *.unv 
     2. cargar entorno: source /home/carlos/OpenFOAM/OpenFOAM-4.0/etc/bashrc
     3. ejecutar  ideasUnvToFoam  *.unv 
-Hacer OpenFoam2kiva4grid
  -Hacer lector de OpenFoam
  -Escribir kiva4gird con nodos
  -Conectividades: Usar "pointfield" para sacar los labels de los ptos y con los esos mismos labels, sacar los puntos
http://www.openfoam.com/documentation/cpp-guide/html/classFoam_1_1cell.html#a77ca2209afb8888c352cf6c8d4b4380c
  -Tipos soportados:  Pirámide, tetrahedro, hexaedro.
  -Tener en cuenta la numeración de los nodos que es diferente
  -¿Cómo sería el caso de mallas no hexagonales? -> Manual de Kiva4 Fig 2 y 3

  -Condiciones de frontera: Asociar caras a las celdas
        https://github.com/OpenFOAM/OpenFOAM-4.x/blob/master/applications/utilities/mesh/conversion/kivaToFoam/kivaToFoam.C
    -Celda
    -Nodos (Con la convención de KIVA4)
    -Caras (en la posición que indique OpenFOAM luego de hacer la conversión de vértices)
    -Cond. frontera -> Caras 
    -Tipo de Celda -> Aparecen los grupos de volúmenes en el formato (Carlos Barrera)
    
  -MALLA ESTRUCTURADA: Carlos Barrera
    -%TODO : Para sacar las celdas vecinas y esas cosas, al parecer toca invocar las librerías de OpenFoam: https://www.cfd-online.com/Forums/openfoam-meshing-technical/61906-how-identify-cellsudop-neighbours.html

########################
    
*-Hacer la ingeniería inversa de las librerías de los formatos de archivo
*  -Hacer lista de posibles funciones :https://www.cfd-online.com/Forums/openfoam-meshing-technical/61906-how-identify-cell-neighbours.html
www.cfd-online.com/Forums/openfoam-meshing-technical/61906-how-identify-cellsudop-neighbours.html

