import math
import numpy as np
import matplotlib.pyplot as plt
import os
#############################################
############we write a .dat archive###########
##############################################

archivo = open("example","w")
n = 10
xs = np.linspace(0,10,n)
ys = np.linspace(0,100,n)
decimal_places = 4
total_width = 8
format_str = f"{{:.{decimal_places}f}}" 
#format_str = f"{{:{total_width}.{decimal_places}f}}"
for i in range(n):
    archivo.write( format_str.format(xs[i])) #guarda el grado del nodo
    #archivo.write(str(xs[i])) 
    archivo.write(" ") 
    #archivo.write(str(ys[i]))
    archivo.write( format_str.format(ys[i]))
    archivo.write("\n")



