import math
import numpy as np
import matplotlib.pyplot as plt
import os


x_data = []
y_data = []


datos = 'example'
fichero = open("example")
for item in [data.split()[0] for data in fichero]: #saca la columna 0
    x_data.append(float(item))
fichero.close()

nuevo = open(datos)
for item in [data.split()[1] for data in nuevo]:  
    y_data.append(float(item))   
nuevo.close()
        

print("X Data:", x_data)
print("Y Data:", y_data)