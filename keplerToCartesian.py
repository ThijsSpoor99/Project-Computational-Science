"""This code file contains code to convert the Ke"""
import numpy as np
import pandas as pd
import pyorb

Centaurs = pd.read_csv("Data/Centaurs.csv", index_col=0)

Centaurs_before_conversion = Centaurs[['a', 'e', 'i', 'w', 'om', 'ma']]
Centaurs_before_conversion.loc[:, 'a'] = Centaurs_before_conversion['a'] * 1.496e+11
Centaurs_before_conversion.loc[:, 'w'] = np.radians(Centaurs_before_conversion['w']) 
Centaurs_before_conversion.loc[:, 'om'] = np.radians(Centaurs_before_conversion['om'])
Centaurs_before_conversion.loc[:, 'ma'] = np.radians(Centaurs_before_conversion['ma'])

Centaurs_cartesian = []
for i in range(Centaurs_before_conversion.shape[0]):
    Centaurs_cartesian.append(pyorb.kep_to_cart(Centaurs_before_conversion.iloc[i].values))

Centaurs_cartesian = pd.DataFrame(Centaurs_cartesian, columns=['x','y','z','Vx','Vy','Vz'])
Centaurs_cartesian.to_csv("Data/CentaursCartesian.csv")