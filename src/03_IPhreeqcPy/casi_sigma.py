# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt
#%%
#casi = np.array([1.67, 1.25, 0.83])
#sigma = np.array([3.18485736443e-06, 1.25655425766e-07, 1.63726878765e-08])

casi = np.array([1.67, 1.65, 1.60, 1.55, 1.50, 1.45, 1.40, 1.35, 1.30, 1.25,
                  1.20, 1.15, 1.10, 1.05, 1.00, 0.95, 0.90, 0.85, 0.83])
d = np.array([-8.40, -8.47, -8.53, -8.53, -8.56, -8.67, -8.89, -9.20, -9.58,
              -9.98, -10.36, -10.67, -10.88, -10.97, -10.95, -10.86, -10.79, 
              -10.87, -10.99])
s = np.array([865, 880, 925,  970, 1015, 1060, 1105, 1150, 1195, 1240, 1285, 1330,
              1375, 1420, 1465, 1510, 1555, 1600, 1615])
sigma = (10**d) * s
plt.figure()
plt.plot(casi, sigma)
plt.yscale("log")
plt.xlabel("Ca/Si")
plt.ylabel(r"$\sigma$")

#%% predict
from sklearn.linear_model import LinearRegression

X = casi.reshape(-1,1) # put your dates in here
y = sigma # put your kwh in here

model = LinearRegression()
model.fit(X, y)

X_predict = np.array([1.52654, 0.99]).reshape(-1,1)  # put the dates of which you want to predict kwh here
y_predict = model.predict(X_predict)
print(y_predict)