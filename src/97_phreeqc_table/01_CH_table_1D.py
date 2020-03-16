# -*- coding: utf-8 -*-

import IPhreeqcPy
import numpy as np
import time

#%% FUNCTIONS
def phrqc_string1(c, ca, ch, cc, water, n):
    phrqc_input = []    
    for i in np.arange(1, n+1):
        phrqc_input.append('SOLUTION\t10000' + str(i))
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t'+ str(water[i-1]))
        phrqc_input.append('\tpH\t7\tcharge')
        phrqc_input.append('\tC\t' + str(c[i-1]) )
        phrqc_input.append('\tCa\t' + str(ca[i-1]) + '\n')
        phrqc_input.append('EQUILIBRIUM_PHASES\t10000'+ str(i))
        phrqc_input.append('\tcalcite\t0\t' + str(cc[i-1]) )
        phrqc_input.append('\tportlandite\t0\t' + str(ch[i-1]) + '\n')
    
   
    phrqc_input.append('SELECTED_OUTPUT')
    phrqc_input.append('\t-reset false')
    phrqc_input.append('\t-time false')
    phrqc_input.append('\t-high_precision true')
    phrqc_input.append('\t-solution true')
    phrqc_input.append('\t-pH true')
    phrqc_input.append('\t-pe true')
    phrqc_input.append('\t-charge_balance false')
    phrqc_input.append('\t-alkalinity true')
    phrqc_input.append('\t-ionic_strength false')
    phrqc_input.append('\t-percent_error false')
    
    phrqc_input.append('USER_PUNCH')
    phrqc_input.append('\t10\tpunch\ttot("C")')
    phrqc_input.append('\t20\tpunch\ttot("Ca")')
    phrqc_input.append('\t30\tpunch\ttot("O")')
    phrqc_input.append('\t40\tpunch\ttot("H")')
    phrqc_input.append('\t50\tpunch\tequi("portlandite")')
    phrqc_input.append('\t60\tpunch\tSI("portlandite")')
    phrqc_input.append('\t70\tpunch\tequi("calcite")')
    phrqc_input.append('\t80\tpunch\tSI("calcite")')
    phrqc_input.append('\t90\tpunch\ttot("water")')
    phrqc_input.append('\t100\tpunch\ttot("H2O")')
    phrqc_input.append('\t110\tpunch')
    phrqc_input.append('\t-headings\tC\tCa\tO\tH\tportlandite\tSI_portlandite' +
                       '\tcalcite\tSI_calcite\twater\tH2O')
    phrqc_input.append('\t-start')
    phrqc_input.append('\t-end')
    phrqc_input.append('end')
    return '\n'.join(phrqc_input)


def phrqc_string2(n):
    phrqc_input = []   
    for i in np.arange(1, n+1):
        phrqc_input.append('copy cell 10000'+str(i) + ' ' + str(i))
    phrqc_input.append('end')
    phrqc_input.append('RUN_CELLS')
    phrqc_input.append('\t-cells 1-' + str(n+1))
    phrqc_input.append('\t-start_time 0')
    phrqc_input.append('\t-time_step 0')
    #phrqc_input.append('end')    
   
    phrqc_input.append('SELECTED_OUTPUT')
    phrqc_input.append('\t-reset false')
    phrqc_input.append('\t-time false')
    phrqc_input.append('\t-high_precision true')
    phrqc_input.append('\t-solution true')
    phrqc_input.append('\t-pH true')
    phrqc_input.append('\t-pe true')
    phrqc_input.append('\t-charge_balance false')
    phrqc_input.append('\t-alkalinity true')
    phrqc_input.append('\t-ionic_strength false')
    phrqc_input.append('\t-percent_error false')
    
    phrqc_input.append('USER_PUNCH')
    phrqc_input.append('\t10\tpunch\ttot("C")')
    phrqc_input.append('\t20\tpunch\ttot("Ca")')
    phrqc_input.append('\t30\tpunch\ttot("O")')
    phrqc_input.append('\t40\tpunch\ttot("H")')
    phrqc_input.append('\t50\tpunch\tequi("portlandite")')
    phrqc_input.append('\t60\tpunch\tSI("portlandite")')
    phrqc_input.append('\t70\tpunch\tequi("calcite")')
    phrqc_input.append('\t80\tpunch\tSI("calcite")')
    phrqc_input.append('\t90\tpunch\ttot("water")')
    phrqc_input.append('\t100\tpunch\ttot("H2O")')
    phrqc_input.append('\t110\tpunch')
    phrqc_input.append('\t-headings\tC\tCa\tO\tH\tportlandite\tSI_portlandite' +
                       '\tcalcite\tSI_calcite\twater\tH2O')
    phrqc_input.append('\t-start')
    phrqc_input.append('\t-end')
    phrqc_input.append('end')
    return '\n'.join(phrqc_input)
#%% MESH CONCENTRATIONS 1D
n =100
x1 = np.linspace(0.0, 0.1, n)# C
c= x1#map(np.ravel, np.meshgrid(x1))
#x2 = np.linspace(0.01, 0.1, n)#Ca
#x3 = np.linspace(0, 1, n)#CH
#x4 =  np.linspace(0, 1, n) #CC
#x5 = np.linspace(0.5, 1, n) #water=porosity
#c,ca,water = map(np.ravel, np.meshgrid(x1,x2,x5))

ca = np.ones(np.shape(c))*0.01
water = np.ones(np.shape(c))
ch = np.zeros(np.shape(c))
cc = np.zeros(np.shape(c))
#%%
it=time.time() 

n_cases = len(ca)
ps1 = phrqc_string1(c, ca, ch, cc, water, n_cases)
ps2 = phrqc_string2(n_cases)
IPhreeqc = IPhreeqcPy.IPhreeqc()
IPhreeqc.LoadDatabase('C:\Anaconda2\lib\site-packages\databases\cemdata07.dat')
IPhreeqc.RunString(ps1)
so1 = IPhreeqc.GetSelectedOutputArray()
IPhreeqc.RunString(ps2)
so2 = IPhreeqc.GetSelectedOutputArray()
simulation_time = time.time()-it
print(simulation_time)
#print(so1)
#print(so2)
#print([row[0] for row in so]) #C
#print([row[1] for row in so]) #Ca
import pandas as pd
res = {'i_sol':[row[0] for row in so1][1:-1],
       'i_pH': [row[1] for row in so1][1:-1],
       'i_pe': [row[2] for row in so1][1:-1],
       'i_Alk': [row[3] for row in so1][1:-1],
       'i_C': [row[4] for row in so1][1:-1],
       'i_Ca': [row[5] for row in so1][1:-1],
       'i_O': [row[6] for row in so1][1:-1],
       'i_H': [row[7] for row in so1][1:-1],
       'i_CH': ch,#[row[8] for row in so1][1:-1],
       'i_SI_CH': [row[9] for row in so1][1:-1],
       'i_CC': cc,#[row[10] for row in so1][1:-1],
       'i_SI_CC': [row[11] for row in so1][1:-1],
       'i_poros': [row[12] for row in so1][1:-1],
       'i_H20': [row[13] for row in so1][1:-1],
       'o_sol':[row[0] for row in so2][1::],
       'o_pH': [row[1] for row in so2][1::],
       'o_pe': [row[2] for row in so2][1::],
       'o_Alk': [row[3] for row in so2][1::],
       'o_C': [row[4] for row in so2][1::],
       'o_Ca': [row[5] for row in so2][1::],
       'o_O': [row[6] for row in so2][1::],
       'o_H': [row[7] for row in so2][1::],
       'o_CH': [row[8] for row in so2][1::],
       'o_SI_CH': [row[9] for row in so2][1::],
       'o_CC': [row[10] for row in so2][1::],
       'o_SI_CC': [row[11] for row in so2][1::],
       'o_poros': [row[12] for row in so2][1::],
       'o_H20': [row[13] for row in so2][1::]}
df = pd.DataFrame(data=res)
#%% Plot
pr = 'o_pH'#'o_CC'#
import matplotlib.pyplot as plt
x = df['i_C'].tolist()
#plt.imshow(np.array(df['i_pH'].tolist()))
#plt.show()
#plt.plot(df['i_C'].tolist(),df['i_pH'].tolist())
plt.plot(df['i_C'].tolist(),df[pr].tolist())
plt.show()
#%% LEARN
#ndf = df[['i_Ca','i_C','i_CH','i_CC','i_poros','i_pH']]

feature_cols = ['i_C']#,'i_poros''i_Ca',
X = df.loc[:,feature_cols]
print(X.shape)
#y = df.i_pH
y = df[pr]
print(y.shape)

from sklearn.model_selection import train_test_split
xTrain, xTest, yTrain, yTest = train_test_split(X, y, test_size = 0.2, random_state = 0)
#%% KNN
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler(feature_range=(0, 1))

t0 = time.time()
x_train_scaled = scaler.fit_transform(xTrain)
xTrain2 = pd.DataFrame(x_train_scaled)
x_test_scaled = scaler.fit_transform(xTest)
xTest2 = pd.DataFrame(x_test_scaled)
knn_scaling = time.time() - t0
print("Scaling in %.3f s" % knn_scaling)

from sklearn import neighbors
K = 4
model = neighbors.KNeighborsRegressor(n_neighbors = K, weights='distance')


t0 = time.time()
model.fit(xTrain2, yTrain)  #fit the model
knn_fit = time.time() - t0
print("Fitting in %.3f s" % knn_fit)
t0 = time.time()
pred=model.predict(xTest2) #make prediction on test set
knn_pr = time.time() - t0
print("Fitting in %.3f s" % knn_pr)

print(model.score(xTest2, yTest))
n2 = 500
temp = {'i_C': np.linspace(0.01, 0.1, n2)}
Xnew= pd.DataFrame(scaler.fit_transform(pd.DataFrame(data=temp)))
Ynew = model.predict(Xnew)
#model2 = neighbors.KNeighborsRegressor(n_neighbors = K, weights= 'uniform')
#model2.fit(xTrain2, yTrain)  #fit the model
#Ynew2 = model2.predict(Xnew)

plt.plot(Xnew.loc[:,0], Ynew, c='g', label='kNN' )
#plt.plot(Xnew.loc[:,0], Ynew, c='r', label='kNN2' )
plt.plot((df['i_C']*10).tolist(),df[pr].tolist(), '.', label = "o")
plt.xlabel('data')
plt.ylabel('target')
plt.legend()
plt.show()
#%% Kernel Ridge
'''
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import GridSearchCV


kr = GridSearchCV(KernelRidge(kernel='poly', gamma=0.1, ),
                  param_grid={"alpha": [1e1, 1e0, 0.1, 1e-2, 1e-3],
                              "gamma": np.logspace(-2, 2, 10)})

t0 = time.time()
kr.fit(xTrain, yTrain)
kr_fit = time.time() - t0
print("KRR complexity and bandwidth selected and model fitted in %.3f s"
      % kr_fit)

t0 = time.time()
y_kr = kr.predict(xTest)
kr_predict = time.time() - t0
print("KRR prediction for %d inputs in %.3f s"
      % (xTest.shape[0], kr_predict))
print(kr.score(xTest, yTest))
a = len(y_kr)
#plt.scatter(X.i_Ca, y, c='k', label='data', zorder=1,edgecolors=(0, 0, 0))
plt.scatter(xTest.i_C[:a], yTest[:a], c='r',
         label='Data')
plt.scatter(xTest.i_C[:a], y_kr[:a], c='g',
         label='KRR (fit: %.3fs, predict: %.3fs)' % (kr_fit, kr_predict))
plt.xlabel('data')
plt.ylabel('target')
plt.title('Kernel Ridge')
plt.legend()

# Visualize training and prediction time
plt.figure()

n2 = 1000
temp = {'i_C': np.linspace(0.01, 0.1, n2)}
Xnew= pd.DataFrame(data=temp)
Ynew = kr.predict(Xnew)

plt.plot(Xnew.i_C, Ynew, c='g',
         label='Gauss' )
plt.xlabel('data')
plt.ylabel('target')
plt.title('GPR')
plt.legend()
plt.show()
'''

#%% gauss

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import ConstantKernel
#gp_kernel = RBF(10, (5, 13))
gp_kernel = ConstantKernel(1000, (1, 5e5))*RBF(1000, (1, 5e5))
gpr = GaussianProcessRegressor(kernel=gp_kernel, 
                               n_restarts_optimizer=20, 
                               alpha=1e-7,  #noise
                               normalize_y=True)
stime = time.time()
gpr.fit(xTrain, yTrain)
print("Time for GPR fitting: %.3f" % (time.time() - stime))

stime = time.time()
y_gpr = gpr.predict(xTest, return_std=False)
print("Time for GPR prediction: %.3f" % (time.time() - stime))

print(gpr.score(xTest, yTest))

a = 100
#plt.scatter(X.i_Ca, y, c='k', label='data', zorder=1,edgecolors=(0, 0, 0))
plt.scatter(xTest.i_C[:a], yTest[:a], c='r',
         label='Data')
plt.scatter(xTest.i_C[:a], y_gpr[:a], c='g',
         label='Gauss' )
plt.xlabel('data')
plt.ylabel('target')
plt.title('GPR')
plt.legend()
plt.show()

# generate new test
n2 = 1000
temp = {'i_C': np.linspace(0.00, 0.1, n2)}
Xnew= pd.DataFrame(data=temp)
Ynew = gpr.predict(Xnew, return_std=False)

plt.plot(Xnew.i_C, Ynew, c='g',
         label='Gauss' )
plt.xlabel('data')
plt.ylabel('target')
plt.title('GPR')
plt.legend()
plt.show()

#%% POLYNOMIAL
'''
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import PolynomialFeatures


polynomial_features= PolynomialFeatures(degree=17)
x_poly = polynomial_features.fit_transform(xTrain)

model = LinearRegression()
model.fit(x_poly, yTrain)
y_poly_pred = model.predict(x_poly)


rmse = np.sqrt(mean_squared_error(yTrain,y_poly_pred))
r2 = r2_score(yTrain,y_poly_pred)
print(rmse)
print(r2)

n2 = 10000
temp = {'i_C': np.linspace(0.01, 0.1, n2)}
Xnew= polynomial_features.fit_transform(pd.DataFrame(data=temp))
Ynew = model.predict(Xnew)

plt.plot( Ynew, c='g',
         label='Poly' )
plt.xlabel('data')
plt.ylabel('target')
plt.title('Poly')
plt.legend()
plt.show()
'''
