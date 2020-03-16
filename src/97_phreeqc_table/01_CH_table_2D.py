# -*- coding: utf-8 -*-

import IPhreeqcPy
import numpy as np
import matplotlib.pyplot as plt
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
x1 = np.linspace(0.0, 0.02, 100)# C
x2 = np.linspace(0.0, 0.02, 100)#Ca
c,ca= map(np.ravel, np.meshgrid(x1,x2))
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
import pandas as pd
#'''
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
'''
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
       'o_sol':[row[0] for row in so2][1:-1],
       'o_pH': [row[1] for row in so2][1:-1],
       'o_pe': [row[2] for row in so2][1:-1],
       'o_Alk': [row[3] for row in so2][1:-1],
       'o_C': [row[4] for row in so2][1:-1],
       'o_Ca': [row[5] for row in so2][1:-1],
       'o_O': [row[6] for row in so2][1:-1],
       'o_H': [row[7] for row in so2][1:-1],
       'o_CH': [row[8] for row in so2][1:-1],
       'o_SI_CH': [row[9] for row in so2][1:-1],
       'o_CC': [row[10] for row in so2][1:-1],
       'o_SI_CC': [row[11] for row in so2][1:-1],
       'o_poros': [row[12] for row in so2][1:-1],
       'o_H20': [row[13] for row in so2][1:-1]}
'''
df = pd.DataFrame(data=res)
#%% Plot
pr = 'o_CC'#'o_CC'#'o_pH'
feature_cols = ['i_C', 'i_Ca']#,'i_poros''i_Ca',
X = df.loc[:,feature_cols]
print(X.shape)
y = df[pr]
print(y.shape)

from sklearn.model_selection import train_test_split
xTrain, xTest, yTrain, yTest = train_test_split(X, y, test_size = 0.2, random_state = 0)
#%% KNN
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler(feature_range=(0, 1))

t0 = time.time()
x_train_scaled = scaler.fit_transform(xTrain)
xTrain_s = pd.DataFrame(x_train_scaled)
x_test_scaled = scaler.fit_transform(xTest)
xTest_s = pd.DataFrame(x_test_scaled)
knn_scaling = time.time() - t0
print("Scaling in %.3f s" % knn_scaling)

from sklearn import neighbors
K = 6
model = neighbors.KNeighborsRegressor(n_neighbors = K, weights='distance')


t0 = time.time()
model.fit(xTrain_s, yTrain)  #fit the model
knn_fit = time.time() - t0
print("Fitting in %.3f s" % knn_fit)
t0 = time.time()
pred=model.predict(xTest_s) #make prediction on test set
knn_pr = time.time() - t0
print("Fitting in %.3f s" % knn_pr)

print(model.score(xTest_s, yTest))

a = 10
plt.scatter(xTest.i_C[1:a], yTest[1:a], c='r',
         label='Data')
plt.scatter(xTest.i_C[1:a], pred[1:a], c='g',
         label='kNN' )
plt.xlabel('data')
plt.ylabel('target')
plt.title('kNN')
plt.legend()
plt.show()


#%% Add more
n =100
x1 = np.linspace(0.02, 0.04, 100)# C
x2 = np.linspace(0.02, 0.04, 100)#Ca
c,ca= map(np.ravel, np.meshgrid(x1,x2))
water = np.ones(np.shape(c))
ch = np.zeros(np.shape(c))
cc = np.zeros(np.shape(c))
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
import pandas as pd
#'''
res2 = {'i_sol':[row[0] for row in so1][1:-1],
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

df2 = pd.DataFrame(data=res2)
X2 = df2.loc[:,feature_cols]
y2 = df2[pr]

from sklearn.model_selection import train_test_split
xTrain2, xTest2, yTrain2, yTest2 = train_test_split(X2, y2, test_size = 0.2, random_state = 0)

x_train_scaled = scaler.fit_transform(xTrain2)
xTrain2_s = pd.DataFrame(x_train_scaled)
x_test_scaled = scaler.fit_transform(xTest2)
xTest2_s = pd.DataFrame(x_test_scaled)
model.fit(xTrain2_s, yTrain2)  #fit the model
pred2=model.predict(xTest2_s) #make prediction on test set
print(model.score(xTest_s, yTest))
print(model.score(xTest2_s, yTest2))


a = 10
plt.scatter(xTest2.i_C[1:a], yTest2[1:a], c='r',
         label='Data')
plt.scatter(xTest2.i_C[1:a], pred2[1:a], c='g',
         label='kNN' )
plt.xlabel('data')
plt.ylabel('target')
plt.title('kNN')
plt.legend()
plt.show()