import subprocess
import time
import pandas as pd
import numpy as np
from contextlib import redirect_stdout

parse = lambda x: "".join(x.split()) # Remove Whitespaces.

def wam(datapoints,n):
    w = np.empty((n,n))
    for i in range(0,n):
        w[i][i] = 0
        for j in range(i+1, n):
            sum = np.sum((datapoints[i] -datapoints[j])**2)
            w[i][j] = np.exp(-np.sqrt(sum)/2)
            w[j][i] = w[i][j]
    return w

def ddg(W):
    d = np.sum(W,axis=1)
    d = np.diag(d)
    return d

def Lnorm(W,D,n):
    for i in range(n):
         D[i][i] = 1/np.sqrt(D[i][i])
    l = np.eye(n)-(D@W@D)
    return l

def print_matrix(matrix,n):
    for i in range(0,n):
        for j in range(0,n-1):
            num = matrix[i][j]
            if 0>num and -0.00005<num :
                num = 0
            print(f'{num:.4f},',end="")
        num = matrix[i][n-1]
        if 0>num and -0.00005<num :
            num = 0
        if (i != n-1):
            print(f'{num:.4f}')
        else:
            print(f'{num:.4f}',end="")



PYTHON_IMPL = "python3 spkmeans.py "
C_IMPL = "./spkmeans "
csv = ".csv"
input = "inputs/input"
list = list(map(str, range(2, 51)))
tests = [input + index + csv for index in list]
count = 2
df_wam = pd.DataFrame()
df_ddg = pd.DataFrame()
df_lnorm = pd.DataFrame()

for test in tests:

            datapoints = pd.read_csv(test,header=None).to_numpy()
            n,d = datapoints.shape

            ####wam
            c_start = time.time()
            c_output = subprocess.check_output(f"{C_IMPL} 0 wam {test}", shell=True).decode()
            c_end = time.time()
            p_start = time.time()
            p_output = subprocess.check_output(f"{PYTHON_IMPL} 0 wam {test}", shell=True).decode()
            p_end = time.time()
            print("\033c", end="")
            w = wam(datapoints, n)

            with open('my_output.txt', 'w') as f:
                with redirect_stdout(f):
                    print_matrix(w, n)

            bool1 = ' '
            bool2 = ' '
            with open('my_output.txt', "r") as f:
                my_output = f.read()

            if parse(c_output) != parse(my_output):
                bool1 = 'F'
            if parse(p_output) != parse(my_output):
                bool2 = 'F'

            dict = {'file name':f'input{count}','C result':bool1,'python result':bool2,'C run time':c_end - c_start,
                    'python run time':p_end - p_start}
            df_wam = df_wam.append(dict,ignore_index=True)

            ####ddg
            c_start = time.time()
            c_output = subprocess.check_output(f"{C_IMPL} 0 ddg {test}", shell=True).decode()
            c_end = time.time()
            p_start = time.time()
            p_output = subprocess.check_output(f"{PYTHON_IMPL} 0 ddg {test}", shell=True).decode()
            p_end = time.time()
            print("\033c", end="")
            d = ddg(w)

            with open('my_output.txt', 'w') as f:
                with redirect_stdout(f):
                    print_matrix(d, n)

            bool1 = ' '
            bool2 = ' '
            with open('my_output.txt', "r") as f:
                my_output = f.read()

            if parse(c_output) != parse(my_output):
                bool1 = 'F'
            if parse(p_output) != parse(my_output):
                bool2 = 'F'

            dict = {'file name':f'input{count}','C result':bool1,'python result':bool2,'C run time':c_end - c_start,
                    'python run time':p_end - p_start}
            df_ddg = df_ddg.append(dict,ignore_index=True)

            ####Lnorm
            c_start = time.time()
            c_output = subprocess.check_output(f"{C_IMPL} 0 lnorm {test}", shell=True).decode()
            c_end = time.time()
            p_start = time.time()
            p_output = subprocess.check_output(f"{PYTHON_IMPL} 0 lnorm {test}", shell=True).decode()
            p_end = time.time()
            print("\033c", end="")
            lnorm = Lnorm(w,d,n)


            with open('my_output.txt', 'w') as f:
                with redirect_stdout(f):
                    print_matrix(lnorm,n)

            bool1 = ' '
            bool2 = ' '
            with open('my_output.txt', "r") as f:
                my_output = f.read()

            if parse(c_output) != parse(my_output):
                bool1 = 'F'
            if parse(p_output) != parse(my_output):
                bool2 = 'F'

            dict = {'file name':f'input{count}','C result':bool1,'python result':bool2,'C run time':c_end - c_start,
                    'python run time':p_end - p_start}
            df_lnorm = df_lnorm.append(dict,ignore_index=True)

            count = count+1

df_wam = df_wam.reindex(columns=['file name','C result','python result','C run time','python run time'])
pd.DataFrame(df_wam).to_csv('wam tester result.csv', index=False)
df_ddg = df_ddg.reindex(columns=['file name','C result','python result','C run time','python run time'])
pd.DataFrame(df_ddg).to_csv('ddg tester result.csv', index=False)
df_lnorm = df_lnorm.reindex(columns=['file name','C result','python result','C run time','python run time'])
pd.DataFrame(df_lnorm).to_csv('Lnorm tester result.csv', index=False)
print('End')

