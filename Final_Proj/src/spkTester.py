import subprocess
import time
import pandas as pd

parse = lambda x: "".join(x.split()) # Remove Whitespaces.
PYTHON_IMPL = "python3 spkmeans.py "
C_IMPL = "./spkmeans "
txt= ".txt"
csv = ".csv"
input = "0 spk inputs/input"
outputc = "inputs/outputc"
outputp = "inputs/outputp"
list = list(map(str, range(2,51)))
tests = [input + index + csv for index in list]
outputs_c = [outputc + index + txt for index in list]
outputs_p = [outputp + index + txt for index in list]
count = 2
df = pd.DataFrame()

# c
for output,test in zip(outputs_c, tests):

	c_start = time.time()
	c_output = subprocess.check_output(f"{C_IMPL}{test}", shell=True).decode()
	c_end = time.time()
	bool = ' '
	with open(output, "r") as f:
		my_output = f.read()

	if parse(c_output) != parse(my_output):
		bool = 'F'

	dict = {'C run time':c_end - c_start,'C result':bool,'file name':f'input{count}'}
	df = df.append(dict,ignore_index=True)
	count = count+1

p_runtime = []
p_result = []
count = 2

# python
for output,test in zip(outputs_p, tests):

	python_start = time.time()
	python_output = subprocess.check_output(f"{PYTHON_IMPL}{test}", shell=True).decode()
	python_end = time.time()
	bool = ' '
	with open(output, "r") as f:
		my_output = f.read()

	if parse(python_output) != parse(my_output):
		bool = 'F'

	p_runtime.append(python_end - python_start)
	p_result.append(bool)
	count = count+1

df['python result'] = p_result
df['python run time'] = p_runtime
df = df.reindex(columns=['file name','C result','python result','C run time','python run time'])
pd.DataFrame(df).to_csv('spk tester result.csv', index=False)
print('End')
