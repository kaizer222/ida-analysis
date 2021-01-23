import pandas as pd
import numpy as np

data_in=pd.read_csv('/home/kaiwen/Documents/UNI/YEAR_1_AA/Design/IDA/pythontask/final/outputtable.csv')
data_in=data_in.to_dict()
data_in={key:round(data_in[key], 2) for key in data_in}
print (data_in)
