import pandas as pd
import numpy as np

bus = pd.read_excel("ac_case25.xlsx", sheet_name='bus')
branch = pd.read_excel("ac_case25.xlsx", sheet_name='branch')
transformer = pd.read_excel("ac_case25.xlsx", sheet_name='transformer')
generator = pd.read_excel("ac_case25.xlsx", sheet_name='generator')
param = pd.read_excel("ac_case25.xlsx", sheet_name='param', )

# Cal unit to pu
print(param)