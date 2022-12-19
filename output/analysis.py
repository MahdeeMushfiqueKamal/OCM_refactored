import pandas as pd
import math
import sklearn.metrics
df = pd.read_csv('query_result.csv')
print(df.head())

print("Checking if there exists any row for which truecount > estimated count")
print(df.loc[df['true_count'] > df['estimated_count']])

print("sparman correlation coffecient:")
sp_corr = df.corr(method='spearman')['true_count'].drop('true_count')
print(sp_corr)


rms_error = math.sqrt(sklearn.metrics.mean_squared_error(df['true_count'], df['estimated_count']))
print("RMS Error: ",rms_error)
