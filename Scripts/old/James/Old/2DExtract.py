import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

#reading csv data
df = pd.read_csv('C:\Users\Lab.DESKTOP-HJRL65I\Documents\James_2D extracts\20190411_1PC40V_2Dextract_python.csv')
print(df.head(20))

#Create array of symbol and change
symbol = ((np.asarray(df['Symbol'])).reshape(5,4))
perchange = ((np.asarray(df['Change'])).reshape(5,4))

print(symbol)
print(perchange)

#Create a pivot table
result = df.pivot(index='Yrows',columns='Xcols',values='Change')
print(result)

#Array to label heatmap
labels = (np.asarray(["{0} \n {1:.2f}".format(symb,value) for symb, value in zip(symbol.flatten(),perchange.flatten())])).reshape(5,4)

#Define the plot
fig, ax = plt.subplots(figsize=(12,7))

#Add title to heat map
title = "PC1 Batch 2 40V 2D Grid"

#Set font size and distance of title
plt.title(title,fontsize=18)
ttl = ax.title
ttl.set_position([0.5,1.05])

#Hide ticks for X and Y axis
ax.set_xticks([])
ax.set_yticks([])

ax.axis('off')

#Heatmap function from seaborn
sns.heatmap(result,annot=labels,fmt="",cmap='viridis',linewidths=0.30,ax=ax)

plt.show()


