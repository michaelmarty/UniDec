import pandas as pd
import numpy as np
from sklearn.cross_decomposition import PLSRegression
import matplotlib.pyplot as plt

filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\GradStudentDatabase3.xlsx"
df = pd.read_excel(filepath)
df = df[df["FinalDegree"]=="PhD"]

X = df.drop("Total", axis="columns")
#X = X.drop("Start", axis="columns")
#X = X.drop("End", axis="columns")
X = pd.get_dummies(X).fillna("missing")

Y = df['Total']

pls2 = PLSRegression(n_components=2, scale=False)
scores= pls2.fit_transform(X,Y)

xw = pls2.x_weights_

xpos=np.arange(0, len(xw))
plt.subplot(121)
plt.bar(xpos, xw[:, 1])
plt.xticks(xpos, X.keys())
plt.subplot(122)
plt.scatter(scores[0][:,0], scores[0][:,1], c=X["Gender_Male"]+2*X["FinalDegree_PhD"])
plt.show()

"""
famd = FAMD(n_components=2, n_iter=3, random_state=101)
famd.fit(df2)


fig = plt.figure()
ax = fig.add_axes()

ax = famd.plot_row_coordinates(df, ax=ax, x_component=0, y_component=1, color_labels=df2["Gender"], ellipse_outline=False, ellipse_fill=True, show_points=True)

plt.show()
"""