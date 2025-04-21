import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

def string_to_year(date):
    term, year = date.split(" ")
    year = float(year)
    if term.lower() == "fall":
        year += 11/12.
    if term.lower() == "spring":
        year += 5/12.
    if term.lower() == "summer":
        year +=8/12.
    return year

filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\GradStudentDatabase2.xlsx"
df = pd.read_excel(filepath)


b1 = df["TTD"]>3
df2 = df[b1]
'''
x=df2["GPA"]
y=df2["TTD"]
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
print(slope, intercept, r_value**2, p_value)
xsim = np.arange(3, 4.1, 0.1)
ysim = xsim * slope+intercept

plt.scatter(df2["GPA"], df2["TTD"])
plt.plot(xsim, ysim)
plt.ylabel("Time To Degree (Years)")
plt.xlabel("GPA")
plt.show()

exit()'''

x=df2["TTD"]
y=df2["Total"]
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
print(slope, intercept, r_value**2, p_value)
fit = np.polyfit(df2["TTD"], df2["Total"],1)
xsim = np.arange(4, 10)
ysim = xsim * fit[0]+fit[1]

plt.scatter(df2["TTD"], df2["Total"])
plt.plot(xsim, ysim)
plt.xlabel("Time To Degree (Years)")
plt.ylabel("Total Publications")
plt.show()

exit()
for index, row in df.iterrows():
    try:
        start = float(string_to_year(row["Start_Sem"]))
        end = float(string_to_year(row["End_Sem"]))

        ttd = end-start
    except Exception as e:
        ttd = ""
    if row["FinalDegree"] != "PhD":
        ttd = ""
    df.loc[index, "TTD"] = ttd

df.to_excel(filepath, index=False)