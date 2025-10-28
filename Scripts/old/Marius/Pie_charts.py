import matplotlib.pyplot as plt

labels = 'POPC', 'SOPE','POPS',# 'SSM'
colors = [ "crimson", "turquoise","palegreen"]#"lavender" ]
sizes = [45,40,15]
explode = (0,0,0)  # only "explode" the 2nd slice (i.e. 'Hogs')


fig1, ax1 = plt.subplots()
ax1.pie(sizes, colors= colors, explode=explode, labels=labels, autopct='%1.1f%%',
        shadow=True, startangle=90,radius= 5 , textprops={'fontsize': 14})

ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

plt.show()