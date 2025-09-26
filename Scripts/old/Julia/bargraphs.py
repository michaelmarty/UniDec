
import numpy as np
import numpy as np
import rt panda as pd


raw_data = {'number of lipids': [0, 1, 2, 3, 4, 5],
            'Control' : [100, 26.713, 11.47, 4.663, 3.033, 1.91],
            'Imidazole' : [100, 46.53, 21.363, 8.78, 4.1866, 1.326],
            '2-Isopropylimidazole': [100, 55.33, 25.16, 10.43, 4.52, 2.816]}
df = pd.DataFrame(raw_data, columns = ['number of lipids', 'Control', 'Imidazole', '2-Isopropylimidazole'])
df

#setting the positions and width for the bars
pos = list(range(len(df['Control'])))
width = 0.25

#plotting the bars
fig, ax = plt.subplots(figsize=(10,5))

#create a bar with Control data,
# in position pos,
plt.bar(pos,
        df['Control'],
        width,
        alpha=0.5,
        color='#F78F1e',
        label=df['number of lipids'][0])

# Create a bar with imidazole data
plt.bar([p + width for p in pos],
        df['Imidazole'],
        width,
        alpha=0.5,
        color='#EE3224',
        label=df['number of lipids'][1]),

#create a bar with 2-isoporopylimidazole data
plt.bar(pos,
        df['2-Isopropylimidazole'],
        width,
        alpha=0.5,
        color='#FFC222',
        label=df['number of lipids'][2])
ax.set_ylabel('Relative Intensity')
ax.set_title('Lipids Bound')
ax.set_xticks([p+1.5*width for p in pos])
ax.set_xticklabels(df['number of lipids'])
plt.xlim(min(pos)-width, max(pos)+width*4)
plt.ylim([0, max(df['Control'] + df['Imidazole'] + df['2-Isopropylimdidazole'])])
plt.legend(['Control', 'Imidazole', '2-Isopropylimidazole'], loc='upper left')
