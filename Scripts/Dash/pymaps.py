import pandas as pd

df = df.readcsv("C:\Users\admin\University of Arizona\martylab - Documents\Presentations\Dash\Independent Proposal\Assets\IHME_AMR_BURDEN_2019_DATA_NUMBER\IHME_AMR_BURDEN_2019_NUMBER_Y2022M01D20.CSV")
values = df.loc[df['Staphylococcus aureus'] == 'row_value'].values.tolist()values = df.loc[df['column_name'] == 'row_value'].values.tolist()