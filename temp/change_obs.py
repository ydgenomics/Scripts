import pandas as pd

df= pd.DataFrame(columns=['old', 'new'])

df = pd.DataFrame({
    'old': ['D1', 'D2', 'D3', 'D5', 'D10', 'D30'],
    'new': [1, 2, 3, 5, 10, 30]
})

df.to_csv('old_new.csv', index=False)

df = pd.read_csv('old_new.csv')

time_map = dict(zip(df['old'], df['new']))
print(time_map)

# 建立映射
time_map = {'D1': 1, 'D2': 2, 'D3': 3, 'D5': 5, 'D10': 10, 'D30': 30}

# 一次性替换
adata.obs['time_num'] = adata.obs['time'].map(time_map)