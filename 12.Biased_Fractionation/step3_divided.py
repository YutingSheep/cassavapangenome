import pandas as pd

# 读取数据（新增第一列）
df = pd.read_csv('step3_input.txt', sep='\s+', header=None, 
                names=['new_col', 'col1', 'col2', 'col3', 'block'])  # 新增列名

# 修改分组聚合逻辑（新增列处理）
agg_df = df.groupby('block').agg({
    'new_col': ['first', 'last'],  # 新增列的聚合
    'col1': ['first', 'last'],
    'col2': ['first', 'last'],
    'col3': ['first', 'last']
}).reset_index()

# 调整展平后的列名（增加新列）
agg_df.columns = ['block', 
                 'first_new', 'last_new',  # 新增列首尾
                 'first_col1', 'last_col1',
                 'first_col2', 'last_col2',
                 'first_col3', 'last_col3']

# 重组字段顺序（新列排最前）
final_df = agg_df[[
    'first_new', 'first_col1', 'first_col2', 'first_col3',
    'last_new', 'last_col1', 'last_col2', 'last_col3',
    'block'  # 保持block在末尾
]]

# 输出文件（格式不变）
final_df.to_csv('step3_output.txt', sep='\t', index=False, header=False)

