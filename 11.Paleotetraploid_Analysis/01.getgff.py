import sys
import pandas as pd

# 尝试跳过注释行和处理数据
try:
    data = pd.read_csv(sys.argv[1], sep="\t", header=None, skiprows=3, comment='#')
except pd.errors.ParserError as e:
    print(f"解析错误：{e}")
    sys.exit(1)

# 筛选 'mRNA' 类型的行
data = data[data[2] == 'mRNA']

# 选择需要的列
data = data.loc[:, [0, 8, 3, 4, 6]]

# 确保第0列是字符串类型，并处理NaN
data[0] = data[0].fillna('').astype(str)  # 用空字符串替换NaN，然后转换为字符串

# 处理第8列，提取 '=' 后的值
data[8] = data[8].str.split(':|=|;', expand=True)[1]

# 处理第0列，去掉 'Chr_' 前缀
data[0] = data[0].str.replace('Chr_0?', '', regex=True)

# 保存为新的文件
data.to_csv(sys.argv[2], sep="\t", header=None, index=False)
