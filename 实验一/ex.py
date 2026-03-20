import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# 设置中文字体（Windows用SimHei，Mac用Arial Unicode MS）
plt.rcParams['font.sans-serif'] = ['SimHei']  # 替换为你的系统支持的中文字体
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# 读取Excel文件
file_path = "data2.xls"
df = pd.read_excel(file_path, engine='xlrd', sheet_name='Sheet1')

# 绘制折线图
plt.figure(figsize=(10, 6))

# 绘制每条算法曲线
plt.plot(df["阶数"], df["平凡算法"], marker='o', linestyle='-', label="平凡算法")
plt.plot(df["阶数"], df["二重循环"], marker='s', linestyle='--', label="二重循环")
plt.plot(df["阶数"], df["多路链路"], marker='^', linestyle=':', label="多路链路")
plt.plot(df["阶数"], df["递归"], marker='d', linestyle='-.', label="递归")

# 设置中文标签和标题
plt.title("不同矩阵阶数下的算法执行时间", fontsize=14)
plt.xlabel("矩阵阶数", fontsize=12)
plt.ylabel("执行时间 (秒)", fontsize=12)
plt.legend(fontsize=10)  # 显示图例
plt.grid(True, linestyle='--', alpha=0.6)  # 虚线网格

# 优化刻度标签（避免拥挤）
plt.xticks(rotation=45)
plt.tight_layout()  # 自动调整布局
plt.show()