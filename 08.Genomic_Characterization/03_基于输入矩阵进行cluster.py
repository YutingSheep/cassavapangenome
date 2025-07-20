import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import silhouette_score
import argparse

def load_transposed_data(file_path):
    """加载转置后的变异矩阵数据"""
    df = pd.read_csv(file_path, index_col=0)
    print("\n加载的转置数据矩阵 (变异位点×样本):")
    print(df.head())
    return df.T  # 转置回样本×变异位点格式

def calculate_distance(df):
    """计算样本间距离矩阵"""
    dist_matrix = pdist(df.values, metric='hamming')
    dist_matrix = squareform(dist_matrix)
    print("\n样本间距离矩阵 (前5×5):")
    print(pd.DataFrame(dist_matrix, index=df.index, columns=df.index).round(3).iloc[:5,:5])
    return dist_matrix

def hierarchical_clustering(dist_matrix, sample_names):
    """执行层次聚类"""
    Z = linkage(dist_matrix, method='ward')
    
    #plt.figure(figsize=(12, 6))
    dendrogram(Z, labels=sample_names, leaf_rotation=90)
    #plt.title('样本层次聚类树状图')
    #plt.xlabel('样本ID')
    #plt.ylabel('距离')
    #plt.tight_layout()
    #plt.show()
    
    return Z

def determine_clusters(Z, df):
    """通过轮廓系数确定最佳聚类数"""
    max_clusters = min(10, len(df)-1)
    best_score = -1
    best_n = 2
    
    for n in range(2, max_clusters+1):
        clusters = fcluster(Z, t=n, criterion='maxclust')
        score = silhouette_score(df, clusters, metric='hamming')
        print(f"聚类数={n} 轮廓系数={score:.3f}")
        if score > best_score:
            best_score = score
            best_n = n
    
    print(f"\n建议聚类数: {best_n} (轮廓系数={best_score:.3f})")
    return best_n

def renumber_clusters_by_size(clusters):
    """按聚类大小重新编号（从大到小）"""
    # 统计每个原始簇的样本数
    cluster_counts = pd.Series(clusters).value_counts()
    
    # 创建从原始簇号到新编号的映射（按样本数降序）
    size_rank = cluster_counts.rank(method='first', ascending=False).astype(int)
    mapping = {k: size_rank[k] for k in cluster_counts.index}
    
    # 应用新编号
    new_clusters = [mapping[c] for c in clusters]
    return new_clusters

def main(input_file, output_file):
    # 加载数据
    df = load_transposed_data(input_file)
    
    # 计算距离矩阵
    dist_matrix = calculate_distance(df)
    
    # 层次聚类
    Z = hierarchical_clustering(dist_matrix, df.index)
    
    # 确定最佳聚类数
    best_n = determine_clusters(Z, df)
    
    # 获取初始聚类结果
    initial_clusters = fcluster(Z, t=best_n, criterion='maxclust')
    
    # 按聚类大小重新编号
    final_clusters = renumber_clusters_by_size(initial_clusters)
    
    # 创建结果DataFrame
    result = pd.DataFrame({'sample': df.index, 
                          'initial_cluster': initial_clusters,
                          'cluster': final_clusters})
    
    # 打印聚类统计信息
    print("\n聚类大小统计:")
    print(result['cluster'].value_counts().sort_index())
    
    print("\n最终聚类结果 (已按聚类大小重新编号):")
    print(result[['sample', 'cluster']])
    
    # 保存结果
    result.to_csv(output_file, index=False)
    print(f"\n聚类结果已保存到 {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform hierarchical clustering on genetic variant data.')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file path')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file path')
    
    args = parser.parse_args()
    
    main(args.input, args.output)