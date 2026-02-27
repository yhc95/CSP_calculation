import numpy as np
import math

# ==================== 硬编码数据（源自文档 data.xlsx Sheet1）====================
AMINO_ACID_DATA = [
    ["Ala", "β", 1.353, 0.276, 19.028, 2.911],
    ["Ile", "δ", 0.674, 0.326, 13.489, 3.318],
    ["Ile", "γ", 0.77, 0.302, 17.601, 3.15],
    ["Leu", "δ1", 0.747, 0.327, 24.654, 2.01],
    ["Leu", "δ2", 0.729, 0.383, 24.13, 2.085],
    ["Thr", "γ", 1.139, 0.273, 21.592, 1.855],
    ["Val", "γ1", 0.819, 0.328, 21.534, 2.344],
    ["Val", "γ2", 0.802, 0.419, 21.346, 2.44],
    ["Met", "ε", 1.787, 1.469, 17.238, 3.992],
    ["Phe", "δ1", 7.04, 0.393, 131.207, 5.749],
    ["Phe", "δ2", 7.041, 0.405, 131.344, 4.401],
    ["Phe", "ε1", 7.065, 0.444, 130.351, 5.705],
    ["Phe", "ε2", 7.064, 0.439, 130.544, 4.096],
    ["Phe", "ζ", 6.995, 0.702, 129.035, 4.103],
    ["Trp", "δ", 7.128, 0.36, 126.35, 4.321],
    ["Trp", "ε3", 7.302, 0.514, 120.216, 5.345],
    ["Trp", "η2", 6.958, 0.447, 123.566, 4.837],
    ["Trp", "ζ2", 7.27, 0.404, 114.076, 4.446],
    ["Trp", "ζ3", 6.854, 0.464, 121.186, 4.493],
    ["Tyr", "δ1", 6.921, 0.367, 132.388, 5.213],
    ["Tyr", "δ2", 6.918, 0.371, 132.395, 5.177],
    ["Tyr", "ε1", 6.691, 0.303, 117.748, 3.941],
    ["Tyr", "ε2", 6.692, 0.313, 117.79, 3.221],
    ["His", "ε2", 7.138, 3.154, 119.91, 5.5],
    ["His", "ε", 7.841, 2.454, 137.258, 5.512],
]

# ==================== 核心计算函数 ====================
def probability_density(a, b, mu_H, sigma_H, mu_C, sigma_C):
    """计算二维正态分布的概率密度值 f(a,b)"""
    exponent = -0.5 * (((a - mu_H) / sigma_H) ** 2 + ((b - mu_C) / sigma_C) ** 2)
    coefficient = 1.0 / (2.0 * math.pi * sigma_H * sigma_C)
    return coefficient * math.exp(exponent)

def calculate_amino_acid_probabilities(a, b, verbose=True):
    """
    计算给定化学位移 (a, b) 的热点残基属于每种氨基酸类型的概率。
    参数 verbose: 是否打印详细结果
    """
    density_dict = {}
    amino_acids = list(set([entry[0] for entry in AMINO_ACID_DATA]))
    
    for aa in amino_acids:
        aa_entries = [entry for entry in AMINO_ACID_DATA if entry[0] == aa]
        densities = []
        for entry in aa_entries:
            _, _, mu_H, sigma_H, mu_C, sigma_C = entry
            density = probability_density(a, b, mu_H, sigma_H, mu_C, sigma_C)
            densities.append(density)
        f_x = max(densities) if densities else 0.0
        density_dict[aa] = f_x
    
    total_density = sum(density_dict.values())
    prob_dict = {}
    for aa, f_x in density_dict.items():
        prob = f_x / total_density if total_density > 0 else 0.0
        prob_dict[aa] = prob
    
    sorted_probs = sorted(prob_dict.items(), key=lambda item: item[1], reverse=True)
    
    if verbose:
        print(f"\n对于化学位移 (氢={a}, 碳={b}) 的热点残基：")
        print("-" * 70)
        print(f"{'氨基酸类型':<8} | {'概率 (P)':<12} | {'概率密度 (f)':<20} | {'主要参考位置':<10}")
        print("-" * 70)
        
        for aa, prob in sorted_probs:
            max_density_entry = None
            max_density = -1
            for entry in AMINO_ACID_DATA:
                if entry[0] == aa:
                    _, pos, mu_H, sigma_H, mu_C, sigma_C = entry
                    density = probability_density(a, b, mu_H, sigma_H, mu_C, sigma_C)
                    if density > max_density:
                        max_density = density
                        max_density_entry = entry
            
            main_position = max_density_entry[1] if max_density_entry else "N/A"
            print(f"{aa:<10} | {prob:<12.6f} | {density_dict[aa]:<20.6e} | {main_position:<10}")
        
        print("-" * 70)
        most_likely_aa, highest_prob = sorted_probs[0]
        print(f"最可能的氨基酸类型是: {most_likely_aa} (概率 = {highest_prob:.4f})\n")
    
    return prob_dict, sorted_probs

def interactive_mode():
    """交互模式：用户输入化学位移值进行分析"""
    print("=" * 60)
    print("氨基酸类型概率分析工具")
    print("=" * 60)
    
    while True:
        print("\n请选择操作:")
        print("1. 分析单个热点残基")
        print("2. 批量分析多个热点残基")
        print("3. 退出")
        
        choice = input("\n请输入选项 (1/2/3): ").strip()
        
        if choice == '1':
            print("\n" + "=" * 50)
            print("单个热点残基分析")
            print("=" * 50)
            try:
                a = float(input("请输入氢维度化学位移 (例如 7.0): "))
                b = float(input("请输入碳维度化学位移 (例如 130.0): "))
                calculate_amino_acid_probabilities(a, b)
            except ValueError:
                print("错误：请输入有效的数值！")
        
        elif choice == '2':
            print("\n" + "=" * 50)
            print("批量分析多个热点残基")
            print("=" * 50)
            print("输入格式：每行输入一对化学位移，格式为 '氢位移 碳位移'")
            print("例如: 7.0 130.0")
            print("输入 'done' 结束输入")
            
            hotspot_list = []
            while True:
                user_input = input("输入化学位移对 (或 'done'): ").strip()
                if user_input.lower() == 'done':
                    break
                try:
                    parts = user_input.split()
                    if len(parts) == 2:
                        a_val = float(parts[0])
                        b_val = float(parts[1])
                        hotspot_list.append((a_val, b_val))
                        print(f"已添加: ({a_val}, {b_val})")
                    else:
                        print("错误：请输入两个数值，用空格分隔")
                except ValueError:
                    print("错误：请输入有效的数值！")
            
            if hotspot_list:
                print(f"\n开始分析 {len(hotspot_list)} 个热点残基...")
                for i, (a_val, b_val) in enumerate(hotspot_list):
                    print(f"\n>>> 热点残基 #{i+1} (氢={a_val}, 碳={b_val}):")
                    prob_dict_i, sorted_probs_i = calculate_amino_acid_probabilities(a_val, b_val, verbose=True)
            else:
                print("没有输入任何热点残基数据。")
        
        elif choice == '3':
            print("感谢使用！程序退出。")
            break
        
        else:
            print("无效的选项，请重新输入。")

# ==================== 主程序 ====================
if __name__ == "__main__":
    interactive_mode()