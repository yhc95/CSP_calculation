import math
from datetime import datetime

def calculate_delta_comb(H1, C1, H2, C2, residue_type):
    """
    Calculate combined chemical shift perturbation (∆δ_comb) from two chemical shift states.
    
    Parameters:
    H1, C1 (float): Chemical shifts in 1H and 13C dimensions for state 1 (in ppm).
    H2, C2 (float): Chemical shifts in 1H and 13C dimensions for state 2 (in ppm).
    residue_type (str): 'aliphatic' or 'aromatic' to determine the normalization factor ω_C.
    
    Returns:
    tuple: (∆δ_comb, ∆δ_H, ∆δ_C) all in ppm.
    
    Formula:
    ∆δ_H = |H2 - H1|, ∆δ_C = |C2 - C1|
    ∆δ_comb = √(ω_H * (∆δ_H)^2 + ω_C * (∆δ_C)^2)
    where ω_H = 1.00, ω_C = 0.34 for aliphatic or 0.07 for aromatic.
    """
    # 计算化学位移变化
    delta_H = abs(H2 - H1)
    delta_C = abs(C2 - C1)
    
    # 设置权重因子
    omega_H = 1.00
    residue_type = residue_type.lower()
    if residue_type == 'aliphatic':
        omega_C = 0.34
    elif residue_type == 'aromatic':
        omega_C = 0.07
    else:
        raise ValueError("residue_type must be either 'aliphatic' or 'aromatic'")
    
    # 计算综合化学位移扰动
    delta_comb = math.sqrt(omega_H * (delta_H ** 2) + omega_C * (delta_C ** 2))
    return delta_comb, delta_H, delta_C

def save_results_to_file(H1, C1, H2, C2, residue_type, delta_H, delta_C, delta_comb, filename="nmr_results.txt"):
    """保存计算结果到文件"""
    with open(filename, 'a', encoding='utf-8') as f:
        f.write("="*50 + "\n")
        f.write(f"计算时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"状态1: 1H={H1} ppm, 13C={C1} ppm\n")
        f.write(f"状态2: 1H={H2} ppm, 13C={C2} ppm\n")
        f.write(f"残基类型: {residue_type}\n")
        f.write(f"Δδ_H = {delta_H:.4f} ppm\n")
        f.write(f"Δδ_C = {delta_C:.4f} ppm\n")
        f.write(f"Δδ_comb = {delta_comb:.4f} ppm\n")
        f.write("="*50 + "\n\n")
    return filename

def batch_calculate_from_file(filename, residue_type, output_file=None):
    """
    从文件批量计算化学位移扰动并保存结果
    
    文件格式应为（每行）：
    残基编号 H1 C1 H2 C2
    例如：
    A1 1.2 110 1.25 109.5
    A2 0.9 115 0.92 114.8
    """
    results = []
    with open(filename, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split()
                if len(parts) >= 5:
                    residue_id = parts[0]
                    H1, C1, H2, C2 = map(float, parts[1:5])
                    delta_comb, delta_H, delta_C = calculate_delta_comb(H1, C1, H2, C2, residue_type)
                    results.append({
                        'residue': residue_id,
                        'delta_H': delta_H,
                        'delta_C': delta_C,
                        'delta_comb': delta_comb
                    })
    
    # 如果提供了输出文件名，则保存结果
    if output_file:
        with open(output_file, 'w') as f:
            f.write("Residue\tΔδ_H(ppm)\tΔδ_C(ppm)\tΔδ_comb(ppm)\n")
            for res in results:
                f.write(f"{res['residue']}\t{res['delta_H']:.4f}\t{res['delta_C']:.4f}\t{res['delta_comb']:.4f}\n")
        print(f"Results saved to {output_file}")
    
    return results

# 主程序
if __name__ == "__main__":
    print("化学位移扰动计算")
    print("=" * 40)
    
    try:
        # 获取用户输入
        H1 = float(input("Enter 1H chemical shift for state 1 (ppm): "))
        C1 = float(input("Enter 13C chemical shift for state 1 (ppm): "))
        H2 = float(input("Enter 1H chemical shift for state 2 (ppm): "))
        C2 = float(input("Enter 13C chemical shift for state 2 (ppm): "))
        residue_type = input("Enter residue type (aliphatic/aromatic): ").strip()
        
        # 计算化学位移扰动
        delta_comb, delta_H, delta_C = calculate_delta_comb(H1, C1, H2, C2, residue_type)
        
        # 显示结果
        print("\n" + "="*40)
        print("计算结果:")
        print(f"Δδ_H = {delta_H:.4f} ppm")
        print(f"Δδ_C = {delta_C:.4f} ppm")
        print(f"Δδ_comb = {delta_comb:.4f} ppm")
        
        # 保存结果到文件
        filename = save_results_to_file(H1, C1, H2, C2, residue_type, delta_H, delta_C, delta_comb)
        print(f"\n结果已保存到文件: {filename}")
        
        # 使用示例2的数据进行验证（根据您的图片输入）
        print("\n" + "="*40)
        print("根据您的输入验证:")
        print(f"状态1: 1H = {H1} ppm, 13C = {C1} ppm")
        print(f"状态2: 1H = {H2} ppm, 13C = {C2} ppm")
        print(f"残基类型: {residue_type}")
        
        if residue_type.lower() == 'aliphatic':
            print("\n验证计算:")
            print(f"Δδ_H = |{H2} - {H1}| = {delta_H:.4f} ppm")
            print(f"Δδ_C = |{C2} - {C1}| = {delta_C:.4f} ppm")
            print(f"Δδ_comb = √(1.00 × ({delta_H:.4f})² + 0.34 × ({delta_C:.4f})²)")
            print(f"       = √(1.00 × {delta_H**2:.4f} + 0.34 × {delta_C**2:.4f})")
            print(f"       = √({delta_H**2:.4f} + {0.34*delta_C**2:.4f})")
            print(f"       = √({delta_H**2 + 0.34*delta_C**2:.4f})")
            print(f"       = {delta_comb:.4f} ppm")
        
        # 等待用户按键，防止窗口关闭
        input("\n按Enter键退出程序...")
        
    except ValueError as e:
        print(f"\n错误: {e}")
        print("请确保输入正确的数值和残基类型 (aliphatic/aromatic)")
        input("\n按Enter键退出程序...")