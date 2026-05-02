import subprocess
import os
import argparse

def find_center_atom(file='step5_1.gro', bottom=None):
    '''
    将指定的.gro文件中的距离蛋白几何中心最近的原子序号提取出来
    '''
    print(f"传递参数 file={file}; bottom={bottom}, type={type(bottom)}")
    aminoacids = [
        'ALA','ASN','CYS','GLN','GLY','ILE','LEU','MET',
        'PHE','PRO','SER','THR','TRP','TYR','VAL','ASP',
        'ARG','HIS','LYS','GLU',
        'PYL','SEC',    # 吡咯赖氨酸、 硒代半胱氨酸
        'HID','HSD', 	# δ-氮 (ND1) 质子化 HIS
        'HIE','HSE',     # ε-氮 (NE2) 质子化 HIS
        'HIP','HSP',    # 双质子化，带正电荷 HIS
        'ASH', # ASP 羧基未离子化，即带氢，电中性
        'GLH', # 谷氨酸羧基未离子化
        'CYX', # 形成二硫键（CYS-CYS）
        'CYN', # 去质子化CYS，带负电
        'LYN', # 去质子化赖氨酸，通常是负电
        'TYM', # 酪氨酸去质子化，负电状态
        'ARN', # 去质子化的精氨酸，负电状态
        ]
    # 定义几何中心的变量
    center_X, center_Y, center_Z = 0, 0, 0 
    i = 0
    position = []
    with open(file, "r", encoding='utf-8') as fo:
        try:
            next(fo)  # 跳过第1行
            next(fo)  # 跳过第2行
        except StopIteration:
            pass  # 文件不足两行时直接跳过
        if bottom:
            for line in fo:
                if line[5:8] == bottom:
                    break
                else:
                    # 蛋白原子的坐标在第4、5、6列
                    x, y, z, atom = float(line[20:28]), float(line[28:36]), float(line[36:44]), int(line[15:20])
                    center_X += x
                    center_Y += y
                    center_Z += z
                    i += 1
                    position.append([atom, x, y, z])
        else:
            for line in fo:
                if line[5:8] in aminoacids:
                    # 蛋白原子的坐标在第4、5、6列
                    x, y, z, atom = float(line[20:28]), float(line[28:36]), float(line[36:44]), int(line[15:20])
                    center_X += x
                    center_Y += y
                    center_Z += z
                    i += 1
                    position.append([atom, x, y, z])
                else:
                    break            
    # 计算几何中心的坐标
    center_X = center_X / i
    center_Y = center_Y / i
    center_Z = center_Z / i
    print(f'几何中心: {center_X, center_Y, center_Z}')
    dist = 10
    atom_info = ''
    # 再遍历一遍原子，寻找离几何中心最近的原子
    for entry in position:
        # 计算距离，小于dist就刷新atom_info和dist
        atom_center_dist = ((entry[1]-center_X)**2 + (entry[2]-center_Y)**2 + (entry[3]-center_Z)**2 )**0.5
        if atom_center_dist < dist:
            dist = atom_center_dist
            # 记录此原子信息
            atom_info = entry[0]
    print(f'最靠近几何中心的原子: {entry}')
    return atom_info


def add_first_line(file_path, new_line):
    '''
    在文件首行添加内容
    '''
    with open(file_path, 'r+', encoding='utf-8') as f:
        content = f.read()
        if new_line in content:
            print(f'[Read {file_path}] 无需更新 {file_path} 文件的[center]内容')
        else:
            f.seek(0)  # 将指针移动到文件开头
            f.write(new_line + '\n' + content)
            print(f'[Read {file_path}] 已更新 {file_path} 文件的[center]内容')
        f.close()


def gmx_trjconv_trans(start, end, prod_step, number):
    '''
    平移轨迹文件几何位置
    '''
    print('[gmx_trjconv_noline] 运行xtc轨迹平移')
    for i in range(start, end + 1):
        istee = f"{prod_step}_{i}"
        isteq = f"{istee}_{number}_center"        
        # 删除冗余文件（添加异常处理）
        try:
            for file in [f"{isteq}.xtc"]:
                if os.path.exists(file):
                    os.remove(file)
                    print(f'> 冗余文件 {file} 已删除')
        except Exception as e:
            print(f"文件删除失败: {e}")
            continue  # 根据需求决定是否终止循环
        # 执行命令（严格串行）
        try:
            # 使用 subprocess.run 替代 Popen 简化同步逻辑
            result = subprocess.run(
                [
                    "gmx", "trjconv",
                    "-s", f"{istee}.tpr",
                    "-f", f"{istee}.xtc",
                    "-n", "index.ndx",
                    "-o", f"{isteq}.xtc",
                    "-pbc", "atom", "-center",
                ],
                input="0\n3\n".encode(),  # 一次性发送输入
                stdout=subprocess.PIPE,     # 可选：捕获输出
                stderr=subprocess.PIPE,     # 可选：捕获错误
                check=True,               # <=== 加上这个！
            )
            print(result.stdout.decode())
            print(f'[gmx_trjconv] {isteq}.xtc 已写入\n'+'-'*10)
        
        except subprocess.CalledProcessError as e:
            print(f"命令执行失败，返回码 {e.returncode}")
            print("错误输出:", e.stderr.decode())
            continue  # 或终止循环（raise）
        except subprocess.TimeoutExpired:
            print("命令执行超时！")
            continue
        except Exception as e:
            print(f"未知异常: {e}")
            continue


def gmx_trjconv_noline(start, end, prod_step, number):
    '''
    消除 VMD 可视化的线框问题
    '''
    print('[gmx_trjconv_noline] 运行可视化线框消除')
    for i in range(start, end + 1):
        istee = f"{prod_step}_{i}"
        isteq = f"{istee}_{number}_center"
        istep = f"{isteq}_noline"
        # 删除冗余文件（添加异常处理）
        try:
            for file in [f"{istep}.xtc"]:
                if os.path.exists(file):
                    os.remove(file)
                    print(f'> 冗余文件 {file} 已删除')
        except Exception as e:
            print(f"文件删除失败: {e}")
            continue  # 根据需求决定是否终止循环
        try:
            result2 = subprocess.run(
            [
                "gmx", "trjconv", 
                "-f", f"{isteq}.xtc",
                "-o", f"{istep}.xtc",
                "-s", f"{istee}.tpr",
                "-pbc", "mol",
            ],
            input="0\n".encode(),  # 一次性发送输入
            stdout=subprocess.PIPE,     # 可选：捕获输出
            stderr=subprocess.PIPE,     # 可选：捕获错误
            check=True,               # <=== 加上这个！
            )
            print(result2.stdout.decode())
            print(f'[gmx_trjconv] {istep}.xtc 已写入\n'+'-'*10)
        
        except subprocess.CalledProcessError as e:
            print(f"命令执行失败，返回码 {e.returncode}")
            print("错误输出:", e.stderr.decode())
            continue  # 或终止循环（raise）
        except subprocess.TimeoutExpired:
            print("命令执行超时！")
            continue
        except Exception as e:
            print(f"未知异常: {e}")
            continue

def gmx_trjcat(start, end, prod_step, number):
    '''
    合并xtc轨迹文件
    '''
    out_xtc = f"{prod_step}_all.xtc"
    if os.path.exists(out_xtc):
        os.remove(out_xtc)
        print(f'> 冗余文件 {out_xtc} 已删除')
    entf = [f"{prod_step}_{j}_{number}_center_noline.xtc" 
            for j in range(start, end + 1)]
    print(f'[gmx_trjcat] 正在合并 \n - {'\n - '.join(entf)}')
    try:
        result3 = subprocess.run(
        ["gmx", "trjcat", "-f",] + entf + 
        ["-o", out_xtc, "-sort", "yes", "-cat", "-settime"],
        input=str("c\n"*int(1+end-start)).encode(),  # 一次性发送输入
        stdout=subprocess.PIPE,     # 可选：捕获输出
        stderr=subprocess.PIPE,     # 可选：捕获错误
        check=True,
        )
        print(result3.stdout.decode())
        print(f'[gmx_trjcat] {out_xtc} 已合并')
    
    except subprocess.CalledProcessError as e:
        print(f"命令执行失败，返回码 {e.returncode}")
        print("错误输出:", e.stderr.decode())
    except subprocess.TimeoutExpired:
        print("命令执行超时！")
    except Exception as e:
        print(f"未知异常: {e}")

gro_file = 'step5_1.gro'
step_start = 1
step_end = 5
xtc_step = 'step5'


if __name__ == '__main__':
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser(description="Run deal_xtc.py")
    # 添加命令行参数，并设置默认值
    parser.add_argument('--gro', type=str, default=gro_file, help=f".gro file (default: {gro_file})")
    parser.add_argument('--gro_until', type=str, default=None, help=f"End step5_*.xtc index (default: None)")
    parser.add_argument('--start', type=int, default=step_start, help=f"Start step5_*.xtc index (default: {step_start})")
    parser.add_argument('--end', type=int, default=step_end, help=f"End step5_*.xtc index (default: {step_end})")
    parser.add_argument('--prod_step', type=str, default=xtc_step, help=f"Production step (default: {xtc_step})")

    # 解析命令行参数
    args = parser.parse_args()
    
    # 输出最靠近几何中心的原子序号
    center_id = find_center_atom(file=args.gro, bottom=args.gro_until)
    # 在 index.ndx 文件首行添加新内容
    add_first_line("index.ndx", f'[ center ]\n{str(center_id)}\n')
    # 轨迹文件前处理
    gmx_trjconv_trans(start=args.start, end=args.end, prod_step=args.prod_step, number=center_id)
    gmx_trjconv_noline(start=args.start, end=args.end, prod_step=args.prod_step, number=center_id)
    gmx_trjcat(start=args.start, end=args.end, prod_step=args.prod_step, number=center_id)
    
    
    
    