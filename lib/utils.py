"""
工具函数
"""
import os
import sys
from datetime import datetime


def ensure_dir(path):
    """确保目录存在"""
    os.makedirs(path, exist_ok=True)
    return path


def log(msg):
    """日志输出"""
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{ts}] {msg}")


def check_file(filepath, desc="文件"):
    """检查文件是否存在"""
    if not os.path.exists(filepath):
        print(f"错误: {desc}不存在: {filepath}")
        sys.exit(1)


def write_tsv(filepath, header, rows):
    """写入TSV文件"""
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write('\t'.join(str(x) for x in header) + '\n')
        for row in rows:
            f.write('\t'.join(str(x) for x in row) + '\n')


def read_lines(filepath):
    """读取文件所有行"""
    with open(filepath, 'r', encoding='utf-8') as f:
        return [line.strip() for line in f if line.strip()]
