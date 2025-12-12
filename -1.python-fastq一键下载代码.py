#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# ----------------- 配置 -----------------
filereport_file = r"C:\Users\XZQ\Desktop\NCBI\Plant communications_RNA_seq\filereport_read_run_PRJNA1043393_tsv.txt"
download_dir = r"C:\Users\XZQ\Desktop\NCBI\Plant communications_RNA_seq\rna-seq-project\sra"
max_workers = 8  # 并发线程数，16核可以用 8~16 线程
chunk_size = 1024 * 1024  # 每次下载块大小 1MB

# 创建下载目录
os.makedirs(download_dir, exist_ok=True)

# ---------- 读取 filereport 文件 ----------
df = pd.read_csv(filereport_file, sep='\t', dtype=str)

# ---------- 定义下载函数 ----------
def download_file(url, out_file, base_name):
    """下载单个文件，带进度条"""
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            total = int(r.headers.get('content-length', 0))
            with open(out_file, 'wb') as f, tqdm(
                total=total, unit='B', unit_scale=True, desc=os.path.basename(out_file)
            ) as pbar:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))
    except Exception as e:
        print(f"下载失败: {url}")
        print(e)

# ---------- 准备下载任务 ----------
tasks = []
for idx, row in df.iterrows():
    urls = row['fastq_aspera'].split(';')
    base_name = row['sample_title'].replace(',', '_').replace(' ', '')
    for i, url in enumerate(urls, 1):
        if pd.isna(url) or url.strip() == "":
            continue
        url = url.replace("fasp.sra.ebi.ac.uk:", "https://ftp.sra.ebi.ac.uk")
        out_file = os.path.join(download_dir, f"{base_name}_{i}.fastq.gz")
        tasks.append((url, out_file, base_name))

# ---------- 并发下载 ----------
with ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(download_file, url, out_file, base_name) for url, out_file, base_name in tasks]
    for future in as_completed(futures):
        pass  # tqdm 会显示进度条

print("所有下载完成！")
