#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import re
import sys
from pathlib import Path

# 输入 GFF3 文件完整路径
gff_file = sys.argv[1]
gff_path = Path(gff_file)

# 输出文件：和输入同目录，命名为 TAIR10.62_genes.csv
output_file = gff_path.parent / "TAIR10.62_genes.csv"

with open(gff_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    writer = csv.writer(outfile)
    # 写表头
    writer.writerow([
        'feature_type', 'featureID', 'ParentID', 'symbol', 'biotype',
        'region', 'region_type', 'start', 'end', 'strand', 'length', 'source', 'description'
    ])

    region_types = {}  # 记录染色体还是 fragment

    for line in infile:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        seqid, source, ftype, start, end, score, strand, phase, attr = fields

        # 只提取 gene 类型
        if ftype.lower() != "gene":
            continue

        # region_type 判断
        if seqid not in region_types:
            if ftype.lower() == 'chromosome' or seqid.lower().startswith('chr'):
                region_types[seqid] = 'chromosome'
            else:
                region_types[seqid] = 'fragment'

        region_type = region_types[seqid]

        # 提取属性字段
        def get_attr(pattern):
            m = re.search(pattern, attr)
            return m.group(1) if m else '-'

        featureID = get_attr(r'ID=([^;]+)')
        ParentID = get_attr(r'Parent=([^;]+)')
        symbol = get_attr(r'Name=([^;]+)')
        biotype = get_attr(r'(biotype|gene_biotype)=([^;]+)')  # ensembl 或 ncbi
        description = get_attr(r'description=([^;]+)')

        start, end = int(start), int(end)
        length = end - start + 1

        writer.writerow([
            ftype, featureID, ParentID, symbol, biotype,
            seqid, region_type, start, end, strand, length, source, description
        ])

print(f"CSV file generated: {output_file}")
