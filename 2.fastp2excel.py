import json
import os
import pandas as pd

input_dir = "qc/cleanFq"  # JSON 文件所在目录
output_file = "qc/cleanFq_basic_stats.xlsx"

samples = []
raw_reads = []
raw_bases = []
clean_reads = []
clean_bases = []
avg_q20 = []
avg_q30 = []
avg_gc = []
passed_filter_reads = []
low_quality_reads = []
too_many_N_reads = []
duplication = []

for fname in os.listdir(input_dir):
    if fname.endswith(".filtered.fastp.json"):
        path = os.path.join(input_dir, fname)
        sample_name = fname.replace(".filtered.fastp.json", "")
        with open(path) as f:
            j = json.load(f)

            # 处理 before_filtering
            before = j['summary']['before_filtering']
            if 'read1' in before:  # 双端
                r1 = before['read1']['total_reads']
                r2 = before.get('read2', {}).get('total_reads', 0)
                raw_reads.append(r1 + r2 if r2 else r1)
                raw_bases.append(before['read1']['total_bases'])
            else:  # 单端
                raw_reads.append(before['total_reads'])
                raw_bases.append(before['total_bases'])

            # 处理 after_filtering
            after = j['summary']['after_filtering']
            if 'read1' in after:  # 双端
                r1_clean = after['read1']['total_reads']
                r2_clean = after.get('read2', {}).get('total_reads', 0)
                clean_reads.append(r1_clean + r2_clean if r2_clean else r1_clean)
                clean_bases.append(after['read1']['total_bases'])
                avg_q20.append(round(after['read1']['q20_rate']*100, 2))
                avg_q30.append(round(after['read1']['q30_rate']*100, 2))
                avg_gc.append(round(after['read1']['gc_content'], 2))
            else:  # 单端
                clean_reads.append(after['total_reads'])
                clean_bases.append(after['total_bases'])
                avg_q20.append(round(after['q20_rate']*100, 2))
                avg_q30.append(round(after['q30_rate']*100, 2))
                avg_gc.append(round(after['gc_content'], 2))

            # 可选字段，尽量使用 get 避免 KeyError
            passed_filter_reads.append(after.get('passed_filter_reads', 0))
            low_quality_reads.append(after.get('low_quality_reads', 0))
            too_many_N_reads.append(after.get('too_many_N_reads', 0))
            duplication.append(after.get('duplication', 0))

            samples.append(sample_name)

# 构建 DataFrame
df = pd.DataFrame({
    "Sample": samples,
    "Raw Reads.No": raw_reads,
    "Raw Bases(bp)": raw_bases,
    "Clean Reads.No": clean_reads,
    "Clean Bases(bp)": clean_bases,
    "Avg_Q20(%)": avg_q20,
    "Avg_Q30(%)": avg_q30,
    "Avg_GC(%)": avg_gc,
    "Passed_filter_reads": passed_filter_reads,
    "Low_quality_reads": low_quality_reads,
    "Too_many_N_reads": too_many_N_reads,
    "Duplication(%)": duplication
})

# 保存 Excel
df.to_excel(output_file, index=False)
print(f"成功生成 Excel: {output_file}")
