import argparse
import os
import re

import pandas as pd


def read_featurecounts(file):
    df = pd.read_csv(file, sep='\t', comment='#')
    df.columns = [re.sub(r'.*/|\.sorted\.bam$', '', col) if col not in ['Geneid', 'Chr', 'Start', 'End', 'Strand',
                                                                        'Length'] else col for col in df.columns]
    return df


def calculate_fpkm(counts_df, lengths):
    counts_per_kb = counts_df.div(lengths, axis=0) * 1e3
    total_counts_per_million = counts_df.sum(axis=0) / 1e6
    fpkm = counts_per_kb.div(total_counts_per_million, axis=1)
    fpkm = fpkm.round(2)  # 保留两位小数
    return fpkm


def calculate_tpm(counts_df, lengths):
    counts_per_kb = counts_df.div(lengths, axis=0) * 1e3
    sum_counts_per_kb = counts_per_kb.sum(axis=0)
    tpm = counts_per_kb.div(sum_counts_per_kb, axis=1) * 1e6
    tpm = tpm.round(2)  # 保留两位小数
    return tpm


def save_to_excel(readcount_df, fpkm_df, tpm_df, output_file):
    merged_df = pd.DataFrame()
    merged_df['Geneid'] = readcount_df['Geneid']

    for sample in readcount_df.columns[1:]:
        merged_df[f'{sample}:FPKM'] = fpkm_df[sample]
        merged_df[f'{sample}:TPM'] = tpm_df[sample]
        merged_df[f'{sample}:read_count'] = readcount_df[sample]

    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        workbook = writer.book
        worksheet = workbook.add_worksheet('Expression')

        # Write the header with the desired format
        header_format = workbook.add_format({
            'font_name': 'Arial',
            'font_size': 10,
            'bold': True,
            'bg_color': '#4682b4',
            'font_color': 'white',
            'align': 'center',
            'valign': 'vcenter'
        })

        for col_num, value in enumerate(merged_df.columns.values):
            worksheet.write(0, col_num, value, header_format)

        # Set the row height and freeze panes
        worksheet.set_row(0, 20)
        worksheet.freeze_panes(1, 0)

        # Set all cells format
        cell_format = workbook.add_format({
            'font_name': 'Arial',
            'font_size': 10,
            'align': 'center',
            'valign': 'vcenter'
        })

        # Add the dataframe to the worksheet
        for r_idx, row in enumerate(merged_df.itertuples(index=False)):
            for c_idx, value in enumerate(row):
                worksheet.write(r_idx + 1, c_idx, value, cell_format)

        # Auto-adjust column width
        for i, col in enumerate(merged_df.columns):
            max_len = max(merged_df[col].astype(str).map(len).max(), len(col)) + 2
            worksheet.set_column(i, i, max_len, cell_format)


def main():
    parser = argparse.ArgumentParser(description="Process featurecounts.tsv to generate readcount, FPKM, and TPM files")
    parser.add_argument("input_file", help="Input featurecounts.tsv file")
    parser.add_argument("-o", "--output_dir", default=".", help="Output directory")
    parser.add_argument("-p", "--prefix", default="gene", help="Output file prefix")

    args = parser.parse_args()

    input_file = args.input_file
    output_dir = args.output_dir
    prefix = args.prefix

    # Read the featurecounts file
    df = read_featurecounts(input_file)

    gene_ids = df['Geneid']
    lengths = df['Length']
    sample_counts = df.drop(columns=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'])

    # Calculate read counts
    readcount_df = pd.concat([gene_ids, sample_counts], axis=1)
    readcount_df.to_csv(os.path.join(output_dir, f"{prefix}.readcount.csv"), index=False)

    # Calculate FPKM
    fpkm_df = calculate_fpkm(sample_counts, lengths)
    fpkm_df.insert(0, 'Geneid', gene_ids)
    fpkm_df.to_csv(os.path.join(output_dir, f"{prefix}.count.FPKM.csv"), index=False)

    # Calculate TPM
    tpm_df = calculate_tpm(sample_counts, lengths)
    tpm_df.insert(0, 'Geneid', gene_ids)
    tpm_df.to_csv(os.path.join(output_dir, f"{prefix}.count.TPM.csv"), index=False)

    # Save to Excel
    save_to_excel(readcount_df, fpkm_df, tpm_df, os.path.join(output_dir, "Sample.Expression.count.xlsx"))


if __name__ == "__main__":
    main()