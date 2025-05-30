#!/usr/bin/env python3


###################################
####增加只不包含活性位点得表格制作OneSpeciesChecker，只适用于流程读入，目前不能单脚本执行
###################################

import pandas as pd
import argparse
import yaml
from log_config import setup_logging

def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

class SpeciesChecker:
    def __init__(self, fileA_path, fileB_path, gcf_annotation_path, output_path, log):
        self.fileA_path = fileA_path
        self.fileB_path = fileB_path
        self.gcf_annotation_path = gcf_annotation_path
        self.output_path = output_path
        self.log = log

    def extract_species_and_proteins(self, file_path):
        species_proteins = {}
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith('>'):
                        parts = line.split('_**_')
                        gcf_id = '_'.join(parts[0].split('>')[1].split('_')[:2])
                        protein_id = parts[1].split()[0]
                        species_proteins.setdefault(gcf_id, set()).add(protein_id)
            self.log.info(f"Successfully extracted species and proteins from {file_path}")
        except Exception as e:
            self.log.error(f"Error extracting from {file_path}: {e}")
        return species_proteins

    def extract_gcf_annotations(self):
        try:
            df = pd.read_csv(self.gcf_annotation_path)
            gcf_annotations = df.set_index('GCF_Prefix').to_dict('index')
            self.log.info("Successfully extracted GCF annotations")
            return gcf_annotations
        except Exception as e:
            self.log.error(f"Failed to extract GCF annotations: {e}")
            return {}

    def check_presence(self, species_proteins, file_path):
        presence_info = {}
        try:
            with open(file_path, 'r') as file:
                file_content = file.read()
                for gcf_id, proteins in species_proteins.items():
                    matching_proteins = [protein for protein in proteins if protein in file_content]
                    presence_info[gcf_id] = 'yes:' + '|'.join(matching_proteins) if matching_proteins else 'no'
            self.log.info(f"Checked presence for {file_path}")
        except Exception as e:
            self.log.error(f"Error checking presence in {file_path}: {e}")
        return presence_info

    def generate_report(self):
        try:
            species_proteins_A = self.extract_species_and_proteins(self.fileA_path)
            species_proteins_B = self.extract_species_and_proteins(self.fileB_path)
            gcf_annotations = self.extract_gcf_annotations()

            presence_A = self.check_presence(species_proteins_A, self.fileA_path)
            presence_B = self.check_presence(species_proteins_B, self.fileB_path)
            print(presence_B)

            all_gcf_ids = set(presence_A.keys()) | set(presence_B.keys())

            columns_order = ['GCF_ID', 'strain_Name', 'taxon','Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'homology', 'homology_with_site']

            data = []
            for gcf_id in all_gcf_ids:
                annotation = gcf_annotations.get(gcf_id, {key: 'N/A' for key in columns_order[1:-2]})
                row = {
                    'GCF_ID': gcf_id,
                    'homology': presence_A.get(gcf_id, 'no'),
                    'homology_with_site': presence_B.get(gcf_id, 'no'),
                    **annotation
                }
                ordered_row = {col: str(row.get(col, 'N/A')) if col in ['homology', 'homology_with_site'] else row.get(col, 'N/A') for col in columns_order}
                data.append(ordered_row)

            # 创建DataFrame
            df = pd.DataFrame(data)
            df = df[columns_order]
            df.to_csv(self.output_path, index=False, sep='\t')
            self.log.info(f"Output saved to {self.output_path}")
        except Exception as e:
            self.log.error(f"Failed to generate report: {e}")

            # 计算含有"yes"的行数
            fileA_yes_count = df['homology'].apply(lambda x: 'yes' in x).sum()
            fileB_yes_count = df['homology_with_site'].apply(lambda x: 'yes' in x).sum()

            # 计算差值
            difference = abs(fileA_yes_count - fileB_yes_count)

            # 日志输出结果
            self.log.info(f"homology fileA has {fileA_yes_count} rows with 'yes'.")
            self.log.info(f"homology_with_site fileB has {fileB_yes_count} rows with 'yes'.")
            self.log.info(f"The difference in 'yes' rows between FileA and FileB is {difference}.")


        #     df = pd.DataFrame(data)
        #     df = df[columns_order]
        #     df.to_csv(self.output_path, index=False, sep='\t')
        #     self.log.info(f"Output saved to {self.output_path}")
        # except Exception as e:
        #     self.log.error(f"Failed to generate report: {e}")

class OneSpeciesChecker:
    def __init__(self, fileA_path, gcf_annotation_path, output_path, log):
        self.fileA_path = fileA_path
        self.gcf_annotation_path = gcf_annotation_path
        self.output_path = output_path
        self.log = log

    def extract_species_and_proteins(self, file_path):
        species_proteins = {}
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith('>'):
                        parts = line.split('_**_')
                        gcf_id = '_'.join(parts[0].split('>')[1].split('_')[:2])
                        protein_id = parts[1].split()[0]
                        species_proteins.setdefault(gcf_id, set()).add(protein_id)
            self.log.info(f"Successfully extracted species and proteins from {file_path}")
        except Exception as e:
            self.log.error(f"Error extracting from {file_path}: {e}")
        return species_proteins

    def extract_gcf_annotations(self):
        try:
            df = pd.read_csv(self.gcf_annotation_path)
            gcf_annotations = df.set_index('GCF Prefix').to_dict('index')
            self.log.info("Successfully extracted GCF annotations")
            return gcf_annotations
        except Exception as e:
            self.log.error(f"Failed to extract GCF annotations: {e}")
            return {}

    def check_presence(self, species_proteins, file_path):
        presence_info = {}
        try:
            with open(file_path, 'r') as file:
                file_content = file.read()
                for gcf_id, proteins in species_proteins.items():
                    matching_proteins = [protein for protein in proteins if protein in file_content]
                    presence_info[gcf_id] = 'yes:' + '|'.join(matching_proteins) if matching_proteins else 'no'
            self.log.info(f"Checked presence for {file_path}")
        except Exception as e:
            self.log.error(f"Error checking presence in {file_path}: {e}")
        return presence_info

    def generate_report(self):
        try:
            species_proteins_A = self.extract_species_and_proteins(self.fileA_path)
            gcf_annotations = self.extract_gcf_annotations()

            presence_A = self.check_presence(species_proteins_A, self.fileA_path)

            all_gcf_ids = set(presence_A.keys())

            columns_order = ['GCF_ID', 'strain_Name', 'taxon','Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum','homology']

            data = []
            for gcf_id in all_gcf_ids:
                annotation = gcf_annotations.get(gcf_id, {key: 'N/A' for key in columns_order[1:-2]})
                row = {
                    'GCF_ID': gcf_id,
                    'homology': presence_A.get(gcf_id, 'no'),
                    **annotation
                }
                ordered_row = {col: str(row.get(col, 'N/A')) if col in ['homology'] else row.get(col, 'N/A') for col in columns_order}
                data.append(ordered_row)

            # 创建DataFrame
            df = pd.DataFrame(data, columns=columns_order)
            df.to_csv(self.output_path, index=False, sep='\t')
            self.log.info(f"Output saved to {self.output_path}")
        except Exception as e:
            self.log.error(f"Failed to generate report: {e}")

            # 计算含有"yes"的行数
            #fileA_yes_count = df['homology'].apply(lambda x: 'yes' in x).sum()

            # 日志输出结果
            #self.log.info(f"homology fileA has {fileA_yes_count} rows with 'yes'.")

def main():
    parser = argparse.ArgumentParser(description="Check species presence in files, record associated proteins, and include detailed GCF annotations.")
    parser.add_argument('--extractprotein', required=True, type=str, help='Output file for complete protein sequences.', metavar='PROTEIN_FILE')
    parser.add_argument("--site_protein_output", required=True,  type=str, help="Output file for matched protein sequences.", metavar="PROTEIN_OUT")
    parser.add_argument('--config', type=str, required=True, help='Path to configuration file')
    parser.add_argument('--log_file', '-l', help='Log file path (optional)')
    parser.add_argument('--formatoutput', type=str, required=True, help='Path to output file')
    args = parser.parse_args()

    config = load_config(args.config)
    logger = setup_logging(args.log_file) if args.log_file else setup_logging()

    format_path = config['medusa-annotation']

    checker = SpeciesChecker(args.extractprotein, args.site_protein_output, format_path, args.formatoutput, logger)
    checker.generate_report()

if __name__ == "__main__":
    main()

