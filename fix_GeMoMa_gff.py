#The output Gff from GeMoMa for genome annotation can be an absolute mess. I created this script to try and clean those annotations up. 
#This assumes that you are using a single species as a reference (or using the flag in GeMoMa to generate separate GFFs for each reference)


import sys
import os
from collections import defaultdict
import re

def parse_attributes(attr_str):
    """Parses a GFF attribute column into a dictionary."""
    return dict(item.split("=", 1) for item in attr_str.strip(";").split(";") if "=" in item)

def main():
    if len(sys.argv) < 2:
        print("Usage: python fix_gemoma_gff.py <input_gff> [output_gff]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "fixed_GeMoMa.gff"

    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)

    raw_lines = []
    with open(input_file) as f:
        raw_lines = [line for line in f if line.strip() and not line.startswith("#")]

    entries_by_transcript = defaultdict(lambda: {"CDS": [], "line": None, "ref_gene": ""})

    for line in raw_lines:
        parts = line.strip().split("\t")
        if len(parts) != 9:
            continue

        source = parts[1]
        feature_type = parts[2]
        attributes = parse_attributes(parts[8])

        if source == "GAF":
            continue  

        if feature_type == "mRNA":
            transcript_id = attributes.get("ID", "")
            ref_gene = attributes.get("ref-gene", "")
            entries_by_transcript[transcript_id]["line"] = parts
            entries_by_transcript[transcript_id]["ref_gene"] = ref_gene

        elif feature_type == "CDS":
            parent_id = attributes.get("Parent", "")
            entries_by_transcript[parent_id]["CDS"].append(parts)

    gene_counter = 1
    output_lines = []

    for transcript_id, data in entries_by_transcript.items():
        mRNA_line = data["line"]
        CDS_lines = data["CDS"]
        ref_gene = data["ref_gene"]

        if mRNA_line is None or not CDS_lines:
            continue  
        ref_match = re.search(r'_(.*?)_', ref_gene)
        gene_name = ref_match.group(1) if ref_match else "Unknown"
        gene_id = f"{gene_name}_gene{gene_counter:05d}"

        chrom, source, _, start, end, _, strand, _, _ = mRNA_line

        
        output_lines.append("\t".join([
            chrom, source, "gene", start, end, ".", strand, ".",
            f"ID={gene_id};Name={gene_name};gene={gene_name}"
        ]))

       
        output_lines.append("\t".join([
            chrom, source, "mRNA", start, end, ".", strand, ".",
            f"ID=rna_{gene_id};Parent={gene_id};gene={gene_name}"
        ]))

       
        for cds_parts in CDS_lines:
            cds_start, cds_end, cds_phase = cds_parts[3], cds_parts[4], cds_parts[7]
            output_lines.append("\t".join([
                cds_parts[0], source, "CDS", cds_start, cds_end, ".", cds_parts[6], cds_phase,
                f"ID=cds_{gene_id};Parent=rna_{gene_id};gene={gene_name}"
            ]))

        gene_counter += 1  


    with open(output_file, "w") as out:
        out.write("\n".join(output_lines) + "\n")

    print(f"✅ Fixed GFF written to: {output_file}")

if __name__ == "__main__":
    main()
