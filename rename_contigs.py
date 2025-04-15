#The name of contigs from hifiasm aren't particularly useful. If I don't have a chromsome-level assembly, I still would like the contig names to be meaningful.
#So I will map my reference genome from hifiasm to a chromsome-level reference genome of a closely related species (if it exists) using RagTag. 
#I then take the output AGP file from Ragtag to rename my assembly chromsomes to the mapping from Ragtag.

import argparse
from collections import defaultdict
from Bio import SeqIO

def parse_agp(agp_file):
    """
    Parses the AGP file to extract contig renaming information.
    """
    contig_mapping = {}
    chromosome_counts = defaultdict(int)
    unplaced_counts = 0

    with open(agp_file, 'r') as agp:
        for line in agp:
            if line.startswith("#") or not line.strip():
                continue  # Skip comments or empty lines

            fields = line.strip().split("\t")
            chrom, contig_name = fields[0], fields[5]

            if contig_name == "scaffold":  # Skip gap lines
                continue

            if "RagTag" in chrom:
                if chrom == "Mitochondrial_RagTag":
                    new_name = "Mitochondria"
                else:
                    chrom_num = chrom.split("_")[1]  # Extract chromosome number
                    chromosome_counts[chrom_num] += 1
                    new_name = f"Chr_{chrom_num.zfill(2)}_{str(chromosome_counts[chrom_num]).zfill(2)}"
            else:
                unplaced_counts += 1
                new_name = f"Unplaced_{str(unplaced_counts).zfill(2)}"

            contig_mapping[contig_name] = new_name

    return contig_mapping

def rename_fasta(fasta_file, agp_file, output_file="Renamed.fasta"):
    """
    Renames contigs in the FASTA file based on AGP file mapping.
    """
    contig_mapping = parse_agp(agp_file)

    with open(output_file, "w") as out_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            new_name = contig_mapping.get(record.id, record.id)  # Default to original if not found
            out_fasta.write(f">{new_name}\n{str(record.seq)}\n")

    print(f"Renamed FASTA file saved as: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rename HiFiASM contigs based on RagTag AGP output.")
    parser.add_argument("-q", "--query", required=True, help="Input FASTA file from HiFiASM assembly.")
    parser.add_argument("-a", "--agp", required=True, help="RagTag AGP file.")
    parser.add_argument("-o", "--output", default="Renamed.fasta", help="Output renamed FASTA file (default: Renamed.fasta).")

    args = parser.parse_args()

    rename_fasta(args.query, args.agp, args.output)
