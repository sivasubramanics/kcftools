#!/usr/bin/env python3
# script: GTF_utils.py
# description: A utility to parse and manipulate GTF files (tired of AGAT and bedtools and combinations)
# author: c.s.sivasubramani@gmail.com
# date: 2025-10-07

import argparse
import sys, re, os
from collections import defaultdict
from functools import total_ordering
import networkx as nx
import json
import logging
from argparse import RawTextHelpFormatter

FTYPE_GENE = {"gene", "pseudogene"}
FTYPE_TR = {
    "mRNA",
    "transcript",
    "RNA",
    "lnc_RNA",
    "rRNA",
    "snRNA",
    "snoRNA",
    "tRNA",
    "primary_transcript",
    "antisense_RNA",
    "miRNA",
    "ncRNA",
    "SRP_RNA",
    "pseudogenic_tRNA",
}
FTYPE_EXON = {"exon"}
FTYPE_CDS = {"CDS", "start_codon", "stop_codon"}
FTYPE_UTR = {
    "five_prime_UTR",
    "three_prime_UTR",
    "5UTR",
    "3UTR",
    "five-UTR",
    "three-UTR",
}
FTYPE = FTYPE_GENE | FTYPE_TR | FTYPE_EXON | FTYPE_CDS | FTYPE_UTR
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
CODON_TABLE = {
    # A-starting codons
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    # C-starting codons
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    # G-starting codons
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    # T-starting codons
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGC": "C",
    "TGT": "C",
    "TGA": "*",
    "TGG": "W",
}


class FASTA(object):
    def __init__(self, name, sequence, description=None):
        self.name = name
        self.sequence = sequence
        self.description = description

    def __str__(self):
        return f">{self.name}\n{self.sequence}"

    def __len__(self):
        return len(self.sequence)

    def rev_complement(self):
        return FASTA(self.name, reverse_complement(self.sequence))

    def write_seq(self, handle, fold_len=60):
        if self.description:
            header = f">{self.name} {self.description}\n"
        else:
            header = f">{self.name}\n"
        if fold_len < 0:
            logging.error("fold_len must be non-negative")
            sys.exit(1)
        elif fold_len == 0:
            handle.write(f"{header}")
            handle.write(f"{self.sequence}\n")
        else:
            handle.write(f"{header}")
            for i in range(0, len(self.sequence), fold_len):
                handle.write(self.sequence[i : i + fold_len] + "\n")

    def extract_subseq(self, start: int, end: int):
        if start < 0 or end < 0:
            logging.error("start and end must be positive")
            sys.exit(1)
        if start >= end:
            logging.error("start must be smaller than end")
            sys.exit(1)
        if end > len(self.sequence):
            logging.error("end exceeds sequence length")
            sys.exit(1)
        return self.sequence[start:end]


@total_ordering
class Loci:
    def __init__(self, chrom: str, start: int, end: int, strand: str = "+"):
        if start < 0 or end < 0:
            logging.error("start and end must be positive")
            sys.exit(1)
        if start >= end:
            print(start, end, chrom, strand)
            logging.error("start must be smaller than end")
            sys.exit(1)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand == "+"

    def __repr__(self):
        strand_str = "+" if self.strand else "-"
        return f"Loci(chrom='{self.chrom}', start={self.start}, end={self.end}, strand='{strand_str}')"

    def __str__(self):
        strand_str = "+" if self.strand else "-"
        return f"{self.chrom}:{self.start}-{self.end}:{strand_str}"

    def to_tsv(self):
        strand_str = "+" if self.strand else "-"
        return f"{self.chrom}\t{self.start}\t{self.end}\t{strand_str}"

    def __eq__(self, other):
        return (self.chrom, self.start, self.end, self.strand) == (
            other.chrom,
            other.start,
            other.end,
            other.strand,
        )

    def __lt__(self, other):
        return (self.chrom, self.start, self.end) < (
            other.chrom,
            other.start,
            other.end,
        )

    @property
    def length(self):
        return self.end - self.start

    def overlaps_with(self, other: "Loci") -> bool:
        return (
            self.chrom == other.chrom
            and self.start < other.end
            and other.start < self.end
        )

    def shift(self, offset: int) -> "Loci":
        return Loci(
            self.chrom,
            self.start + offset,
            self.end + offset,
            "+" if self.strand else "-",
        )

    def to_dict(self):
        return {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "strand": "+" if self.strand else "-",
        }


class Feature:
    def __init__(
        self,
        id: str,
        type: str,
        loci: Loci = None,
        score=None,
        frame=None,
        source: str = None,
        attributes: dict = None,
    ):
        self.id = id
        self.type = type
        self.loci = loci
        if score:
            self.score = score
        else:
            self.score = "."
        if frame in {"0", "1", "2"}:
            self.frame = frame
        else:
            self.frame = "."
        if source:
            self.source = source
        else:
            self.source = "gtf_util"
        self.attributes = attributes or {}
        self.parent = None
        self.children = []

    def __len__(self):
        if self.type in FTYPE_GENE:
            return self.loci.length
        if self.type in FTYPE_TR:
            return sum(
                child.loci.length for child in self.children if child.type == "exon"
            )
        if self.type in FTYPE_EXON | FTYPE_CDS | FTYPE_UTR:
            return self.loci.length
        return 0

    def __repr__(self):
        return f"Feature(type='{self.type}', id='{self.id}')"

    def __str__(self):
        attr_str = "; ".join(f"{k}={v}" for k, v in self.attributes.items())
        return f"{self.type}: {self.id} [{attr_str}]"

    def add_child(self, child: "Feature"):
        self.children.append(child)
        child.parent = self

    def get_children(self, type_filter=None):
        return (
            [c for c in self.children if c.type == type_filter]
            if type_filter
            else self.children
        )

    def get_ancestors(self):
        lineage = []
        current = self
        while current.parent:
            lineage.append(current.parent)
            current = current.parent
        return lineage

    def to_dict(self):
        return {
            "id": self.id,
            "type": self.type,
            "source": self.source,
            "score": self.score,
            "frame": self.frame,
            "loci": self.loci.to_dict() if self.loci else None,
            "attributes": self.attributes,
            "parent": self.parent.id if self.parent else None,
            "children": [child.id for child in self.children],
        }

    def to_gtf(self):
        if "protein_id" in self.attributes and self.attributes["protein_id"] == "-":
            del self.attributes["protein_id"]
        attr_str = "; ".join(f'{k} "{v}"' for k, v in self.attributes.items()) + ";"
        gtf_line = [
            self.loci.chrom,
            self.source,
            self.type,
            str(self.loci.start + 1),  # GTF is 1-based
            str(self.loci.end),
            self.score,
            "+" if self.loci.strand else "-",
            self.frame,
            attr_str,
        ]
        return "\t".join(gtf_line)


class GTF:
    def __init__(self, window_size=1_000_000):
        self.graph = nx.DiGraph()
        self.genes = defaultdict(lambda: defaultdict(list))
        self.window_size = window_size

    def add_feature(self, feature: Feature):
        self.graph.add_node(feature.id, feature=feature)
        if feature.type in FTYPE_GENE:
            self.genes[feature.loci.chrom][
                int(feature.loci.start / self.window_size)
            ].append(feature.id)

    def add_relationship(self, parent_id, child_id):
        if parent_id:
            self.graph.add_edge(parent_id, child_id)

    def get_feature(self, feature_id):
        return self.graph.nodes[feature_id]["feature"]

    def get_children(self, feature_id):
        return [self.get_feature(cid) for cid in self.graph.successors(feature_id)]

    def get_parents(self, feature_id):
        return [self.get_feature(pid) for pid in self.graph.predecessors(feature_id)]

    def get_descendants(self, feature_id):
        return [self.get_feature(did) for did in nx.descendants(self.graph, feature_id)]

    # based on the edges, add the children to it's corresponding parent in the node
    def fix_children(self):
        for parent_id, child_id in self.graph.edges:
            parent = self.get_feature(parent_id)
            child = self.get_feature(child_id)
            if child not in parent.children:
                parent.add_child(child)

    def to_json(self, filepath: str = None):
        data = {
            "nodes": {
                node_id: self.graph.nodes[node_id]["feature"].to_dict()
                for node_id in self.graph.nodes
            },
            "edges": list(self.graph.edges),
        }
        json_str = json.dumps(data, indent=2)

        if filepath:
            with open(filepath, "w") as f:
                f.write(json_str)
        return json_str

    def to_gtf(self, filepath: str):
        with open(filepath, "w") as f:
            for node_id in self.graph.nodes:
                feature = self.get_feature(node_id)
                if feature.loci:
                    f.write(feature.to_gtf() + "\n")
        logging.info(f"Written GTF to {filepath}")

    def print_summary(self):
        num_genes = sum(
            1
            for n in self.graph.nodes
            if self.graph.in_degree(n) == 0 and self.get_feature(n).type in FTYPE_GENE
        )
        num_transcripts = sum(
            1
            for n in self.graph.nodes
            if self.graph.in_degree(n) == 1 and self.get_feature(n).type in FTYPE_TR
        )
        num_exons = sum(
            1
            for n in self.graph.nodes
            if self.graph.in_degree(n) == 1 and self.get_feature(n).type == "exon"
        )
        num_cds = sum(
            1
            for n in self.graph.nodes
            if self.graph.in_degree(n) == 1 and self.get_feature(n).type == "CDS"
        )
        num_chrom = len(set(self.genes.keys()))
        logging.info("=" * 50)
        logging.info(f"{'Feature Summary':^50}")
        logging.info("=" * 50)
        logging.info(f"{'Category':<25} | {'Count':>20}")
        logging.info("-" * 50)
        logging.info(f"{'Total features':<25} | {self.graph.number_of_nodes():>20,}")
        logging.info(
            f"{'Total relationships':<25} | {self.graph.number_of_edges():>20,}"
        )
        logging.info(f"{'Chromosomes':<25} | {num_chrom:>20,}")
        logging.info(f"{'Genes':<25} | {num_genes:>20,}")
        logging.info(f"{'Transcripts':<25} | {num_transcripts:>20,}")
        logging.info(f"{'Exons':<25} | {num_exons:>20,}")
        logging.info(f"{'CDS':<25} | {num_cds:>20,}")
        logging.info("=" * 50)

    def gene2tr(self, output_tsv):
        with open(output_tsv, "w") as out_f:
            for parent_id in [
                n for n in self.graph.nodes if self.graph.in_degree(n) == 0
            ]:
                parent = self.get_feature(parent_id)

                transcript_ids = [
                    child.id
                    for child in self.get_children(parent_id)
                    if child.type in FTYPE_TR
                ]
                for tr_id in transcript_ids:
                    out_f.write(f"{parent.id}\t{tr_id}\n")

    def extract_genes(self, chrom: str = None, start: int = None, end: int = None):
        genes = set()
        if chrom in self.genes:
            start_mb = int(start / 1_000_000) - 1
            end_mb = int(end / 1_000_000) + 1
            for mb in range(start_mb, end_mb):
                for gene_id in self.genes[chrom].get(mb, []):
                    gene = self.get_feature(gene_id)
                    if gene.loci.overlaps_with(Loci(chrom, start, end, "+")):
                        genes.add(gene_id)
        return genes

    def filter(self, genes):
        logging.info(f"Filtering GTF to retain {len(genes)} genes")
        for gene_id in [n for n in self.graph.nodes if self.graph.in_degree(n) == 0]:
            feature = self.get_feature(gene_id)
            if feature.type in FTYPE_GENE and feature.id not in genes:
                # remove it's children before removing the node
                descendants = list(nx.descendants(self.graph, gene_id))
                for desc_id in descendants:
                    self.graph.remove_node(desc_id)
                self.graph.remove_node(gene_id)
        self.genes = {
            chrom: {k: [g for g in v if g in genes] for k, v in megas.items()}
            for chrom, megas in self.genes.items()
        }
        self.fix_children()
        self.print_summary()
        logging.info(f"Filtered GTF has {self.graph.number_of_nodes()} features")

    def to_tsv(self, filepath: str):
        with open(filepath, "w") as out_f:
            out_f.write(
                "gene_id\ttranscript_id\tsource\tfeature_type\tfeature_length\tis_longest\tnum_exon\tnum_cds\tprotein_id\tprotein_len\tchrom\tstart\tend\tstrand\n"
            )
            # gtf.to_json(output_tsv + ".json")
            for parent_id in [
                n for n in self.graph.nodes if self.graph.in_degree(n) == 0
            ]:
                parent = self.get_feature(parent_id)

                for child in self.get_children(parent_id):
                    gene_id = parent.id
                    transcript_id = child.id
                    if child.source.lower() not in {"agat", "gtf_util"}:
                        source = child.source
                    else:
                        source = parent.source
                    feature_type = child.type
                    feature_length = len(child)
                    protein_id = child.attributes.get("protein_id", "NA")
                    num_exon = 0
                    num_cds = 0
                    cds_len = 0
                    for b_child in self.get_children(child.id):
                        if b_child.type in {"exon", "EXON"}:
                            num_exon += 1
                        if b_child.type in {"CDS", "cds"}:
                            num_cds += 1
                            cds_len += len(b_child)
                    if cds_len % 3 != 0:
                        logging.warning(
                            f"CDS length {cds_len} for transcript {transcript_id} is not a multiple of 3. Could be truncated."
                        )
                    if cds_len == 0:
                        prot_len = "-"
                    else:
                        prot_len = int(cds_len / 3) - 1
                    out_f.write(
                        f"{gene_id}\t{transcript_id}\t{source}\t{feature_type}\t{feature_length}\t{child.attributes.get('is_longest', '-')}\t"
                    )
                    out_f.write(
                        f"{num_exon}\t{num_cds}\t{protein_id}\t{prot_len}\t{child.loci.to_tsv()}\n"
                    )
        logging.info(f"Written TSV to {filepath}")

    def add_protein_id(self):
        for feature_id in self.graph.nodes:
            feature = self.get_feature(feature_id)
            children = self.get_children(feature_id)

            cds_children = [
                child
                for child in children
                if child.type == "CDS" or child.type == "exon"
            ]
            if not cds_children:
                continue

            # Try to extract protein_id from any CDS child
            protein_id = "-"
            is_protein_coding = False
            for child in cds_children:
                protein_id = child.attributes.get("protein_id")
                if child.type in FTYPE_CDS:
                    is_protein_coding = True
                if protein_id:
                    break
            if not protein_id:
                if is_protein_coding:
                    protein_id = feature.attributes.get("transcript_id")
                else:
                    protein_id = "-"
            feature.attributes["protein_id"] = protein_id
            for child in cds_children:
                if child.type in FTYPE_CDS:
                    child.attributes["protein_id"] = protein_id

    def extract_fna(self, fasta_dict, feature_type, output_fasta, features=None):
        if feature_type == "tr":
            feature_type = FTYPE_TR
        elif feature_type == "ge":
            feature_type = FTYPE_GENE
        elif feature_type == "utr":
            feature_type = FTYPE_UTR
        else:
            feature_type = {feature_type}

        seq_count = 0
        with open(output_fasta, "w") as out_f:
            for node in self.graph.nodes:
                feature = self.get_feature(node)
                seq = ""
                if feature.type not in feature_type:
                    continue
                if feature.type in FTYPE_GENE:
                    locis = []
                    for exon in self.get_descendants(feature.id):
                        if exon.type == "exon":
                            locis.append(exon.loci)
                    # merge overlapping loci
                    if not locis:
                        continue
                    locis.sort(key=lambda x: x.start)
                    merged_locis = [locis[0]]
                    for loci in locis[1:]:
                        last = merged_locis[-1]
                        if loci.start <= last.end:
                            merged_locis[-1] = Loci(
                                last.chrom,
                                last.start,
                                max(last.end, loci.end),
                                "+" if last.strand else "-",
                            )
                        else:
                            merged_locis.append(loci)

                    for loci in merged_locis:
                        if loci.chrom not in fasta_dict:
                            logging.warning(
                                f"Chromosome {loci.chrom} not found in FASTA"
                            )
                            continue
                        seq += fasta_dict[loci.chrom].extract_subseq(
                            loci.start, loci.end
                        )
                    locis = merged_locis
                    if not feature.loci.strand:
                        seq = reverse_complement(seq)
                elif feature.type in FTYPE_TR:
                    locis = []
                    for exon in self.get_children(feature.id):
                        if exon.type == "exon":
                            locis.append(exon.loci)
                    if not locis:
                        continue
                    locis.sort(key=lambda x: x.start)
                    for loci in locis:
                        if loci.chrom not in fasta_dict:
                            logging.warning(
                                f"Chromosome {loci.chrom} not found in FASTA"
                            )
                            continue
                        seq += fasta_dict[loci.chrom].extract_subseq(
                            loci.start, loci.end
                        )
                    if not feature.loci.strand:
                        seq = reverse_complement(seq)
                elif feature.type in FTYPE_EXON | FTYPE_CDS | FTYPE_UTR:
                    loci = feature.loci
                    if loci.chrom not in fasta_dict:
                        logging.warning(f"Chromosome {loci.chrom} not found in FASTA")
                        continue
                    seq = fasta_dict[loci.chrom].extract_subseq(loci.start, loci.end)
                    locis = [loci]
                    if not feature.loci.strand:
                        seq = reverse_complement(seq)
                if seq:
                    seq_count += 1
                    fasta_record = FASTA(
                        feature.id,
                        seq,
                        f"type={feature.type} len={len(seq)} loci={feature.loci.__str__()} segments=[{','.join(map(str, locis))}]",
                    )
                    fasta_record.write_seq(out_f)
        logging.info(
            f"Extracted sequences for feature type(s) {feature_type} to {output_fasta} ({seq_count} sequences)"
        )

    def extract_faa(self, fasta_dict, output_fasta):
        seq_count = 0
        with open(output_fasta, "w") as out_f:
            for node in self.graph.nodes:
                feature = self.get_feature(node)
                if feature.type not in FTYPE_TR:
                    continue
                protein_id = feature.attributes.get("protein_id")
                if not protein_id or protein_id == "-":
                    protein_id = feature.id
                cds_list = [
                    child
                    for child in self.get_children(feature.id)
                    if child.type == "CDS"
                ]
                if not cds_list:
                    continue
                cds_list.sort(key=lambda x: x.loci.start)
                seq = ""
                for cds in cds_list:
                    loci = cds.loci
                    if loci.chrom not in fasta_dict:
                        logging.warning(f"Chromosome {loci.chrom} not found in FASTA")
                        continue
                    seq += fasta_dict[loci.chrom].extract_subseq(loci.start, loci.end)
                if not feature.loci.strand:
                    seq = reverse_complement(seq)
                # translate the sequence to protein
                protein_seq = tranlate(seq)
                if protein_seq:
                    seq_count += 1
                    fasta_record = FASTA(
                        protein_id,
                        protein_seq,
                        f"gene_id={feature.get_ancestors()[0].id if feature.get_ancestors() else 'NA'} transcript_id={feature.id} len={len(protein_seq)} loci={feature.loci.__str__()}",
                    )
                    fasta_record.write_seq(out_f)
        logging.info(
            f"Extracted protein sequences to {output_fasta} ({seq_count} sequences)"
        )

    def strip(self):
        # remove all the attributes, leaving gene_id, transcript_id, and protein_id
        for node in self.graph.nodes:
            feature = self.get_feature(node)
            new_attributes = {}
            if "gene_id" in feature.attributes:
                new_attributes["gene_id"] = feature.attributes["gene_id"]
            if "transcript_id" in feature.attributes:
                new_attributes["transcript_id"] = feature.attributes["transcript_id"]
            if "protein_id" in feature.attributes:
                new_attributes["protein_id"] = feature.attributes["protein_id"]
            if "is_longest" in feature.attributes:
                new_attributes["is_longest"] = feature.attributes["is_longest"]
            feature.attributes = new_attributes

    def is_longest(self):
        """
        For each gene, label transcripts as:
            - is_longest = "yes" for the longest transcript (based on CDS length, then exon length)
            - is_longest = "no"  for all other transcripts
        Criteria:
            1. The transcript with the longest CDS length
            2. If multiple transcripts share the same CDS length, keep the one with the longest exon length
            3. If no CDS exist, use exon length only
        """
        logging.info("Marking longest transcripts per gene...")

        num_genes = 0
        num_transcripts = 0

        # Iterate through genes (nodes with no parents)
        for gene_id in [n for n in self.graph.nodes if self.graph.in_degree(n) == 0]:
            gene = self.get_feature(gene_id)
            if gene.type not in FTYPE_GENE:
                continue

            transcripts = [
                tr for tr in self.get_children(gene_id) if tr.type in FTYPE_TR
            ]
            if not transcripts:
                continue

            num_genes += 1
            num_transcripts += len(transcripts)

            best_tr = None
            best_cds_len = -1
            best_exon_len = -1

            # Find the "best" transcript according to CDS/exon length rules
            for tr in transcripts:
                cds_len = 0
                exon_len = 0
                for child in self.get_children(tr.id):
                    if child.type.lower() == "cds":
                        cds_len += len(child)
                    elif child.type.lower() == "exon":
                        exon_len += len(child)

                if cds_len > best_cds_len:
                    best_tr = tr
                    best_cds_len = cds_len
                    best_exon_len = exon_len
                elif cds_len == best_cds_len and exon_len > best_exon_len:
                    best_tr = tr
                    best_exon_len = exon_len

            # If no CDS at all, use exon length only
            if best_cds_len == 0:
                best_tr = max(
                    transcripts,
                    key=lambda tr: sum(
                        len(child)
                        for child in self.get_children(tr.id)
                        if child.type.lower() == "exon"
                    ),
                )

            # Mark transcripts
            for tr in transcripts:
                if tr.id == best_tr.id:
                    tr.attributes["is_longest"] = "yes"
                else:
                    tr.attributes["is_longest"] = "no"

    def longest(self):
        # remove all the transcripts (and it's descendants) where is_longest is no
        logging.info("Retaining only longest transcripts per gene...")
        for gene_id in [n for n in self.graph.nodes if self.graph.in_degree(n) == 0]:
            gene = self.get_feature(gene_id)
            if gene.type not in FTYPE_GENE:
                continue

            transcripts = [
                tr for tr in self.get_children(gene_id) if tr.type in FTYPE_TR
            ]
            if not transcripts:
                continue

            for tr in transcripts:
                if tr.attributes.get("is_longest") == "no":
                    descendants = list(nx.descendants(self.graph, tr.id))
                    for desc_id in descendants:
                        self.graph.remove_node(desc_id)
                    self.graph.remove_node(tr.id)
        self.fix_children()
        self.print_summary()

    def sort(self):
        # sort the descendants of each Transcripts based on loci, and order the nodes
        dag = nx.DiGraph()
        logging.info("Sorting Features based on loci...")
        for chrom in self.genes:
            for mb in self.genes[chrom]:
                for gene_id in self.genes[chrom][mb]:
                    gene = self.get_feature(gene_id)
                    dag.add_node(gene.id, feature=gene)
                    transcripts = [
                        tr for tr in self.get_children(gene.id) if tr.type in FTYPE_TR
                    ]
                    transcripts.sort(key=lambda x: (x.loci.start, x.loci.end))
                    for tr in transcripts:
                        dag.add_node(tr.id, feature=tr)
                        dag.add_edge(gene.id, tr.id)
                        # get all children and sort them by loci and add to dag
                        children = self.get_children(tr.id)
                        children.sort(key=lambda x: (x.loci.start))
                        for child in children:
                            dag.add_node(child.id, feature=child)
                            dag.add_edge(tr.id, child.id)
                            # get all grandchildren and sort them by loci and add to dag
                            grandchildren = self.get_children(child.id)
                            grandchildren.sort(key=lambda x: (x.loci.start))
                            for grandchild in grandchildren:
                                dag.add_node(grandchild.id, feature=grandchild)
                                dag.add_edge(child.id, grandchild.id)
        self.graph = dag

    def gene_density(self, filepath: str):
        with open(filepath, "w") as out_f:
            out_f.write("chrom\tstart\tend\tgene_count\n")
            for chrom in sorted(self.genes.keys()):
                for mb in sorted(self.genes[chrom].keys()):
                    start = mb * self.window_size
                    end = start + self.window_size
                    gene_count = len(self.genes[chrom][mb])
                    out_f.write(f"{chrom}\t{start}\t{end}\t{gene_count}\n")
        logging.info(f"Written gene density to {filepath}")

    def rename_features(self, prefix: str, chrom: bool = False):
        """
        Rename gene / transcript related features in self.graph.

        - Genes (features with in_degree == 0) that match FTYPE_GENE are enumerated.
        - Transcripts children of each gene that match FTYPE_TR are enumerated per-gene.
        - Child features (CDS, exon, UTRs, start/stop codons) receive updated attributes
          linking them to the new gene_id and transcript_id, while preserving 'old_*' attributes.
        - Builds a new directed graph (new_dag) containing renamed nodes and relationships,
          then replaces self.graph and fixes relationships via self.fix_children().
        """

        logging.info("Renaming features with prefix '%s'", prefix)

        # helper: safely set old/new attributes on a feature object (feature.attributes is a dict)
        def _set_attrs(feat, **attrs):
            for k, v in attrs.items():
                if v is None:
                    continue
                # preserve existing old_* if already present
                if k.startswith("old_"):
                    # only set old_ if not already present
                    feat.attributes.setdefault(k, v)
                else:
                    feat.attributes[k] = v

        # normalized type matching sets
        CDS_TYPES = {"cds", "CDS"}
        EXON_TYPES = {"exon", "EXON"}
        FIVE_UTR_TYPES = {"five_prime_UTR", "5UTR", "five-UTR"}
        THREE_UTR_TYPES = {"three_prime_UTR", "3UTR", "three-UTR"}
        START_CODON_TYPES = {"start_codon", "start-codon", "startCodon"}
        STOP_CODON_TYPES = {"stop_codon", "stop-codon", "stopCodon"}

        new_dag = nx.DiGraph()
        gene_index = 0

        # iterate over genes: nodes with in_degree == 0
        gene_nodes = [n for n in self.graph.nodes if self.graph.in_degree(n) == 0]
        for gene_id in gene_nodes:
            gene = self.get_feature(gene_id)
            # skip if not a gene type
            if gene.type not in FTYPE_GENE:
                continue

            gene_index += 1
            if chrom:
                new_gene_id = f"{prefix}_{gene.loci.chrom}_g{gene_index:06d}"
            else:
                new_gene_id = f"{prefix}_g{gene_index:06d}"

            # update gene attributes
            _set_attrs(gene, gene_id=new_gene_id, old_gene_id=gene.id)

            # add gene node only if not present
            if new_gene_id not in new_dag:
                new_dag.add_node(new_gene_id, feature=gene)

            # transcripts should be numbered per gene
            transcript_index = 0
            transcripts = [
                tr for tr in self.get_children(gene.id) if tr.type in FTYPE_TR
            ]
            for tr in transcripts:
                transcript_index += 1
                new_tr_id = f"{new_gene_id}.t{transcript_index}"

                # update transcript attributes and preserve old ids
                _set_attrs(
                    tr,
                    gene_id=new_gene_id,
                    old_gene_id=gene.id,
                    transcript_id=new_tr_id,
                    old_transcript_id=tr.id,
                )

                if new_tr_id not in new_dag:
                    new_dag.add_node(new_tr_id, feature=tr)
                # connect gene -> transcript
                new_dag.add_edge(new_gene_id, new_tr_id)

                # process children of the transcript
                children = self.get_children(tr.id)
                for child in children:
                    # prepare common replacements for all relevant child types
                    child_replacements = {
                        "gene_id": new_gene_id,
                        "transcript_id": new_tr_id,
                        "old_gene_id": gene.id,
                        "old_transcript_id": tr.id,
                    }

                    # handle CDS specially (also set protein_id and preserve old_protein_id)
                    child_type = child.type
                    if child_type in CDS_TYPES:
                        # preserve existing protein_id into old_protein_id if present
                        old_prot = child.attributes.get("protein_id", "-")
                        _set_attrs(child, protein_id=new_tr_id, old_protein_id=old_prot)
                        _set_attrs(child, **child_replacements)
                    elif child_type in EXON_TYPES:
                        _set_attrs(child, **child_replacements)
                    elif child_type in FIVE_UTR_TYPES:
                        _set_attrs(child, **child_replacements)
                    elif child_type in THREE_UTR_TYPES:
                        _set_attrs(child, **child_replacements)
                    elif child_type in START_CODON_TYPES:
                        _set_attrs(child, **child_replacements)
                    elif child_type in STOP_CODON_TYPES:
                        _set_attrs(child, **child_replacements)
                    else:
                        # For any other child types that might be relevant, set the basic linkage
                        # but avoid overwriting unrelated attributes.
                        # This keeps behavior conservative and extendable.
                        _set_attrs(child, **child_replacements)

                    # add child node (use original id as node key to preserve uniqueness)
                    if child.id not in new_dag:
                        new_dag.add_node(child.id, feature=child)
                    new_dag.add_edge(new_tr_id, child.id)

                    # add grandchildren (one level deeper), preserving their IDs/feature objects
                    grandchildren = self.get_children(child.id)
                    for grandchild in grandchildren:
                        if grandchild.id not in new_dag:
                            new_dag.add_node(grandchild.id, feature=grandchild)
                        new_dag.add_edge(child.id, grandchild.id)

        # replace the graph and repair relationships
        self.graph = new_dag
        self.fix_children()
        self.print_summary()


def parse_gtf(gtf_file, do_sort=False, window_size=1_000_000):
    dag = GTF(window_size=window_size)
    cds = 0

    logging.info(f'Parsing GTF file "{gtf_file}"')
    with open(gtf_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue

            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                continue

            chrom, source, feature_type, start, end, score, strand, phase, attr_str = (
                fields
            )
            if feature_type in {
                "repeat_region",
                "dispersed_repeat",
                "sequence_feature",
                "region",
                "inverted_repeat",
                "intron",
            }:
                logging.warning(f"Skipping GTF entry {line}")
                continue

            if feature_type not in FTYPE:
                logging.error(f"Unknown feature type {feature_type} in line: {line}")
                sys.exit(1)

            start, end = int(start) - 1, int(end)  # GFF is 1-based inclusive
            attr_dict = parse_attributes(attr_str)
            loci = Loci(chrom, start, end, strand)

            if feature_type in FTYPE_GENE:
                exon = 0
                cds = 0
                five_utr = 0
                three_utr = 0
                id = attr_dict.get("gene_id")
                parent_id = ""
                # dag.genes[chrom][int(start / 1_000_000)].append(id)
            if feature_type in FTYPE_TR:
                id = attr_dict.get("transcript_id")
                parent_id = attr_dict.get("gene_id")
            if feature_type in FTYPE_EXON:
                exon += 1
                id = f"{attr_dict.get('transcript_id')}:exon:{exon}"
                parent_id = attr_dict.get("transcript_id")
            if feature_type in {"five_prime_UTR", "5UTR", "five-UTR"}:
                five_utr += 1
                id = f"{attr_dict.get('transcript_id')}:5UTR:{five_utr}"
                parent_id = attr_dict.get("transcript_id")
            if feature_type in {"three_prime_UTR", "3UTR", "three-UTR"}:
                three_utr += 1
                id = f"{attr_dict.get('transcript_id')}:3UTR:{three_utr}"
                parent_id = attr_dict.get("transcript_id")
            if feature_type == "CDS":
                cds += 1
                id = f"{attr_dict.get('transcript_id')}:CDS:{cds}"
                parent_id = attr_dict.get("transcript_id")
            if feature_type == "start_codon":
                id = f"{attr_dict.get('transcript_id')}:start_codon"
                parent_id = attr_dict.get("transcript_id")
            if feature_type == "stop_codon":
                id = f"{attr_dict.get('transcript_id')}:stop_codon"
                parent_id = attr_dict.get("transcript_id")

            feature = Feature(
                id=id,
                type=feature_type,
                loci=loci,
                score=score,
                frame=phase,
                source=source,
                attributes=attr_dict,
            )
            dag.add_feature(feature)
            if parent_id:
                if parent_id == id:
                    logging.error(
                        f"Parent ID and Child ID are same. Fix the GTF. Error at: {line}"
                    )
                    sys.exit(1)
                if parent_id in dag.graph:
                    dag.add_relationship(parent_id, id)
                elif parent_id:
                    logging.warning(f"Parent {parent_id} not in graph")

    dag.add_protein_id()
    dag.fix_children()
    dag.is_longest()
    if do_sort:
        dag.sort()
    dag.print_summary()
    logging.info(f'Finished paring GTF file "{gtf_file}"')
    return dag


def parse_attributes(attr_string):
    attr_pairs = [
        field.strip() for field in attr_string.strip().split(";") if field.strip()
    ]
    attr_dict = {}
    for pair in attr_pairs:
        if pair:
            attr = pair.split()
            key = attr[0]
            value = " ".join(attr[1:]).strip('"')
            # attr_dict[key] = strip_tag(value)
            attr_dict[key] = value
    return attr_dict


def print_tree(dag, feature_id, level=0):
    feature = dag.get_feature(feature_id)
    indent = "  " * level
    print(f"{indent}- {feature.id} ({feature.type})")

    for child in dag.get_children(feature_id):
        print_tree(dag, child.id, level + 1)


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def parse_fasta(fasta_file):
    logging.info(f'Parsing FASTA file "{fasta_file}"')
    with open(fasta_file) as f:
        sequences = {}
        seq_id = None
        seq_desc = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    sequences[seq_id] = FASTA(seq_id, "".join(seq_lines), seq_desc)
                seq_id = line[1:].split()[0]
                seq_desc = (
                    " ".join(line[1:].split()[1:]) if len(line.split()) > 1 else None
                )
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_id:
            sequences[seq_id] = FASTA(seq_id, "".join(seq_lines), seq_desc)
    logging.info(f'Parsed {len(sequences)} sequences from FASTA file "{fasta_file}"')
    return sequences


def tranlate(seq):
    seq = seq.upper().replace("U", "T")
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        amino_acid = CODON_TABLE.get(codon, "X")  # Use 'X' for unknown codons
        if amino_acid == "*":  # Stop codon
            continue
        protein.append(amino_acid)
    # remove the last amino acid if it's a stop codon
    if protein and (protein[-1] == "*" or protein[-1] == "-"):
        protein = protein[:-1]
    return "".join(protein)


def main():
    parser = argparse.ArgumentParser(
        description="GTF utils (Just another more simpler one)"
    )
    subparsers = parser.add_subparsers(dest="command", help="sub-command help")

    gtf2tsv_parser = subparsers.add_parser(
        "gtf2tsv", help="Convert GTF to TSV (each transcript/mRNA)"
    )
    gtf2tsv_parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input GTF file [AGAT fixed GTF file is recommended]",
    )
    gtf2tsv_parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    gtf2tsv_parser.add_argument(
        "-s", "--sort", action="store_true", help="Sort features based on loci"
    )

    gene2tr_parser = subparsers.add_parser(
        "gene2tr", help="Map genes to transcripts (Like Trinity gene_to_tr_map)"
    )
    gene2tr_parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input GTF file [AGAT fixed GTF file is recommended]",
    )
    gene2tr_parser.add_argument("-o", "--output", required=True, help="Output TSV file")

    extract_genes_parser = subparsers.add_parser(
        "extract_genes", help="Extract genes from GTF based on coordinates"
    )
    extract_genes_parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input GTF file [AGAT fixed GTF file is recommended]",
    )
    extract_genes_parser.add_argument(
        "-o", "--output", required=True, help="Output GTF file"
    )
    extract_genes_parser.add_argument("-c", "--chrom", help="Chromosome")
    extract_genes_parser.add_argument("-s", "--start", type=int, help="Start position")
    extract_genes_parser.add_argument("-e", "--end", type=int, help="End position")
    extract_genes_parser.add_argument(
        "-b", "--bed", help="Input BED file with regions to extract genes"
    )
    extract_genes_parser.add_argument(
        "-O", "--format", help="Output format ['gtf', 'tsv']", default="gtf"
    )
    extract_genes_parser.add_argument(
        "-g", "--gene", help="Gene ID(s) to extract", nargs="+"
    )
    extract_genes_parser.add_argument(
        "-G", "--gene_list", help="File with Gene IDs to extract"
    )

    extract_seq_parser = subparsers.add_parser(
        "extract_seq",
        help="Extract sequences from FASTA based on GTF",
        formatter_class=RawTextHelpFormatter,
    )
    extract_seq_parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input GTF file [AGAT fixed GTF file is recommended]",
    )
    extract_seq_parser.add_argument(
        "-f", "--fasta", required=True, help="Input FASTA file"
    )
    extract_seq_parser.add_argument(
        "-o", "--output", required=True, help="Output FASTA file"
    )
    extract_seq_parser.add_argument(
        "-t",
        "--type",
        help=(
            "Feature type(s) to extract sequences for.\n\n"
            "Predefined groups:\n"
            "  tr  →  all transcript-like features:\n"
            "         {'mRNA', 'transcript', 'RNA', 'lnc_RNA', 'rRNA', 'snRNA', 'snoRNA', 'tRNA'}\n"
            "  ge  →  all gene-like features:\n"
            "         {'gene', 'pseudogene'}\n\n"
            "You may also provide any specific feature type directly from GTF, e.g. 'gene', 'mRNA', 'exon', 'CDS', 'UTR', etc.\n\n"
        ),
    )
    extract_seq_parser.add_argument(
        "-p",
        "--protein",
        help="Extract protein sequences (FAA) to this file (Only for the transcripts with valid CDS features).",
    )

    gtf_strip_parser = subparsers.add_parser(
        "gtf_strip",
        help="Strip GTF attributes, retaining only gene_id, transcript_id, and protein_id (Kind of cleaning)",
    )
    gtf_strip_parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input GTF file [AGAT fixed GTF file is recommended]",
    )
    gtf_strip_parser.add_argument(
        "-o", "--output", required=True, help="Output GTF file"
    )
    gtf_strip_parser.add_argument(
        "-s", "--sort", action="store_true", help="Sort features based on loci"
    )

    longest_tr_parser = subparsers.add_parser(
        "longest_tr",
        help="For each gene, retain only the transcript with the longest CDS (or exon if no CDS)",
    )
    longest_tr_parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input GTF file [AGAT fixed GTF file is recommended]",
    )
    longest_tr_parser.add_argument(
        "-o", "--output", required=True, help="Output GTF file"
    )
    longest_tr_parser.add_argument(
        "-O", "--format", help="Output format ['gtf', 'tsv']", default="gtf"
    )

    gene_density_parser = subparsers.add_parser(
        "gene_density",
        help="Calculate gene density per chromosome for a given window size",
    )
    gene_density_parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input GTF file [AGAT fixed GTF file is recommended]",
    )
    gene_density_parser.add_argument(
        "-o", "--output", required=True, help="Output TSV file"
    )
    gene_density_parser.add_argument(
        "-w",
        "--window",
        type=int,
        default=1_000_000,
        help="Window size (default: 1,000,000)",
    )

    rename_features_parser = subparsers.add_parser(
        "rename_features",
        formatter_class=RawTextHelpFormatter,
        help="Rename and reformat GTF feature IDs",  # short summary for parent help
        description=(
            "Rename features in GTF to the format: [PREFIX]_[CHROM]_g{n}_t{m}\n"
            "(affects gene, transcript, protein IDs only)\n\n"
            "WARNING:\n"
            "This will replace the 'old_gene_id', 'old_transcript_id', and 'old_protein_id' "
            "attributes if they already exist.\n"
        ),
    )
    rename_features_parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input GTF file [AGAT fixed GTF file is recommended]",
    )
    rename_features_parser.add_argument(
        "-o", "--output", required=True, help="Output GTF file"
    )
    rename_features_parser.add_argument(
        "-p",
        "--prefix",
        default="",
        help="Prefix for gene IDs (default: '')",
    )
    rename_features_parser.add_argument(
        "-c",
        "--chr",
        action="store_true",
        help="Add Chromsome Name to the prefix (e.g. prefix_chr1_g000001)",
    )
    rename_features_parser.add_argument(
        "-x",
        "--strip",
        action="store_true",
        help="Also strip attributes, retaining only gene_id, transcript_id, and protein_id",
    )

    extract_genes_ibs_parser = subparsers.add_parser(
        "extract_genes_ibs",
        help="Extract genes from GTF based on coordinates from kcftools findIBS output",
    )
    extract_genes_ibs_parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input GTF file [AGAT fixed GTF file is recommended]",
    )
    extract_genes_ibs_parser.add_argument(
        "-o", "--output", required=True, help="Output TSV file"
    )
    extract_genes_ibs_parser.add_argument(
        "-f", "--findibs", required=True, help="Input kcftools findIBS output file"
    )

    print_tree_parser = subparsers.add_parser(
        "print_tree", help="Print the feature tree for a given feature ID"
    )
    print_tree_parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input GTF file [AGAT fixed GTF file is recommended]",
    )
    print_tree_parser.add_argument(
        "-f", "--feature_id", required=True, help="Feature ID to print the tree for"
    )

    if len(sys.argv) == 1:
        sys.argv.append("--help")

    elif len(sys.argv) == 2 and sys.argv[1] not in ["-h", "--help"]:
        sys.argv.append("--help")
    args = parser.parse_args(args=None if sys.argv[1:] else ["--help"])

    # set the logging level to INFO in format 'YYYY-MM-DD HH:MM:SS - LEVEL - MESSAGE'
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # log the command line arguments provided, with the message "cmdline: ", and print on the stdout
    logging.info("cmdline: %s", " ".join(sys.argv))

    if os.path.exists(args.input) is False:
        logging.error(f"Input file {args.input} does not exist")
        sys.exit(1)

    if args.command == "gtf2tsv":
        if args.sort:
            gtf = parse_gtf(args.input, do_sort=True)
        else:
            gtf = parse_gtf(args.input)
        gtf.to_tsv(args.output)
    elif args.command == "gene2tr":
        gtf = parse_gtf(args.input)
        gtf.gene2tr(args.output)
    elif args.command == "extract_genes":
        if args.bed:
            if args.chrom or args.start or args.end or args.gene or args.gene_list:
                logging.error(
                    "Cannot use --bed with --chrom, --start, --end, --gene, or --gene_list"
                )
                sys.exit(1)
            genes = set()
            gtf = parse_gtf(args.input)
            with open(args.bed) as bed_f:
                for line in bed_f:
                    line = line.strip()
                    if line.startswith("#") or not line:
                        continue
                    fields = line.split("\t")
                    if len(fields) < 3:
                        logging.error(f"Invalid BED entry: {line}")
                        sys.exit(1)
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    genes.update(gtf.extract_genes(chrom=chrom, start=start, end=end))
            logging.info(f"Found {len(genes)} genes for the provided BED regions")
        elif args.chrom and args.start is not None and args.end is not None:
            if args.bed or args.gene or args.gene_list:
                logging.error(
                    "Must provide --chrom, --start, and --end without --bed, --gene, or --gene_list"
                )
                sys.exit(1)
            gtf = parse_gtf(args.input)
            genes = gtf.extract_genes(chrom=args.chrom, start=args.start, end=args.end)
            logging.info(
                f"Found {len(genes)} genes for the region {args.chrom}:{args.start}-{args.end}"
            )
        elif args.gene:
            # warn if any of chrom, start, end, bed is provided with gene
            if args.chrom or args.start or args.end or args.bed or args.gene_list:
                logging.error(
                    "Cannot use --gene with --chrom, --start, --end, --bed, or --gene_list"
                )
                sys.exit(1)
            gtf = parse_gtf(args.input)
            genes = set(args.gene or [])
        elif args.gene_list:
            if args.chrom or args.start or args.end or args.bed or args.gene:
                logging.error(
                    "Cannot use --gene_list with --chrom, --start, --end, or --bed"
                )
                sys.exit(1)
            gtf = parse_gtf(args.input)
            genes = set()
            with open(args.gene_list) as gl_f:
                for line in gl_f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        genes.add(line)
            logging.info(f"Loaded {len(genes)} genes from {args.gene_list}")
        else:
            logging.error(
                "Must provide either --bed, or --chrom, --start, and --end, or --gene, or --gene_list"
            )
            sys.exit(1)
        gtf.filter(genes)
        if args.format == "tsv":
            gtf.to_tsv(args.output)
        else:
            gtf.to_gtf(args.output)
    elif args.command == "print_tree":
        gtf = parse_gtf(args.input)
        if args.feature_id not in gtf.graph:
            logging.error(f"Feature ID {args.feature_id} not found in the GTF")
            sys.exit(1)
        print_tree(gtf, args.feature_id)
    elif args.command == "extract_seq":
        if args.type not in {"tr", "ge"}:
            logging.warning(
                f"Only the spcified type (STRICTLY) '{args.type}' will be extracted from GTF."
            )
            if args.type in FTYPE_GENE:
                missing = FTYPE_GENE - {args.type}
                logging.warning(f"{missing} features will NOT be extracted.")
                logging.warning(
                    f"consider using 'ge' to extract all gene-like features."
                )
            elif args.type in FTYPE_TR:
                missing = FTYPE_TR - {args.type}
                logging.warning(f"{missing} features will NOT be extracted.")
                logging.warning(
                    f"consider using 'tr' to extract all transcript-like features."
                )

        gtf = parse_gtf(args.input)
        fasta = parse_fasta(args.fasta)
        gtf.extract_fna(fasta, args.type, args.output)
        if args.protein:
            gtf.extract_faa(fasta, args.protein)
    elif args.command == "gtf_strip":
        if args.sort:
            gtf = parse_gtf(args.input, do_sort=True)
        else:
            gtf = parse_gtf(args.input)
        gtf.strip()
        gtf.to_gtf(args.output)
    elif args.command == "longest_tr":
        gtf = parse_gtf(args.input)
        gtf.longest()
        if args.format == "tsv":
            gtf.to_tsv(args.output)
        else:
            gtf.to_gtf(args.output)
    elif args.command == "gene_density":
        gtf = parse_gtf(args.input, window_size=args.window)
        gtf.gene_density(args.output)
    elif args.command == "rename_features":
        gtf = parse_gtf(args.input, do_sort=True)
        gtf.rename_features(prefix=args.prefix, chrom=args.chr)
        if args.strip:
            gtf.strip()
        gtf.to_gtf(args.output)
    elif args.command == "extract_genes_ibs":
        gtf = parse_gtf(args.input)
        line_num = 0
        logging.info(f'Parsing findIBS file "{args.findibs}"')
        with open(args.findibs) as fibs_f, open(args.output, "w") as out_f:
            for line in fibs_f:
                line = line.strip()
                if line.startswith("#") or not line:
                    out_f.write(line + "\n")
                    continue
                line_num += 1
                fields = line.split("\t")
                if line_num == 1:
                    # header line
                    out_f.write(line + "\tgenes\n")
                    continue
                chrom = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                genes = gtf.extract_genes(chrom=chrom, start=start, end=end)
                out_f.write(line + "\t" + ",".join(sorted(genes)) + "\n")
        logging.info(f'Finished parsing findIBS file "{args.findibs}"')


if __name__ == "__main__":
    main()
