"""Codon optimization for expression in E. coli.

When you express a gene from one organism in another, the host may use
different codons preferentially. Recoding the gene to the host's favored
codons (without changing the protein) can improve expression. This module
does the simplest, most legible version: swap each residue for E. coli's
most frequently used codon.

This is a teaching-grade model. Real optimization also balances GC content,
avoids unwanted restriction sites and mRNA secondary structure, and tunes
codon adaptation index rather than always taking the single top codon.
"""

from __future__ import annotations

from .sequence import CODON_TABLE, translate, validate

# Most frequently used codon per amino acid in highly expressed E. coli genes.
ECOLI_PREFERRED = {
    "F": "TTC", "L": "CTG", "I": "ATC", "M": "ATG", "V": "GTG",
    "S": "AGC", "P": "CCG", "T": "ACC", "A": "GCG", "Y": "TAC",
    "H": "CAC", "Q": "CAG", "N": "AAC", "K": "AAA", "D": "GAT",
    "E": "GAA", "C": "TGC", "W": "TGG", "R": "CGT", "G": "GGC",
    "*": "TAA",
}


def optimize(seq: str) -> str:
    """Recode a DNA coding sequence to E. coli preferred codons.

    The input is translated to protein, then each residue is written back
    using E. coli's favored codon. The output encodes the same protein.
    """
    seq = validate(seq)
    protein = translate(seq, to_stop=False)
    return "".join(ECOLI_PREFERRED[aa] for aa in protein)


def codon_usage(seq: str) -> dict[str, float]:
    """Return the fraction of each codon used across the sequence (frame 1)."""
    seq = validate(seq)
    counts: dict[str, int] = {}
    total = 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        counts[codon] = counts.get(codon, 0) + 1
        total += 1
    if total == 0:
        return {}
    return {c: round(n / total, 3) for c, n in sorted(counts.items())}
