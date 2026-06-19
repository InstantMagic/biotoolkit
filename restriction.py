"""Core DNA/RNA sequence operations.

Everything here is implemented from scratch using only the standard library,
so the logic is visible rather than hidden behind a dependency.
"""

from __future__ import annotations

# Standard genetic code (DNA codons, T not U). "*" marks a stop codon.
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

VALID_BASES = set("ATGCN")


def clean(seq: str) -> str:
    """Uppercase a sequence and strip whitespace. Does not validate."""
    return "".join(seq.split()).upper()


def validate(seq: str) -> str:
    """Return a cleaned sequence, raising ValueError on any non-DNA base."""
    seq = clean(seq)
    bad = set(seq) - VALID_BASES
    if bad:
        raise ValueError(f"Sequence contains non-DNA characters: {sorted(bad)}")
    return seq


def gc_content(seq: str) -> float:
    """Fraction of G and C bases, returned as a percentage (0-100)."""
    seq = validate(seq)
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in ("G", "C"))
    return round(100 * gc / len(seq), 2)


def reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA sequence."""
    seq = validate(seq)
    return "".join(_COMPLEMENT[b] for b in reversed(seq))


def transcribe(seq: str) -> str:
    """Transcribe the coding strand of DNA to mRNA (T becomes U)."""
    return validate(seq).replace("T", "U")


def translate(seq: str, to_stop: bool = True) -> str:
    """Translate DNA to a protein string using the standard genetic code.

    Reads frame 1 (starting at index 0). Trailing bases that do not form a
    full codon are ignored. If to_stop is True, translation halts at the
    first stop codon and the stop is not included.
    """
    seq = validate(seq)
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        residue = CODON_TABLE[codon]
        if residue == "*":
            if to_stop:
                break
            protein.append("*")
        else:
            protein.append(residue)
    return "".join(protein)
