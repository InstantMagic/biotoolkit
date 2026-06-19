"""biotoolkit: a small, dependency-free molecular biology toolkit.

Implements the everyday sequence operations used in synthetic biology and
molecular cloning, written in plain Python so the logic stays readable.
"""

from __future__ import annotations

from .codon import codon_usage, optimize
from .orf import ORF, find_orfs
from .primers import tm, tm_gc, tm_wallace
from .restriction import ENZYMES, Hit, find_sites, single_cutters
from .sequence import (
    gc_content,
    reverse_complement,
    transcribe,
    translate,
    validate,
)

__version__ = "0.1.0"

__all__ = [
    "gc_content", "reverse_complement", "transcribe", "translate", "validate",
    "find_orfs", "ORF",
    "tm", "tm_wallace", "tm_gc",
    "find_sites", "single_cutters", "ENZYMES", "Hit",
    "optimize", "codon_usage",
    "read_fasta",
]


def read_fasta(path: str) -> dict[str, str]:
    """Read a FASTA file into a dict of {header: sequence}.

    Headers are taken from the text after '>' up to the first whitespace.
    """
    records: dict[str, str] = {}
    header = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks).upper()
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records[header] = "".join(chunks).upper()
    return records
