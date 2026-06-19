"""Restriction site mapping.

Given a sequence, find where common restriction enzymes cut. This is the core
of choosing cloning strategies and reading a plasmid map.
"""

from __future__ import annotations

from dataclasses import dataclass

from .sequence import reverse_complement, validate

# Recognition sequences for a set of common enzymes. All listed here are
# palindromic, so they cut the same site on both strands.
ENZYMES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC",
    "XhoI": "CTCGAG",
    "PstI": "CTGCAG",
    "SacI": "GAGCTC",
    "KpnI": "GGTACC",
    "SmaI": "CCCGGG",
    "SpeI": "ACTAGT",
    "NdeI": "CATATG",
    "SalI": "GTCGAC",
}


@dataclass
class Hit:
    enzyme: str
    site: str
    position: int  # 0-based index where the recognition site starts


def find_sites(seq: str, enzymes: dict[str, str] | None = None) -> list[Hit]:
    """Return every recognition site found, sorted by position."""
    seq = validate(seq)
    enzymes = enzymes or ENZYMES
    hits: list[Hit] = []
    for name, site in enzymes.items():
        start = seq.find(site)
        while start != -1:
            hits.append(Hit(name, site, start))
            start = seq.find(site, start + 1)
    return sorted(hits, key=lambda h: h.position)


def single_cutters(seq: str, enzymes: dict[str, str] | None = None) -> list[str]:
    """Enzymes that cut the sequence exactly once.

    Single cutters are the useful ones for linearizing a plasmid or dropping
    in an insert, so this is a question you ask constantly when cloning.
    """
    hits = find_sites(seq, enzymes)
    counts: dict[str, int] = {}
    for h in hits:
        counts[h.enzyme] = counts.get(h.enzyme, 0) + 1
    return sorted(name for name, c in counts.items() if c == 1)
