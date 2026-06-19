"""Primer melting temperature (Tm) estimates.

Two standard approximations are provided. Nearest-neighbor thermodynamics
are more accurate; these formulas are the quick estimates used at the bench
for routine primer design.
"""

from __future__ import annotations

from .sequence import validate


def tm_wallace(seq: str) -> float:
    """Wallace rule: Tm = 2*(A+T) + 4*(G+C).

    A reasonable estimate for short oligos (under ~14 bases).
    """
    seq = validate(seq)
    at = sum(1 for b in seq if b in ("A", "T"))
    gc = sum(1 for b in seq if b in ("G", "C"))
    return float(2 * at + 4 * gc)


def tm_gc(seq: str) -> float:
    """Salt-independent GC formula for longer oligos.

    Tm = 64.9 + 41 * (GC - 16.4) / length
    """
    seq = validate(seq)
    n = len(seq)
    if n == 0:
        return 0.0
    gc = sum(1 for b in seq if b in ("G", "C"))
    return round(64.9 + 41 * (gc - 16.4) / n, 1)


def tm(seq: str) -> float:
    """Pick a sensible estimate based on primer length."""
    seq = validate(seq)
    if len(seq) < 14:
        return tm_wallace(seq)
    return tm_gc(seq)
