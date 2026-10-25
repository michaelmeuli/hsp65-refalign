"""In-silico hsp65 PCR on an assembly (vendored from mlsa-kansasii mlsa/loci.py and mlsa/protocols.py; keep in sync by hand)."""
from __future__ import annotations

import numpy as np
import numpy.typing as npt

from .refalign import revcomp

# TB-11 forward; TB-12w tried before TB-12 as reverse (lab's routine-then-fallback order); 439-441 bp amplicon
HSP65_FWD = "ACCAACGATGGTGTGTCCAT"
HSP65_REV = ["CTTGTCGAASCGCATACCCT", "CTTGTCGAACCGCATACCCT"]  # TB-12w (S = C/G), TB-12

_BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
_IUPAC_EXPAND = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}


def _seq_to_codes(seq: str) -> npt.NDArray[np.int8]:
    arr = np.frombuffer(seq.upper().encode("ascii", errors="replace"), dtype=np.uint8)
    codes = np.full(arr.shape, -1, dtype=np.int8)
    for base, code in _BASE_CODE.items():
        codes[arr == ord(base)] = code
    return codes


def find_primer_sites(seq: str, primer: str, max_mismatches: int) -> list[int]:
    """0-based start of every window matching `primer` (IUPAC honoured) with at most max_mismatches substitutions."""
    codes = _seq_to_codes(seq)
    L = len(primer)
    if len(codes) - L + 1 <= 0:
        return []
    allowed = np.zeros((L, 4), dtype=bool)
    for i, base in enumerate(primer.upper()):
        for allowed_base in _IUPAC_EXPAND.get(base, "ACGT"):
            allowed[i, _BASE_CODE[allowed_base]] = True
    windows = np.lib.stride_tricks.sliding_window_view(codes, L)
    matches = np.zeros(windows.shape, dtype=bool)
    for base_code in range(4):  # per-base broadcast keeps memory flat on multi-Mb genomes
        matches |= (windows == base_code) & allowed[:, base_code][None, :]
    mismatches = L - matches.sum(axis=1)
    hit_starts: list[int] = np.where(mismatches <= max_mismatches)[0].tolist()
    return hit_starts


def in_silico_pcr(contigs: dict[str, str], fwd_primer: str = HSP65_FWD, rev_primer_options: list[str] = HSP65_REV,
                  max_amplicon: int = 2000, max_mismatches: int = 3) -> str | None:
    """Amplicon (inclusive of both primers) of the first reverse-primer option that pairs with the forward primer."""
    for rev_primer in rev_primer_options:
        rev_rc = revcomp(rev_primer)
        for contig_seq in contigs.values():
            for strand_seq in (contig_seq, revcomp(contig_seq)):
                fwd_sites = find_primer_sites(strand_seq, fwd_primer, max_mismatches)
                if not fwd_sites:
                    continue
                rev_sites = find_primer_sites(strand_seq, rev_rc, max_mismatches)
                if not rev_sites:
                    continue
                for fpos in sorted(fwd_sites):
                    for rpos in sorted(rev_sites):
                        rend = rpos + len(rev_rc)
                        if fpos < rend <= fpos + max_amplicon:
                            return strand_seq[fpos:rend]
    return None
