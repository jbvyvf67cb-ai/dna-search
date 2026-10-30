"""Minimal streaming VCF reader (stdlib only).

Handles plain or bgzipped VCFs, splits multi-allelic records into one
Variant per ALT allele, and parses VEP CSQ annotations per allele.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from typing import Iterator

HOM_REF, HET, HOM_ALT = 0, 1, 2


def open_text(path: str):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


@dataclass
class Call:
    """One sample's genotype at one ALT allele (after splitting)."""

    gt: int | None  # HOM_REF / HET / HOM_ALT / None (missing)
    dp: int | None = None
    gq: int | None = None
    ad_ref: int | None = None
    ad_alt: int | None = None
    raw_gt: str = "."

    @property
    def has_alt(self) -> bool:
        return self.gt in (HET, HOM_ALT)

    @property
    def vaf(self) -> float | None:
        if self.ad_ref is None or self.ad_alt is None:
            return None
        tot = self.ad_ref + self.ad_alt
        return self.ad_alt / tot if tot else None


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    filter: str
    info: dict[str, str]
    calls: dict[str, Call]
    csq: list[dict[str, str]] = field(default_factory=list)

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def passed(self) -> bool:
        return self.filter in ("PASS", ".", "")


class VcfReader:
    def __init__(self, path: str):
        self.path = path
        self.header_lines: list[str] = []
        self.samples: list[str] = []
        self.csq_fields: list[str] = []
        self._fh = open_text(path)
        self._read_header()

    def _read_header(self) -> None:
        for line in self._fh:
            line = line.rstrip("\n")
            if line.startswith("##"):
                self.header_lines.append(line)
                if line.startswith("##INFO=<ID=CSQ,"):
                    m = re.search(r"Format: ([^\"]+)", line)
                    if m:
                        self.csq_fields = m.group(1).strip().split("|")
            elif line.startswith("#CHROM"):
                self.samples = line.split("\t")[9:]
                return
        raise ValueError(f"{self.path}: no #CHROM header line")

    def __iter__(self) -> Iterator[Variant]:
        for line in self._fh:
            if not line.strip() or line.startswith("#"):
                continue
            yield from self._parse(line.rstrip("\n"))

    def close(self) -> None:
        self._fh.close()

    def _parse(self, line: str) -> Iterator[Variant]:
        f = line.split("\t")
        chrom, pos, ref, alts = f[0], int(f[1]), f[3], f[4].split(",")
        flt = f[6]
        info = parse_info(f[7])
        fmt = f[8].split(":") if len(f) > 8 else []
        raw_samples = [s.split(":") for s in f[9:]]
        csq_by_allele = self._csq_by_allele(info.get("CSQ"), ref, alts)
        for i, alt in enumerate(alts, start=1):
            if alt == "*":  # spanning deletion placeholder
                continue
            calls = {
                name: _make_call(fmt, vals, i)
                for name, vals in zip(self.samples, raw_samples)
            }
            yield Variant(chrom, pos, ref, alt, flt, info, calls,
                          csq_by_allele.get(i, []))

    def _csq_by_allele(self, raw: str | None, ref: str, alts: list[str]):
        out: dict[int, list[dict[str, str]]] = {}
        if not raw or not self.csq_fields:
            return out
        vep_alleles = {vep_allele(ref, a, alts): i
                       for i, a in enumerate(alts, start=1)}
        for entry in raw.split(","):
            vals = entry.split("|")
            d = dict(zip(self.csq_fields, vals))
            if d.get("ALLELE_NUM"):
                idx = int(d["ALLELE_NUM"])
            else:
                idx = vep_alleles.get(d.get("Allele", ""), 1 if len(alts) == 1 else None)
            if idx is not None:
                out.setdefault(idx, []).append(d)
        return out


def parse_info(raw: str) -> dict[str, str]:
    info: dict[str, str] = {}
    if raw in (".", ""):
        return info
    for item in raw.split(";"):
        k, _, v = item.partition("=")
        info[k] = v if _ else "1"
    return info


def vep_allele(ref: str, alt: str, all_alts: list[str]) -> str:
    """Allele string VEP writes into CSQ for this ALT (VEP trims a shared
    first base only when every ALT of the record shares it with REF)."""
    if all(len(a) != len(ref) and a[:1] == ref[:1] for a in all_alts if a != "*"):
        trimmed = alt[1:]
        return trimmed or "-"
    return alt


def _int(x: str) -> int | None:
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def _make_call(fmt: list[str], vals: list[str], allele: int) -> Call:
    d = dict(zip(fmt, vals))
    raw_gt = d.get("GT", ".")
    alleles = re.split(r"[/|]", raw_gt)
    gt: int | None
    if any(a in (".", "") for a in alleles):
        gt = None
    else:
        idx = [int(a) for a in alleles]
        n_alt = sum(1 for a in idx if a == allele)
        if len(idx) == 1:  # haploid call (e.g. male chrX from some callers)
            gt = HOM_ALT if n_alt else HOM_REF
        elif n_alt == 0:
            gt = HOM_REF
        elif n_alt == len(idx):
            gt = HOM_ALT
        else:
            gt = HET
    ad = d.get("AD", "").split(",")
    ad_ref = _int(ad[0]) if ad and ad[0] not in ("", ".") else None
    ad_alt = _int(ad[allele]) if len(ad) > allele else None
    return Call(gt=gt, dp=_int(d.get("DP")), gq=_int(d.get("GQ")),
                ad_ref=ad_ref, ad_alt=ad_alt, raw_gt=raw_gt)
