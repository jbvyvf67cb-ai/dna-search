"""Pedigree handling (PED-style: sex 1=male 2=female, phenotype 1=unaffected 2=affected)."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Person:
    id: str
    sex: int  # 1 male, 2 female, 0 unknown
    affected: bool | None  # None = unknown
    father: str | None = None
    mother: str | None = None


@dataclass
class Pedigree:
    people: dict[str, Person]
    proband: str

    @property
    def proband_person(self) -> Person:
        return self.people[self.proband]

    @property
    def father(self) -> Person | None:
        f = self.proband_person.father
        return self.people.get(f) if f else None

    @property
    def mother(self) -> Person | None:
        m = self.proband_person.mother
        return self.people.get(m) if m else None

    @property
    def is_trio(self) -> bool:
        return self.father is not None and self.mother is not None

    def validate(self, vcf_samples: list[str]) -> None:
        missing = [p for p in self.people if p not in vcf_samples]
        if missing:
            raise ValueError(f"pedigree samples not in VCF: {missing}; VCF has {vcf_samples}")
        if self.proband not in self.people:
            raise ValueError(f"proband {self.proband!r} not in pedigree")


_SEX = {"male": 1, "m": 1, "1": 1, "female": 2, "f": 2, "2": 2}


def from_config(cfg: dict) -> Pedigree:
    people = {}
    for s in cfg["samples"]:
        aff = s.get("affected")
        people[s["id"]] = Person(
            id=s["id"],
            sex=_SEX.get(str(s.get("sex", "")).lower(), 0),
            affected=aff if isinstance(aff, bool) else None,
            father=s.get("father") or None,
            mother=s.get("mother") or None,
        )
    proband = cfg.get("proband") or next(
        (p.id for p in people.values() if p.affected), next(iter(people)))
    return Pedigree(people, proband)
