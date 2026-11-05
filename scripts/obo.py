"""Minimal, dependency-free OBO reader: enough for validation and hierarchy queries.

For production work use OAK (oaklib) or pronto, which handle the full OBO/OWL
spec, relationships beyond is_a, and imports. This reader is deliberately small
so every step of the validation logic is visible.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field


def normalize_curie(raw: str) -> str:
    """EFO's OBO export writes 'efo:EFO_0009922'; everyone else writes 'EFO:0009922'."""
    raw = raw.strip()
    m = re.match(r"^[a-z]+:([A-Za-z]+)_(\d+)$", raw)
    return f"{m.group(1)}:{m.group(2)}" if m else raw


@dataclass
class Term:
    id: str
    name: str = ""
    is_a: list[str] = field(default_factory=list)
    exact_synonyms: list[str] = field(default_factory=list)
    obsolete: bool = False
    replaced_by: list[str] = field(default_factory=list)
    consider: list[str] = field(default_factory=list)


class Ontology:
    def __init__(self, path: str, prefix: str):
        self.prefix = prefix
        self.terms: dict[str, Term] = {}
        self.data_version = ""
        self._ancestor_cache: dict[str, set[str]] = {}
        self._parse(path)

    def _parse(self, path: str) -> None:
        text = open(path, encoding="utf-8").read()
        header, *stanzas = re.split(r"\n\[(?:Term|Typedef|Instance)\]\n", "\n" + text)
        m = re.search(r"^data-version: (.+)$", header, re.M)
        self.data_version = m.group(1).strip() if m else ""
        # Typedef stanzas are split out too; only keep blocks that look like terms.
        for block in stanzas:
            tid = re.search(r"^id: (\S+)", block, re.M)
            if not tid:
                continue
            t = Term(id=normalize_curie(tid.group(1)))
            name = re.search(r"^name: (.+)$", block, re.M)
            t.name = name.group(1).strip() if name else ""
            t.is_a = [normalize_curie(x) for x in re.findall(r"^is_a: (\S+)", block, re.M)]
            t.exact_synonyms = re.findall(r'^synonym: "(.+?)" EXACT', block, re.M)
            t.obsolete = bool(re.search(r"^is_obsolete: true", block, re.M))
            t.replaced_by = [normalize_curie(x) for x in re.findall(r"^replaced_by: (\S+)", block, re.M)]
            t.consider = [normalize_curie(x) for x in re.findall(r"^consider: (\S+)", block, re.M)]
            self.terms.setdefault(t.id, t)

    def get(self, tid: str) -> Term | None:
        return self.terms.get(tid)

    def ancestors(self, tid: str) -> set[str]:
        """Transitive is_a closure (excluding the term itself)."""
        if tid in self._ancestor_cache:
            return self._ancestor_cache[tid]
        seen: set[str] = set()
        stack = list(self.terms[tid].is_a) if tid in self.terms else []
        while stack:
            p = stack.pop()
            if p in seen:
                continue
            seen.add(p)
            if p in self.terms:
                stack.extend(self.terms[p].is_a)
        self._ancestor_cache[tid] = seen
        return seen

    def is_descendant_of(self, tid: str, root: str) -> bool:
        return tid == root or root in self.ancestors(tid)

    def exact_lookup(self, text: str, root: str | None = None) -> list[str]:
        """Terms whose label or EXACT synonym equals `text` (case-insensitive),
        optionally restricted to the branch under `root`. Obsolete terms are skipped."""
        key = text.strip().lower()
        hits = []
        for t in self.terms.values():
            if t.obsolete:
                continue
            if t.name.lower() == key or key in (s.lower() for s in t.exact_synonyms):
                if root is None or self.is_descendant_of(t.id, root):
                    hits.append(t.id)
        return hits
