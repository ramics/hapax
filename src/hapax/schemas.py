"""A schema per stage output, so a present-but-empty file is not a filled section.

WHY SCHEMAS AND NOT JUST FILE EXISTENCE. The failure this is built against is not a
missing file; it is a file that is there and says nothing. `load()` returns None for an
absent file and the page prints "absent", which is at least visible. The worse case is the
one the first production pass over an untouched library produced: a page that rendered,
looked finished, and was a column of dashes, because the four fields the table read live
under a different key. Existence cannot catch that and a schema can.

WHAT EACH SCHEMA REQUIRES. The keys the REPORT reads, and nothing else. A schema that
described a stage's whole output would fail whenever a stage grew a field, which trains
people to loosen schemas; a schema that names what the page needs fails only when the page
would have drawn nothing. So these are contracts between a stage and the report, not
descriptions of the stages.

Each also requires its own non-emptiness where emptiness is the failure: a `counts` of
`{}` and a `subjects` of `[]` are both a section with no content, and both have happened.
An explicit census of nothing found is a different thing and is allowed, because "no
significant hits" has to be recordable - an empty table and a search that never ran look
identical, and conflating them cost a withdrawn claim. Where a stage can legitimately
measure nothing, the schema requires the field that SAYS so.
"""
from __future__ import annotations

NONEMPTY_OBJECT = {"type": "object", "minProperties": 1}
NONEMPTY_ARRAY = {"type": "array", "minItems": 1}

#: schema name (as declared on a registry Output) -> JSON Schema, Draft 2020-12.
SCHEMAS: dict[str, dict] = {
    "tiers": {
        "type": "object",
        "required": ["pairs", "tiers", "curve", "residue_pairs"],
        "properties": {"pairs": {"type": "integer", "minimum": 1},
                       "tiers": NONEMPTY_OBJECT,
                       # The cascade sweep's own measured points. The table that read a
                       # column of dashes was reading tiers[name] for fields that live
                       # here, so the contract names both.
                       "curve": {"type": ["object", "array"]}},
    },
    "cascade": {"type": "object", "required": ["pairs"],
                "properties": {"pairs": {"type": "integer", "minimum": 1}}},
    "residue_rows": {"type": "object"},
    "contigs": {"type": "object"},
    "reads": {"type": "object"},
    "tier_table": {
        "type": "object",
        "required": ["pairs", "operating_point", "breakdown", "curve",
                     "interval_method"],
        "properties": {"breakdown": {"type": ["object", "array"]},
                       # Every number ships with its operating point, so a tier table
                       # without one is not readable and must not render.
                       "operating_point": NONEMPTY_OBJECT,
                       "interval_method": {"type": "string", "minLength": 1}},
    },
    "layers": {
        "type": "object",
        "required": ["counts", "totals", "rule", "against_known_positive"],
        "properties": {"counts": NONEMPTY_OBJECT,
                       # A layer's fields measured only on the residue cannot be compared
                       # against the known pile, and its zero there is an absence rather
                       # than a result. The comparison arm is therefore required, not
                       # optional.
                       "against_known_positive": {"type": ["object", "array"]},
                       "rule": {"type": ["object", "string", "array"]}},
    },
    "layer_rows": {"type": "object"},
    "fields": {
        "type": "object",
        "required": ["pairs_scored", "fields"],
        "properties": {"pairs_scored": {"type": "integer", "minimum": 1},
                       "fields": {"type": ["object", "array"]}},
    },
    "linkage": {
        "type": "object",
        "required": ["result", "parameters", "pairs_used"],
        "properties": {"pairs_used": {"type": "integer", "minimum": 1},
                       "result": NONEMPTY_OBJECT},
    },
    "census": {
        "type": "object",
        # `classes` and `result` are the census itself; `bias` is the direction of the
        # bias, which every number here ships with or is unreadable.
        "required": ["classes", "result", "bias"],
        "properties": {"classes": {"type": ["object", "array"]}},
    },
    "strains": {
        "type": "object",
        # `state` is what makes a zero readable: "searched and found nothing" is a
        # legitimate state and is not the same as not having searched, so `subjects` may
        # be empty and `state` may not be missing.
        "required": ["state", "searched", "subjects", "significance_evalue"],
        "properties": {"state": {"type": "string", "minLength": 1},
                       "subjects": {"type": "array"},
                       "significance_evalue": {"type": ["number", "string"]}},
    },
    "divergent": {"type": "object"},
    "shortlist": {"type": "object"},
    "crossplant": {"type": "object"},
    "profile": {"type": "object", "minProperties": 1},
    "taxa": {"type": "object", "minProperties": 1},
    "ranking": {
        "type": "object",
        "required": ["unit", "formula", "piles", "blind_to"],
        # `blind_to` is required because the ranking is blind to anything present once,
        # measured, and a ranking that does not say what it cannot see reads as complete.
        "properties": {"piles": {"type": ["object", "array"]},
                       "unit": {"type": "string", "minLength": 1},
                       "blind_to": {"type": ["object", "array", "string"]}},
    },
    "calibration": {"type": "object"},
    "panel": {"type": "object"},
    "report": {},            # the page is HTML; its check is the render check
    "render_check": {
        "type": "object",
        "required": ["top"],
        "properties": {"top": NONEMPTY_ARRAY},
    },
}


class SchemaError(Exception):
    """A stage output that does not meet its contract with the report."""


def validate(schema_name: str, data) -> list[str]:
    """Problems with `data` under the named schema. Empty means it meets the contract.

    Returns rather than raises, because the caller collects every problem across a run
    and prints them together: a build that stops at the first absence makes the next one
    invisible, and the point of this pass is one list of what the report cannot draw.
    """
    if schema_name not in SCHEMAS:
        raise SchemaError(
            f"no schema named {schema_name!r}. A stage output with no schema is an output "
            f"nothing checks, which is how a page came to render a column of dashes.")
    schema = SCHEMAS[schema_name]
    if not schema:
        return []
    import jsonschema

    validator = jsonschema.Draft202012Validator(schema)
    return [f"{'/'.join(str(p) for p in e.absolute_path) or '(root)'}: {e.message}"
            for e in sorted(validator.iter_errors(data), key=lambda e: e.path)]
