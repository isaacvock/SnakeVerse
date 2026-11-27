from __future__ import annotations

import shlex
from typing import Any

from samples import sample_ids, sample_metadata


VALID_BASES = frozenset(("A", "C", "G", "T"))
VALID_MUT_TYPES = frozenset(
    f"{r}{a}" for r in "ACGT" for a in "ACGT" if r != a
)
ANOMALY_MODES = frozenset(("auto", "controls", "treated", "all_samples"))
VALID_LIBRARY_TYPES = frozenset(("treated", "negative_control"))


def _module(config: dict[str, Any], name: str) -> dict[str, Any]:
    value = (config.get("modules") or {}).get(name, {})
    return value if isinstance(value, dict) else {}


def normalize_mutprof_library_type(raw: str) -> str:
    value = (raw or "").strip().lower()
    if value == "negative control":
        return "negative_control"
    return value


def parse_configured_mut_types(config: dict[str, Any]) -> list[str]:
    """Return validated uppercase 2-letter RNA mutation types (e.g. ['AG', 'CT'])."""
    raw = _module(config, "mutation_counts").get("mut_types")
    tokens: list[str] = []
    if isinstance(raw, str):
        tokens = [item.strip().upper() for item in raw.split(",") if item.strip()]
    elif isinstance(raw, (list, tuple)):
        for item in raw:
            for piece in str(item).split(","):
                cleaned = piece.strip().upper()
                if cleaned:
                    tokens.append(cleaned)
    elif raw is not None:
        raise ValueError(
            "modules.mutation_counts.mut_types must be a list or comma-separated string"
        )

    deduped: list[str] = []
    for mt in tokens:
        if mt not in VALID_MUT_TYPES:
            raise ValueError(
                f"Invalid mutation type '{mt}' in modules.mutation_counts.mut_types; "
                "expected two distinct ACGT bases such as AG or CT"
            )
        if mt not in deduped:
            deduped.append(mt)

    if not deduped:
        raise ValueError(
            "modules.mutation_counts.mut_types must specify at least one mutation type (e.g. [AG])"
        )
    return deduped


def samples_by_library_type(
    samples: list[dict[str, str]], library_type: str
) -> list[str]:
    target = normalize_mutprof_library_type(library_type)
    matched: list[str] = []
    for sid in sample_ids(samples):
        val = normalize_mutprof_library_type(sample_metadata(samples, sid, "library_type"))
        if val == target:
            matched.append(sid)
    return matched


def anomaly_calling_selection(
    config: dict[str, Any], samples: list[dict[str, str]]
) -> tuple[list[str], bool]:
    """Return `(selected_sample_ids, use_target_thresholds)` for anomaly calling."""
    mode = str(_module(config, "anomaly_calling").get("mode", "auto")).strip().lower()
    controls = samples_by_library_type(samples, "negative_control")
    treated = samples_by_library_type(samples, "treated")
    all_sids = sample_ids(samples)

    if mode == "auto":
        if controls:
            return controls, False
        return treated, True
    if mode == "controls":
        return controls, False
    if mode == "treated":
        return treated, True
    if mode == "all_samples":
        return all_sids, True
    raise ValueError(
        f"Unsupported modules.anomaly_calling.mode '{mode}'; "
        f"expected one of {', '.join(sorted(ANOMALY_MODES))}"
    )


def anomaly_bam_samples(
    config: dict[str, Any], samples: list[dict[str, str]]
) -> list[str]:
    selected, _ = anomaly_calling_selection(config, samples)
    return selected


def render_anomaly_calling_args(
    config: dict[str, Any], samples: list[dict[str, str]]
) -> str:
    settings = _module(config, "anomaly_calling")
    _, use_target_thresholds = anomaly_calling_selection(config, samples)

    min_bq = int(settings.get("min_base_quality", 30))
    min_depth = int(settings.get("min_depth", 5))
    min_alt = int(settings.get("min_alt", 3))
    min_freq = float(settings.get("min_freq", 0.25))
    err_rate = float(settings.get("err_rate", 0.002))
    p_cutoff = float(settings.get("p_cutoff", 1e-4))

    parts = [
        "--min-base-quality",
        str(min_bq),
        "--min-depth",
        str(min_depth),
        "--min-alt",
        str(min_alt),
        "--min-freq",
        str(min_freq),
        "--err-rate",
        str(err_rate),
        "--p-cutoff",
        str(p_cutoff),
    ]

    if use_target_thresholds and bool(_module(config, "mutation_counts").get("enabled", False)):
        mut_types = parse_configured_mut_types(config)
        target_min_freq = float(settings.get("target_min_freq", 0.90))
        parts.extend(
            [
                "--target-mut-types",
                shlex.quote(",".join(mut_types)),
                "--target-min-freq",
                str(target_min_freq),
            ]
        )

    mpileup_extra = str(settings.get("mpileup_extra") or "").strip()
    if mpileup_extra:
        parts.append(f"--mpileup-extra={shlex.quote(mpileup_extra)}")

    return " ".join(parts)


def render_mutation_counts_args(config: dict[str, Any]) -> str:
    settings = _module(config, "mutation_counts")
    min_bq = int(settings.get("min_base_quality", 30))
    min_dist = int(settings.get("min_end_distance", 5))
    indel_flank = int(settings.get("indel_flank", 3))
    min_depth = int(settings.get("min_depth", 1))
    return (
        f"--min-base-quality {min_bq} "
        f"--min-end-distance {min_dist} "
        f"--indel-flank {indel_flank} "
        f"--min-depth {min_depth}"
    )


def render_editing_tracks_args(config: dict[str, Any]) -> str:
    settings = _module(config, "editing_tracks")
    mut_types = parse_configured_mut_types(config)
    min_depth = int(settings.get("min_depth", 5))
    min_edit_count = int(settings.get("min_edit_count", 1))
    parts = [
        "--mut-types",
        shlex.quote(",".join(mut_types)),
        "--min-depth",
        str(min_depth),
        "--min-edit-count",
        str(min_edit_count),
    ]
    if bool(settings.get("emit_bedgraph", False)):
        parts.append("--emit-bedgraph")
    return " ".join(parts)


def anomaly_targets(results_dir: str) -> list[str]:
    return [
        f"{results_dir}/mutations/anomalies/snps.txt",
        f"{results_dir}/mutations/anomalies/anomalies.tsv",
    ]


def per_base_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    return [
        f"{results_dir}/mutations/per_base/{sample}.tsv.gz"
        for sample in sample_ids(samples)
    ]


def editing_track_patterns(config: dict[str, Any], results_dir: str) -> list[str]:
    """Return rule output templates parameterized by `{sample}`."""
    mut_types = parse_configured_mut_types(config)
    emit_bg = bool(_module(config, "editing_tracks").get("emit_bedgraph", False))
    outputs: list[str] = []
    for mt in mut_types:
        for metric in ("edit_rate", "edit_count"):
            for strand in ("plus", "minus"):
                outputs.append(
                    f"{results_dir}/tracks/editing/{{sample}}.{mt}.{metric}.{strand}.bw"
                )
                if emit_bg:
                    outputs.append(
                        f"{results_dir}/tracks/editing/{{sample}}.{mt}.{metric}.{strand}.bedGraph.gz"
                    )
    return outputs


def editing_track_targets(
    config: dict[str, Any], samples: list[dict[str, str]], results_dir: str
) -> list[str]:
    patterns = editing_track_patterns(config, results_dir)
    return [
        pattern.format(sample=sample)
        for sample in sample_ids(samples)
        for pattern in patterns
    ]
