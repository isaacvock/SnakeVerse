#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover - exercised only in missing dependency envs.
    yaml = None


def require_yaml() -> None:
    if yaml is None:
        raise SystemExit(
            "PyYAML is required to read SnakeVerse YAML files. "
            "Install it in the active environment or use the snakeverse-dev conda env."
        )


def script_paths() -> tuple[Path, Path, Path]:
    config_root = Path(__file__).resolve().parents[1]
    project_root = config_root.parent
    shipped_root = config_root / "_ngsflow"
    return project_root, config_root, shipped_root


def load_yaml(path: str | Path) -> dict[str, Any]:
    require_yaml()
    yaml_path = Path(path)
    if not yaml_path.exists():
        raise FileNotFoundError(f"YAML file does not exist: {yaml_path}")
    with yaml_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML file must contain a mapping at top level: {yaml_path}")
    return data


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = deepcopy(base)
    for key, value in override.items():
        if (
            key in merged
            and isinstance(merged[key], dict)
            and isinstance(value, dict)
        ):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def project_path(value: str | Path, project_root: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (project_root / path).resolve()


def display_path(path: Path, project_root: Path) -> str:
    try:
        return path.relative_to(project_root).as_posix()
    except ValueError:
        return path.as_posix()


def manifest() -> dict[str, Any]:
    _, _, shipped_root = script_paths()
    return load_yaml(shipped_root / "manifest.yaml")


def assay_entry(assay: str) -> dict[str, Any]:
    assays = manifest().get("assays", {})
    if assay not in assays:
        raise SystemExit(f"Unknown assay '{assay}'. Run: python config/bin/ngsflow.py list assays")
    return assays[assay]


def preset_entry(assay: str, preset: str) -> dict[str, Any]:
    entry = assay_entry(assay)
    presets = entry.get("presets", {})
    if preset not in presets:
        raise SystemExit(
            f"Unknown preset '{preset}' for assay '{assay}'. "
            f"Run: python config/bin/ngsflow.py list presets --assay {assay}"
        )
    return presets[preset]


def render_template(text: str, replacements: dict[str, str]) -> str:
    rendered = text
    for key, value in replacements.items():
        rendered = rendered.replace("{{" + key + "}}", value)
    return rendered


def planned_copy(src: Path, dst: Path, replacements: dict[str, str] | None = None) -> tuple[Path, Path, dict[str, str] | None]:
    return src, dst, replacements


def write_text_from_template(
    src: Path,
    dst: Path,
    replacements: dict[str, str] | None,
    overwrite: bool,
    skip_existing: bool,
) -> str:
    existed = dst.exists()
    if dst.exists() and not overwrite:
        if skip_existing:
            return f"skip existing {dst}"
        raise FileExistsError(f"Refusing to overwrite existing file: {dst}")

    dst.parent.mkdir(parents=True, exist_ok=True)
    text = src.read_text(encoding="utf-8")
    if replacements:
        text = render_template(text, replacements)
    dst.write_text(text, encoding="utf-8")
    action = "overwrite" if existed and overwrite else "write"
    return f"{action} {dst}"


def write_active_pointer(config_root: Path, run_name: str) -> Path:
    active_config = config_root / "config.yaml"
    active_config.write_text(
        "# Active SnakeVerse run pointer.\n"
        f"run_config: config/runs/{run_name}.yaml\n",
        encoding="utf-8",
    )
    return active_config


def resolve_config(configfile: str | Path) -> dict[str, Any]:
    project_root, _, _ = script_paths()
    config_path = project_path(configfile, project_root)
    pointer = load_yaml(config_path)
    run_config_value = pointer.get("run_config")
    if not run_config_value:
        raise ValueError(f"{display_path(config_path, project_root)} must define run_config")

    run_config_path = project_path(run_config_value, project_root)
    run_config = load_yaml(run_config_path)

    resolved: dict[str, Any] = {}
    loaded_profiles: list[str] = []
    for profile_value in run_config.get("profile_stack", []) or []:
        profile_path = project_path(profile_value, project_root)
        profile = load_yaml(profile_path)
        resolved = deep_merge(resolved, profile)
        loaded_profiles.append(display_path(profile_path, project_root))

    tools_dir = config_path.parent / "profiles" / "tools"
    tools: dict[str, dict[str, Any]] = {}
    if tools_dir.exists():
        for tool_path in sorted(tools_dir.glob("*.yaml")):
            profile = load_yaml(tool_path)
            tools[str(profile.get("tool") or tool_path.stem)] = profile

    resolved = deep_merge(resolved, {"tools": tools})
    resolved = deep_merge(resolved, pointer)
    resolved = deep_merge(resolved, run_config)
    resolved["_ngsflow"] = {
        "project_root": project_root.as_posix(),
        "configfile": display_path(config_path, project_root),
        "run_config": display_path(run_config_path, project_root),
        "loaded_profiles": loaded_profiles,
        "loaded_tool_profiles": sorted(tools),
    }
    return resolved


def required_sample_columns(assay: str) -> list[str]:
    if assay == "rnaseq":
        return [
            "sample_id",
            "unit_id",
            "fastq_1",
            "fastq_2",
            "condition",
            "replicate",
            "strandedness",
        ]
    return ["sample_id", "unit_id", "fastq_1", "fastq_2"]


def read_samples(sample_path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with sample_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
            if any((value or "").strip() for value in row.values())
        ]
    return reader.fieldnames or [], rows


def command_list(args: argparse.Namespace) -> int:
    data = manifest()
    assays = data.get("assays", {})
    if args.kind == "assays":
        for assay, entry in assays.items():
            print(f"{assay}\t{entry.get('description', '')}")
        return 0

    if not args.assay:
        raise SystemExit("--assay is required when listing presets")
    entry = assay_entry(args.assay)
    for preset, preset_data in entry.get("presets", {}).items():
        print(f"{preset}\t{preset_data.get('description', '')}")
    return 0


def command_init_run(args: argparse.Namespace) -> int:
    if args.overwrite and args.skip_existing:
        raise SystemExit("Use either --overwrite or --skip-existing, not both.")

    project_root, config_root, shipped_root = script_paths()
    preset = preset_entry(args.assay, args.preset)
    replacements = {
        "RUN_NAME": args.run_name,
        "ASSAY": args.assay,
        "PRESET": args.preset,
        "GENOME": args.genome,
    }

    templates_root = shipped_root / "templates"
    operations: list[tuple[Path, Path, dict[str, str] | None]] = []
    operations.append(
        planned_copy(
            templates_root / preset["run_template"],
            config_root / "runs" / f"{args.run_name}.yaml",
            replacements,
        )
    )
    operations.append(
        planned_copy(
            templates_root / preset["sample_template"],
            config_root / "samples" / f"{args.run_name}.tsv",
            replacements,
        )
    )

    profile_groups = preset.get("profiles", {})
    for group, filenames in profile_groups.items():
        for filename in filenames:
            src = templates_root / "profiles" / group / filename
            dst = config_root / "profiles" / group / filename
            operations.append(planned_copy(src, dst, replacements))

    missing_sources = [src for src, _, _ in operations if not src.exists()]
    if missing_sources:
        raise SystemExit(
            "Manifest references missing template files:\n  - "
            + "\n  - ".join(str(path) for path in missing_sources)
        )

    conflicts = [
        dst
        for _, dst, _ in operations
        if dst.exists() and not args.overwrite and not args.skip_existing
    ]
    if conflicts:
        raise SystemExit(
            "Refusing to overwrite existing files. Use --overwrite or --skip-existing:\n  - "
            + "\n  - ".join(display_path(path, project_root) for path in conflicts)
        )

    for src, dst, replacements_for_file in operations:
        message = write_text_from_template(
            src,
            dst,
            replacements_for_file,
            overwrite=args.overwrite,
            skip_existing=args.skip_existing,
        )
        print(message.replace(str(project_root) + "\\", "").replace(str(project_root) + "/", ""))

    active_config = write_active_pointer(config_root, args.run_name)
    print(f"activated {display_path(active_config, project_root)}")
    print()
    print("Next steps:")
    print(f"  1. Edit config/samples/{args.run_name}.tsv")
    print(f"  2. Edit profiles under config/profiles/ as needed")
    print("  3. Run: python config/bin/ngsflow.py validate --configfile config/config.yaml")
    print("  4. Run: snakemake --configfile config/config.yaml --use-conda --cores 16")
    return 0


def command_activate_run(args: argparse.Namespace) -> int:
    project_root, config_root, _ = script_paths()
    run_name = args.run_name.removesuffix(".yaml")
    run_path = config_root / "runs" / f"{run_name}.yaml"
    if not run_path.exists():
        raise SystemExit(f"Run config does not exist: {display_path(run_path, project_root)}")
    active_config = write_active_pointer(config_root, run_name)
    print(f"Activated {display_path(run_path, project_root)} via {display_path(active_config, project_root)}")
    return 0


def command_explain(args: argparse.Namespace) -> int:
    config = resolve_config(args.configfile)
    project = config.get("project", {})
    print("SnakeVerse run summary")
    print(f"  active config: {config['_ngsflow']['configfile']}")
    print(f"  run config: {config['_ngsflow']['run_config']}")
    print(f"  project: {project.get('name', '<unset>')}")
    print(f"  run: {project.get('run_name', '<unset>')}")
    print(f"  assay: {config.get('assay', '<unset>')}")
    print(f"  preset: {config.get('preset', '<unset>')}")
    print(f"  samples: {config.get('samples', '<unset>')}")
    print(f"  results: {config.get('results_dir', '<unset>')}")
    print("  loaded profiles:")
    for profile in config["_ngsflow"]["loaded_profiles"]:
        print(f"    - {profile}")
    print("  loaded tool profiles:")
    for tool in config["_ngsflow"]["loaded_tool_profiles"]:
        print(f"    - {tool}")
    print("  enabled outputs:")
    for name, enabled in (config.get("outputs", {}) or {}).items():
        print(f"    - {name}: {enabled}")
    return 0


def validation_messages(config: dict[str, Any]) -> tuple[list[str], list[str]]:
    project_root = Path(config["_ngsflow"]["project_root"])
    errors: list[str] = []
    warnings: list[str] = []

    for profile in config["_ngsflow"].get("loaded_profiles", []):
        if not project_path(profile, project_root).exists():
            errors.append(f"Missing profile: {profile}")

    samples_value = config.get("samples")
    if not samples_value:
        errors.append("Run config is missing samples")
        return errors, warnings

    sample_path = project_path(samples_value, project_root)
    if not sample_path.exists():
        errors.append(f"Samples file does not exist: {samples_value}")
        return errors, warnings

    fieldnames, rows = read_samples(sample_path)
    required = required_sample_columns(str(config.get("assay")))
    missing_columns = [column for column in required if column not in fieldnames]
    if missing_columns:
        errors.append("Samples file is missing columns: " + ", ".join(missing_columns))
    for row in rows:
        sample = row.get("sample_id", "<unknown>")
        for column in required:
            if not row.get(column):
                errors.append(f"Sample {sample} is missing {column}")
        if config.get("assay") == "rnaseq" and row.get("strandedness") not in {
            "unstranded",
            "forward",
            "reverse",
        }:
            errors.append(
                f"Sample {sample} has invalid strandedness '{row.get('strandedness')}'"
            )
        for column in ("fastq_1", "fastq_2"):
            value = row.get(column)
            if value and not project_path(value, project_root).exists():
                warnings.append(f"FASTQ path does not exist yet: {value}")

    aligner = (config.get("alignment") or {}).get("tool")
    required_tools = ["fastqc", "samtools", "multiqc"]
    if (config.get("steps") or {}).get("trimming", False):
        required_tools.append("cutadapt")
    if aligner:
        required_tools.append(aligner)
    if config.get("assay") == "rnaseq":
        required_tools.append("featurecounts")
    if (config.get("steps") or {}).get("coverage", False):
        required_tools.append("deeptools")
    missing_tools = [
        tool for tool in sorted(set(required_tools)) if tool not in config.get("tools", {})
    ]
    if missing_tools:
        errors.append("Missing active tool profiles: " + ", ".join(missing_tools))

    genome = config.get("genome") or {}
    if aligner == "bowtie2" and not genome.get("bowtie2_index"):
        errors.append("Genome profile must define genome.bowtie2_index for Bowtie2 runs")
    if aligner == "star" and not genome.get("star_index"):
        errors.append("Genome profile must define genome.star_index for STAR runs")
    if config.get("assay") == "rnaseq" and not genome.get("gtf"):
        errors.append("Genome profile must define genome.gtf for RNA-seq runs")

    for key in ("fasta", "gtf", "chrom_sizes", "bowtie2_index", "star_index"):
        value = genome.get(key)
        if value and not project_path(str(value), project_root).exists():
            warnings.append(f"Reference path for genome.{key} does not exist yet: {value}")

    return errors, warnings


def command_validate(args: argparse.Namespace) -> int:
    try:
        config = resolve_config(args.configfile)
        errors, warnings = validation_messages(config)
    except Exception as exc:
        print(f"ERROR: {exc}")
        return 1

    for warning in warnings:
        print(f"WARNING: {warning}")
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1

    if warnings:
        print("Validation passed with warnings.")
    else:
        print("Validation passed.")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ngsflow.py",
        description="Local helper for materializing and validating SnakeVerse configs.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    list_parser = subparsers.add_parser("list", help="List available shipped templates.")
    list_parser.add_argument("kind", choices=["assays", "presets"])
    list_parser.add_argument("--assay", help="Assay name when listing presets.")
    list_parser.set_defaults(func=command_list)

    init_parser = subparsers.add_parser("init-run", help="Initialize an active run config.")
    init_parser.add_argument("--assay", required=True)
    init_parser.add_argument("--preset", required=True)
    init_parser.add_argument("--genome", required=True)
    init_parser.add_argument("--run-name", required=True)
    init_parser.add_argument("--overwrite", action="store_true")
    init_parser.add_argument("--skip-existing", action="store_true")
    init_parser.set_defaults(func=command_init_run)

    activate_parser = subparsers.add_parser("activate-run", help="Switch config/config.yaml to an existing run.")
    activate_parser.add_argument("run_name")
    activate_parser.set_defaults(func=command_activate_run)

    explain_parser = subparsers.add_parser("explain", help="Explain the active resolved config.")
    explain_parser.add_argument("--configfile", default="config/config.yaml")
    explain_parser.set_defaults(func=command_explain)

    validate_parser = subparsers.add_parser("validate", help="Run basic active config checks.")
    validate_parser.add_argument("--configfile", default="config/config.yaml")
    validate_parser.set_defaults(func=command_validate)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
