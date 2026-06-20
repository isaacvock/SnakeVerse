#!/usr/bin/env python3
"""Local configuration helper for SnakeVerse.

This script materializes editable files from config/_catalog. It is a helper,
not a workflow runner: Snakemake reads the generated files directly.
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


INDEX_FILES = {
    "star": (
        "Genome",
        "SA",
        "SAindex",
        "chrLength.txt",
        "chrName.txt",
        "chrNameLength.txt",
        "chrStart.txt",
        "genomeParameters.txt",
    ),
    "bowtie2": (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    "bwa_mem2": (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
}
ADAPTER_TOOL_KEYS = {
    "fastp": {"detect_adapter_for_pe", "adapter_sequence", "adapter_sequence_r2"},
    "cutadapt": {"adapter_r1", "adapter_r2"},
}
FEATURECOUNTS_RUN_KEYS = {
    "feature_type",
    "attribute_type",
    "strand",
    "strandedness",
    "paired_end",
    "count_read_pairs",
    "require_both_ends_mapped",
    "count_multimapping_reads",
    "count_overlapping_features",
}


def require_yaml() -> None:
    if yaml is None:
        raise SystemExit(
            "PyYAML is required. Install it in the active environment or use "
            "the snakeverse-dev conda environment."
        )


def script_paths() -> tuple[Path, Path, Path]:
    config_root = Path(__file__).resolve().parents[1]
    return config_root.parent, config_root, config_root / "_catalog"


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


def project_path(value: str | Path, project_root: Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (project_root / path).resolve()


def display_path(path: Path, project_root: Path) -> str:
    try:
        return path.resolve().relative_to(project_root.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def catalog_manifest() -> dict[str, Any]:
    return load_yaml(script_paths()[2] / "manifest.yaml")


def assay_entry(assay: str) -> dict[str, Any]:
    assays = catalog_manifest().get("assays", {})
    if assay not in assays:
        raise ValueError(
            f"Unknown assay '{assay}'. Run: python config/bin/ngsflow.py list assays"
        )
    return assays[assay]


def recipe_entry(assay: str, recipe: str) -> dict[str, Any]:
    recipes = assay_entry(assay).get("recipes", {})
    if recipe not in recipes:
        raise ValueError(
            f"Unknown recipe '{recipe}' for assay '{assay}'. Run: "
            f"python config/bin/ngsflow.py list recipes --assay {assay}"
        )
    return recipes[recipe]


def render_template(text: str, replacements: dict[str, str]) -> str:
    for key, value in replacements.items():
        text = text.replace("{{" + key + "}}", value)
    return text


def materialize(
    source: Path,
    destination: Path,
    replacements: dict[str, str],
    overwrite: bool,
    skip_existing: bool,
) -> str:
    existed = destination.exists()
    if existed and not overwrite:
        if skip_existing:
            return f"skip  {destination}"
        raise FileExistsError(
            f"Refusing to overwrite existing file: {destination}. "
            "Use --skip-existing to keep it or --overwrite to replace it."
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        render_template(source.read_text(encoding="utf-8"), replacements),
        encoding="utf-8",
    )
    return f"{'update' if existed else 'create'} {destination}"


def write_active_pointer(config_root: Path, run_name: str) -> Path:
    pointer = config_root / "config.yaml"
    pointer.write_text(
        "# Active SnakeVerse run pointer.\n"
        "# Switch runs with: python config/bin/ngsflow.py activate-run <run-name>\n"
        f"run_config: config/runs/{run_name}.yaml\n",
        encoding="utf-8",
    )
    return pointer


def resolve_config(configfile: str | Path) -> dict[str, Any]:
    project_root, _, _ = script_paths()
    pointer_path = project_path(configfile, project_root)
    pointer = load_yaml(pointer_path)
    run_value = pointer.get("run_config")
    if not run_value:
        raise ValueError(f"{display_path(pointer_path, project_root)} must define run_config")
    run_path = project_path(str(run_value), project_root)
    run = load_yaml(run_path)
    if run.get("schema_version") != 2:
        raise ValueError(
            f"Run config must use schema_version: 2: {display_path(run_path, project_root)}"
        )

    reference_value = run.get("reference")
    if not reference_value:
        raise ValueError(f"Run config must define reference: {display_path(run_path, project_root)}")
    reference_path = project_path(str(reference_value), project_root)
    reference_document = load_yaml(reference_path)
    reference = reference_document.get("reference")
    if not isinstance(reference, dict):
        raise ValueError(f"Reference file must contain a top-level reference mapping: {reference_path}")

    config_root = pointer_path.parent
    tools: dict[str, dict[str, Any]] = {}
    tool_files: dict[str, str] = {}
    for tool_path in sorted((config_root / "tools").glob("*.yaml")):
        tool_config = load_yaml(tool_path)
        name = str(tool_config.get("tool") or tool_path.stem)
        if name in tools:
            raise ValueError(f"Duplicate active tool configuration for '{name}'")
        tools[name] = tool_config
        tool_files[name] = display_path(tool_path, project_root)

    resolved = dict(run)
    resolved["reference"] = reference
    resolved["tools"] = tools
    resolved["_ngsflow"] = {
        "project_root": project_root.as_posix(),
        "configfile": display_path(pointer_path, project_root),
        "run_config": display_path(run_path, project_root),
        "reference_file": display_path(reference_path, project_root),
        "tool_files": tool_files,
    }
    return resolved


def module(config: dict[str, Any], name: str) -> dict[str, Any]:
    value = (config.get("modules") or {}).get(name, {})
    return value if isinstance(value, dict) else {}


def module_enabled(config: dict[str, Any], name: str) -> bool:
    return bool(module(config, name).get("enabled", False))


def required_tools(config: dict[str, Any], rows: list[dict[str, str]]) -> list[str]:
    tools: list[str] = []
    if module_enabled(config, "fastq_qc"):
        tools.append("fastqc")
    if module_enabled(config, "sra_download") and any(row.get("sra_id") for row in rows):
        tools.append("sra_tools")
    for name in ("trimming", "alignment", "gene_counts", "peak_calling"):
        if module_enabled(config, name):
            tools.append(str(module(config, name).get("tool") or ""))
    if any(module_enabled(config, name) for name in ("alignment", "mark_duplicates", "bam_filter", "bam_qc")):
        tools.append("samtools")
    if module_enabled(config, "coverage"):
        tools.append("deeptools")
    if module_enabled(config, "multiqc"):
        tools.append("multiqc")
    if module_enabled(config, "salmon_quantification"):
        tools.append("salmon")
    if module_enabled(config, "rsem_quantification"):
        tools.append("rsem")
    if module_enabled(config, "frip") or module_enabled(config, "tss_enrichment"):
        tools.append("bedtools")
    return sorted({name for name in tools if name})


def load_samples(config: dict[str, Any]) -> tuple[list[dict[str, str]], list[str]]:
    project_root = Path(config["_ngsflow"]["project_root"])
    sample_path = project_path(str(config.get("samples") or ""), project_root)
    if not sample_path.exists():
        raise FileNotFoundError(f"Samples file does not exist: {sample_path}")
    with sample_path.open("r", encoding="utf-8", newline="") as handle:
        text = "".join(line for line in handle if not line.lstrip().startswith("#"))
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    fields = reader.fieldnames or []
    rows = [
        {key: (value or "").strip() for key, value in row.items()}
        for row in reader
        if any((value or "").strip() for value in row.values())
    ]
    return rows, fields


def index_is_ready(config: dict[str, Any], aligner: str) -> bool:
    index_value = str((config.get("reference", {}).get("indexes", {}) or {}).get(aligner) or "").rstrip("/")
    if not index_value:
        return False
    root = Path(config["_ngsflow"]["project_root"])
    path = project_path(index_value, root)
    if aligner == "star":
        return all((path / name).exists() for name in INDEX_FILES[aligner])
    suffixes = INDEX_FILES[aligner]
    if aligner == "bowtie2":
        return all(Path(str(path) + suffix).exists() for suffix in suffixes) or all(
            Path(str(path) + suffix + "l").exists() for suffix in suffixes
        )
    return all(Path(str(path) + suffix).exists() for suffix in suffixes)


def validate_config(config: dict[str, Any]) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    rows, fields = load_samples(config)
    assay = str(config.get("assay") or "")
    required_columns = {
        "generic": ["sample_id", "unit_id", "fastq_1"],
        "rnaseq": ["sample_id", "unit_id", "fastq_1", "condition", "replicate", "strandedness"],
        "atacseq": ["sample_id", "unit_id", "fastq_1", "condition", "replicate"],
    }.get(assay, [])
    if not required_columns:
        errors.append(f"Unsupported assay: {assay}")
    missing_columns = [name for name in required_columns if name not in fields]
    if missing_columns:
        errors.append("Samples file is missing columns: " + ", ".join(missing_columns))
    if not rows:
        errors.append("Samples file contains no sample rows")
    units = [(row.get("sample_id"), row.get("unit_id")) for row in rows]
    if len(units) != len(set(units)):
        errors.append("Samples file contains a duplicate sample_id/unit_id combination")

    active_tools = config.get("tools", {}) or {}
    missing_tools = [name for name in required_tools(config, rows) if name not in active_tools]
    if missing_tools:
        errors.append(
            "Missing active tool configurations: " + ", ".join(missing_tools)
            + ". Add each with: python config/bin/ngsflow.py add-tool <tool>"
        )

    alignment = module(config, "alignment")
    aligner = str(alignment.get("tool") or "")
    if not module_enabled(config, "alignment") or aligner not in INDEX_FILES:
        errors.append("modules.alignment must select star, bowtie2, or bwa_mem2")
    elif not index_is_ready(config, aligner) and not config.get("reference", {}).get("fasta"):
        errors.append(f"reference.fasta is required because the {aligner} index must be built")

    trimming = module(config, "trimming")
    if module_enabled(config, "trimming"):
        tool = str(trimming.get("tool") or "")
        adapters = trimming.get("adapters", {}) or {}
        mode = adapters.get("mode")
        if tool not in {"fastp", "cutadapt"}:
            errors.append("modules.trimming.tool must be fastp or cutadapt")
        if mode not in {"auto", "sequences", "none"}:
            errors.append("modules.trimming.adapters.mode must be auto, sequences, or none")
        if mode == "auto" and tool != "fastp":
            errors.append("Adapter auto-detection is available only with fastp")
        if mode == "sequences" and not adapters.get("read1"):
            errors.append("Adapter mode sequences requires modules.trimming.adapters.read1")
        paired = any(row.get("fastq_2") or row.get("sra_layout") == "paired" for row in rows)
        if mode == "sequences" and paired and not adapters.get("read2"):
            errors.append("Paired-end adapter mode sequences requires modules.trimming.adapters.read2")
        if mode == "none":
            warnings.append("Trimming is enabled without adapter removal")
        conflicts = sorted(
            ADAPTER_TOOL_KEYS.get(tool, set()).intersection(
                (active_tools.get(tool, {}).get("params", {}) or {}).keys()
            )
        )
        if conflicts:
            errors.append(
                "Adapter settings belong in modules.trimming.adapters; remove from "
                f"config/tools/{tool}.yaml: " + ", ".join(conflicts)
            )

    reference = config.get("reference", {}) or {}
    if assay == "rnaseq" and module_enabled(config, "gene_counts") and not reference.get("gtf"):
        errors.append("reference.gtf is required by modules.gene_counts")
    if assay == "rnaseq" and module_enabled(config, "gene_counts"):
        modes = module(config, "gene_counts").get("modes", {}) or {}
        if not any(bool(value) for value in modes.values()):
            errors.append("modules.gene_counts is enabled but no counting mode is enabled")
        conflicts = sorted(
            FEATURECOUNTS_RUN_KEYS.intersection(
                (active_tools.get("featurecounts", {}).get("params", {}) or {}).keys()
            )
        )
        if conflicts:
            errors.append(
                "Counting semantics belong in modules.gene_counts; remove from "
                "config/tools/featurecounts.yaml: " + ", ".join(conflicts)
            )
        layouts = {
            "paired" if row.get("fastq_2") or row.get("sra_layout") == "paired" else "single"
            for row in rows
        }
        if len(layouts) > 1:
            errors.append("featureCounts requires one read layout per run; samples mix paired and single")
        if module(config, "gene_counts").get("strandedness", "sample") == "sample":
            strand_values = {row.get("strandedness") for row in rows}
            if len(strand_values) > 1:
                errors.append("featureCounts requires one strandedness per run; sample values differ")
    for name in ("salmon_quantification", "rsem_quantification"):
        settings = module(config, name)
        if module_enabled(config, name) and not (
            settings.get("gene_results", True) or settings.get("isoform_results", True)
        ):
            errors.append(f"modules.{name} is enabled but requests no result type")
    if assay == "atacseq" and module_enabled(config, "bam_filter") and not reference.get("blacklist"):
        warnings.append("ATAC-seq blacklist filtering is recommended; reference.blacklist is blank")
    if module_enabled(config, "tss_enrichment"):
        for key in ("tss_bed", "chrom_sizes"):
            if not reference.get(key):
                errors.append(f"reference.{key} is required by modules.tss_enrichment")

    project_root = Path(config["_ngsflow"]["project_root"])
    for key in ("fasta", "gtf", "chrom_sizes", "blacklist", "tss_bed"):
        value = reference.get(key)
        if value and not project_path(str(value), project_root).exists():
            warnings.append(f"Reference path does not exist yet: reference.{key} = {value}")
    for row in rows:
        sample = row.get("sample_id") or "<unknown>"
        if not row.get("fastq_1") and not row.get("sra_id"):
            errors.append(f"Sample {sample} needs fastq_1 or sra_id")
        if assay == "rnaseq" and row.get("strandedness") not in {"unstranded", "forward", "reverse"}:
            errors.append(f"Sample {sample} has invalid strandedness: {row.get('strandedness')}")
        if row.get("sra_id") and not module_enabled(config, "sra_download"):
            errors.append(f"Sample {sample} uses sra_id but modules.sra_download is disabled")
        for key in ("fastq_1", "fastq_2"):
            value = row.get(key)
            if value and not project_path(value, project_root).exists():
                warnings.append(f"FASTQ path does not exist yet: {value}")
    return errors, warnings


def command_list(args: argparse.Namespace) -> None:
    if args.kind == "assays":
        for name, data in catalog_manifest().get("assays", {}).items():
            print(f"{name}\t{data.get('description', '')}")
        return
    if not args.assay:
        raise ValueError("--assay is required when listing recipes")
    for name, data in assay_entry(args.assay).get("recipes", {}).items():
        print(f"{name}\t{data.get('description', '')}")


def command_init_run(args: argparse.Namespace) -> None:
    project_root, config_root, catalog_root = script_paths()
    recipe = recipe_entry(args.assay, args.recipe)
    replacements = {
        "RUN_NAME": args.run_name,
        "GENOME": args.genome,
        "ASSAY": args.assay,
        "RECIPE": args.recipe,
    }
    jobs: list[tuple[Path, Path]] = [
        (catalog_root / recipe["run_template"], config_root / "runs" / f"{args.run_name}.yaml"),
        (catalog_root / recipe["sample_template"], config_root / "samples" / f"{args.run_name}.tsv"),
        (catalog_root / "references" / f"{args.genome}.yaml", config_root / "references" / f"{args.genome}.yaml"),
    ]
    jobs.extend(
        (catalog_root / "tools" / f"{tool}.yaml", config_root / "tools" / f"{tool}.yaml")
        for tool in recipe.get("tools", [])
    )
    missing_sources = [source for source, _ in jobs if not source.exists()]
    if missing_sources:
        raise FileNotFoundError("Catalog files are missing: " + ", ".join(map(str, missing_sources)))
    conflicts = [destination for _, destination in jobs if destination.exists()]
    if conflicts and not (args.overwrite or args.skip_existing):
        raise FileExistsError(
            "Initialization would overwrite existing files:\n  "
            + "\n  ".join(display_path(path, project_root) for path in conflicts)
            + "\nUse --skip-existing to preserve them or --overwrite to replace them."
        )
    for source, destination in jobs:
        print(materialize(source, destination, replacements, args.overwrite, args.skip_existing))
    pointer = write_active_pointer(config_root, args.run_name)
    print(f"activate {display_path(pointer, project_root)}")
    print("\nEdit these files in order:")
    print(f"  1. config/samples/{args.run_name}.tsv - replace example samples and FASTQ paths")
    print(f"  2. config/references/{args.genome}.yaml - set FASTA, annotation, and any ready index")
    print(f"  3. config/runs/{args.run_name}.yaml - review adapters, modules, and resources")
    alignment_tool = load_yaml(config_root / "runs" / f"{args.run_name}.yaml").get("modules", {}).get("alignment", {}).get("tool")
    print(f"  4. config/tools/{alignment_tool}.yaml - review alignment policy")
    print("  5. Other files in config/tools/ - adjust advanced tool behavior only when needed")
    print("\nThen run:")
    print("  python config/bin/ngsflow.py explain")
    print("  python config/bin/ngsflow.py validate")
    print("  snakemake --configfile config/config.yaml --use-conda --cores 16")


def command_add_tool(args: argparse.Namespace) -> None:
    project_root, config_root, catalog_root = script_paths()
    source = catalog_root / "tools" / f"{args.tool}.yaml"
    if not source.exists():
        available = ", ".join(path.stem for path in sorted((catalog_root / "tools").glob("*.yaml")))
        raise ValueError(f"Unknown catalog tool '{args.tool}'. Available: {available}")
    destination = config_root / "tools" / source.name
    print(materialize(source, destination, {}, args.overwrite, args.skip_existing))
    print(f"Review {display_path(destination, project_root)} before running the workflow.")


def command_activate_run(args: argparse.Namespace) -> None:
    project_root, config_root, _ = script_paths()
    run_path = config_root / "runs" / f"{args.run_name}.yaml"
    if not run_path.exists():
        raise FileNotFoundError(f"Run config does not exist: {display_path(run_path, project_root)}")
    pointer = write_active_pointer(config_root, args.run_name)
    print(f"Activated {display_path(run_path, project_root)} via {display_path(pointer, project_root)}")


def print_default_explanation(config: dict[str, Any]) -> None:
    meta = config["_ngsflow"]
    project = config.get("project", {}) or {}
    print(f"Run: {project.get('run_name', '<unnamed>')} ({config.get('assay', '<unset>')})")
    recipe = (config.get("created_from") or {}).get("recipe")
    if recipe:
        print(f"Created from recipe: {recipe} (provenance only; it does not affect runtime)")
    print("\nFiles that control this run:")
    print(f"  run        {meta['run_config']}")
    print(f"  samples    {config.get('samples')}")
    print(f"  reference  {meta['reference_file']}")
    print("\nPipeline:")
    for name, settings in (config.get("modules") or {}).items():
        state = "ON " if isinstance(settings, dict) and settings.get("enabled") else "off"
        tool = f" [{settings.get('tool')}]" if isinstance(settings, dict) and settings.get("tool") else ""
        print(f"  {state:3} {name}{tool}")
    print("\nActive tool files:")
    for name, path in meta.get("tool_files", {}).items():
        print(f"  {name:18} {path}")
    trimming = module(config, "trimming")
    alignment = module(config, "alignment")
    print("\nReview first:")
    print(f"  adapters   mode={(trimming.get('adapters') or {}).get('mode', '<unset>')} (run file)")
    print(f"  alignment  {alignment.get('tool', '<unset>')} (run file + matching tool file)")
    print("  reference  required fields depend on the enabled modules")
    print("\nMore detail: ngsflow.py explain trimming | alignment | reference | tool <name>")


def command_explain(args: argparse.Namespace) -> None:
    config = resolve_config(args.configfile)
    topic = args.topic or []
    if not topic:
        print_default_explanation(config)
        return
    if topic[0] == "trimming":
        settings = module(config, "trimming")
        print("Trimming is configured in the run file under modules.trimming.")
        print(f"  enabled: {settings.get('enabled', False)}")
        print(f"  tool: {settings.get('tool', '<unset>')}")
        print(f"  adapters: {settings.get('adapters', {})}")
        print("Use mode auto with fastp, sequences for explicit adapters, or none for quality-only trimming.")
    elif topic[0] == "alignment":
        settings = module(config, "alignment")
        tool = str(settings.get("tool") or "")
        print("The run file selects the aligner; its active tool file controls advanced CLI behavior.")
        print(f"  selection: {settings}")
        print(f"  tool file: {config['_ngsflow'].get('tool_files', {}).get(tool, '<missing>')}")
        print(f"  index: {(config.get('reference', {}).get('indexes', {}) or {}).get(tool) or '<build from FASTA>'}")
    elif topic[0] == "reference":
        print(f"Reference file: {config['_ngsflow']['reference_file']}")
        print("Blank optional fields are valid; requirements are derived from enabled modules.")
        for key, value in config.get("reference", {}).items():
            print(f"  {key}: {value if value not in (None, '') else '<blank>'}")
    elif topic[0] == "tool" and len(topic) == 2:
        name = topic[1]
        if name not in config.get("tools", {}):
            raise ValueError(f"Tool '{name}' is not active. Add it with: ngsflow.py add-tool {name}")
        print(f"Tool file: {config['_ngsflow']['tool_files'][name]}")
        print(yaml.safe_dump(config["tools"][name], sort_keys=False).rstrip())
    else:
        raise ValueError("Explain topics: trimming, alignment, reference, or tool <name>")


def command_validate(args: argparse.Namespace) -> None:
    config = resolve_config(args.configfile)
    errors, warnings = validate_config(config)
    for warning in warnings:
        print(f"WARNING: {warning}")
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
    print(f"Valid SnakeVerse configuration: {config['_ngsflow']['run_config']}")
    if warnings:
        print(f"Validation passed with {len(warnings)} warning(s).")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Materialize and inspect SnakeVerse configuration.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    list_parser = subparsers.add_parser("list", help="List catalog assays or recipes.")
    list_parser.add_argument("kind", choices=["assays", "recipes"])
    list_parser.add_argument("--assay", help="Assay name when listing recipes.")
    list_parser.set_defaults(func=command_list)

    init_parser = subparsers.add_parser("init-run", help="Create one complete editable run.")
    init_parser.add_argument("--assay", required=True)
    init_parser.add_argument("--recipe", required=True)
    init_parser.add_argument("--genome", required=True)
    init_parser.add_argument("--run-name", required=True)
    init_mode = init_parser.add_mutually_exclusive_group()
    init_mode.add_argument("--overwrite", action="store_true")
    init_mode.add_argument("--skip-existing", action="store_true")
    init_parser.set_defaults(func=command_init_run)

    tool_parser = subparsers.add_parser("add-tool", help="Copy one advanced tool config from the catalog.")
    tool_parser.add_argument("tool")
    tool_mode = tool_parser.add_mutually_exclusive_group()
    tool_mode.add_argument("--overwrite", action="store_true")
    tool_mode.add_argument("--skip-existing", action="store_true")
    tool_parser.set_defaults(func=command_add_tool)

    activate_parser = subparsers.add_parser("activate-run", help="Point config/config.yaml at a run.")
    activate_parser.add_argument("run_name")
    activate_parser.set_defaults(func=command_activate_run)

    explain_parser = subparsers.add_parser("explain", help="Explain the active run or one topic.")
    explain_parser.add_argument("topic", nargs="*")
    explain_parser.add_argument("--configfile", default="config/config.yaml")
    explain_parser.set_defaults(func=command_explain)

    validate_parser = subparsers.add_parser("validate", help="Validate the active run.")
    validate_parser.add_argument("--configfile", default="config/config.yaml")
    validate_parser.set_defaults(func=command_validate)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        args.func(args)
    except (FileNotFoundError, FileExistsError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
