from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml


class ConfigError(RuntimeError):
    """Raised when the active SnakeVerse configuration cannot be resolved."""


def load_yaml(path: str | Path) -> dict[str, Any]:
    """Load a YAML mapping and provide a path-specific error on failure."""
    yaml_path = Path(path)
    if not yaml_path.exists():
        raise ConfigError(f"YAML file does not exist: {yaml_path}")
    with yaml_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ConfigError(f"YAML file must contain a mapping at top level: {yaml_path}")
    return data


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge mappings; retained as a general library utility."""
    merged = deepcopy(base)
    for key, value in override.items():
        if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def _project_root_from_configfile(configfile: str | Path) -> Path:
    configfile_path = Path(configfile)
    if not configfile_path.is_absolute():
        configfile_path = Path.cwd() / configfile_path
    configfile_path = configfile_path.resolve()
    if configfile_path.parent.name == "config":
        return configfile_path.parent.parent
    return Path.cwd().resolve()


def _resolve_project_path(path: str | Path, project_root: Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return (project_root / candidate).resolve()


def _display_project_path(path: Path, project_root: Path) -> str:
    try:
        return path.relative_to(project_root).as_posix()
    except ValueError:
        return path.as_posix()


def _config_root_from_run_config(run_config_path: Path) -> Path:
    if run_config_path.parent.name == "runs":
        return run_config_path.parent.parent
    return run_config_path.parent


def load_reference(
    run_config: dict[str, Any], project_root: str | Path
) -> tuple[Path, dict[str, Any]]:
    """Load the single reference file named by a v2 run config."""
    reference_value = run_config.get("reference")
    if not isinstance(reference_value, str) or not reference_value:
        raise ConfigError("Run config must define reference as a YAML file path")
    root = Path(project_root)
    reference_path = _resolve_project_path(reference_value, root)
    document = load_yaml(reference_path)
    reference = document.get("reference")
    if not isinstance(reference, dict):
        raise ConfigError(
            f"Reference file must contain a top-level 'reference' mapping: {reference_path}"
        )
    return reference_path, reference


def load_tool_profiles(
    config_root: str | Path,
) -> tuple[dict[str, dict[str, Any]], dict[str, Path]]:
    """Load available editable tool profiles from config/tools/*.yaml."""
    tools_dir = Path(config_root) / "tools"
    tools: dict[str, dict[str, Any]] = {}
    paths: dict[str, Path] = {}
    if not tools_dir.exists():
        return tools, paths
    for tool_path in sorted(tools_dir.glob("*.yaml")):
        profile = load_yaml(tool_path)
        tool_name = str(profile.get("tool") or tool_path.stem)
        if tool_name in tools:
            raise ConfigError(f"Duplicate tool profile for '{tool_name}'")
        tools[tool_name] = profile
        paths[tool_name] = tool_path.resolve()
    return tools, paths


def resolve_config(configfile: str | Path) -> dict[str, Any]:
    """Resolve an active pointer into one run, one reference, and available tools."""
    project_root = _project_root_from_configfile(configfile)
    config_path = _resolve_project_path(configfile, project_root)
    pointer_config = load_yaml(config_path)
    return resolve_config_from_mapping(
        pointer_config,
        project_root=project_root,
        configfile_label=_display_project_path(config_path, project_root),
        fallback_config_root=config_path.parent,
    )


def resolve_config_from_mapping(
    pointer_config: dict[str, Any],
    project_root: str | Path | None = None,
    configfile_label: str = "<snakemake config>",
    fallback_config_root: str | Path | None = None,
) -> dict[str, Any]:
    """Resolve Snakemake's merged mapping using its active run pointer."""
    root = Path(project_root).resolve() if project_root else Path.cwd().resolve()
    run_config_value = pointer_config.get("run_config")
    if not run_config_value:
        raise ConfigError(f"{configfile_label} must define run_config")

    run_config_path = _resolve_project_path(str(run_config_value), root)
    run_config = load_yaml(run_config_path)
    if run_config.get("schema_version") != 2:
        raise ConfigError(
            f"Run config must use schema_version: 2: {_display_project_path(run_config_path, root)}"
        )

    config_root = (
        Path(fallback_config_root).resolve()
        if fallback_config_root is not None
        else _config_root_from_run_config(run_config_path)
    )
    reference_path, reference = load_reference(run_config, root)
    tools, tool_paths = load_tool_profiles(config_root)

    resolved = deepcopy(run_config)
    resolved["reference"] = reference
    resolved["tools"] = tools
    resolved["_ngsflow"] = {
        "project_root": root.as_posix(),
        "config_root": _display_project_path(config_root, root),
        "configfile": configfile_label,
        "run_config": _display_project_path(run_config_path, root),
        "reference_file": _display_project_path(reference_path, root),
        "tool_files": {
            tool: _display_project_path(path, root) for tool, path in sorted(tool_paths.items())
        },
    }
    return resolved


def module_config(config: dict[str, Any], name: str) -> dict[str, Any]:
    """Return one module mapping, or an empty mapping when it is absent."""
    value = (config.get("modules") or {}).get(name, {})
    return value if isinstance(value, dict) else {}


def module_enabled(config: dict[str, Any], name: str) -> bool:
    """Return whether a named v2 module is enabled."""
    return bool(module_config(config, name).get("enabled", False))


def module_tool(config: dict[str, Any], name: str, default: str = "") -> str:
    """Return the tool selected by a named module."""
    return str(module_config(config, name).get("tool") or default)


def get_results_dir(config: dict[str, Any]) -> str:
    results_dir = config.get("results_dir")
    if not results_dir:
        raise ConfigError("Resolved config is missing results_dir")
    return str(results_dir).rstrip("/")


def write_resolved_config(resolved_config: dict[str, Any], results_dir: str | Path) -> Path:
    """Write the resolved configuration for provenance and debugging."""
    outdir = Path(results_dir) / "config"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "resolved_config.yaml"
    with outpath.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(resolved_config, handle, sort_keys=False)
    return outpath
