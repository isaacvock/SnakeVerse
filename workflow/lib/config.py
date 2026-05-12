from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml


class ConfigError(RuntimeError):
    """Raised when the active SnakeVerse configuration cannot be resolved."""


def load_yaml(path: str | Path) -> dict[str, Any]:
    """Load a YAML document and return an empty dict for empty files."""
    yaml_path = Path(path)
    if not yaml_path.exists():
        raise ConfigError(f"YAML file does not exist: {yaml_path}")
    with yaml_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ConfigError(f"YAML file must contain a mapping at top level: {yaml_path}")
    return data


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge two dictionaries without mutating either input."""
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


def load_profile_stack(
    run_config: dict[str, Any], project_root: str | Path
) -> list[tuple[Path, dict[str, Any]]]:
    """Load every profile listed in a run config, preserving listed order."""
    root = Path(project_root)
    loaded: list[tuple[Path, dict[str, Any]]] = []
    for raw_path in run_config.get("profile_stack", []) or []:
        profile_path = _resolve_project_path(raw_path, root)
        loaded.append((profile_path, load_yaml(profile_path)))
    return loaded


def load_tool_profiles(config_root: str | Path) -> dict[str, dict[str, Any]]:
    """Load active tool profiles from config/profiles/tools/*.yaml."""
    root = Path(config_root)
    tools_dir = root / "profiles" / "tools"
    tools: dict[str, dict[str, Any]] = {}
    if not tools_dir.exists():
        return tools
    for tool_path in sorted(tools_dir.glob("*.yaml")):
        profile = load_yaml(tool_path)
        tool_name = str(profile.get("tool") or tool_path.stem)
        tools[tool_name] = profile
    return tools


def resolve_config(configfile: str | Path) -> dict[str, Any]:
    """Resolve the pointer config, profile stack, tool profiles, and run config."""
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
    """Resolve config from Snakemake's merged config mapping.

    Snakemake merges config files listed in a Snakefile with config files supplied
    via ``--configfile``. This entry point honors the merged ``run_config`` value
    instead of assuming the Snakefile's default configfile is the active pointer.
    """
    root = Path(project_root).resolve() if project_root else Path.cwd().resolve()

    run_config_value = pointer_config.get("run_config")
    if not run_config_value:
        raise ConfigError(f"{configfile_label} must define run_config")

    run_config_path = _resolve_project_path(run_config_value, root)
    run_config = load_yaml(run_config_path)
    config_root = (
        Path(fallback_config_root)
        if fallback_config_root is not None
        else _config_root_from_run_config(run_config_path)
    )

    resolved: dict[str, Any] = {}
    loaded_profiles: list[str] = []
    for profile_path, profile in load_profile_stack(run_config, root):
        resolved = deep_merge(resolved, profile)
        loaded_profiles.append(_display_project_path(profile_path, root))

    tool_profiles = load_tool_profiles(config_root)
    resolved = deep_merge(resolved, {"tools": tool_profiles})
    resolved = deep_merge(resolved, pointer_config)
    resolved = deep_merge(resolved, run_config)

    resolved["_ngsflow"] = {
        "project_root": root.as_posix(),
        "configfile": configfile_label,
        "run_config": _display_project_path(run_config_path, root),
        "loaded_profiles": loaded_profiles,
        "loaded_tool_profiles": sorted(tool_profiles),
    }
    return resolved


def get_results_dir(config: dict[str, Any]) -> str:
    results_dir = config.get("results_dir")
    if not results_dir:
        raise ConfigError("Resolved config is missing results_dir")
    return str(results_dir).rstrip("/")


def write_resolved_config(resolved_config: dict[str, Any], results_dir: str | Path) -> Path:
    """Write the resolved config to results_dir/config/resolved_config.yaml."""
    outdir = Path(results_dir) / "config"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "resolved_config.yaml"
    with outpath.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(resolved_config, handle, sort_keys=False)
    return outpath
