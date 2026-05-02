# -*- coding: utf-8 -*-
"""
pipeline_utils.py

Tiny utilities shared by scripts in this repository.

Design goals:
- Cross-platform (Windows/macOS/Linux)
- No heavy dependencies
- Clear error messages when required inputs are missing
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class ProjectPaths:
    """Resolved project paths (project-root relative)."""
    root: Path
    inputs: Path
    outputs: Path
    outputs_radclock: Path
    outputs_functional: Path

    def ensure(self) -> None:
        self.inputs.mkdir(parents=True, exist_ok=True)
        self.outputs.mkdir(parents=True, exist_ok=True)
        self.outputs_radclock.mkdir(parents=True, exist_ok=True)
        self.outputs_functional.mkdir(parents=True, exist_ok=True)


def resolve_project_paths(project_root: Optional[str] = None) -> ProjectPaths:
    """
    Resolve project root.

    If project_root is None, assume current working directory is project root.
    """
    root = Path(project_root) if project_root else Path.cwd()
    root = root.resolve()
    paths = ProjectPaths(
        root=root,
        inputs=root / "inputs",
        outputs=root / "outputs",
        outputs_radclock=root / "outputs_radclock",
        outputs_functional=root / "outputs_functional_shielding",
    )
    paths.ensure()
    return paths


def find_first_existing(
    search_dirs: Sequence[Path],
    candidate_names: Sequence[str],
) -> Path:
    """
    Find the first existing file among candidate_names across search_dirs.

    Raises FileNotFoundError with a helpful message if none are found.
    """
    tried: List[str] = []
    for d in search_dirs:
        for name in candidate_names:
            p = (d / name).resolve()
            tried.append(str(p))
            if p.exists():
                return p
    msg = "Required input file not found. Tried:\n" + "\n".join(f"  - {t}" for t in tried)
    raise FileNotFoundError(msg)


def ensure_parent_dir(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
