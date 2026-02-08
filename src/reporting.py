"""
Report generation utilities for the RNA-seq pipeline.

Provides a MarkdownReport builder class and an optional HTML renderer.
"""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

logger = logging.getLogger(__name__)


class MarkdownReport:
    """Incrementally build a Markdown report document."""

    def __init__(self, title: str = "RNA-seq Analysis Report"):
        self._lines: List[str] = []
        self.title = title
        self._lines.append(f"# {title}\n")
        self._lines.append(
            f"*Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}*\n"
        )

    # ------------------------------------------------------------------
    # Section helpers
    # ------------------------------------------------------------------

    def h2(self, text: str) -> "MarkdownReport":
        """Add a level-2 heading."""
        self._lines.append(f"\n## {text}\n")
        return self

    def h3(self, text: str) -> "MarkdownReport":
        """Add a level-3 heading."""
        self._lines.append(f"\n### {text}\n")
        return self

    def h4(self, text: str) -> "MarkdownReport":
        """Add a level-4 heading."""
        self._lines.append(f"\n#### {text}\n")
        return self

    def paragraph(self, text: str) -> "MarkdownReport":
        """Add a paragraph."""
        self._lines.append(f"\n{text}\n")
        return self

    def bullet(self, text: str) -> "MarkdownReport":
        """Add a bullet point."""
        self._lines.append(f"- {text}")
        return self

    def code_block(self, code: str, lang: str = "") -> "MarkdownReport":
        """Add a fenced code block."""
        self._lines.append(f"\n```{lang}")
        self._lines.append(code)
        self._lines.append("```\n")
        return self

    def link(self, text: str, url: str) -> "MarkdownReport":
        """Add a Markdown link on its own line."""
        self._lines.append(f"[{text}]({url})")
        return self

    def horizontal_rule(self) -> "MarkdownReport":
        self._lines.append("\n---\n")
        return self

    # ------------------------------------------------------------------
    # Tables
    # ------------------------------------------------------------------

    def table(
        self,
        headers: Sequence[str],
        rows: Sequence[Sequence[Any]],
        alignments: Optional[Sequence[str]] = None,
    ) -> "MarkdownReport":
        """Add a Markdown table.

        Parameters
        ----------
        headers : sequence of str
            Column headers.
        rows : sequence of sequences
            Table rows (one inner sequence per row).
        alignments : sequence of str, optional
            Per-column alignment: ``"l"``, ``"r"``, or ``"c"``.
            Defaults to left-aligned.
        """
        if not headers:
            return self

        n = len(headers)
        if alignments is None:
            alignments = ["l"] * n

        # Header row
        self._lines.append("\n| " + " | ".join(str(h) for h in headers) + " |")

        # Separator row
        seps: List[str] = []
        for a in alignments:
            if a == "r":
                seps.append("---:")
            elif a == "c":
                seps.append(":---:")
            else:
                seps.append("---")
        self._lines.append("| " + " | ".join(seps) + " |")

        # Data rows
        for row in rows:
            cells = [str(c) for c in row]
            # Pad / truncate to match header count
            while len(cells) < n:
                cells.append("")
            self._lines.append("| " + " | ".join(cells[:n]) + " |")

        self._lines.append("")
        return self

    # ------------------------------------------------------------------
    # File-level table from TSV
    # ------------------------------------------------------------------

    def table_from_tsv(self, tsv_path: Union[str, Path]) -> "MarkdownReport":
        """Read a TSV file and render it as a Markdown table."""
        tsv_path = Path(tsv_path)
        if not tsv_path.exists():
            self.paragraph(f"*(File not found: {tsv_path})*")
            return self

        with open(tsv_path, "r", encoding="utf-8") as fh:
            lines = [l.strip() for l in fh if l.strip()]

        if not lines:
            self.paragraph("*(Empty table)*")
            return self

        headers = lines[0].split("\t")
        rows = [l.split("\t") for l in lines[1:]]
        return self.table(headers, rows)

    # ------------------------------------------------------------------
    # Config provenance
    # ------------------------------------------------------------------

    def add_config_provenance(self, cfg: Dict[str, Any]) -> "MarkdownReport":
        """Add a provenance section with the full config dump."""
        import yaml

        self.h2("Provenance")
        self.paragraph("Full configuration used for this analysis:")
        # Remove internal keys
        clean = {k: v for k, v in cfg.items() if not k.startswith("_")}
        self.code_block(yaml.dump(clean, default_flow_style=False, sort_keys=True), "yaml")

        if "_config_hash" in cfg:
            self.paragraph(f"**Config SHA-256:** `{cfg['_config_hash']}`")

        return self

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def render(self) -> str:
        """Return the full Markdown document as a string."""
        return "\n".join(self._lines) + "\n"

    def write(self, path: Union[str, Path]) -> Path:
        """Write Markdown to a file, return the path."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.render(), encoding="utf-8")
        logger.info("Wrote report -> %s", path)
        return path


# ---------------------------------------------------------------------------
# Optional HTML rendering
# ---------------------------------------------------------------------------

def render_html(md_path: Union[str, Path], html_path: Optional[Union[str, Path]] = None) -> Path:
    """Convert a Markdown file to HTML.

    Requires the ``markdown`` package.  Falls back gracefully if not installed.

    Parameters
    ----------
    md_path : path-like
        Input Markdown file.
    html_path : path-like, optional
        Output HTML file.  Defaults to same stem with ``.html``.

    Returns
    -------
    Path
        The written HTML file path.
    """
    md_path = Path(md_path)
    if html_path is None:
        html_path = md_path.with_suffix(".html")
    else:
        html_path = Path(html_path)

    try:
        import markdown as md_lib
    except ImportError:
        logger.warning(
            "Python 'markdown' package not installed; skipping HTML rendering. "
            "Install with: pip install markdown"
        )
        return md_path  # return the .md path as fallback

    md_text = md_path.read_text(encoding="utf-8")
    html_body = md_lib.markdown(md_text, extensions=["tables", "fenced_code"])

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>RNA-seq Report</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
           max-width: 960px; margin: 2em auto; padding: 0 1em; line-height: 1.6; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1em 0; }}
    th, td {{ border: 1px solid #ddd; padding: 0.5em 0.8em; text-align: left; }}
    th {{ background: #f5f5f5; }}
    code {{ background: #f5f5f5; padding: 0.15em 0.3em; border-radius: 3px; }}
    pre {{ background: #f5f5f5; padding: 1em; overflow-x: auto; border-radius: 5px; }}
    h1, h2, h3 {{ color: #333; }}
  </style>
</head>
<body>
{html_body}
</body>
</html>"""

    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(html_doc, encoding="utf-8")
    logger.info("Wrote HTML report -> %s", html_path)
    return html_path
