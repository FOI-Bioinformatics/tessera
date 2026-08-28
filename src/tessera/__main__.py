"""Entry point for ``python -m tessera`` (equivalent to the ``tessera`` console script)."""

from __future__ import annotations

from .cli.main import app

if __name__ == "__main__":
    app()
