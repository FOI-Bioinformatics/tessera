"""Entry point for ``python -m tessera``.

Delegating here rather than relying on ``python -m tessera.cli.main`` is deliberate.
Running a module with ``-m`` executes it as ``__main__``, so the ``cmd_*`` modules'
``from .main import app`` imports a *second*, separate ``tessera.cli.main`` module
object and registers every subcommand on that copy's ``app`` -- leaving the
``__main__`` copy with no commands at all. Importing the module normally from here
means there is only ever one ``app``.
"""

from __future__ import annotations

from .main import app

if __name__ == "__main__":
    app()
