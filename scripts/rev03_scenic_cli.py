#!/usr/bin/env python3
"""rev03 — pySCENIC 0.12.1 CLI shim (NumPy-2 compatibility) + entry point.

pySCENIC 0.12.1 (2022) still references the removed NumPy 2.0 aliases `np.object`
and `np.float` in `pyscenic/transform.py`, `pyscenic/rss.py` and
`pyscenic/diptest.py`, so every sub-command dies at import time on this
environment (NumPy 2.4). Both aliases map 1:1 onto the Python builtins, so
restoring them before importing pySCENIC is behaviour-preserving.

Sub-commands:
    grn     -> pyscenic.cli.arboreto_with_multiprocessing (GRNBoost2, step 1)
    <other> -> pyscenic.cli.pyscenic (ctx, aucell, ...)

Unlike the sibling wrapper scripts of Paper 1, this shim is invoked through the
interpreter (`python scripts/rev03_scenic_cli.py grn ...`) so it does not depend on
the console-script shebangs, which point at a now-relocated venv path.

The canonical invocation, with all parameters, is recorded in scripts/rev03_run_scenic.sh.
"""
import sys

import numpy as np

for _alias, _target in (("object", object), ("float", float)):
    if not hasattr(np, _alias):
        setattr(np, _alias, _target)


def main(argv):
    if not argv:
        raise SystemExit("usage: rev03_scenic_cli.py {grn|ctx|aucell|...} [args...]")
    if argv[0] == "grn":                                   # GRNBoost2 (arboreto)
        from pyscenic.cli.arboreto_with_multiprocessing import main as grn_main
        # arboreto's entry point parses sys.argv itself and takes no arguments,
        # so drop our own sub-command token before handing over.
        del sys.argv[1]
        return grn_main()
    from pyscenic.cli.pyscenic import main as cli_main      # ctx / aucell / ...
    return cli_main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
