# -*- coding: utf-8 -*-

"""A command line shaped like the one planned for JCAST on Typer (PR 18): a group with `run` and `gui`, `-s` as a list
of splice types, `-q` as two floats, and a help string with a hand-written [default: ...]. Run by scope_evidence.py in
a virtual environment holding one Typer version; prints the arguments it parsed.
"""

from enum import Enum
from pathlib import Path
from typing import Annotated, Optional

import typer


class SpliceType(str, Enum):
    MXE = 'MXE'
    SE = 'SE'
    RI = 'RI'
    A5SS = 'A5SS'
    A3SS = 'A3SS'


app = typer.Typer(no_args_is_help=True, add_completion=False, context_settings={'help_option_names': ['-h', '--help']})


@app.command()
def run(rmats_folder: Path,
        gtf_file: Path,
        genome: Path,
        read: Annotated[int, typer.Option('-r', '--read',
                                          help='the lowest read count [default: 1, which keeps every junction]')] = 1,
        splice_type: Annotated[Optional[list[SpliceType]], typer.Option('-s', '--splice_type')] = None,
        qvalue: Annotated[tuple[float, float], typer.Option('-q', '--qvalue')] = (0, 1),
        cache_dir: Annotated[Optional[Path], typer.Option('--cache_dir')] = None):
    """Translate splice junctions."""
    print('parsed:', f'read={read}', f'splice_type={[s.value for s in splice_type or []]}', f'qvalue={qvalue}',
          f'cache_dir={cache_dir}')


@app.command()
def gui():
    """Open the graphical interface."""
    print('parsed: gui')


if __name__ == '__main__':
    app()
