#!/usr/bin/env python3

import click

add_header = click.option(
        '--header',
        type=str,
        help='Put this in the multiqc report')
