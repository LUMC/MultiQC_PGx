#!/usr/bin/env python3

import click

target_genes = click.option(
        '--target-genes',
        type=click.Path(exists=True),
        help='BED file of genes of interest')

whatshap_blocklist = click.option(
        '--whatshap-blocklist',
        multiple=True,
        type=click.Path(exists=True),
        help='Phased regions produced by WhatHap using --block-list')

whatshap_sample = click.option(
        '--whatshap-sample',
        multiple=True,
        type=str,
        help='Sample names, in the same order as --whatshap-blocklist')
