#!/usr/bin/env python3

import click

target_genes = click.option(
        '--target-genes',
        help='BED file of genes of interest')

whatshap_blocklist = click.option(
        '--whatshap-blocklist',
        help='Phased regions produced by WhatHap using --block-list')
