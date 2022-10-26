#!/usr/bin/env python3

"""
MultiQC_PGx is a plugin for MultiQC developed for the LUMC PGx
project. It adds additional statistics related to the phasing of the target
regions, which are not suitable to add as full MultiQC modules.
"""

from setuptools import setup, find_packages

version = '0.1.2'

setup(
    name = 'multiqc_PGx',
    version = version,
    author = 'Redmar van den Berg',
    description = 'MultiQC plugin for the LUMC PGx project',
    long_description = __doc__,
    keywords = 'bioinformatics',
    license = 'MIT',
    packages = find_packages(),
    include_package_data = True,
    install_requires = [
        'multiqc'
    ],
    entry_points = {
        'multiqc.cli_options.v1': [
            'target_genes = multiqc_pgx.cli:target_genes',
            'whatshap_blocklist = multiqc_pgx.cli:whatshap_blocklist',
            'whatshap_sample = multiqc_pgx.cli:whatshap_sample',
        ],
        'multiqc.modules.v1': [
            'target_phasing = multiqc_pgx.modules.target_phasing:MultiqcModule'
        ],
        'multiqc.hooks.v1': [
            'before_config = multiqc_pgx.modules.target_phasing:add_fake_file_pattern'
        ]

    }
)
