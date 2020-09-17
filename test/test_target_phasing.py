#!/usr/bin/env python3

import pytest
import sys

# This line allows the tests to run if you just naively run this script.
# But the preferred way is to use run_tests.sh
sys.path.insert(0,'../MultiQC_PGx')

from multiqc_pgx.modules.target_phasing import Target

# target, phased_blocks, result
TARGETS = [
        (Target('chr1', 0, 5, 'test'), [], '-----'),
        (Target('chr1', 0, 5, 'test'), [('chr2', 0, 5)], '-----'),
        (Target('chr1', 0, 5, 'test'), [('chr1', 0, 5)], '+++++'),
        (Target('chr1', 0, 5, 'test'), [('chr1', 0, 1)], '+----'),
        (Target('chr1', 0, 5, 'test'), [('chr1', 0, 10)], '+++++'),
        (Target('chr1', 0, 5, 'test'), [('chr1', 5, 10)], '-----'),
        (Target('chr1', 0, 5, 'test'), [('chr1', 0, 2), ('chr1', 2, 5)], '+++++'),
        (Target('chr1', 0, 5, 'test'), [('chr1', 0, 2), ('chr1', 3, 5)], '++-++'),
        (Target('chr1', 0, 5, 'test'), [('chr1', 3, 5), ('chr1', 0, 2)], '++-++'),
        (Target('chr1', 5, 10, 'test'), [('chr1', 0, 20)], '+++++'),
        (Target('chr1', 5, 10, 'test'), [('chr1', 0, 7)], '++---'),
]

PHASED = [
        (Target('chr1', 5, 10, 'test'), '+++++', [(0,5)]),
        #(Target('chr1', 0, 5, 'test'), [('chr1', 3, 5)], [(3, 5)]),
        #(Target('chr1', 0, 5, 'test'), [('chr1', 3, 5),('chr1', 3, 5)], [(3, 5)]),
]

def test_target():
    T = Target('chr1', 10, 20, 'test')
    assert T.name == 'test'
    assert T.chrom == 'chr1'
    assert T.phasing == '-'*10


@pytest.mark.parametrize(['target', 'phased', 'result'], TARGETS)
def test_phased_target(target, phased, result):
    target.update(phased)
    assert str(target) == result
