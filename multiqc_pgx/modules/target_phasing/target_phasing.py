from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.utils import config

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # If the command line parameters were not specified, we are done
        # instantly
        target_genes = config.kwargs['target_genes']
        blocklist = config.kwargs['whatshap_blocklist']
        if not (target_genes and blocklist):
            raise UserWarning

        print(target_genes)
        print(blocklist)

        # Initialse the parent object
        super(MultiqcModule, self).__init__(name='PGx', anchor='pgx',
                info=(
                    'The PGx module is used by the LUMC PGx project to '
                    'visualise phasing data for the target genes.'
                )
        )


class Target():
    """
    Class to store the target region in. At initialisation, the entire region
    is unphased. The phasing can be updated by passing phased regions to
    Target.
    """
    def __init__(self, chromosome, begin, end, name):
        self.chrom = chromosome
        self.name = name
        self.begin = begin
        self.end = end
        self.phasing = '-' * (end-begin)

    def __repr__(self):
        return self.phasing

    def phased(self):
        phased = False
        phase_start = None
        for i in range(len(self.phasing)):
            # If we find a new phased block
            if self.phasing[i] == '+' and not phased:
                phase_start = i
                phased = True
            # If we come to the end of the phased block
            if self.phasing[i] == '-' and phased:
                # We have to offset the start of self
                yield (phase_start+self.begin, i+self.begin)
                phased = False
        # If we are at the end of the Target and we were in a phased block
        else:
            if phased:
                yield(phase_start+self.begin, i+1+self.begin)

    def update(self, regions):
        """ Update the phasing according to the regions """
        old_length = len(self.phasing)
        for region in regions:
            chrom, begin, end = region
            size = end-begin
            # If the region is on a different chromosome, we are done
            if chrom != self.chrom:
                continue

            # If both begin and end are in the target
            if begin >= self.begin and end <= self.end:
                self.phasing = self.phasing[:begin] + '+'*size + self.phasing[end:]
            # If the beginning is in the target, but the end is not
            if begin >= self.begin and end > self.end:
                # The size of the target region counting from begin
                rest = len(self.phasing)-begin
                self.phasing = self.phasing[:begin] + '+'*rest
            # If region completely overlaps the target
            if begin <= self.begin and end >= self.end:
                self.phasing = '+'*len(self.phasing)
            # If the beginning is before the target, but the ending is in the
            # target
            if begin <= self.begin and end <= self.end:
                rest = end-self.begin
                self.phasing = '+' * rest + self.phasing[rest:]
        assert old_length == len(self.phasing)
