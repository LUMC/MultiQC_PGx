from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.utils import config

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # If the command line parameters were not specified, we are done
        # instantly
        self.target_genes = config.kwargs['target_genes']
        self.blocklist = config.kwargs['whatshap_blocklist']
        if not (self.target_genes and self.blocklist):
            raise UserWarning

        # Initialse the parent object
        super(MultiqcModule, self).__init__(name='PGx', anchor='pgx',
                info=(
                    'The PGx module is used by the LUMC PGx project to '
                    'visualise phasing data for the target genes.'
                )
        )

        self.whatshap = dict()
        self.parse_blocklist_files()

    def parse_blocklist_files(self):
        for filename in self.blocklist:
            sample = self.get_sample(filename)
            for target in self.parse_target_genes():
                print(f'{target.name} is {target.chrom}:{target.begin}-{target.end}')
                self.update_phasing(filename, target)
                exit()

    def parse_target_genes(self):
        with open(self.target_genes) as fin:
            for line in fin:
                chrom, begin, end, name = line.strip().split()
                begin = int(begin)
                end = int(end)
                yield Target(chrom, begin, end, name)

    def get_sample(self, filename):
        with open(filename) as fin:
            # skip header
            next(fin)
            return next(fin).split()[0]

    def update_phasing(self, filename, target):
        with open(filename) as fin:
            header =next(fin).strip().split()
            assert header == ['#sample', 'chromosome', 'phase_set', 'from', 'to', 'variants']
            for line in fin:
                spline = line.strip().split()
                chrom = spline[1]
                begin = int(spline[3])
                end = int(spline[4])
                #print(chrom, begin, end)
                target.update([(chrom, begin, end)])



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
                # The phasing of the first part stays the same
                first = self.phasing[:begin-self.begin]
                # Then we add the newly phased part, which is smaller than wat
                # is left
                size = end - begin
                phased = '+' * size
                # Then the rest of the phasing
                rest = len(first) + size
                last = self.phasing[rest:]
                self.phasing = first + phased + last

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
        assert old_length == len(self.phasing), f'Error in {region}'
