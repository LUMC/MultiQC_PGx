from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.utils import config
from multiqc.plots import bargraph

import json

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
        self.genes = list()
        self.parse_blocklist_files()
        self.write_data_files()
        # Try to plot the gene data for a single sample
        sample = 'sample1_pre_opt0_5x'
        data = self.whatshap[sample]
        print(json.dumps(data, indent=True))

        self.add_section(
                name = 'Module section',
                anchor = 'multiqc_pgx_phasing',
                plot = bargraph.plot(data))


    def parse_blocklist_files(self):
        # For each sample (defined in a blocklist)
        for filename in self.blocklist:
            sample = self.get_sample(filename)
            self.whatshap[sample] = dict()
            # For each of the target genes
            for target in self.parse_target_genes():
                # We store a dictionary for every gene
                self.whatshap[sample][target.name] = dict()
                self.genes.append(target.name)
                print(f'{target.name} is {target.chrom}:{target.begin}-{target.end}')
                # We update the phasing of the target based on the blocklist
                self.update_phasing(filename, target)

                # We store the size of each phased block
                phase_counter = 0
                for begin, end in target.phased():
                    phase_counter += 1
                    self.whatshap[sample][target.name][f'phased-{phase_counter}'] = end-begin

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
            header = next(fin).strip().split()
            assert header == ['#sample', 'chromosome', 'phase_set', 'from', 'to', 'variants']
            for line in fin:
                spline = line.strip().split()
                chrom = spline[1]
                begin = int(spline[3])
                end = int(spline[4])
                target.update([(chrom, begin, end)])

    def write_data_files(self):
        self.write_data_file(self.whatshap, 'multiqc_pgx_phasing')


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
                #print(f'{region} is contained within target ({self.chrom}:{self.begin}-{self.end})')
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
            elif begin >= self.begin and end > self.end:
                #print(f'{region} overlaps the end of target ({self.chrom}:{self.begin}-{self.end}))')
                # The size of the target region counting from begin
                rest = len(self.phasing)-begin
                self.phasing = self.phasing[:begin] + '+'*rest

            # If region completely overlaps the target
            elif begin <= self.begin and end >= self.end:
                #print(f'{region} overlaps the target completely ({self.chrom}:{self.begin}-{self.end}))')
                self.phasing = '+'*len(self.phasing)

            # If the beginning is before the target, but the ending is in the
            # target
            elif begin <= self.begin and end <= self.end and end >= self.begin:
                #print(f'{region} overlaps the beginning of target ({self.chrom}:{self.begin}-{self.end}))')
                rest = end-self.begin
                self.phasing = '+' * rest + self.phasing[rest:]
        assert old_length == len(self.phasing), f'Error in {region}'
