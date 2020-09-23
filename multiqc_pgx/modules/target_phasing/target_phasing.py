from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.utils import config
from multiqc.plots import bargraph

import json

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Check if the command line arguments were specified, and specified
        # correctly. Raise an error or warning if this is not the case.
        self.check_command_line()

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
        self.plot_phasing_per_sample()
        self.plot_phasing_per_gene()

    def check_command_line(self):
        """ Make sure the command line arguments are usable """
        # Get the command line arguments for this module
        self.target_genes = config.kwargs['target_genes']
        self.blocklist = config.kwargs['whatshap_blocklist']
        self.samples = config.kwargs['whatshap_sample']

        arguments = [self.target_genes, self.blocklist, self.samples]

        # If there were no target genes specified, we don't have to do anything
        if not self.target_genes:
            raise UserWarning

        # If not all of the module arguments were specified, raise an error
        if not all(arguments):
            msg = ('--target-genes, --whatshap-blocklist and '
                   '--whatshap-sample must all be specified.')
            raise RuntimeError(msg)

        # If we did not get a --sample-name for each --whatshap-blocklist,
        # raise an error
        if len(self.blocklist) != len(self.samples):
            msg = ('Please specify --whatshap-sample for each '
                   '--whatshap-blocklist file')
            raise RuntimeError(msg)

    def parse_blocklist_files(self):
        # For each sample (defined in a blocklist)
        for sample, filename in zip(self.samples, self.blocklist):
            self.whatshap[sample] = dict()
            # For each of the target genes
            for target in self.parse_target_genes():
                # We store a dictionary for every gene
                self.whatshap[sample][target.name] = dict()
                self.genes.append(target.name)
                print(f'{target.name} is {target.chrom}:{target.begin}-{target.end}')
                # We update the phasing of the target based on the blocklist
                self.update_phasing(filename, target)

                # We store the size of each phased and unphased block
                for begin, end, phasing in target.all_regions():
                    self.whatshap[sample][target.name][phasing] = end-begin

    def parse_target_genes(self):
        with open(self.target_genes) as fin:
            for line in fin:
                chrom, begin, end, name = line.strip().split()
                begin = int(begin)
                end = int(end)
                yield Target(chrom, begin, end, name)

    def update_phasing(self, filename, target):
        with open(filename) as fin:
            # If filename is an empty file, we are done
            try:
                header = next(fin).strip().split()
            except StopIteration:
                return

            # Did we get the expected header
            assert header == ['#sample', 'chromosome', 'phase_set', 'from', 'to', 'variants']

            # Start parsing the file
            for line in fin:
                spline = line.strip().split()
                chrom = spline[1]
                begin = int(spline[3])
                end = int(spline[4])
                target.update([(chrom, begin, end)])

    def write_data_files(self):
        self.write_data_file(self.whatshap, 'multiqc_pgx_phasing')

    def plot_phasing_per_sample(self):
        """ Plot the phasing of all genes for each sample """
        pdata = list()
        categories = list()
        for sample in self.whatshap:
            data = self.whatshap[sample]
            print(json.dumps(data, indent=True))
            # We want to have all unphased regions in black
            formatting = dict()
            for gene in data:
                for block in data[gene]:
                    d = dict()
                    d['name'] = block
                    if block.startswith('unphased'):
                        d['color'] = '#000000'
                    if block.startswith('phased'):
                        d['color'] = '#7CB5EC'
                    formatting[block] = d
            pdata.append(data)
            categories.append(formatting)

        configuration = {
            'id': 'multiqc_pgx_by_sample',
            'title': 'Phasing per Sample',
            'anchor': 'multiqc_pgx_phasing_sample',
            'data_labels': [
                {
                    'name': sample,
                    'cpswitch_counts_label': sample
                } for sample in self.whatshap
            ]
        }


        self.add_section(
                name='Phasing per Sample',
                anchor='multiqc_pgx_phasing_sample',
                description=
                """
                    This plot shows the phased and unphased blocks for each
                    gene of interest. You can use the buttons at the top of the
                    plot to switch between different samples. The phased blocks
                    are shown in the order they are on the genome. All unphased
                    blocks are black.
                """,
                helptext=
                """
                    MultiQC interprets the gene names in this plot as sample
                    names, so you can use the **Toolbox** on the right to
                    filter out specific genes. For example, filtering *DPYD*
                    gives a much clearer overview of the phasing of the other
                    genes.
                """,
                plot = bargraph.plot(pdata, categories, configuration))

    def plot_phasing_per_gene(self):
        """ Plot the phasing of samples for each gene """
        pdata = list()
        categories = list()

        # Get the genes of interest
        genes = list(self.whatshap.values())[0]

        # Get the gene data fore each sample
        for gene in genes:
            gene_data = dict()
            for sample in self.whatshap:
                data = self.whatshap[sample][gene]
                gene_data[sample] = data
            print(json.dumps(gene_data, indent=True))

            # We want to have all unphased regions in black
            formatting = dict()
            for sample in gene_data:
                for block in gene_data[sample]:
                    d = dict()
                    d['name'] = block
                    if block.startswith('unphased'):
                        d['color'] = '#000000'
                    if block.startswith('phased'):
                        d['color'] = '#7CB5EC'
                    formatting[block] = d
            pdata.append(gene_data)
            categories.append(formatting)

        configuration = {
            'id': 'multiqc_pgx_by_gene',
            'title': 'Phasing per Gene',
            'anchor': 'multiqc_pgx_phasing_gene',
            'data_labels': [
                {
                    'name': gene,
                    'cpswitch_counts_label': gene
                } for gene in genes
            ]
        }


        self.add_section(
                name='Phasing per Gene',
                anchor='multiqc_pgx_phasing_gene',
                description=
                """
                    This plot shows the phased and unphased blocks for each
                    sample. You can use the buttons at the top of the
                    plot to switch between different genes. The phased blocks
                    are shown in the order they are on the genome. All unphased
                    blocks are black.
                """,
                helptext=
                """
                    You can use the **Toolbox** on the right to
                    filter out specific samples.
                """,
                plot = bargraph.plot(pdata, categories, configuration))


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

    def all_regions(self):
        phased = self.phasing[0] == '+'
        block_start = 0
        phased_count = 0
        unphased_count = 0
        for i in range(len(self.phasing)):
            # If we find a phased block while we were unphased
            if self.phasing[i] == '+' and not phased:
                unphased_count += 1
                yield (block_start + self.begin, i + self.begin, f'unphased-{unphased_count}')
                block_start = i
                phased = True
            # If we find an unphased block while we were phased
            if self.phasing[i] == '-' and phased:
                phased_count += 1
                # We have to offset the start of self
                yield (block_start + self.begin, i + self.begin, f'phased-{phased_count}')
                block_start = i
                phased = False
        # If we are at the end of the Target
        if phased:
            yield(block_start+self.begin, i+1+self.begin, f'phased-{phased_count+1}')
        else:
            yield(block_start+self.begin, i+1+self.begin, f'unphased-{unphased_count+1}')

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
