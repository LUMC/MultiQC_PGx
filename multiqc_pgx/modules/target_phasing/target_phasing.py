from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.utils import config

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        print(config.kwargs['header'])
