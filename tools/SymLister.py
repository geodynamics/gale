import config
from Compiler import Compiler

class SymLister(config.Tool):

    def __init__(self, ctx):
        config.Tool.__init__(self, ctx)

    def _setup_dependencies(self):
        self.com = self.add_dependency(Compiler, required=True, mode=0)

    def _setup_members(self):
        self.add_member(config.tools.nm)

    def setup_trial(self, trial):
        config.Tool.setup_trial(self, trial)
        trial["list_symbols"] = self.list_symbols
        trial["get_symbols"] = self.get_symbols

    def list_symbols(self, cfg, args):
        return self.execute(cfg, args)

    def get_symbols(self, out):
        return {}
