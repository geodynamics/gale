import platform
import config
from Compiler import Compiler

class DLTool(config.Tool):

    def __init__(self, ctx):
        config.Tool.__init__(self, ctx)

    def _setup_dependencies(self):
        self.com = self.add_dependency(Compiler, required=True, mode=0)

    def setup_trial(self, trial):
        config.Tool.setup_trial(self, trial)
        cfg["list_libraries"] = self.list_libraries
        cfg["get_libraries"] = self.get_libraries

    def _setup_members(self):
        self.add_member(config.tools.ldd)
        if platform.system() == "Darwin":
            self.add_member(config.tools.otool)

    def list_libraries(self, cfg, args):
        return self.execute(cfg, args)

    def get_libraries(self, out):
        return {}
