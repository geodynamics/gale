class Exporter(object):

    def __init__(self, ctx):
        self.ctx = ctx
        self.global_env = ctx.option_dict

    def export(self):
        self.modules = list(self.ctx.order)
        while len(self.modules):
            m = self.modules.pop()
            if m.configs and m.mode > 0 and m.enabled:
                cfg = m.configs[0]
                m.visit(cfg, self)
