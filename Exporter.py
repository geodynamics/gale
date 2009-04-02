class Exporter(object):

    def __init__(self, ctx):
        self.ctx = ctx

    def export(self):
        done = []
        for m in self.ctx.order:
            if m.configs and m.mode > 0:
                cfg = m.configs[0]
                m.visit(cfg, self)
            done.append(m)
