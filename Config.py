from utils.Environment import Environment
import utils

class Config(Environment):

    def __init__(self, ctx, id, **kw):
        Environment.__init__(self, **kw)
        self.ctx = ctx
        self.id = id
        self.base_modules = []
        self.cand_modules = []
        self.valid_modules = []
        self.done = False

    def dup(self):
        return Config(self.ctx, self.id)

    def test(self):
        if not self.done:
            self.done = True
            for m in self.cand_modules:
                m.test(self)

    def validate_modules(self, mods):
        mods = utils.conv.to_list(mods)
        for m in mods:
            m = self.ctx.get_module(m)
            if m not in self.valid_modules:
                self.valid_modules.append(m)

#     def validate_modules(self, mods):
#         mods = utils.conv.to_list(mods)
#         for m in mods:

#             m = self.ctx.get_module(m)
#             if m in self.valid_modules:
#                 raise "Validating for overlapping modules."
#             self.valid_modules[m] = m
#             if m not in self._all_valid_mods:
#                 self._all_valid_mods.append(m)

#             cur = m._member_of
#             while len(cur):
#                 new_set = []
#                 for vm in cur:
#                     if vm not in self.valid_modules:
#                         self.valid_modules[vm] = m
#                     if vm not in self._all_valid_mods:
#                         self._all_valid_mods.append(vm)
#                     for mem in vm._member_of:
#                         if mem not in new_set:
#                             new_set.append(mem)
#                 cur = new_set

#         import pdb
#         pdb.set_trace()
#         print "Huzzah"

    def is_valid(self, mod):
        return self.are_all_valid(mod)

    def are_any_valid(self, mods):
        mods = utils.conv.to_list(mods)
        for m in mods:
            m = self.ctx.get_module(m)
            if m in self.valid_modules:
                return True
        return False

    def are_all_valid(self, mods):
        mods = utils.conv.to_list(mods)
        for m in mods:
            m = self.ctx.get_module(m)
            if m not in self.valid_modules:
                return False
        return True

    def add_base_modules(self, mods):
        mods = utils.conv.to_list(mods)
        for m in mods:
            if m not in self.base_modules:
                self.base_modules.append(m)

    def add_cand_modules(self, mods):
        mods = utils.conv.to_list(mods)
        for m in mods:
            if m not in self.cand_modules:
                self.cand_modules.append(m)

    def same_candidate(self, cand):
        for m in self.base_modules:
            if not m.same_candidates(self, cand):
                return False
        return True

    def print_summary(self):
        pass

    def print_fail_reason(self):
        pass

    def __repr__(self):
        return self.id
