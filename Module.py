import sys
from Config import Config
import utils

class Module(object):

    def __init__(self, ctx, **kw):
        self.name = self.__module__
        self.options = utils.options.OptionSet()
        self.option_dict = {}
        self.mode = 0
        self.ctx = ctx
        self.required = kw.get("required", False)
        self.candidates = []
        self.configs = []
        self.depends = []
        self.dependees = []
        self.members = []
        self.member_of = []
        self.help = {}

        self.do_search = True
        self.enabled = True

        self._done = False
        self._to_scan = []

    @staticmethod
    def new_candidate(ctx, id, *args, **kw):
        return Config(ctx, id, *args, **kw)

    def copy_candidate(self, cand):
        return cand.clone()

    @staticmethod
    def merge_candidates(c0, c1):
        c0.validate_modules(c1.valid_modules)
        c0.add_base_modules(c1.base_modules)
        c0.add_cand_modules(c1.cand_modules)
        for k, v in c1.iteritems():
            if k not in c0:
                c0[k] = v

    def merge_trial(self, cand, trial):
        self.merge_candidates(cand, trial)
        to_del = []
        for k, v in cand.iteritems():
            if k not in trial:
                to_del.append(k)
        for k in to_del:
            del cand[k]

    @staticmethod
    def same_candidates(c0, c1):
        if len(c0.base_modules) != len(c1.base_modules):
            return False
        for m in c0.base_modules:
            if m in c1.base_modules:
                continue
            found = False
            for m2 in c1.base_modules:
                if m2.has_member(m):
                    found = True
                    break
            if not found:
                return False
        for k in c0.iterkeys():
            if k not in c1:
                return False
        for k in c1.iterkeys():
            if k not in c0:
                return False
#             if k not in exclude and v0 != c1[k]:
#                 return False
        return True

    def setup_dependencies(self):
        utils.log("Setting up dependencies for %s."%self.name, post_indent=1)
        self._setup_dependencies()
        utils.log.unindent()

    def _setup_dependencies(self):
        pass

    def setup_members(self):
        utils.log("Setting up members for %s."%self.name, post_indent=1)
        self._setup_members()
        utils.log.unindent()

    def _setup_members(self):
        pass

    def setup(self):
        pass

    def setup_trial(self, trial):
        trial["apply_env"] = self.apply_env

    def search(self):
        if not self.do_search:
            utils.log("Not searching %s."%(repr(self.name)))
            return
        utils.log("Searching for candidates for %s."%self.name, post_indent=1)
        self._search()
        utils.log.unindent()

    def _search(self):
        pass

    def scan(self):
        utils.log("Scanning module %s."%self.name, post_indent=1)
        f = False
        while len(self._to_scan):
            c = self._to_scan.pop(0)
            if self.scan_candidate(c):
                f = True
        utils.log.unindent()
        return f

    def scan_candidate(self, cand):
        utils.log("Scanning candidate: %s"%cand.id, post_indent=1)
        fnd_more = self._scan_candidate(cand)
        utils.log.unindent()
        return fnd_more

    def _scan_candidate(self, cand):
        return False

    def test(self, cand):
        utils.log("Running tests for module %s."%self.name, post_indent=1)
        trial = self.copy_candidate(cand)
        self.setup_trial(trial)
        if self._test(trial):
            utils.log("Passed.")
            self.merge_trial(cand, trial)
        else:
            utils.log("Failed.")
        utils.log.unindent()

    def combine_depends(self):
        req_deps = [d.configs for d, r, c in self.depends if r and c and d.configs]
        opt_deps = [d.configs for d, r, c in self.depends if not r and c and d.configs]
        deps = req_deps + opt_deps
        for i in range(len(req_deps), len(deps) + 1):
            for p in utils.perm.combine2(deps, i):
                yield p

#     def count_candidates(self):
#         n = len(self.candidates)
#         for m in self.members:
#             n += len(m.configs)
#         return n

#     def _gen_cfgs(self):
#         for m in self.members:
#             for c in m.configs:
#                 if not self.has_config(c):
#                     yield c
#         for c in self.candidates:
#             if not self.has_config(c):
#                 yield c

    #
    # Handle dependencies.
    #

    def add_dependency(self, mod, required=False, combine=False, mode=1):
        m = self.ctx.get_module(mod, mode)
        d = self._find_dep(mod)
        if d is not None:
            d[1] = required
            d[2] = combine
        else:
            self.depends.append([m, required, combine])
            m.dependees.append(self)
        return m

    def add_member(self, mod):
        m = self.ctx.get_module(mod)
        if not self.has_member(m):
            self.members.append(m)
            m.member_of.append(self)
#         self.add_dependency(m, False, False)
        return m

    def _find_dep(self, mod):
        for d in self.depends:
            if d[0] == mod:
                return d
        return None

    def set_mode(self, mode):
        if mode > self.mode:
            self.mode = mode

    def apply_env(self, env):
        pass

    #
    # Handle configuration of candidates/uses.
    #

    def configure(self):
        utils.log("Configuring module %s."%self.name, post_indent=1)
        if self._done:
            utils.log("Already done, skipping.", post_indent=-1)
            return len(self.configs) > 0
        self._done = True
        if not self.enabled:
            sys.stdout.write("disabled.\n")
            utils.log("Disabled, skipping.", post_indent=-1)
            return False
        res = self._configure()
        utils.log.unindent()
        return res

    def _configure(self):
        res = self.config_depends()
        if res:
            self.config_members()
            self.config_cands()
        return len(self.configs) > 0

    def config_depends(self):
        utils.log("Configuring %d dependencies."%len(self.depends), post_indent=1)
        res = self._config_depends()
        utils.log.unindent()
        return res

    def _config_depends(self):
        res = True
        for d in self.depends:
            r = d[0].configure()
            if not r and not d[1]:
                utils.log("Module dependency failed, but not required.");
            if not r and d[1]:
                utils.log("Required module dependency failed.")
                res = False
                break
        return res

    def config_members(self):
        utils.log("Configuring %d members."%len(self.members), post_indent=1)
        self._config_members()
        utils.log.unindent()

    def _config_members(self):
        for m in self.members:
            if m.mode > 0:
                m.configure()
            else:
                m.pass_members()

    def pass_members(self):
        if not self.enabled:
            return
        for m in self.members:
            m.pass_members()
            for cfg in m.configs:
                if cfg.is_valid(self):
                    utils.log("Adding member configuration %s."%repr(cfg))
                    self.add_config(cfg)

    def generate_candidates(self):
        if not self.enabled:
            return
        for cand in self.candidates:
            yield cand
        for m in self.members:
            for cand in m.generate_candidates():
                yield cand

    def config_cands(self):
        utils.log("Configuring candidates.", post_indent=1)
        self._config_cands()
        utils.log.unindent()

    def _config_cands(self):
        self.pass_members()
        if not self.configs:
            cands = list(self.generate_candidates())
            if not cands:
                utils.log("No candidates found!")
                return False

            for cand in cands:
                if cand in self.configs:
                    continue
                utils.log("Configuring candidate %s."%cand.id, post_indent=1)
                cand.test()
                if cand.is_valid(self):
                    utils.log("Valid.", pre_indent=-1)
                    self.add_config(cand)
                    break
                else:
                    utils.log("Invalid.", pre_indent=-1)

    #
    # Configuration manipulation.
    #

    def add_candidate(self, cand):
        if self.members:
            for m in self.members:
                m.add_candidate(cand)
        elif not self.has_candidate(cand):
            self.candidates.append(cand)
            self._to_scan.append(cand)
            cand.add_cand_modules(self)
        else:
            utils.log("Not adding '%s' to '%s' (already exists)."%(repr(cand), self.name))

    def has_candidate(self, cand):
        for c in self.candidates:
            if c.same_candidate(cand):
                return True
        return False

    def add_config(self, cfg):
        if not self.has_config(cfg):
            self.configs.append(cfg)

    def has_config(self, cfg):
        for c in self.configs:
            if c.same_candidate(cfg):
                return True
        return False

    def has_member(self, mod):
        if mod in self.members:
            return True
        for m in self.members:
            if m.has_member(mod):
                return True
        return False

    def get_summary_items(self, cfg):
        return []

    #
    # Help and options.
    #

    def setup_options(self):
        # Create the option name.
        i = self.name.rfind(".")
        if i != -1:
            self.opt_name = self.name[i + 1:].lower()
        else:
            self.opt_name = self.name.lower()

        self.make_help_dict()

        # Add options.
        if "with" not in self.options:
            self.options.add_option(
                utils.options.Option("with", "--with-" + self.opt_name, type="bool",
                                     pref_true_flag="yes", pref_false_flag="no", pref_sep="=",
                                     help=self.help["enable"]))

    def parse_options(self):
        opts = self.options.parse_option_list(sys.argv[1:])
        self.option_dict = self.options.gather(opts)

    def process_options(self):
        w = self.get_option("with", None)
        if w is not None:
            if w:
                self.set_mode(1)
                self.required = True
                self.enabled = True
            else:
                self.enabled = False

    def get_option(self, option, default):
        if option in self.option_dict:
            return self.option_dict[option]
        for m in self.member_of:
            if option in m.option_dict:
                return m.option_dict[option]
        return default

    def make_help_dict(self):
        self.help["enable"] = "Enable %s."%self.opt_name

    def print_summary(self):
        if self.configs:
            cfg = self.configs[0]
            sys.stdout.write(utils.format.box(cfg.id, 80, indent_size=2) + "\n")

        else:
            pass
