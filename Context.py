import os, sys, shutil, platform, tempfile, platform
from Module import Module
import utils

class Context(object):

    def __init__(self, keep_tmp_dir=False, log_filename="config.log", log_level=0):
        utils.log.level = log_level
        utils.log.set_filename(log_filename)
        utils.log("Launched with command line:")
        utils.log.indent()
        utils.log(" ".join(sys.argv))
        utils.log.unindent()

        self.modules = []
        self.module_dict = {}
        self.candidates = []
        self.order = []

        self.add_options_funcs = []
        self.process_options_funcs = []
        self.add_modules_funcs = []

        self.options = utils.options.OptionSet()
        self.option_dict = {}

        self._old_wd = os.getcwd()
        self._tmp_dir = ""
        self.keep_tmp_dir = keep_tmp_dir

        self.setup_platform()

    def scan_subdirs(self):
        dir_list = os.listdir(".")
        for dir in dir_list:
            if not os.path.isdir(dir):
                continue
            for path, dirs, files in os.walk(dir):
                if "configure.proj" in files:
                    sys.stdout.write("Found sub-project in %s\n"%repr(path))
                    loc = locals()
                    glo = globals()
                    glo["config"] = __import__("config")
                    execfile(os.path.join(path, "configure.proj"), glo, loc)
                    for k, l in [("add_options", self.add_options_funcs),
                                 ("process_options", self.process_options_funcs),
                                 ("add_modules", self.add_modules_funcs)]:
                        if k in loc.keys():
                            l.append(loc[k])

    def setup_options(self):
        for c in self.add_options_funcs:
            c(self)
        for m in self.order:
            m.setup_options()

    def parse_options(self):
        opts = self.options.parse_option_list(sys.argv[1:])
        self.option_dict = self.options.gather(opts)
        for c in self.process_options_funcs:
            c(self)
        for m in self.order:
            m.parse_options()
            m.process_options()

    def print_summary(self):
        sys.stdout.write("\n")
        for m in self.order:
            if not m.enabled or m.mode == 0:
                continue
            if m.configs:
                sys.stdout.write(m.configs[0].get_summary(m) + "\n")
#             else:
#                 sys.stdout.write("Module %s has no configuration.\n"%repr(m.name))

    #
    # Configuring and access to results.
    #

    def configure(self, *args):

        self.scan_subdirs()

        # Add modules.
        for c in self.add_modules_funcs:
            c(self)

        # Need to parse options for the context ahead of schedule.
        self.options.add_option(
            utils.options.Option("build_report", "--build-report", type="bool",
                                 help="Create a tarball with bug reports."))
        opts = self.options.parse_option_list(sys.argv[1:])
        self.option_dict = self.options.gather(opts)

        # If the user wants to build a bug report, make sure we configure a
        # compressor.
        if "build_report" in self.option_dict:
            import tools
            self.get_module(tools.Compressor)

        utils.log("Beginning configuration.")
        utils.log("Was given %d modules as arguments to configure."%len(args), post_indent=1)
        for m in args:
            tmp = self.new_module(m)
            utils.log("%s"%tmp.name)
        utils.log.unindent()

        self.setup_depends()
        self.order = self.order_modules()

        # Parse any options.
        self.setup_options()
        self.parse_options()

        # Need to setup requirement flags so the plans can
        # determine which modules can have a '--with-*' option.
#         self.setup_requirements()
        if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
            self.print_help()

        else:
            self.setup_modes()
            self.setup()
            self.search()
            self.scan()

            # Now setup the requirement flags again to finalise
            # any additional modules that have been created
            # in 'search' and 'scan'.
#             self.setup_requirements()
#             self.order = self.order_modules()

            utils.log("Configuring %d modules."%len(self.order), post_indent=1)
            for m in self.order:

                # A mode of >0 indicates this module was requested by the user.
                if m.mode == 0:
                    continue

                sys.stdout.write("Configuring %s ... "%m.name)
                sys.stdout.flush()
                r = m.configure()
                if m.enabled:
                    sys.stdout.write("%s.\n"%(r and "success" or "failed"))
                if r:
                    utils.log("Module %s successfully configured." % m.name)
                if not r and m.required:
                    utils.log("Required module %s failed." % m.name)
                    sys.stdout.write(utils.format.box("\nFailed to find a valid configuration of " \
                                                          "module %s. Detailed logs can be found in " \
                                                          "%s.\n"%(m.name, repr(utils.log.filename)), 80))
                    if m.options.options:
                        sys.stdout.write("You can influence search parameters using the following:\n\n")
                        sys.stdout.write(m.options.make_help_string())
                        self.exit()
            utils.log.unindent()

        # Free temporary stuff.
        self.__del__()

    def print_help(self):
        help = ""
        done_flags = []
        for n in self.options._order:
            o = self.options.options[n]
            f = o._vars[0].flgs[0]
            if f in done_flags:
                continue
            done_flags.append(f)
            help += o.make_help_string() + "\n"
        for m in self.modules:
            if m.mode == 0:
                continue
            for n in m.options._order:
                o = m.options.options[n]
                f = o._vars[0].flgs[0]
                if f in done_flags:
                    continue
                done_flags.append(f)
                help += o.make_help_string() + "\n"
        sys.stdout.write(help)

    def setup_depends(self):
        utils.log("Setting up dependencies.", post_indent=1)
        d = []
        c = list(self.modules)
        utils.log("Beginning with %d modules to setup."%len(c))
        while c:
            for m in c:
                m.setup_dependencies()
                m.setup_members()
            d += c
            c = [m for m in self.modules if m not in d]
            utils.log("Found %d new modules to setup dependencies for:"%len(c), post_indent=1)
            for m in c:
                utils.log("%s"%m.name)
            utils.log.unindent()
        utils.log.unindent()

    def setup(self):
        for m in self.modules:
            m.setup()

    def order_modules(self):
        """Here we want to analyse the dependency graph to
        determine what order to configure the modules so to
        avoid any nasties."""
        rem = []
        for m in self.modules:
            deps = [d[0] for d in m.depends]
            rem.append((m, deps))
        tiers = utils.graph.build_tiers(rem)
        return sum(tiers, [])

    def setup_modes(self):
        utils.log("Setting up modes.")
#         r = list(self.modules)
#         while len(r):
#             m = r.pop()
#             if m.mode > 0:
#                 for mem in m.members:
#                     if m.mode > 1 and mem.mode == 0:
#                         mem.mode = 1
#                         if mem not in r:
#                             r.append(mem)
        utils.log.indent()
        for m in self.modules:
            utils.log("%s mode %d"%(m.name, m.mode))
        utils.log.unindent()

    def search(self):
        utils.log("Beginning module search.", post_indent=1)
        for m in self.modules:
            m.search()
        utils.log.unindent()

    def scan(self):
        utils.log("Scanning all candidates.", post_indent=1)
        while 1:
            d = True
            for m in self.modules:
                if m.scan():
                    d = False
            if d:
                break
        utils.log.unindent()

#     def setup_requirements(self):
#         rem = list(self.modules)
#         while len(rem):
#             m = rem.pop()
#             if m.required:
#                 for d, r, c in m.depends:
#                     if r and not d.required:
#                         d.required = True
#                         if d not in rem:
#                             rem.append(d)

    def __getitem__(self, k):
        return self.module_dict[k]

    #
    # Configuration manipulation.
    #

    def new_candidate(self, mods, *args, **kw):
        # Create each module and candidate.
        mods = utils.conv.to_list(mods)
        ncs = []
        insts = []
        for m in mods:
            m = self.get_module(m)
            insts.append(m)
            nc = m.new_candidate(self, *args, **kw)
            nc.add_base_modules(m)
            ncs.append(nc)

        # Merge all the candidates.
        nc = ncs[0]
        for i in xrange(1, len(ncs)):
            mods[0].merge_candidates(nc, ncs[i])

        # See if we can find the same candidate.
        c = self.find_candidate(nc)
        if c is None:
            utils.log("Created new candidate: %s" % nc.id)
            self.candidates.append(nc)
        else:
            utils.log("Found pre-existing candidate: %s" % c.id)
            nc = c
        return nc

    def find_candidate(self, cand):
        for c in self.candidates:
            if c.same_candidate(cand):
                return c
        return None

#     def iter_valid_configs(self):
#         done = []
#         for m in self.modules:
#             for c in m.configs:
#                 if c not in done:
#                     done.append(c)
#                     yield c

    #
    # Make sure modules are singletons.
    #

    def new_module(self, mod, **kw):
        return self.get_module(mod, mode=2, **kw)

    def get_module(self, mod, mode=0, **kw):
        if isinstance(mod, Module):
            return mod
        m = self.module_dict.get(mod, None)
        if m:
            if mode > m.mode:
                m.mode = mode
            return m
        m = mod(self, **kw)
        self.module_dict[mod] = m
        self.modules.append(m)
        m.mode = mode
        return m

    #
    # Misc. stuff.
    #

    def setup_platform(self):
        if platform.system() in ["Linux", "Darwin"]:
            self.tool_dirs = ["/bin", "/usr/bin", "/usr/local/bin"]
            self.shared_lib_prefixes = ["lib"]
            self.shared_lib_exts = [".so"]
            self.static_lib_prefixes = ["lib"]
            self.static_lib_exts = [".a"]
            self.pkg_base_dirs = ["/usr", "/usr/local"]
            self.pkg_inc_dirs = [["include"]]
            self.pkg_lib_dirs = [["lib"]]
            self.pkg_bin_dirs = ["bin"]

            if platform.architecture()[0] == "64bit":
                self.pkg_lib_dirs.insert(0, ["lib64"])

            if platform.system() == "Darwin":
                self.shared_lib_prefixes.append("")
                self.shared_lib_exts.append("")
                self.pkg_inc_dirs.append(["Headers"])
                self.tool_dirs.append("/sw/bin")
                self.shared_lib_exts.insert(0, ".dylib")
                self.pkg_base_dirs.append("/sw")
        else:
            # TODO: Handle other platforms.
            pass

        # Setup standard language extensions.
        self.lang_exts = {"c": [".c"],
                          "cpp": [".cxx", ".cpp", ".cc"],
                          "f77": [".f"],
                          "f90": [".f90"],
                          "f95": [".f95"]}

        # Create library prefixes/suffixes.
        self.lib_prefixes = list(self.static_lib_prefixes)
        for p in self.shared_lib_prefixes:
            if p not in self.lib_prefixes:
                self.lib_prefixes.append(p)
        self.lib_exts = list(self.static_lib_exts)
        for p in self.shared_lib_exts:
            if p not in self.lib_exts:
                self.lib_exts.append(p)

    def make_lib_names(self, lib):
        names = utils.path.make_file_names(lib, self.lib_prefixes, self.lib_exts)
	# If we're using Darwin include the lib name without
        # prefix or suffix to catch any framework libraries.
        if platform.system() == "Darwin":
            names.append(lib)
        return names

    def make_static_lib_names(self, lib):
        return utils.path.make_file_names(lib, self.static_lib_prefixes, self.static_lib_exts)

    def make_shared_lib_names(self, lib):
        return utils.path.make_file_names(lib, self.shared_lib_prefixes, self.shared_lib_exts)
	# If we're using Darwin include the lib name without
        # prefix or suffix to catch any framework libraries.
        if platform.system() == "Darwin":
            names.append(lib)
        return names

    def make_temp_name(self, prefix="tmp", suffix=None, ext=None, len=8):
        return utils.path.make_temp_name(prefix, suffix, ext, len)

    def make_temp_file(self, contents=None, prefix="tmp", suffix=None, ext=None):
        if not self._tmp_dir:
            self._tmp_dir = tempfile.mkdtemp(dir=os.getcwd())
            os.chdir(self._tmp_dir)
        return utils.path.make_temp_file(contents, prefix=prefix, suffix=suffix, ext=ext)
#         kw = {"dir": self._tmp_dir, "text": True}
#         if prefix is not None: kw["prefix"] = prefix
#         if suffix is not None: kw["suffix"] = suffix
#         fd, fn = tempfile.mkstemp(**kw)
#         f = os.fdopen(fd, "w")
#         if contents is not None: f.write(contents)
#         f.close()
#         if contents is None: os.remove(fn)
#         return os.path.basename(fn)

    #
    # Destroy stuff.
    #

    def exit(self):
        self.__del__()
        sys.exit()

    def __del__(self):
        # Change back to project root.
        if self._tmp_dir:
            os.chdir(self._old_wd)

        if "build_report" in self.option_dict:
            report_fn = "report.tgz"
            import tools
            comp = self[tools.Compressor]
            if not comp.configs:
                print "Couldn't find any compression tools, submit report manually."
            else:
                comp = comp.configs[0]
                comp["compress"](comp, source=[utils.log.filename, self._tmp_dir],
                                 target=os.path.join("report" + comp["extension"]))

        if self._tmp_dir and not self.keep_tmp_dir:
            shutil.rmtree(self._tmp_dir)
