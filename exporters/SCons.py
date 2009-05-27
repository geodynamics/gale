import platform
import config
import config.utils as utils

class File(object):

    def __init__(self, f):
        self.f = f

    def __repr__(self):
        return "File('" + self.f + "')"

class Command(object):

    def __init__(self, s):
	self.s = s

    def __repr__(self):
        return self.s

class SCons(config.Exporter):

    def __init__(self, ctx):
        config.Exporter.__init__(self, ctx)
        self.env = config.Environment()

    def export(self, fn):
        self.env.clear()
        self.export_global_env()
        config.Exporter.export(self)

        if "LINKFLAGS" in self.env:
            self.env.append("LINKFLAGS", "$__RPATH")

        out_f = open(fn, "w")
        for k, v in self.env.iteritems():
            out_f.write(k + " = " + repr(v) + "\n")
        out_f.close()

    def export_global_env(self):
        for k, v in self.global_env.iteritems():
            self.env.append_unique(k, v)

    def ccompiler(self, mod, cfg):
        # Store the compiler so we can figure out whether to include
        # run-time libraries for packages.
        self._com = cfg

        self.env["CC"] = cfg["command"]
        self.env.append_unique("CPPDEFINES", cfg.subst("$pp_defs").split())
        self.env.append_unique("CFLAGS", cfg.subst("$bool($comflags)").split())
        self.env.append_unique("CFLAGS", cfg.subst("$dcomflags").split())

        # Need to handle Darwin.
        if platform.system() == "Darwin":
            self.env["_RPATH"] = ""
            self.env.append("SHLINKFLAGS", "-dynamiclib", "-flat_namespace", "-single_module",
                            "-undefined", "suppress", "-install_name", "${_abspath(TARGET)}")
            self.env.append("_abspath", Command("lambda x: File(x).abspath"))

    def linker(self, mod, cfg):
        self.env["LINKER"] = cfg["command"]
        arch = self.global_env.get("lib_arch", None)
        if arch == "32":
            self._add_opt(cfg["options"], "LINKFLAGS", ("32bit", True))
        elif arch == "64":
            self._add_opt(cfg["options"], "LINKFLAGS", ("64bit", True))

        type = self.global_env.get("build_type", None)
        if type == "debug":
            self._add_opt(cfg["options"], "LINKFLAGS", ("debug", True))

        self.env.append_unique("LINKFLAGS", cfg.subst("$dlnkflags").split())

    def package(self, mod, cfg):
        # Should really iterate over package compilers and languages.
        if hasattr(self, "_com"):
            com = self._com
        else:
            com = cfg["com"].keys()[0]
        com_env = cfg["com"][com]
        lang = com["languages"][0]

        # If we have the compiler use it, otherwise pick any of them.
        if com in cfg["lnk"]:
            lnk = com
        else:
            lnk = cfg["lnk"].keys()[0]

        # Make the environement.
        local = com.clone()
        cfg["apply_compile"](cfg, local, com, lang)

        # Produce the compile command line.
        cmd_str = local.subst("$compile_cmd_" + lang)
        opt_dict = com["options"].gather(com["options"].parse_option_string(cmd_str))

        self.env.append_unique("CPPDEFINES", make_list=True,
                               *utils.conv.to_list(opt_dict.get("pp_def", [])))
        self.env.prepend_unique("CPPPATH", make_list=True,
                                *utils.conv.to_list(opt_dict.get("inc_dir", [])))

        # Produce program link command line.
        local = lnk.clone()
        cfg["apply_link"](cfg, local, lnk)

        cmd_str = local.subst("$link_prog_cmd")
        opt_dict = com["options"].gather(com["options"].parse_option_string(cmd_str))

        self.env.prepend_unique("LIBPATH", make_list=True,
                                *utils.conv.to_list(opt_dict.get("lib_dir", [])))
        self.env.prepend_unique("LIBS", make_list=True,
                                *utils.conv.to_list(opt_dict.get("lib", [])))
        self.env.prepend_unique("RPATH", make_list=True,
                                *utils.conv.to_list(opt_dict.get("rpath", [])))
        lib_paths = [File(d) for d in cfg["lib_paths"]]
        self.env.prepend_unique("LIBS", make_list=True,
                                *utils.conv.to_list(lib_paths))

        # Add any special link options.
        self.env.append_unique("LINKFLAGS", cfg.subst("$lnkprogflags").split())

        # Add some names to identify which packages have
        # been found. Not sure I'll keep this one.
        self.env.append_unique("found_packages", mod.name[mod.name.rfind(".") + 1:],
                               make_list=True)

    def _add_opt(self, opts, scons_name, opt_args):
        if opt_args is None:
            return
        if not isinstance(opt_args, list):
            opt_args = [opt_args]
        flg = opts.make_option_string(opt_args)
        self.env.prepend_unique(scons_name, *flg.split())
