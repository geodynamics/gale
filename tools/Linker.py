import os, platform
import config
import config.utils as utils
from config.utils.options import Option

class Linker(config.Tool):

    def __init__(self, ctx, **kw):
        config.Tool.__init__(self, ctx, **kw)
        self.default_lib_dirs = []
        self.help["location"] = "Specifies the command to use " + \
            "for linking object files with the plan <name>."
        self.setup_defaults()
        self.link_env = config.Environment()

    def setup_defaults(self):
        if os.name in ["posix", "mac"]:
            self.default_lib_dirs = ["/lib", "/usr/lib", "/usr/lib64"]
        else:
            # TODO
            pass

    def _setup_dependencies(self):
        self.com = self.add_dependency(config.tools.Compiler, required=True, mode=0)

    def _setup_members(self):
        if self.__module__ == "config.tools.Linker":
            self.add_member(config.tools.Compiler)
            self.add_member(config.tools.ld)

    def setup_trial(self, cfg):
        config.Tool.setup_trial(self, cfg)
        cfg["link_prog"] = self.link_prog
        cfg["link_shared"] = self.link_shared
        cfg["link_success"] = self.success
        cfg["default_lib_dirs"] = self.default_lib_dirs
        if platform.system() == "Darwin":
            cfg["link_prog_cmd"] = "$command $lnkprogflags $bool($lnkflags) $dlnkflags " \
                "$string(output,$target) " \
                "$__lib_dirs $_rt_begin " \
                "$source $lib_paths $__libs " \
                "$_rt_end"
            cfg["link_shared_cmd"] = "$command $bool(shared) "\
                "$string(install_name,$target_abs) " \
                "$string(undefined,suppress) $bool(flat_namespace) " \
                "$bool($lnkflags) $dlnkflags $string(output,$target) " \
                "$__lib_dirs $source $lib_paths $__libs"
        else:
            cfg["link_prog_cmd"] = "$command $lnkprogflags $bool($lnkflags) $dlnkflags " \
                "$string(output,$target) $__rpaths " \
                "$__lib_dirs $_rt_begin " \
                "$source $lib_paths $__libs " \
                "$_rt_end"
            cfg["link_shared_cmd"] = "$command $bool(shared) $bool($lnkflags) $dlnkflags $string(output,$target) " \
                "$__rpaths $__lib_dirs $source $lib_paths $__libs"

        cfg["_lib_dirs"] = "$extend($strip($default_lib_dirs,$lib_dirs)," \
            "$strip($default_lib_dirs,$rt_lib_dirs))"
        cfg["__lib_dirs"] = "$string(lib_dir,$_lib_dirs)"

        cfg["_libs"] = "$extend($libs,$_rt_libs)"
        cfg["__libs"] = "$string(lib,$_libs)"

        cfg["_rpaths"] = "$extend($rpaths,$rt_rpaths)"
        cfg["__rpaths"] = "$string(rpath,$strip($default_lib_dirs,$_rpaths))"

        cfg["_rt_libs"] = "$strip($rt.static.rt_libs,$rt_libs)"
        cfg["_rt_begin"] = "$strip($rt.static.rt_begin,$rt_begin)"
        cfg["_rt_end"] = "$strip($rt.static.rt_end,$rt_end)"
        cfg.merge(self.link_env)

    def _test(self, cfg):
        if not self.find_obj(cfg):
            return False
        if not self.test_core_link_options(cfg):
            return False
        if not self.test_link_options(cfg):
            return False
        cfg.validate_modules([self, Linker])
        self.parse_version_string(cfg)
        return True

    def test_core_link_options(self, cfg):
        # in order to use the extra run-time libraries for unnatural linkers we
        # need to have -L, -rpath and -l all ready.
        cfg["options"].add_option(Option("lib_dir", "-L", sep=["", " "]))
        cfg["options"].add_option(Option("rpath", "-Wl,-rpath", sep=['=', ',', ' ']))
        cfg["options"].add_option(Option("lib", "-l", sep=""))

        # output filename
        out_fn = self.ctx.make_temp_name()
        opts = Option("output", "-o")
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("output", out_fn)] + cfg.get("lnkprogflags", []) + \
                                    cfg.get("_static_start_objs", []) + \
                                    [cfg["_simple_obj"]] + cfg.get("_static_end_objs", []) + \
                                    cfg.get("_static_rt_args", [])):
            return False
        cfg["_simple_prog"] = out_fn

        return True

    def test_link_options(self, cfg):
        # shared library
        out_fn = self.ctx.make_temp_file(prefix=self.ctx.shared_lib_prefixes[0],
                                         ext=self.ctx.shared_lib_exts[0])
        opts = [Option("shared", "-shared", type="bool"),
                Option("shared", "-dynamiclib", type="bool")]
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("shared", True), ("output", out_fn)] + \
                                    cfg.get("_shared_start_objs", []) + \
                                    [cfg["_pic_obj"]] + cfg.get("_shared_rt_args", []) + \
                                    cfg.get("_shared_end_objs", [])):
            return False
        cfg["_shared_lib"] = out_fn

        # install name (darwin)
        out_fn = self.ctx.make_temp_file(prefix=self.ctx.shared_lib_prefixes[0],
                                         ext=self.ctx.shared_lib_exts[0])
        opts = Option("install_name", "-install_name", sep=[' '])
        self.test_option(cfg, opts, products=out_fn,
                         args=[("shared", True), ("install_name", out_fn), ("output", out_fn)] + \
                             cfg.get("_shared_start_objs", []) + [cfg["_pic_obj"]] + \
                             cfg.get("_shared_end_objs", []) + \
                             cfg.get("_shared_rt_args", []))

        # rpath
        out_fn = self.ctx.make_temp_name()
        opts = Option("rpath", "-Wl,-rpath", sep=['=', ',', ' '])
        self.test_option(cfg, opts, products=out_fn,
                         args=[("output", out_fn), ("rpath", os.getcwd())] + \
                             cfg.get("_static_start_objs", []) + \
                             [cfg["_simple_obj"]] + cfg.get("_static_end_objs", []) + \
                             cfg.get("_static_rt_args", []))

        # library path
        out_fn = self.ctx.make_temp_file(prefix=self.ctx.shared_lib_prefixes[0],
                                         ext=self.ctx.shared_lib_exts[0])
        opts = Option("lib_dir", "-L", sep=["", " "])
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("output", out_fn), ("lib_dir", os.getcwd())] + \
                                    cfg.get("_static_start_objs", []) + [cfg["_simple_obj"]] + \
                                    cfg.get("_static_end_objs", []) + \
                                    cfg.get("_static_rt_args", [])):
            return False

        # library
        # NOTE: Can't check this just yet, have to assume it works.
        if "lib" not in cfg["options"]:
            cfg["options"].add_option(Option("lib", "-l", sep=""))

        return True

    def apply_env(self, cfg, env):
        self.apply_link_env(cfg, env)

    def apply_link_env(self, cfg, env):
        env.merge(self.compile_env)
        env["lnk"] = cfg["command"]

    def link_prog(self, cfg, env={}, **kw):
        local = self.link_env.clone()
        local.merge(env)
        return self.execute(cfg, "link_prog_cmd", local, **kw)

    def link_shared(self, cfg, env, **kw):
        local = self.link_env.clone()
        local.merge(env)
        return self.execute(cfg, "link_shared_cmd", local, **kw)

    def find_obj(self, cfg):
        for com in self.com.configs:
            if "_simple_obj" in com:
                cfg["_simple_obj"] = com["_simple_obj"]
                cfg["_pic_obj"] = com["_pic_obj"]
                rt = com["rt"]
                cfg["_static_start_objs"] = rt["static"]["rt_begin"]
                cfg["_static_end_objs"] = rt["static"]["rt_end"]
                cfg["_static_rt_args"] = [("rpath", p) for p in rt["static"]["rt_rpaths"]] + \
                    [("lib_dir", d) for d in rt["static"]["rt_lib_dirs"]] + \
                    [("lib", l) for l in rt["static"]["rt_libs"]]
                cfg["_shared_start_objs"] = rt["shared"]["rt_begin"]
                cfg["_shared_end_objs"] = rt["shared"]["rt_end"]
                cfg["_shared_rt_args"] = [("rpath", p) for p in rt["shared"]["rt_rpaths"]] + \
                    [("lib_dir", d) for d in rt["shared"]["rt_lib_dirs"]] + \
                    [("lib", l) for l in rt["shared"]["rt_libs"]]
                return True
        return False

    def setup_options(self):
        config.Tool.setup_options(self)
        self.options.add_option(
            utils.options.Option("lflags", "--link-flags", pref_sep="=",
                                 help="Specify additional linker flags."))

    def process_options(self):
        config.Tool.process_options(self)
        flags = self.get_option("lflags", [])
        if flags:
            flags = [flags]
        self.link_env.append_unique("dlnkflags", make_list=True, *flags)

    def parse_version_string(self, cfg):
        pass

    def visit(self, cfg, visitor):
        visitor.linker(self, cfg)
