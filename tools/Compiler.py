import os
import config
import config.utils as utils
from config.utils.options import Option

class Compiler(config.Tool):

    def __init__(self, ctx, **kw):
        config.Tool.__init__(self, ctx, **kw)
        self.default_inc_dirs = []
        self.language = ""
        self.source_empty = {}
        self.source_lib = {}

    def setup_defaults(self):
        if os.name in ["posix", "mac"]:
            self.default_inc_dirs = ["/usr/include"]
        else:
            # TODO
            pass

    def _setup_members(self):
        self.cc = self.add_member(config.tools.CCompiler)
        self.f77 = self.add_member(config.tools.F77Compiler)

    def setup_trial(self, cfg):
        config.Tool.setup_trial(self, cfg)
        cfg["compile"] = self.compile
        cfg["compile_cmd"] = "$command $bool(compile) $bool($comflags) $string(output,$target) " \
            "$string(pp_def,$pp_defs) " \
            "$string(inc_dir,$strip($default_inc_dirs,$inc_dirs)) $source"
        cfg["compile_success"] = self.success
        cfg.append_unique("default_inc_dirs", *self.default_inc_dirs)
        cfg.append_unique("languages", self.language, make_list=True)
        if "exts" not in cfg:
            cfg["exts"] = {}
        cfg["exts"][self.language] = self.ctx.lang_exts[self.language][0]

    def _test(self, cfg):
        if not self.test_core_compile_options(cfg):
            return False
        if not self.test_compile_options(cfg):
            return False
        cfg.validate_modules([self, Compiler])
        return True

    def test_core_compile_options(self, cfg):
        lang = self.language
        src_empty = self.source_empty[lang]

        # compile only
        src_fn = self.ctx.make_temp_file(src_empty, ext=cfg["exts"][lang])
        opts = Option("compile", "-c", type="bool")
        if not self.test_option(cfg, opts,
                                args=[("compile", True), src_fn]):
            return False

        # output filename
        src_fn = self.ctx.make_temp_file(src_empty, ext=cfg["exts"][lang])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("output", "-o")
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("compile", True), ("output", out_fn),
                                      src_fn]):
            return False
        cfg["_simple_obj"] = out_fn

        # header directory
        src_fn = self.ctx.make_temp_file(src_empty, ext=cfg["exts"][lang])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("inc_dir", "-I", sep=["", " "])
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("compile", True), ("output", out_fn), ("inc_dir", "."),
                                      src_fn]):
            return False

        # CPP definition
        src_fn = self.ctx.make_temp_file(src_empty, ext=cfg["exts"][lang])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("pp_def", "-D", sep="")
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("compile", True), ("output", out_fn), ("pp_def", "DUMMY"),
                                      src_fn]):
            return False

        # debug only
        src_fn = self.ctx.make_temp_file(src_empty, ext=cfg["exts"][lang])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("debug", "-g", type="bool")
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("compile", True), ("output", out_fn), ("debug", True),
                                      src_fn]):
            return False

        return True

    def test_options(cfg):
        return True

    def compile(self, cfg, lang="", env={}, **kw):
        cmd = "compile_cmd"
        if lang:
            cmd += "_" + lang
        return self.execute(cfg, cmd, env, **kw)

    def parse_rt_libs(self, cfg):
        # Do this for both static and shared linking.
        utils.log("Stripping static run-time libraries.", post_indent=1)
        res, out, err = self.execute(cfg, args=[("dump_cmds", True), cfg["_simple_obj"]])
        self.parse_rt_output(cfg, out, err, "static")
        utils.log.unindent()

        utils.log("Stripping shared run-time libraries.", post_indent=1)
        res, out, err = self.execute(cfg, args=[("shared", True), ("dump_cmds", True),
                                                cfg["_simple_obj"]])
        self.parse_rt_output(cfg, out, err, "shared")
        utils.log.unindent()

    def parse_rt_output(self, cfg, out, err, prefix):
        # Need to run compiler specific code to extract the
        # relevant arguments.
        cmd_args = self.get_link_rt_args(cfg, out, err)

        # Parse into options and store a dict version.
        opts = cfg["options"].parse_option_list(cmd_args)

        # We need to include every library after the object file and
        # every lib dir.
        libs = []
        lib_dirs = []
        rpaths = []
        for o in opts:
            if isinstance(o, tuple):
                if o[0] == "lib_dir":
                    d = os.path.realpath(o[1])
                    if d not in lib_dirs:
                        lib_dirs.append(d)
                elif o[0] == "rpath":
                    d = os.path.realpath(o[1])
                    if d not in rpaths:
                        rpaths.append(d)
        opts.reverse()
        for o in opts:
            if isinstance(o, tuple):
                if o[0] == "lib":
                    if o[1] not in libs:
                        libs.insert(0, o[1])
            else:
                if o == cfg["_simple_obj"]:
                    break
        opts.reverse()

        # Turns out we also need to get hold of the *.o files. We'll need to group
        # them into start and end, that being those before the object file and those
        # after the object file.
        start_objs = []
        end_objs = []
        dst = start_objs
        for o in opts:
            if isinstance(o, tuple):
                if o[0] == "lib":
                    dst = end_objs
            elif o[0] == "-":
                continue
            elif o[-2:] == ".o":
                dst.append(os.path.realpath(o))

        # Store on the compiler configuration.
        if "rt" not in cfg:
            cfg["rt"] = {}
        cfg["rt"][prefix] = {}
        cfg["rt"][prefix]["rt_rpaths"] = rpaths
        cfg["rt"][prefix]["rt_lib_dirs"] = lib_dirs
        cfg["rt"][prefix]["rt_libs"] = libs
        cfg["rt"][prefix]["rt_begin"] = start_objs
        cfg["rt"][prefix]["rt_end"] = end_objs

        # Log details.
        if rpaths:
            utils.log("Found rpaths: %s"%repr(rpaths))
        if lib_dirs:
            utils.log("Found library directories: %s"%repr(lib_dirs))
        if libs:
            utils.log("Found libraries: %s"%repr(libs))
        if start_objs:
            utils.log("Found start objects: %s"%repr(start_objs))
        if end_objs:
            utils.log("Found end objects: %s"%repr(end_objs))

    def get_link_rt_args(self, cfg, out, err):
        return err.splitlines()[-1].split()
