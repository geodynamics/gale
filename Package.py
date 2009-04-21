import os, glob, platform
import config
from Module import Module
import utils

class Package(Module):

    def __init__(self, ctx, **kw):
        Module.__init__(self, ctx, **kw)
        self.patterns = []
        self.base_dirs = []
        self.forced_base_dirs = []
        self.inc_exts = []
        self.lib_exts = []
        self.forced_inc_dirs = []
        self.forced_lib_dirs = []
        self.forced_libs = []
        self.forced_lib_paths = []
        self._lib_sets = []
        self._aux_lib_cands = {}
        self._aux_inc_dirs = []
        self._aux_lib_dirs = []
        self._aux_rpaths = []
        self._rt_libs = {}
        self._rt_lib_dirs = []
        self._rt_rpaths = []

    @staticmethod
    def new_candidate(ctx, base_dir, inc_dirs, lib_dirs, incs, lib_paths):
        inc_dirs = utils.conv.to_list(inc_dirs)
        lib_dirs = utils.conv.to_list(lib_dirs)
        incs = utils.conv.to_list(incs)
        lib_paths = utils.conv.to_list(lib_paths)
        id = str([base_dir, inc_dirs, lib_dirs, incs, lib_paths])
        cfg = config.Config(ctx, id, base_dir=base_dir,
                            inc_dirs=inc_dirs, lib_dirs=lib_dirs,
                            incs=incs, lib_paths=lib_paths)
        return cfg

    def copy_candidate(self, cand):
        cpy = Module.copy_candidate(self, cand)
        return cpy

    @staticmethod
    def merge_candidates(cand, trial):
        # TODO: For now we just overwrite. Should be more
        # complicated than this.
        if "valid_deps" in cand:
            del cand["valid_deps"]
        if "com" in cand:
            del cand["com"]
        if "lnk" in cand:
            del cand["lnk"]
        Module.merge_candidates(cand, trial)

    @staticmethod
    def same_candidates(c0, c1):
        if not Module.same_candidates(c0, c1):
            return False
        if os.path.realpath(c0["base_dir"]) != os.path.realpath(c1["base_dir"]) or \
                len(c0["inc_dirs"]) != len(c1["inc_dirs"]) or \
                len(c0["lib_dirs"]) != len(c1["lib_dirs"]) or \
                c0["incs"] != c1["incs"] or \
                len(c0["lib_paths"]) != len(c1["lib_paths"]):
            return False
        mine = set([os.path.realpath(d) for d in c0["inc_dirs"]])
        for d in c1["inc_dirs"]:
            if os.path.realpath(d) not in mine:
                return False
        mine = set([os.path.realpath(d) for d in c0["lib_dirs"]])
        for d in c1["lib_dirs"]:
            if os.path.realpath(d) not in mine:
                return False
        mine = set([os.path.realpath(d) for d in c0["lib_paths"]])
        for d in c1["lib_paths"]:
            if os.path.realpath(d) not in mine:
                return False
        return True

    def _setup_dependencies(self):
        self.com = self.add_dependency(config.tools.Compiler, required=True, mode=0)
        self.lnk = self.add_dependency(config.tools.Linker, required=True, mode=0)
        self.setup_libraries()

    def setup_libraries(self):
        pass

    def combine_depends(self):
        for deps in Module.combine_depends(self):
            sets = []
            for d in deps:
                vd = d.get("valid_deps", [])
                if vd:
                    sets.append(vd)
            if sets:
                for more_deps in utils.perm.combine2(sets):
                    all_deps = list(deps)
                    for ds in more_deps:
                        for d in ds:
                            if d not in all_deps:
                                all_deps.append(d)
                    yield all_deps
            else:
                yield list(deps)


    def setup_trial(self, trial):
        Module.setup_trial(self, trial)
        trial["apply_compile"] = self.apply_compile
        trial["apply_link"] = self.apply_link
        trial["make_compile_env"] = self.make_compile_env
        trial["make_link_env"] = self.make_link_env


    def _test(self, cfg):
        # Store the configuration temporaryily.
        self._cfg = cfg

        # Setup the build environment.
        env = config.Environment()
        env["run_test_prog"] = self.test_run
        env.update(cfg)
        env["rpaths"] = cfg["lib_dirs"]

        # Try to build and run tests.
        utils.log("Testing libraries: %s"%repr(cfg["lib_paths"]), post_indent=1)
        output = {}
        res = config.checks.check_libraries(self.ctx, env, output,
                                            self.com.configs, self.lnk.configs,
                                            self.apply_compile_dependency, self.apply_link_dependency,
                                            self.gen_dependencies, self.gen_compilers,
                                            self.gen_linkers, self.gen_sources,
                                            self._aux_lib_cands)
        utils.log.unindent()

        # Eliminate the configuration.
        del self._cfg

        if res:

            # Copy the compiler results to a general location.
            cfg.update(output)

            cfg.validate_modules(self)
            return True
        return False


    def gen_dependencies(self):
        dep_cands = list(self.combine_depends())
        for deps in dep_cands:
            yield deps


    def gen_compilers(self):
        for com in self.com.configs:
            yield com


    def gen_linkers(self):
        for lnk in self.lnk.configs:
            yield lnk


    def gen_sources(self, comp, lang):
        try:
            srcs = self.source_code[lang]
            srcs = utils.conv.to_list(srcs)
            for s in srcs:
                yield s
        except:
            pass


    def apply_compile_dependency(self, dep, com, lang, env):
        # If this dependency was unable to be built using the current compiler,
        # skip the compiler altogether.
        if com not in dep["com"]:
            return False
        dep["make_compile_env"](dep, env, com, lang)
        return True


    def apply_link_dependency(self, dep, lnk, env):
        # If this dependency was unable to be built using the current compiler,
        # skip the compiler altogether.
        dep["make_link_env"](dep, env, lnk)
        return True


    def test_run(self, cfg, fn):
        utils.log("Trying to run test code:")
        utils.log.indent()
        cmd_line = self.make_run_command_line(cfg, fn, cfg["_deps"])
        utils.log(cmd_line)
        res, out, err = utils.command.run(cmd_line)
        if res: # TODO: bloody fortran, can't check 'err'
            utils.log("Failed, error code %d, output and error follows:"%res)
            utils.log.indent()
            if out:
                utils.log("")
                for l in out.split("\n"):
                    utils.log(l)
            if err:
                utils.log("")
                for l in err.split("\n"):
                    utils.log(l)
            utils.log.unindent()
        else:
            utils.log("Success.")
        utils.log.unindent()
        return not res, {"err_code": res, "stdout": out, "stderr": err}

    def make_run_command_line(self, cfg, fn, deps):
        return os.path.join(".", fn)

    def find_incs(self, dirs, incs):
        return utils.path.find_all(incs, dirs)

    def find_static_libs(self, dirs, libs):
        return utils.path.find_all(libs, dirs,
                                   self.ctx.static_lib_prefixes,
                                   self.ctx.static_lib_exts)

    def find_shared_libs(self, dirs, libs):
        return utils.path.find_all(libs, dirs,
                                   self.ctx.shared_lib_prefixes,
                                   self.ctx.shared_lib_exts)

    def add_library_set(self, incs, libs, aux_libs=None):
        env = config.Environment(incs=utils.conv.to_list(incs),
                                 libs=utils.conv.to_list(libs))
        self._lib_sets.append(env)
        if aux_libs is not None:
            env["aux_libs"] = utils.conv.to_list(aux_libs)

    def add_auxilliary_libs(self, langs, libs):
        langs = utils.conv.to_list(langs)
        libs = utils.conv.to_list(libs)
        for l in langs:
            if l not in self._aux_lib_cands:
                self._aux_lib_cands[l] = [[]]
            self._aux_lib_cands[l].append(libs)

    def _search(self):
        for cand in self.gen_candidates():
            self.add_candidate(cand)

    def gen_candidates(self):
        utils.log("Generating candidates.", post_indent=1)
        for c in self._gen_candidates():
            yield c
        utils.log.unindent()

    def _gen_candidates(self):
        for ds in self.gen_dir_sets():
            for cfg in self.gen_cands_from_dir_set(ds):
                yield cfg

    def gen_cands_from_dir_set(self, ds):
        for info in self.gen_lib_cands(ds):
            cfg = self.ctx.new_candidate(self, ds[0], ds[1], ds[2],
                                         info["incs"], info["lib_paths"])
            self.finalise_candidate(ds, info, cfg)
            yield cfg

    def finalise_candidate(self, dirs, info, cand):
        pass

    def gen_dir_sets(self):
        for ds in self._gen_dir_sets():
            yield ds

    def _gen_dir_sets(self, **kw):
        for base_dir in self.gen_base_dirs():
            for ds in self.gen_from_base_dir(base_dir):
                yield ds

    def candidate_from_dirs(self, dirs):
        cfg = self.ctx.new_candidate(self, dirs[0], dirs[1], dirs[2])
        return cfg

    def is_candidate(self, cand):
        utils.log("Checking %s:"%repr(cand), post_indent=1)
        res = self._is_candidate(cand)
        if res:
            utils.log("Success.")
        else:
            utils.log("Failure.")
        utils.log.unindent()
        return res

    def _is_candidate(self, cand):
        return self.search_dirs(cand)

    def gen_lib_cands(self, dirs):
        for env in self._lib_sets:

            inc_paths = self.find_incs(dirs[1], env["incs"])
            if inc_paths is None:
                continue

            if self.forced_libs:
                libs = self.forced_libs
            else:
                libs = env["libs"]

            # This could be faster.
            st_lib_paths = []
            sh_lib_paths = []
            for l in libs:
                f = self.find_shared_libs(dirs[2], [l])
                if f is not None:
                    sh_lib_paths.extend(f)
                else:
                    sh_lib_paths.append(None)
                f = self.find_static_libs(dirs[2], [l])
                if f is not None:
                    st_lib_paths.extend(f)
                else:
                    st_lib_paths.append(None)

            # Try the natrual combination of primarily shared libararies and
            # filling in any gaps with static libs.
            lib_paths = []
            for i in xrange(len(libs)):
                if sh_lib_paths[i] is not None:
                    lib_paths.append(sh_lib_paths[i])
                elif st_lib_paths[i] is not None:
                    lib_paths.append(st_lib_paths[i])
                else:
                    break
            if len(lib_paths) == len(libs):
                yield {"incs": env["incs"], "lib_paths": lib_paths}

            # Try all static libs.
            st_lib_paths = [p for p in st_lib_paths if p is not None]
            if len(st_lib_paths) == len(libs) and st_lib_paths != lib_paths:
                yield {"incs": env["incs"], "lib_paths": st_lib_paths}


    def gen_base_dirs(self):
        if self.forced_base_dirs:
            dirs = self.forced_base_dirs
        else:
            dirs = self.base_dirs + self.ctx.pkg_base_dirs
        for d in utils.path.gen_dirs(dirs, self.patterns):
            yield d


    def combine_inc_dirs(self, base_dir, **kw):
        if self.forced_inc_dirs:
            yield self.forced_inc_dirs
            return

        for dir_set in self.gen_inc_dir_sets(base_dir, **kw):
            okay = True
            inc_dirs = [os.path.join(base_dir, h) for h in dir_set]
            for h in inc_dirs:
                if not os.path.exists(h):
                    okay = False
                    break
            if not okay:
                continue
            yield inc_dirs
            for e in self.gen_inc_exts(inc_dirs, **kw):
                found = False
                ext_dirs = []
                for d in inc_dirs:
                    p = os.path.join(d, e)
                    if os.path.exists(p):
                        ext_dirs.append(p)
                        found = True
                    else:
                        ext_dirs.append(d)
                if found:
                    yield inc_dirs + ext_dirs

    def gen_inc_dir_sets(self, base_dir):
        for d in self.ctx.pkg_inc_dirs:
            yield d

    def gen_inc_exts(self, inc_dirs, **kw):
        for e in self.inc_exts:
            yield e

    def combine_lib_dirs(self, base_dir, **kw):
        if self.forced_lib_dirs:
            yield self.forced_lib_dirs
            return

        for dir_set in self.gen_lib_dir_sets(base_dir, **kw):
            okay = True
            lib_dirs = [os.path.join(base_dir, l) for l in dir_set]
            for l in lib_dirs:
                if not os.path.exists(l):
                    okay = False
                    break
            if not okay:
                continue
            yield lib_dirs
            for e in self.lib_exts:
                found = False
                ext_dirs = []
                for d in lib_dirs:
                    p = os.path.join(d, e)
                    if os.path.exists(p):
                        ext_dirs.append(p)
                        found = True
                    else:
                        ext_dirs.append(d)
                if found:
                    yield lib_dirs + ext_dirs


    def gen_lib_dir_sets(self, base_dir):
        yield [''] # Darwin needs the base directory on it's own
                   # to find frameworks.
        for d in self.ctx.pkg_lib_dirs:
            yield d


    def gen_from_base_dir(self, base_dir, **kw):
        lib_dirs_set = [d for d in self.combine_lib_dirs(base_dir, **kw)]
        for inc_dirs in self.combine_inc_dirs(base_dir, **kw):
            for lib_dirs in lib_dirs_set:
                yield [base_dir, inc_dirs, lib_dirs]

    def apply_compile(self, cfg, env, com, lang):
        for d in cfg["valid_deps"][0]:
            self.apply_compile_dependency(d, com, lang, env)
        cfg["make_compile_env"](cfg, env, com, lang)

    def apply_link(self, cfg, env, lnk):
        for d in cfg["valid_deps"][0]:
            self.apply_link_dependency(d, lnk, env)
        cfg["make_link_env"](cfg, env, lnk)

    def make_compile_env(self, cfg, env, com, lang):
        com_env = cfg["com"][com][lang]
        env.append_unique("inc_dirs", make_list=True, *cfg["inc_dirs"])
        env.append_unique(com_env)


    def make_link_env(self, cfg, env, lnk):
        env.append_unique("lib_dirs", make_list=True, *cfg["lib_dirs"])
        env.append_unique("rpaths", *cfg["lib_dirs"])
        env.append_unique("lib_paths", *cfg["lib_paths"])
        env.append_unique(cfg["lnk"][lnk])

    def visit(self, cfg, visitor):
        visitor.package(self, cfg)

    def setup_options(self):
        Module.setup_options(self)

    def setup_options(self):
        Module.setup_options(self)
        self.options.add_option(
            utils.options.Option("base_dir", "--" + self.opt_name + "-dir",
                                 pref_sep="=", help=self.help["base_dir"]),
            utils.options.Option("inc_dir", "--" + self.opt_name + "-inc-dir",
                                 pref_sep="=", delim=","),
            utils.options.Option("lib_dir", "--" + self.opt_name + "-lib-dir",
                                 pref_sep="=", delim=","),
            utils.options.Option("lib", "--" + self.opt_name + "-lib",
                                 pref_sep="="),
            utils.options.Option("lib_path", "--" + self.opt_name + "-lib-path",
                                 pref_sep="="))

    def process_options(self):
        Module.process_options(self)

        base_dir = self.get_option("base_dir", None)
        inc_dirs = self.get_option("inc_dir", None)
        lib_dirs = self.get_option("lib_dir", None)
        libs = self.get_option("lib", None)
        lib_paths = self.get_option("lib_path", None)

        if inc_dirs is not None:
            self.forced_inc_dirs = inc_dirs.split(",")

        if lib_dirs is not None:
            self.forced_lib_dirs = lib_dirs.split(",")

        if self.forced_lib_dirs and self.forced_inc_dirs:
            base_dir = os.path.commonprefix(inc_dirs + lib_dirs)
            ds = [base_dir, inc_dirs, lib_dirs]
            for cand in self.gen_cands_from_dir_set(ds):
                self.add_candidate(cand)
                break
            self.do_search = False

        if libs is not None:
            self.forced_libs = libs.split(",")

        if lib_paths is not None:
            self.forced_lib_paths = lib_paths.split(",")

        if base_dir is not None:
            self.forced_base_dirs = [base_dir]

    def make_help_dict(self):
        Module.make_help_dict(self)
        self.help["base_dir"] ="Specify package location to be used " + \
            "for %s."%repr(self.name)
        self.help["enable"] = "Enable or disable the %s."%repr(self.name)
