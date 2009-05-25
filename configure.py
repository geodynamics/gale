#!/usr/bin/env python
import platform, sys
import config

#
# Functions for handling global options.
#

def add_options(ctx):
    options = ctx.options
    options.add_option(
        options.Option("prefix", "--prefix", pref_sep="=",
                       help="Specify installation directory.",
                       default="/usr/local"),
        options.Option("build_dir", "--build-dir", pref_sep="=",
                       help="Specify the temporary build directory.",
                       default="build"),
        options.Option("build_type", "--build-type", pref_sep="=",
                       type="enum", enum=["debug", "optimised"],
                       help="Build either debug or optimised libraries.",
                       default="debug"),
        options.Option("shared_libraries", "--with-shared-libs", type="bool",
                       pref_true_flag="yes", pref_false_flag="no", pref_sep="=",
                       help="Enable/disable shared libraries.",
                       default="yes"),
        options.Option("static_libraries", "--with-static-libs", type="bool",
                       pref_true_flag="yes", pref_false_flag="no", pref_sep="=",
                       help="Enable/disable static libraries.",
                       default="yes"),
        options.Option("lib_arch", "--lib-arch", type="enum", enum=["32", "64"],
                       help="Link using 32 or 64 bit libraries."),
        options.Option("with_glucifer", "--with-glucifer", type="bool",
                       pref_true_flag="yes", pref_false_flag="no", pref_sep="=",
                       help="Build the gLucifer visualisation package.",
                       default="yes"),
        options.Option("mem_stats", "--with-memory-stats", type="bool",
                       pref_true_flag="yes", pref_false_flag="no", pref_sep="=",
                       help="Enable/disable memory statistics.",
                       default="no"),
        options.Option("cautious", "--with-cautious-mode", type="bool",
                       pref_true_flag="yes", pref_false_flag="no", pref_sep="=",
                       help="Enable extra input checks, even in heavily-used functions.",
                       default="no"))

def process_options(ctx):
    opt_dict = ctx.option_dict
    if opt_dict.get("build_type", None) == "optimised":
        cc = ctx[config.tools.CCompiler]
        cc.compile_env.append_unique("comflags", "opt3", make_list=True)
        cc.compile_env.append_unique("pp_defs", "NDEBUG", make_list=True)
    else:
        cc = ctx[config.tools.CCompiler]
        cc.compile_env.append_unique("comflags", "opt0", "debug", make_list=True)
        cc.compile_env.append_unique("pp_defs", "DEBUG", make_list=True)

    if opt_dict.get("mem_stats", None) == True:
        cc = ctx[config.tools.CCompiler]
        cc.compile_env.append_unique("pp_defs", "MEMORY_STATS", make_list=True)
    
    if opt_dict.get("cautious", None) == True:
        cc = ctx[config.tools.CCompiler]
        cc.compile_env.append_unique("pp_defs", "CAUTIOUS", make_list=True)

def add_modules(ctx):
    cc = ctx.new_module(config.tools.CCompiler, required=True)
    link = ctx.new_module(config.tools.Linker, required=True)
    libxml2 = ctx.new_module(config.packages.libxml2, required=True)
    mpi = ctx.new_module(config.packages.MPI, required=True)
    bl = ctx.new_module(config.packages.BlasLapack, required=True)
    petsc = ctx.new_module(config.packages.PETSc, required=True)
    petscext = ctx.new_module(config.packages.PETScExt)
    gl = ctx.new_module(config.packages.OpenGL)
    osm = ctx.new_module(config.packages.OSMesa)
    png = ctx.new_module(config.packages.libpng)
    jpeg = ctx.new_module(config.packages.libjpeg)
    fame = ctx.new_module(config.packages.libfame)
    sdl = ctx.new_module(config.packages.SDL)
    py = ctx.new_module(config.packages.Python)
    hdf5 = ctx.new_module(config.packages.HDF5)
    if platform.system() == "Darwin":
        carbon = ctx.new_module(config.packages.Carbon)

#
# Create the context.
#

ctx = config.Context()
ctx.add_options_funcs.append(add_options)
ctx.process_options_funcs.append(process_options)
ctx.add_modules_funcs.append(add_modules)

#
# Configure and dump results.
#

# Configure modules.
ctx.configure()
ctx.print_summary()

#print cc.configs[0].subst("$comflags $pp_defs")

# Dump SCons configuration.
config.exporters.SCons(ctx).export("output.cfg")
