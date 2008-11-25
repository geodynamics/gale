import os, shutil
import SConfig
import check_headers

def get_symbols_prototypes(cfg):
    src = ''
    for sym, proto in zip(cfg.inst.syms, cfg.inst.pkg.symbol_prototypes):
        src += (proto % sym) + '\n'
    return src

def get_symbols_calls(cfg):
    # Setup any code required for the symbol check.
    src = ''
    if cfg.inst.pkg.symbol_setup:
        src += cfg.inst.pkg.symbol_setup + '\n'

    # Unpack all the symbols.
    for sym, call in zip(cfg.inst.syms, cfg.inst.pkg.symbol_calls):
        src += (call % sym) + '\n'

    # Include any teardown code.
    if cfg.inst.pkg.symbol_teardown:
        src += cfg.inst.pkg.symbol_teardown + '\n'

    return src

def get_symbols_source(cfg):
    """Build the source code required to check that a set of symbols is
    present in a package installation. We produce two sets of source code:
    the first is the code to build a library with a 'main-like' function
    call. The second is a 'runner' that calls the library's function. We
    do this because we need to check that we can build a library out of
    the symbols, which ensures any cross-dependencies are compatible."""

    # Alias the instance and package.
    inst = cfg.inst
    pkg = inst.pkg

    # Begin with header inclusions and symbol prototypes.
    src = check_headers.get_header_source(cfg)
    src += get_symbols_prototypes(cfg)

    # Setup the main function call.
    src += 'int main(int argc, char* argv[]) {\n'
    src += get_symbols_calls(cfg)
    src += 'return 0;\n}\n'

    return src

def generate_library_paths(cfg, lib):
    lib_name = cfg.inst.pkg.env.subst('${SHLIBPREFIX}' + lib + '${SHLIBSUFFIX}')
    if cfg.inst.lib_dirs:
        for d in cfg.inst.lib_dirs:
            path = os.path.join(cfg.inst.base_dir, d, lib_name)
            yield os.path.abspath(path)
    else:
        yield lib_name

def new_check_shared_exist(cfg):
    """Check for the existance of shared libraries for the given
    configuration only. We want to do this by not caring whether
    dependencies have shared libraries or not. The method we use is
    to build a shared library and link directly against this packages
    shared libraries."""

    # Alias the instance and package.
    inst = cfg.inst
    pkg = inst.pkg

    pkg.ctx.Log('  Checking for existence of shared libraries:\n')

    # If we're configuring the dynamic linker, just return True.
    if pkg.name == 'dl':
        pkg.ctx.Log('    No need for the \'dl\' package.\n')
        return True

    # If we don't have any libraries to check, return True.
    if not cfg.libs:
        pkg.ctx.Log('    No libraries to check!\n')
        return True

    # If this package is shared, it should have a configuration of the
    # 'dl' package as one of it's dependencies.
    dl = None
    for dep in cfg.deps:
        if dep.inst.pkg.name == 'dl':
            dl = dep
            break
    if not dl:
        raise 'Error: No dynamic linker as a dependency!'

    # Begin setting up the shared library source by getting the symbols
    # function.
    src = check_headers.get_header_source(cfg) + '\n'
    src += get_symbols_prototypes(cfg) + '\n'
    src += 'void lib_main( int argc, char* argv[] ) {\n'
    src += get_symbols_calls(cfg)
    src += '}\n'

    # Try building a library, but make sure when we enable the system we
    # don't use the regular '-l' library options, we need to actually
    # specify the shared library names directly.
    if pkg.name == 'PETScExt':
        import pdb
        pdb.set_trace()
    old_state = {}
    dep_cfg.enable(pkg.env, old_state, abs_path=True, shared=True)

    return True

def check_shared_exist(cfg):
    """Run a sanity check on shared libraries to see if they exist."""

    # Alias the instance and package.
    inst = cfg.inst
    pkg = inst.pkg

    pkg.ctx.Log('  Checking for existence of shared libraries:\n')

    # If we're configuring the dynamic linker, just return True.
    if pkg.name == 'dl':
        pkg.ctx.Log('    No need for the \'dl\' package.\n')
        return True

    # If we don't have any libraries to check, return True.
    if not cfg.libs:
        pkg.ctx.Log('    No libraries to check!\n')
        return True

    # If this package is shared, it should have a configuration of the
    # 'dl' package as one of it's dependencies.
    dl = None
    for dep in cfg.deps:
        if dep.inst.pkg.name == 'dl':
            dl = dep
            break
    if not dl:
        raise 'Error: No dynamic linker as a dependency!'

    # Build a binary to try and dynamically open all the libraries that are
    # indicated to be shared.
    result = [1, '', '']
    src = check_headers.get_header_source(dl)
    src += """
int main(int argc, char* argv[]) {
  void* lib[%d];
""" % len(cfg.libs)

    # Need to reverse the list of libraries to account for potential global
    # run-time loading bindings.
    libs = pkg.env.reverse_list(cfg.libs)
    for l in libs:
        if pkg.shared_libraries and l not in pkg.shared_libraries:
            continue
        if l in inst.extra_libs:
            continue

        offs = ''
        for p in generate_library_paths(cfg, l):
            offs += '  '
            src += '%slib[%d] = dlopen("%s", RTLD_LAZY | RTLD_GLOBAL);\n' % (offs, cfg.libs.index(l), p)
            src += '%sif( !lib[%d] ) {\n' % (offs, cfg.libs.index(l))
        src += offs + '  printf( "%s", dlerror() );\n'
        src += offs + '  return 1;\n'
        src += offs + '}\n'
        while len(offs) > 2:
            offs = offs[:-2]
            src += offs + '}\n'
    src += '  return 0;\n}\n'

    old_state = {}
    dl.enable(pkg.env, old_state)
    result = cfg.inst.pkg.run_source(src)
    pkg.env.Replace(**old_state)

    if not result[0]:
        pkg.ctx.Log('    Failed.\n')
        return False

    pkg.ctx.Log('    Success.\n')
    return True

def run_dependency_check(cfg, dep_cfg, dl, use_dep=False):
    # Alias some stuff.
    pkg = cfg.inst.pkg
    dep_pkg = dep_cfg.inst.pkg

    # Setup the source code for the initialisation library.
    lib1_src = check_headers.get_header_source(cfg)
    lib1_src += '\nvoid init( int argc, char* argv[] ) {\n'
    lib1_src += dep_pkg.init_code + '\n'
    lib1_src += '}\n\nvoid fina() {\n'
    lib1_src += dep_pkg.fina_code + '\n'
    lib1_src += '}\n'

    # Enable the main package for building the initialisation library. We need to
    # skip enabling the dependency package; that's the point of this check, to see
    # if the main package already has a connection to a shared library of this
    # dependency.
    old_state = {}
    if use_dep:
        cfg.enable(pkg.env, old_state, abs_path=True)
    else:
        cfg.enable(pkg.env, old_state, abs_path=True, lib_exclude=dep_pkg)

    # Build the initialisation library and grab it's path.
    result = pkg.library_source(lib1_src)
    if not result[0]:
        pkg.ctx.Log('      No shared libraries.\n')
        return False
    init_lib = pkg.ctx.sconf.lastTarget.abspath

    # Disable the main package for building the check library.
    pkg.env.Replace(**old_state)

    # Setup the source code to check if this dependency has been
    # initialised.
    lib2_src = check_headers.get_header_source(dep_cfg)
    lib2_src += '\nint check() {\n'
    lib2_src += dep_pkg.check_code + '\n}\n'

    # Enable the dependency package for building the check library.
    old_state = {}
    dep_cfg.enable(pkg.env, old_state)

    # Build the check library and grab it's path.
    result = pkg.library_source(lib2_src)
    if not result[0]:
        raise 'Broken'
    check_lib = pkg.ctx.sconf.lastTarget.abspath

    # Disable the secondary package for building the loader binary.
    pkg.env.Replace(**old_state)

    # Setup the code for the loader.
    ldr_src = """#include <stdlib.h>
#include <stdio.h>
#include <dlfcn.h>

int main( int argc, char* argv[] ) {
   void* lib1;
   void* lib2;
   void (*init)( int, char** );
   void (*fina)();
   int (*check)();

   lib1 = dlopen( \"""" + init_lib + """\", RTLD_NOW );
   if( !lib1 ) {
      printf( "lib1 open failed\\n" );
      printf( "Error: %s\\n", dlerror() );
      return EXIT_SUCCESS;
   }
   init = dlsym( lib1, "init" );
   if( !init ) {
      printf( "init sym failed\\n" );
      printf( "Error: %s\\n", dlerror() );
      return EXIT_SUCCESS;
   }
   fina = dlsym( lib1, "fina" );
   if( !fina ) {
      printf( "fina sym failed\\n" );
      printf( "Error: %s\\n", dlerror() );
      return EXIT_SUCCESS;
   }

   lib2 = dlopen( \"""" + check_lib + """\", RTLD_NOW );
   if( !lib2 ) {
      printf( "lib2 open failed\\n" );
      printf( "Error: %s\\n", dlerror() );
      return EXIT_SUCCESS;
   }
   check = dlsym( lib2, "check" );
   if( !check ) {
      printf( "check sym failed\\n" );
      printf( "Error: %s\\n", dlerror() );
      return EXIT_SUCCESS;
   }

   init( argc, argv );
   if( check() )
      printf( "Connected.\\n" );
   else
      printf( "Disconnected.\\n" );
   fina();

/*
   if( lib1 )
      dlclose( lib1 );
   if( lib2 )
      dlclose( lib2 );
*/

   return EXIT_SUCCESS;
}
"""

    # Enable the dynamic linker package for the loader.
    old_state = {}
    dl.enable(pkg.env, old_state)

    # Build and run the loader. Unfortunately, because of the way
    # SCons works, I think we need to rebuild this everytime.
    result = pkg.run_source(ldr_src)

    # Disable the main package.
    pkg.env.Replace(**old_state)

    return result

def check_shared_dependencies(cfg):
    # Alias the instance and package.
    inst = cfg.inst
    pkg = inst.pkg

    pkg.ctx.Log('  Check the shared library dependencies:\n')

    # Don't check for shared libraries if we're currently checking the
    # dynamic linker package.
    if pkg.name == 'dl':
        pkg.ctx.Log('    No need for the \'dl\' package.\n')
        return True

    # If this package is shared, it should have a configuration of the
    # 'dl' package as one of it's dependencies.
    dl = None
    for dep in cfg.deps:
        if dep.inst.pkg.name == 'dl':
            dl = dep
            break
    if not dl:
        raise 'No dynamic linker as a dependency!'

    # Check every dependency to see if it works.
    for dep_cfg in cfg.deps:

        # Alias the dependency's package.
        dep_pkg = dep_cfg.inst.pkg

        pkg.ctx.Log('    Checking dependency \'%s\'\n' % dep_pkg.name)

        # Make sure the dependency needs to be checked.
        if not (dep_pkg.init_code and dep_pkg.fina_code and dep_pkg.check_code):
            pkg.ctx.Log('      No code for dependency check, skipping.\n')
            continue

        # Run the test without the dependent package.
        result = run_dependency_check(cfg, dep_cfg, dl)

        # If the link failed it means we have no shared libraries for this
        # package.
        if not result[0]:
            return False

        # Check if we were able to open the initialisation library.
        if result[1].find('lib1 open failed') != -1:

            # If not, we need to know why.
            error = result[1].split('\n')[1]

            # If we can find any of the symbols we were trying to link in the
            # error output, this means the main library is not connected to
            # any library that satsifies the dependency. Theoretically, this
            # means we can link against any dependency installation, however
            # this is not really true. So, the answer is to just try it with
            # the current dependency to see if it works.
            if error.find('undefined') != -1: # TODO: search for symbols
                pkg.ctx.Log('      Symbols not present, trying again with dependency enabled.\n')

                # Yep, bogus. Try it with the dependency thrown in.
                result = run_dependency_check(cfg, dep_cfg, dl, use_dep=True)

                # If it failed for any reason, we can't use this combination.
                if not result[0] or \
                        result[1].find('Error') != -1 or \
                        result[1].find('Disconnected') != -1:
                    pkg.ctx.Log('      Failed: Incompatible libraries.\n')
                    return False

        # If the libraries are disconnected, it means this configuration is
        # invalid. We need to use one with a different dependency.
        if result[1].find('Disconnected') != -1:
            pkg.ctx.Log('      Failed: Shared libraries are not connected.\n')
            return False

        pkg.ctx.Log('      Success: Shared libraries are connected.\n')

    return True

def check_files(cfg):
    """Check for the existence of static and shared library files. Return results
    as a tuple of two booleans, the first for static and the second for shared."""

    # Alias the installation and package for easy access.
    inst = cfg.inst
    pkg = inst.pkg

    pkg.ctx.Log('  Checking for existence of files:\n')

    # Set the result flag to false just in case none of the libraries are
    # in our standard set.
    static_found = False

    # Check each library in the current set.
    for lib in cfg.libs:

        pkg.ctx.Log('    Looking for static %s:\n' % lib)

        # Don't check for the existance of a library if it's considered
        # an extra library. This is because extra libraries can often exist
        # in only compiler known locations.
        if lib in pkg.extra_libraries + inst.extra_libs:
            pkg.ctx.Log('      Is an auxilliary library, skipping.\n')
            continue

        # Check each library directory for the files.
        static_found = False
        for lib_dir in inst.lib_dirs:

            # Combine the library sub-directories with the base directory
            # and both the static and shared library names.
            name = pkg.env.subst('${LIBPREFIX}' + lib + '${LIBSUFFIX}')
            path = os.path.join(inst.base_dir, lib_dir, name)

            pkg.ctx.Log('      Trying %s ... ' % path)

            # Check if either of the paths exist.
            static_found = os.path.exists(path)

            # If so, we can break out of the directory loop and move
            # on to the next library.
            if static_found:
                pkg.ctx.Log('found.\n')
                break

            pkg.ctx.Log('not found.\n')

        # If we couldn't find this library, the configuration is invalid.
        # We can return negative straight away.
        if not static_found: break

    # Set the result flag to false just in case none of the libraries are
    # in our standard set.
    shared_found = False

    # Check each library in the current set.
    for lib in cfg.libs:

        pkg.ctx.Log('    Looking for shared %s:\n' % lib)

        # Don't check for the existance of a library if it's considered
        # an extra library. This is because extra libraries can often exist
        # in only compiler known locations.
        if lib in pkg.extra_libraries + inst.extra_libs:
            pkg.ctx.Log('      Is an auxilliary library, skipping.\n')
            continue

        # Don't check for libraries that are known to not be shared.
        if pkg.shared_libraries is not None and lib not in pkg.shared_libraries:
            pkg.ctx.Log('      Is not required to be a shared library.\n')
            continue

        # Check each library directory for the files.
        shared_found = False
        for lib_dir in inst.lib_dirs:

            # Combine the library sub-directories with the base directory
            # and both the static and shared library names.
            name = pkg.env.subst('${SHLIBPREFIX}' + lib + '${SHLIBSUFFIX}')
            path = os.path.join(inst.base_dir, lib_dir, name)

            pkg.ctx.Log('      Trying %s ... ' % path)

            # Check if either of the paths exist.
            shared_found = os.path.exists(path)

            # If so, we can break out of the directory loop and move
            # on to the next library.
            if shared_found:
                pkg.ctx.Log('found.\n')
                break

            pkg.ctx.Log('not found.\n')

        # If we couldn't find this library, the configuration is invalid.
        # We can return negative straight away.
        if not shared_found: break

    # Return results.
    return (static_found, shared_found)

def check_libraries(cfg):
    # Alias the installation and package for easy access.
    inst = cfg.inst
    pkg = inst.pkg

    pkg.ctx.Log('Checking libraries: %s\n' % str(cfg.libs))

    # If there are no libraries to check, automatically pass, but store
    # unknown values for results.
    if not cfg.libs and not inst.fwork:
        pkg.ctx.Log('  No libraries to check!\n')
        cfg.has_static = None
        cfg.has_shared = None
        return True

    # Try and locate files that look like the library files we need.
    if not inst.fwork:
        cfg.has_static, cfg.has_shared = check_files(cfg)

    # If we couldn't find either, we're broken and should return immediately.
    if not inst.fwork and not (cfg.has_static or cfg.has_shared):
        pkg.ctx.Log('  Could not find files for either static or shared libraries.\n')
        return False

    # Now we need to find a valid set of symbols by trying to build a test
    # program. The symbol type required for any configuration is dependent on
    # the installation, not in it's dependencies. So, if we've already
    # found the right set of symbols for this installation, just use those,
    # don't bother searching.
    pkg.ctx.Log('  Checking for symbols:\n')
    if inst.syms is not None:
        pkg.ctx.Log('    Installation already has symbols defined.\n')
        pre_exist = True
        symbols = [(inst.syms, inst.sym_def)]
    else:
        pre_exist = False

        #If there aren't any symbols defined in the package,
        # use an empty set.
        if not pkg.symbols:
            pkg.ctx.Log('    No symbols given, using an empty set.\n')
            symbols = [([], '')]
        else:
            pkg.ctx.Log('    Searching package defined symbols.\n')
            symbols = pkg.symbols

    # Try out each set of symbols.
    for inst.syms, inst.sym_def in symbols:

        pkg.ctx.Log('    Trying symbol set: %s\n' % str(inst.syms))

        # Build the source required for the symbols check.
        src = get_symbols_source(cfg)

        # Enable the environment and build the source.
        old_state = {}
        cfg.enable(pkg.env, old_state)
        result = pkg.run_source(src)
        pkg.env.Replace(**old_state)

        # Check the results.
        if result[0]:

            # In addition to the results indicating success, we need to check
            # the output for any warnings that may indicate failure.
            if result[1].find('skipping incompatible') != -1 or \
                    result[2].find('skipping incompatible') != -1:

                # The library path specified wasn't actually used, instead the
                # libraries in there were incompatible and the compiler was
                # able to make use of libraries in default search paths.
                pkg.ctx.Log('      Failed: compiler reported incompatible libraries.\n')
                result[0] = 0
                continue

            # We're done, break out of the symbols search.
            pkg.ctx.Log('      Success.\n')
            break

        pkg.ctx.Log('      Failed.\n')

    # If we couldn't find a valid set of symbols we need to clear the symbols
    # stored on the installation only if they had not been found before.
    if not result[0]:
        if not pre_exist:
            inst.syms = None
            inst.sym_def = ''

        # TODO: We need to set the correct error here.
        return False

    # Check for the existence of shared libraries if we were able to find shared
    # library-like files.
#    if not inst.fwork and cfg.has_shared:
#        cfg.has_shared = new_check_shared_exist(cfg)

    # Also check if we know of a dependency for which the package was linked
    # against for each dependency listed.
    if not inst.fwork and cfg.has_shared:
        cfg.has_shared = check_shared_dependencies(cfg)

    # If we've gotten this far and we're using a framework, we know we have
    # the right stuff.
    if inst.fwork:
        cfg.has_static = True
        cfg.has_shared = True
        cfg.libs = None

    # If we don't have all the necessary libraries, return negative.
    if pkg.require_shared and not cfg.has_shared:
        pkg.ctx.Log('  Failed: shared libraries are required.\n')
        return False

    pkg.ctx.Log('  Success.\n')
    return True
