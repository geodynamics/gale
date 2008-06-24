import os

class Installation:
    def __init__(self, pkg, base_dir='', hdr_dirs=[], lib_dirs=[], fwork=''):

        # Keep a reference to the package this installation is of.
        self.pkg = pkg

        # The base directory of this installation. Can be empty.
        self.base_dir = base_dir

        # A list of header/library directory extensions to be appended
        # to the base directory. Can be relative and absolute paths.
        self.hdr_dirs = list(hdr_dirs)
        self.lib_dirs = list(lib_dirs)

        # Lists of headers and libraries to be used in addition to the
        # package's description.
        self.hdrs = []
        self.libs = []

        # A list of library names that, if appearing in the above list
        # of libraries, should not be considered a core part of this
        # installation, i.e. should not be checked for.
        self.extra_libs = []

        # The name of a framework to use instead of everything else.
        self.fwork = fwork

        # A list of symbols that are required by this installation and any
        # preprocessor definitions required to identify them.
        self.syms = None
        self.sym_def = ''

        # A list of pre-processor definitions to be set for this installation.
        self.cpp_defines = []

        # We need this flag to indicate whether this installation has
        # been processed by it's owning package.
        self.is_processed = False

    def __eq__(self, inst):
        return (self.pkg is inst.pkg and
                self.base_dir == inst.base_dir and
                self.hdr_dirs == inst.hdr_dirs and
                self.lib_dirs == inst.lib_dirs)

    def add_hdr_dirs(self, hdr_dirs):
        dir_list = self.pkg.env.make_list(hdr_dirs)
        self.hdr_dirs += [d for d in dir_list if d not in self.hdr_dirs]

    def add_lib_dirs(self, lib_dirs, prepend=False):
        dir_list = self.pkg.env.make_list(lib_dirs)
        for d in dir_list:
            d = os.path.normpath(d)
            if d not in self.lib_dirs:
                if prepend:
                    self.lib_dirs = [d] + self.lib_dirs
                else:
                    self.lib_dirs += [d]

    def add_hdrs(self, hdrs):
        hdr_list = self.pkg.env.make_list(hdrs)
        self.hdrs += [h for h in hdr_list if h not in self.hdrs]

    def add_libs(self, libs):
        lib_list = self.pkg.env.make_list(libs)
        self.libs += [l for l in lib_list if l not in self.libs]

    def add_extra_libs(self, libs):
        lib_list = self.pkg.env.make_list(libs)
        self.extra_libs += [l for l in lib_list if l not in self.extra_libs]

    def add_cpp_defines(self, defs):
        def_list = self.pkg.env.make_list(defs)
        self.cpp_defines += [d for d in def_list if d not in self.cpp_defines + ['']]

    def enable(self, scons_env, old_state=None, libs=[], has_shared=True, lib_exclude=[]):
        """Inserts this installation's information into the provided SCons
        environment. Any environment values that are modified are backed
        up into 'old_state' providing there is not already a backup there."""

        # Make sure we have lists.
        lib_exclude = self.pkg.env.make_list(lib_exclude)

        # Insert any pre-processor definitions.
        if self.cpp_defines or self.sym_def:
            cpp_def = list(self.cpp_defines)
            if self.sym_def:
                cpp_def += [self.sym_def]
            self.pkg.backup_variable(scons_env, 'CPPDEFINES', old_state)
            scons_env.AppendUnique(CPPDEFINES=cpp_def)

        # Insert any header file search paths we may have.
        if self.hdr_dirs:
            self.pkg.backup_variable(scons_env, 'CPPPATH', old_state)

            # The ordering of the header search path list is important.
            # Because we insert them into the environment one at a time we
            # need to reverse the list to make sure the order is maintained.
            rev_hdr_dirs = self.pkg.env.reverse_list(self.hdr_dirs)

            # Process each path in turn.
            for d in rev_hdr_dirs:

                # Combine sub-directories to form the complete search path.
                full_dir = os.path.join(self.base_dir, d)

                # If this path is in a predefined list of system specific search
                # paths, then we need to place this path at the end of the
                # list. This way we can be sure that any default installation
                # will not interfere with custom installations.
                if full_dir in self.pkg.system_header_dirs:

                    # If the path is relative, make sure SCons knows it needs
                    # to treat it as relative to the project root.
                    if not os.path.isabs(full_dir): full_dir = '#' + full_dir
                    scons_env.AppendUnique(CPPPATH=[full_dir])
                else:

                    # If the path is relative, make sure SCons knows it needs
                    # to treat it as relative to the project root.
                    if not os.path.isabs(full_dir): full_dir = '#' + full_dir
                    scons_env.PrependUnique(CPPPATH=[full_dir])

        # Insert any library search paths we may have, but not if this package is in
        # out list of exclusions.
        if self.lib_dirs and self.pkg not in lib_exclude and has_shared:
            self.pkg.backup_variable(scons_env, ['LIBPATH', 'RPATH'], old_state)

            # The ordering of the library search path list is important.
            # Because we insert them into the environment one at a time we
            # need to reverse the list to make sure the order is maintained.
            rev_lib_dirs = self.pkg.env.reverse_list(self.lib_dirs)

            # Process each path in turn.
            for d in rev_lib_dirs:

                # Combine sub-directories to form the complete search path.
                full_dir = os.path.join(self.base_dir, d)

                # We need the absolute path for adding rpaths.
                abs_dir = os.path.abspath(full_dir)

                # If this path is in a predefined list of system specific search
                # paths, then we need to place this path at the end of the
                # list. This way we can be sure that any default installation
                # will not interfere with custom installations.
                if full_dir in self.pkg.system_library_dirs:

                    # If the path is relative, make sure SCons knows it needs
                    # to treat it as relative to the project root.
                    if not os.path.isabs(full_dir): full_dir = '#' + full_dir
                    scons_env.AppendUnique(LIBPATH=[full_dir])
                    scons_env.AppendUnique(RPATH=[abs_dir])
                else:

                    # If the path is relative, make sure SCons knows it needs
                    # to treat it as relative to the project root.
                    if not os.path.isabs(full_dir): full_dir = '#' + full_dir
                    scons_env.PrependUnique(LIBPATH=[full_dir])
                    scons_env.PrependUnique(RPATH=[abs_dir])

        # Add libraries if there are any, unless this package is part of our
        # library exclusions.
        if libs and self.pkg not in lib_exclude:

            # If this package is configured to be using static libraries, then
            # we need to specify the path to the library itself.
            if not has_shared:
                self.pkg.backup_variable(scons_env, 'STATICLIBS', old_state)
                libs = self.find_library(libs)
                scons_env.PrependUnique(STATICLIBS=libs)
            else:
                self.pkg.backup_variable(scons_env, 'LIBS', old_state)
                scons_env.PrependUnique(LIBS=libs)

        # If we have a framework, add it now.
        if self.fwork:
            self.pkg.backup_variable(scons_env, 'FRAMEWORKS', old_state)
            scons_env.PrependUnique(FRAMEWORKS=[self.fwork])

    def find_library(self, lib, static=True, shared=False):
        """Using the search paths we know about, try and locate the files corresponding
        to the library name(s) given."""
        
        libs = self.pkg.env.make_list(lib)
        found_libs = []
        if static:
            for l in libs:
                if l in self.pkg.extra_libraries + self.extra_libs:
                    continue
                for lib_dir in self.lib_dirs:
                    name = self.pkg.env.subst('${LIBPREFIX}' + l + '${LIBSUFFIX}')
                    path = os.path.join(self.base_dir, lib_dir, name)
                    if os.path.exists(path):
                        found_libs += [path]
        if shared:
            for l in libs:
                if l in self.pkg.extra_libraries + self.extra_libs:
                    continue
                if self.pkg.shared_libraries is not None and l not in sef.pkg.shared_libraries:
                    continue
                for lib_dir in self.lib_dirs:
                    name = self.pkg.env.subst('${SHLIBPREFIX}' + l + '${SHLIBSUFFIX}')
                    path = os.path.join(self.base_dir, lib_dir, name)
                    if os.path.exists(path):
                        found_libs += [path]
        return found_libs

    def __str__(self):
        """Convert to printable string."""

        txt = 'Package: %s\n' % self.pkg.name
        if self.base_dir:
            txt += '  Base directory: %s\n' % self.base_dir
        if self.hdr_dirs:
            txt += '  Header extensions: %s\n' % self.hdr_dirs
        if self.lib_dirs:
            txt += '  Library extensions: %s\n' % self.lib_dirs
        if self.libs:
            txt += '  Libraries: %s\n' % self.libs
        if self.fwork:
            txt += '  Framework: %s\n' % self.fwork
        if self.cpp_defines or self.sym_def:
            cpp_def = list(self.cpp_defines)
            if self.sym_def:
                cpp_def += [self.sym_def]
            txt += '  Exporting: %s\n' % (cpp_def)

        return txt
