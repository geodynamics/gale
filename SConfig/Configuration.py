class Configuration:
    def __init__(self, inst):

        # The installation this configuration uses.
        self.inst = inst

        # A list of configurations that are this configuration's
        # dependencies.
        self.deps = []

        # The set of headers/libraries that were found to work for this
        # configuration, or 'None' if not available.
        self.hdrs = None
        self.libs = None

        # Flags whether this installation/configuration combo has
        # static/shared libraries available.
        self.has_static = None
        self.has_shared = None

        # Preprocessor definitions.
        self.cpp_defines = []

        # This list provides a way of modifying the printed output of
        # a configuration without needing to replace the print method.
        self.outputs = [('Static libraries: %s', 'has_static'),
                        ('Shared libraries: %s', 'has_shared')]

    def __eq__(self, cfg):
        if not (self.inst == cfg.inst):
            return False
        if len(self.deps) != len(cfg.deps):
            return False
        for d1, d2 in zip(self.deps, cfg.deps):
            if not (d1 == d2):
                return False
        return True

    def add_dependency(self, cfg):
        """Add a configuration to the list of dependencies."""

        if cfg not in self.deps:
            self.deps += [cfg]

    def enable(self, scons_env, old_state=None, lib_exclude=[]):
        """Inserts this configuration's information into the provided SCons
        environment."""

        # First we have to enable all of our dependencies.
        self.enable_dependencies(scons_env, old_state, lib_exclude=lib_exclude)

        # Then call the installation's enable routine.
        self.inst.enable(scons_env, old_state,
                         libs=self.libs,
                         has_shared=self.has_shared,
                         lib_exclude=lib_exclude)

        # Throw in any preprocessor stuff.
        if self.cpp_defines:
            self.inst.pkg.backup_variable(scons_env, 'CPPDEFINES', old_state)
            scons_env.AppendUnique(CPPDEFINES=self.cpp_defines)

    def enable_dependencies(self, scons_env, old_state={}, lib_exclude=[]):
        """Enables all available dependencies."""

        for d in self.deps:
            d.enable(scons_env, old_state, lib_exclude=lib_exclude)

    def flatten_dependencies(self, deps=[], pkgs={}):
        """Return a list of configurations representing all dependencies of this
        configuration and it's child configurations. Uniqueness is determined by
        the package each configuration belongs to."""

        for d in self.deps:
            pkg = d.inst.pkg

            # If the package alrady exists in the dictionary of packages we've
            # already got a configuration for, then the configurations must be
            # matching, or this is not a valid set.
            if pkg in pkgs and pkgs[pkg] != d:
                return None

            deps += [d]
            pkgs[pkg] = d

            # Have to iterate over every child, not just children that are
            # unique. This is because configurations for the same package can,
            # and will, be different.
            result = d.flatten_dependencies(deps, pkgs)
            if result is None:
                return None
            deps += result

    def value(self):
        """Evaluates a 'value' figure for this configuration. The value
        corresponds to how it should be weighted when considering which
        configuration to use. Includes the values of it's children."""

        val = len(self.deps)
        for d in self.deps:
            val += d.value()
        return val

    def __str__(self, brief=True):
        """Convert to printable string."""

        # Get the installation's string first.
        txt = str(self.inst)

        # Print out our outputs.
        for out in self.outputs:
            line = out[0]
            vals = tuple([getattr(self, v) for v in out[1:]])
            txt += ('  ' + line + '\n') % vals

        # Now produce the dependency text.
        dep_txt = ''
        for dep in self.deps:

            # If we only want brief dependency information only print the
            # dependency's name.
            if brief:
                dep_txt += '    %s\n' % dep.inst.pkg.name
            else:
                dep_txt += '    ' + str(dep)[:-1].replace('\n', '\n    ') + '\n'

        # If there were any dependencies, add them to the installation's
        # text.
        if dep_txt:
            txt += '  Dependencies:\n'
            txt += dep_txt

        return txt
