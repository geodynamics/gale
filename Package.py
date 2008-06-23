import sys, os, platform, glob
import SCons.Script
import SConfig
import checks

class Package(SConfig.Node):
    """An object to describe how to search for a package
    that has been installed on a system. There are a lot of
    options that can be modified to refine how the search
    proceeds and also how results are generated."""

    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Node.__init__(self, scons_env, scons_opts, required)

        # Try and gain access to a dynamic loader in order to check for
        # consistent shared libraries. Make sure we don't try this for
        # the dl package.
        if self.name != 'dl':
            self.dl = self.dependency(SConfig.packages.dl)

        # This will be set in the preprocessor.
        self.have_define         = kw.get('have_define', '')

        # Search options.
        self.base_dirs           = [] #['']
        self.base_patterns       = [] #['']
        self.sub_dirs            = [] #[[[''], ['']]]
        self.header_sub_dir      = [] #['']
        self.system_header_dirs  = []
        self.system_library_dirs = []

        # Which headers do we require?
        self.headers             = [] #[['']]

        # Library options.
        self.libraries           = [] #[['']]
        self.shared_libraries    = None # Only libraries listed here will be considered
                                        # when checking for shared libraries.
        self.extra_libraries     = [] # These libraries are not considered locatable.
        self.require_shared      = kw.get('require_shared', False)
        self.frameworks          = []
        self.symbols             = [([], '')] #[([''], '')]
        self.symbol_setup        = ''
        self.symbol_teardown     = ''
        self.symbol_prototypes   = [] #['']
        self.symbol_calls        = [] #['']
        self.init_code           = ''
        self.fina_code           = ''
        self.check_code          = ''

        # Framework options.
        self.frameworks          = [] #[['']]

        # Used to determine whether we've executed our default method to
        # populate the set of candidate installations.
        self.candidates_built = False

        # We need this to flag whether we've been given one or more candidate
        # installations from options.
        self.given_candidates = False

        # These is our set of candidate installations and configurations.
        self.candidates = []
        self.cand_cfgs = []

        # Lists of all valid installations and configurations. These are
        # populated during the configuration of the package.
        self.installations = []
        self.configurations = []

        # This is set to the configuration that was selected to be used.
        self.selected = None


        # Need this so we can get access to information about the
        # platform we're running on.
        self.platform = self.dependency(SConfig.Platform, True) 

        # We have one configuration check.
        self.checks = [self.check_candidates]

        # Setup search defaults for the platform we're on.
        self.setup_search_defaults()

    def setup_options(self):
        """Two things need to happen here. The first is to tell SCons about
        any options we're going to use. Do that by calling
        'self.opts.AddOptions'. The second is to add entries to
        'self.option_map', which maps a command-line option to it's
        envrionment equivalent."""

        SConfig.Node.setup_options(self)
        self.opts.AddOptions(
            SCons.Script.BoolOption('with_' + self.command_name,
                                    'Turn on/off package %s' % self.name, 1),
            SCons.Script.PathOption(self.command_name + '_dir',
                                    '%s installation path' % self.name,
                                    None, SCons.Script.PathOption.PathIsDir),
            SCons.Script.PathOption(self.command_name + '_inc_dir',
                                    '%s header installation path' % self.name,
                                    None, SCons.Script.PathOption.PathIsDir),
            SCons.Script.PathOption(self.command_name + '_lib_dir',
                                    '%s library installation path' % self.name,
                                    None, SCons.Script.PathOption.PathIsDir),
            (self.command_name + '_lib',
             '%s libraries' % self.name,
             None, None),
            (self.command_name + '_framework',
             '%s framework' % self.name,
             None, None))
        self.option_map = {'with_' + self.command_name: None,
                           self.command_name + '_dir': self.environ_name + '_DIR'}

    def setup_search_defaults(self):
        """Setup the usual search paths for packages depending on the kind of system
        we're on."""

        if self.platform.system in ['Darwin', '*ix']:
            self.base_dirs = ['/usr', '/usr/local']
            self.sub_dirs = [[['include'], ['lib']]]
            if self.platform.bits == 64:
                self.sub_dirs = [[['include'], ['lib64']],
                                 [['include'], [os.path.join('lib', '64')]]] + self.sub_dirs

        # Set Darwin specific defaults.
        if self.platform.system == 'Darwin':
            self.base_dirs += ['/sw']

        # Set Window specific defaults.
        if self.platform.system == 'Windows':
            pass # TODO

        # Combine these guys to build default system paths. We need these to ensure specific
        # include paths are used before generic ones.
        for base_dir in ['/usr', '/sw']:
            for hdr_dirs, lib_dirs in self.combine_sub_dirs(base_dir):
                hdr_dirs = [os.path.join(base_dir, h) for h in hdr_dirs]
                hdr_dirs = [h for h in hdr_dirs if h not in self.system_header_dirs]
                self.system_header_dirs += hdr_dirs
                lib_dirs = [os.path.join(base_dir, l) for l in lib_dirs]
                lib_dirs = [l for l in lib_dirs if l not in self.system_library_dirs]
                self.system_library_dirs += lib_dirs

    def setup(self):
        """Finalise everything before running the configuration checks."""

        # Run the parent setup method.
        SConfig.Node.setup(self)

        # If we havn't already done so, build the list of candidates.
        if not self.candidates_built:
            self.candidates_built = True
            self.setup_candidates()

        # Now we need to traverse our set of candidate installations, looking
        # for any that havn't been processed.
        for inst in self.candidates:
            if not inst.is_processed:
                self.process_installation(inst)

    def setup_candidates(self):
        """Using the system specific search information currently set on
        the package class, build a list of candidate installations in
        the 'candidates' member."""

        # If we've been given options directly specifying the location of this
        # package we need to use those in place of searching for locations.
        base_dir = self.env.get(self.command_name + '_dir', '')
        inc_dir = self.env.get(self.command_name + '_inc_dir', '')
        lib_dir = self.env.get(self.command_name + '_lib_dir', '')
        if inc_dir or lib_dir:
            if not (inc_dir and lib_dir):
                print '   Error: must specify both of'
                print '      ' + self.command_name + '_inc_dir'
                print '      ' + self.command_name + '_lib_dir'
                env.Exit()
            self.add_candidate(SConfig.Installation(self, '', [inc_dir], [lib_dir]))
            self.given_candidates = True
            return

        # Combine all possible base directories.
        if not base_dir:
            base_dirs = list(self.base_dirs)
            for dir in self.base_dirs:
                for ptrn in self.base_patterns:
                    base_dirs += glob.glob(os.path.join(dir, ptrn))
        else:
            base_dirs = [base_dir]

        # Make sure there are no symbolic links in the base directories. If there are
        # any in there, expand them and make sure there are no duplicates.
        new_base_dirs = []
        for d in base_dirs:
            d = os.path.realpath(d)
            if d not in new_base_dirs:
                new_base_dirs += [d]
        base_dirs = new_base_dirs

        # Traverse the list of base directories and form each
        # installation.
        for dir in base_dirs:
            for hdr, lib in self.combine_base_dir(dir):
                self.add_candidate(SConfig.Installation(self, dir, hdr, lib))

        # If we have any frameworks to try, create candidates for them now.
        for fw in self.frameworks:
            self.add_candidate(SConfig.Installation(self, fwork=fw))

        # If we were given a base directory we need to set the 'given_candidates' flag
        # so that other candidate entries are ignored.
        if base_dir:
            self.given_candidates = True

    def combine_base_dir(self, base_dir):
        """Yields combinations of the provided base directory and possible sub-
        directories in the form of a [header sub-directory, library sub-
        directory] list. We use a list so that it can be modified in place."""

        # If the path doesn't exist or isn't a directory, don't yield anything.
        if not (os.path.exists(base_dir) and os.path.isdir(base_dir)):
            return

        # Combine the sub-directories.
        for hdr, lib in self.combine_sub_dirs(base_dir):

            # Yield the immediate results.
            yield [hdr, lib]

            # Also try all combinations of header sub-directories, if they
            # were given.
            for sub in self.combine_header_sub_dir(base_dir, hdr):
                yield [sub, lib]

    def combine_sub_dirs(self, base_dir):
        """Take a base directory and combine it with the set of header and library
        subdirectories. Yields (['header_dirs'], ['lib_dirs'])."""

        for hdr, lib in self.sub_dirs:
            loc_okay = True
            hdr_dirs = []
            lib_dirs = []

            # Combine header subdirectories.
            for h in hdr:
                dir = os.path.join(base_dir, h)
                if not (os.path.exists(dir) and os.path.isdir(dir)):
                    loc_okay = False
                    break
                hdr_dirs += [h]
            if not loc_okay:
                continue

            # Combine library subdirectories.
            for l in lib:
                dir = os.path.join(base_dir, l)
                if not (os.path.exists(dir) and os.path.isdir(dir)):
                    loc_okay = False
                    break
                lib_dirs += [l]
            if not loc_okay:
                continue

            yield (hdr_dirs, lib_dirs)

    def combine_header_sub_dir(self, base_dir, hdr_dirs):
        if not self.header_sub_dir or not hdr_dirs:
            return
        for sub_dir in self.header_sub_dir:
            cand = [os.path.join(h, sub_dir) for h in hdr_dirs if h]
            for d in cand:
                path = os.path.join(base_dir, d)
                if not (os.path.exists(path) and os.path.isdir(path)):
                    return
            yield cand

    def add_candidate(self, inst):
        """Add a unique candidate installation. Note that if candidate(s) were
        given via options, this method will ignore additional candidates."""

        if not self.given_candidates and inst not in self.candidates:
            self.candidates += [inst]
            self.is_setup = False

    def process_installation(self, inst):
        """This method gives us a chance to modify any of the details of this
        installation before moving on to checking it's validity. Here we can
        also determine which, if any, other dependant installations are given
        by this installation's package config (or any other means)."""

        inst.is_processed = True

    def check_candidates(self):
        """Runs a sequence of tests to confirm the validity of candidate
        installations. All that pass are moved to the 'installations'
        member."""

        # Combine all of the available candidate installations with the
        # available dependencies to build a set of candidate configurations.
        self.setup_configurations()

        # If we have no candidates, report the problem.
        if not self.candidates:
            self.ctx.Display('  No candidate installations found.\n')
            return False
        if not self.cand_cfgs:
            self.ctx.Display('  No candidate configurations found.\n')
            return False

        # Try out all candidates.
        cur = 1 # Keep track of which one we're currently trying.
        for cfg in self.cand_cfgs:
            # Print current status.
            self.ctx.Log('  Trialing candidate %d of %d ...\n' % (cur, len(self.cand_cfgs)))
            self.ctx.Log('    ' + cfg.__str__(False)[:-1].replace('\n', '\n    ') + '\n')
            print '  Trialing candidate %d of %d ...\r' % (cur, len(self.cand_cfgs)),
            cur = cur + 1

            # Check for the headers and libraries.
            if self.check_headers(cfg) and self.check_libraries(cfg):

                # If the checks passed, include the 'have_define' if there and add
                # this installation to the list.
                if self.have_define:
                    cfg.inst.add_cpp_defines(self.have_define)
                self.configurations += [cfg]

        # Print results.
        if len(self.configurations) == 1:
            self.ctx.Display('\n  Found %d valid configuration.\n' % len(self.configurations))
        else:
            self.ctx.Display('\n  Found %d valid configurations.\n' % len(self.configurations))

        # Log the valid configurations.
        for cfg in self.configurations:
            self.ctx.Log(cfg.__str__(False))

        # If we couldn't find a valid installation return negative.
        if len(self.configurations) == 0:
            return False
        return True

    def setup_configurations(self):
        """Copies each of the candidate installations and includes all
        permutations of dependent packages."""

        for deps in self.combine_dependencies():
            for inst in self.candidates:
                cfg = SConfig.Configuration(inst)
                cfg.deps = list(deps)
                self.cand_cfgs += [cfg]

    def combine_dependencies(self, deps=[], pkgs={}, cur_index=0):
        """Each combination of dependent installations represents a
        unique installation of this package. This method generates sets
        of unique dependency combinations."""

        # The dictionary 'pkgs' is to map packages to installations. This is
        # needed to prevent multiple installations being used for the same
        # packages.

        if cur_index == len(self.deps):
            yield deps # Complete permutation.
        else:
            cur_dep, required = self.deps[cur_index]

            # If the dependency isn't actually required by this package, include
            # a combination that doesn't use it.
            if not required:
                for d in self.combine_dependencies(deps, pkgs, cur_index + 1):
                    yield d

            # We can only iterate over installations if we're dealing with a package.
            if isinstance(cur_dep, Package):

                # Check if we already have this package selected.
                if cur_dep in pkgs:
                    deps += [pkgs[cur_dep]]
                    for d in self.combine_dependencies(deps, pkgs, cur_index + 1):
                        yield d
                    del deps[-1]

                elif len(cur_dep.configurations):

                    # Try each installation.
                    for dep_cfg in cur_dep.configurations:

                        # Traverse the dependency and collect any sub-dependencies
                        # into a copied package list.
                        new_pkgs = dict(pkgs)
                        rem = [dep_cfg]
                        while len(rem):
                            cur = rem.pop()
                            if cur.inst.pkg not in new_pkgs:
                                new_pkgs[cur.inst.pkg] = cur
                            rem += cur.deps

                        # Set dependency and recurse.
                        deps += [dep_cfg]
                        for d in self.combine_dependencies(deps, new_pkgs, cur_index + 1):
                            yield d
                        del deps[-1]

                elif required:

                    # There are no installations for this dependency. If it's
                    # a required dependency then something has gone very wrong.
                    # Throw an error here.
                    # TODO
                    sys.exit()
            else:
                for d in self.combine_dependencies(deps, pkgs, cur_index + 1):
                    yield d

    def get_check_headers_fail_reason(self, fail_logs):
        return ''

    def get_check_symbols_fail_reason(self, fail_logs):
        return ''

    def configure(self, scons_ctx):
        # If this package is deslected just return now.
        if not self.required and not self.env['with_' + self.command_name]:
            return

        # Run the configuration.
        SConfig.Node.configure(self, scons_ctx)

    def check_headers(self, cfg):
        """Determine if the required headers are available with the current construction
        environment settings."""

        # If there are no headers to check, automatically pass.
        if not self.headers:
            return True

        # Try and find a functional set of headers.
        for cfg.hdrs in self.headers:

            # Add any installation specific headers.
            cfg.hdrs += [h for h in cfg.inst.hdrs if h not in cfg.hdrs]

            # Run the check.
            if checks.check_headers(cfg):
                return True

        # If we failed to find anything, set the headers entry to 'None'.
        cfg.hdrs = None
        return False

    def check_libraries(self, cfg):
        """Check if the currently selected location is a valid installation of the
        required package. At this stage we know that the paths given in the location
        actually exist and we need to confirm that the libraries in 'libs' exist."""

        # If there are no libraries or frameworks to check, automatically pass.
        if not self.libraries:
            return True

        # Try and find a functional set of libraries.
        for cfg.libs in self.libraries:

            # Add any installation specific libraries.
            cfg.libs += [l for l in cfg.inst.libs if l not in cfg.libs]

            # Run the check.
            if checks.check_libraries(cfg):
                return True

        # If we failed to find anything, set the libraries entry to 'None'.
        cfg.libs = None
        return False

    def enable(self, scons_env, old_state=None):
        SConfig.Node.enable(self, scons_env, old_state)
        if self.selected:
            self.selected.enable(scons_env, old_state)
