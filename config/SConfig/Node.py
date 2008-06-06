import os
import SConfig

class Node(object):
    def __init__(self, scons_env, scons_opts, required=False):
        # Override to give a different name.
        self.name = self.__module__.split('.')[-1]

        # Option variables.
        self.command_name = self.name.lower()
        self.environ_name = self.name.upper()
        self.option_map = {} # Maps command-line options to their environment alternatives.
        self.options_processed = False
        self.options_result = '' # String to store the results of options processing for later
                                 # display.

        # Override with a list of methods to run during configuration.
        self.checks = []

        # Will be set after configuration.
        self.configured = False
        self.result = False

        # Private stuff.
        self.env = scons_env
        self.opts = scons_opts
        self.required = required
        self.deps = []

        # Setup our option database.
        self.setup_options()
        self.opts.Update(self.env)

    def setup_options(self):
        """Setup all the options for this package."""

        pass

    def dependency(self, package_module, required=True, **kw):
        """Add another package as a dependency of this package. If required is False, the
        dependent package is not required, and thus will not cause this package to fail if
        it cannot be found."""

        if self.configured:
            print 'Error: Cannot add a dependency during configuration.'
            self.env.Exit()
        pkg = self.env.Package(package_module, required, **kw)
        if pkg not in [d[0] for d in self.deps]:
            self.deps += [(pkg, required)]
        return pkg

    def setup(self):
        """Anything that needs to be finalised before continuing with the configuration needs
        to go here."""

        # Process options which will print out what options were found.
        if not self.options_processed:
            self.options_processed = True
            self.process_options()

    def configure(self, scons_ctx):
        """Perform the configuration of this package."""

        # Will need to color some stuff here.
        import TerminalController
        term = TerminalController.TerminalController()

        # Basic setup.
        if self.configured:
            return
        self.ctx = scons_ctx
        self.configured = True
        self.process_dependencies()

        # Print opening message.
        self.ctx.Message('Configuring package %s ... ' % self.name)
        self.ctx.Display('\n')

        # Display the result of options processing.
        self.ctx.Log(self.options_result)
        print term.render('${GREEN}%s${NORMAL}\r' % self.options_result),

        # Check we have all dependencies.
        result = True
        for pkg, req in self.deps:
            if req and not pkg.result:
                self.ctx.Log('  Missing dependency: %s\n' % pkg.name)
                print term.render('  ${RED}Missing dependency: %s${NORMAL}\n' % pkg.name),
                result = False

        # Perform as many checks as we can without failing.
        if result:
            for check in self.checks:
                result = check()
                if not result:
                    break

        # If everything succeeded, display configuration results.
        if result:
            self.display_configuration()

        # Handle results.
        self.result = result
        self.ctx.Display('  ')
        if not self.result:
            print term.render('${RED}'),
        self.ctx.Result(result)
        if not self.result:
            print term.render('${NORMAL}\r'),

        # If this was a critical fail, try and help the user.
        if self.required and not result:
            self.ctx.Display('\nThe required package ' + self.name + ' could not be found.\n')
            self.ctx.Display('The printouts above should provide some information on what went\n')
            self.ctx.Display('wrong. To see further details, please read the \'config.log\' file.\n')
            if len(self.option_map.keys()):
                self.ctx.Display('You can directly specify search parameters for this package via\n')
                self.ctx.Display('the following command line options:\n\n')
                for opt in self.option_map.iterkeys():
                    self.ctx.Display('  ' + opt + '\n')
                self.ctx.Display('\nRun \'scons help\' for more details on these options.\n\n')
            self.env.Exit()

    def enable(self, scons_env, old_state=None):
        """Modify the SCons environment to have this package enabled. Begin by inserting
        all options on this node into the environment."""

        for pkg, req in self.deps: # Enable dependencies first.
            if pkg.result:
                pkg.enable(scons_env, old_state)
        for opt in self.option_map.iterkeys(): # Now add options.
            if opt in self.env._dict:
                scons_env[opt] = self.env[opt]

    def backup_variable(self, scons_env, var_name, old_state):
        if old_state is None:
            return
        if var_name not in old_state:
            if var_name in scons_env._dict:
                old_state[var_name] = scons_env[var_name]
            else:
                old_state[var_name] = None

    def restore_state(self, scons_env, old_state):
        for var_name, state in old_state.iteritems():
            if state is None:
                del scons_env[var_name]
            else:
                scons_env[var_name] = state

    def process_options(self):
        """Do any initial option processing, including importing any values from
        the environment and validating that all options are consistent."""

        # Search command line options.
        cmd_opts = False
        for opt in self.option_map.iterkeys():
            if opt in self.opts.args:
                if not cmd_opts:
                    self.options_result += '  Found command line options:\n'
                    cmd_opts = True
                self.options_result += '    %s = %s\n' % (opt, self.opts.args[opt])
                break

        # We don't want to mix command line and evironment options.
        if cmd_opts:
            return

        # Now go for environment options.
        env_opts = False
        for cmd, env in self.option_map.iteritems():
            if cmd not in self.opts.args and env in self.env['ENV']:
                if not env_opts:
                    self.options_result += '  Found environment options:\n'
                    env_opts = True
                self.env[cmd] = self.env['ENV'][env]
                self.options_result += '    %s = %s\n' % (env, self.env[cmd])

    def process_dependencies(self):
        """Ensure all dependencies have been configured before this package."""

#        old_state = {}
        for pkg, req in self.deps:
            pkg.configure(self.ctx)
#            if pkg.result:
#                pkg.enable(self.env, old_state)
#        return old_state

    def compile_source(self, source):
        """At this point we know all our construction environment has been set up,
        so we should be able to compile some source code."""

        result = self.run_scons_cmd(self.ctx.TryCompile, source, '.c')
        return [result[0], result[1]]

    def link_source(self, source):
        """At this point we know all our construction environment has been set up,
        so we should be able to build and run the application."""

        result = self.run_scons_cmd(self.ctx.TryLink, source, '.c')
        return [result[0], result[1]]

    def library_source(self, source):
        """Build a library out of some source code."""

        result = self.run_scons_cmd(self.ctx.TryBuild, self.env.SharedLibrary,
                                    source, '.c')
        return [result[0], result[1]]

    def run_source(self, source):
        """At this point we know all our construction environment has been set up,
        so we should be able to build and run the application."""

        result = self.run_scons_cmd(self.ctx.TryRun, source, '.c')
        return [result[0][0], result[0][1], result[1]]

    def run_scons_cmd(self, cmd, *args, **kw):
        old_log = self.ctx.sconf.logstream
        self.ctx.sconf.logstream = open('sconfig.log', 'w') # Capture the log.
        res = cmd(*args, **kw) # Execute the command.
        try:
            self.ctx.sconf.logstream.close() # Make sure the file is closed.
        finally:
            pass
        self.ctx.sconf.logstream = old_log # Replace the old log.

        # Return results.
        log_file = open('sconfig.log', 'r')
        log = log_file.read()
        log_file.close()
        os.remove('sconfig.log')
        old_log.write(log)
        return [res, log]

    def display_configuration(self):
        """Print out a brief summary of what options we found/used."""
        pass

    def push_state(self, state, append=False):
        old = {}
        copy = dict(state)
        for k, v in copy.iteritems():
            if not v:
                continue
            if not isinstance(v, list):
                copy[k] = [v]
            else:
                copy[k] = v
            old[k] = self.env.get(k, [])
        if append:
            self.env.AppendUnique(**copy)
        else:
            self.env.PrependUnique(**copy)
        return old

    def pop_state(self, old):
        self.env.Replace(**old)
