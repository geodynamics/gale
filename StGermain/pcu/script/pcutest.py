from SCons.Script import *
#from SCons.Builder import Builder

class ToolPCUTestWarning(SCons.Warnings.Warning):
    pass

SCons.Warnings.enableWarningClass(ToolPCUTestWarning)

def multiget(dicts, key, default=None):
    for d in dicts:
        if d.has_key(key):
            return d[key]
    else:
        return default

def build_suite_runner(env, target, hdrs, objs, **kw):
    hdr_txt = ""
    suite_txt = ""
    init = multiget([kw, env], "PCU_INIT", "")

    bld_tests_dir = os.path.join( env['build_dir'], "tests" )

    libheaders = multiget([kw, env], "PCU_LIBHEADERS", "")

    project_name = multiget([kw, env], "PROJECT", "")

    setup = multiget([kw, env], "PCU_SETUP", "")
    if setup:
        setup = '\n   ' + setup

    teardown = multiget([kw, env], "PCU_TEARDOWN", "")
    if teardown:
        teardown = '\n   ' + teardown

    for h in hdrs:
        name = os.path.splitext(os.path.basename(h.path))[0]
        moduleDir = os.path.split( os.path.dirname( h.path ) )[0]
        suite_txt += "   pcu_runner_addSuite( %s, %s, %s );\n"%(name, name + init, moduleDir )
        hdr_txt += "#include \"%s\"\n"%str(h.abspath)

    src = """#include <stdlib.h>
#include <mpi.h>
#include <pcu/pcu.h>
#include <unistd.h>
%s
%s

int main( int argc, char* argv[] ) {
   pcu_listener_t*   lsnr;
   PCU_Runner_Status result;

   chdir( "%s" );

   MPI_Init( &argc, &argv );
   pcu_runner_init( argc, argv );%s

%s
   lsnr = pcu_textoutput_create( PCU_PRINT_DOCS );
   result = pcu_runner_run( lsnr );
   pcu_textoutput_destroy( lsnr );
%s
   pcu_runner_finalise();
   MPI_Finalize();
   if ( result == PCU_RUNNER_ALLPASS ) {
      return EXIT_SUCCESS;
   }
   else {
      return EXIT_FAILURE;
   }
}
"""%(libheaders, hdr_txt, bld_tests_dir, setup, suite_txt, teardown)

    dir_path = os.path.dirname(target.abspath)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    f = open(target.abspath, "w")
    f.write(src)
    f.close()

    return File(target.abspath)

def generate(env, **kw):
    env.SetDefault(PCUALL_TARGET="check-complete")
    env.SetDefault(REGTEST_TARGET="check")
    env.SetDefault(UNIT_TARGET="check-unit")
    env.SetDefault(INTEGRATION_TARGET="check-integration")
    env.SetDefault(CONVERGENCE_TARGET="check-convergence")
    env.SetDefault(LOWRES_TARGET="check-lowres")
    # Make the 'all' target point to everything.

    def PCUSuite(env, target, source, **kw):
        """Create an object/header pair out of a
        *Suite.cxx/*Suite.h pair. The target should just
        be the name of the suite. So, if target were
        "Happy", the sources would be "HappySuite.cxx" and
        "HappySuite.h\""""
        obj = env.Object(target[0], source[0])
        return [obj + [File(os.path.splitext(source[0].abspath)[0] + ".h")]]

    def PCUTest(env, target, source, **kw):
        # Generate a list of headers, one for each suite source.
        hdrs = []
        objs = []
        for s in source:
            objs.append(s)
            hdrs.append(File(os.path.splitext(s.srcnode().abspath)[0] + '.h'))

        # Generate the program source.
        prog_src = build_suite_runner(env, File(str(target[0]) + ".cxx"), hdrs, objs, **kw)

        # Build everything.
        exps = multiget([kw, env], 'PCU_EXP', [])
        inputs = multiget([kw, env], 'PCU_INPUT', [])
        objs = env.StaticObject(os.path.splitext(prog_src.abspath)[0], prog_src) + objs
        libs = multiget([kw, env], 'LIBS', []) + ["pcu"]
        test = env.Program(target[0], objs, LIBS=libs)
        runner = env.Action('-' + test[0].abspath)
        env.Alias(env["UNIT_TARGET"], [exps, inputs, test], runner)
        env.AlwaysBuild(env["UNIT_TARGET"])
        env.Alias(env['PCUALL_TARGET'], env['UNIT_TARGET'])
        env.Alias(env['REGTEST_TARGET'], env['UNIT_TARGET'])
        return test

    def IntegrationTest(env, target, source, **kw):
        script = File(source[0].split()[0]).srcnode().abspath
        script_dir = os.path.dirname(script)
        script = os.path.basename(script)
        args = source[0].split()[1:]

        runner = env.Action('-./' + script + ' ' + ' '.join(args), chdir=script_dir)
        env.Alias(env["INTEGRATION_TARGET"], [], runner)
        env.AlwaysBuild(env["INTEGRATION_TARGET"])
        env.Alias(env['PCUALL_TARGET'], env['INTEGRATION_TARGET'])
        return None


    def LowResTest(env, target, source, **kw):
        script = File(source[0].split()[0]).srcnode().abspath
        script_dir = os.path.dirname(script)
        script = os.path.basename(script)
        args = source[0].split()[1:]

        runner = env.Action('-./' + script + ' ' + ' '.join(args), chdir=script_dir)
        env.Alias(env["LOWRES_TARGET"], [], runner)
        env.AlwaysBuild(env["LOWRES_TARGET"])
        env.Alias(env['PCUALL_TARGET'], env['LOWRES_TARGET'])
        env.Alias(env['REGTEST_TARGET'], env['LOWRES_TARGET'])
        return None

    def ConvergenceTest(env, target, source, **kw):
        script = File(source[0].split()[0]).srcnode().abspath
        script_dir = os.path.dirname(script)
        script = os.path.basename(script)
        args = source[0].split()[1:]

        runner = env.Action('-./' + script + ' ' + ' '.join(args), chdir=script_dir)
        env.Alias(env["CONVERGENCE_TARGET"], [], runner)
        env.AlwaysBuild(env["CONVERGENCE_TARGET"])
        env.Alias(env['PCUALL_TARGET'], env['CONVERGENCE_TARGET'])
        return None

    env.Append(BUILDERS={"PCUSuite": PCUSuite, "PCUTest": PCUTest,
                         "LowResTest": LowResTest, "IntegrationTest": IntegrationTest, "ConvergenceTest": ConvergenceTest})

def exists(env):
    # Should probably have this search for the pcu
    # libraries/source or something.
    return True
