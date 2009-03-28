import os

def build_suite_runner(target, source, env):
    hdrs = ""
    suites = ""
    init = env.get('PCURUNNERINIT', '_init')
    setup = env.get('PCURUNNERSETUP', '')
    if setup: setup = '\n   ' + setup
    teardown = env.get('PCURUNNERTEARDOWN', '')
    if teardown: teardown = '\n   ' + teardown
    for node in source:
        name = os.path.basename(node.path[:node.path.rfind('.')])
        suites += "   pcu_runner_addSuite( %s, %s );\n" % (name, name + init)
        name = node.abspath[:node.abspath.rfind('.')]
        hdrs += "#include \"%s.h\"\n" % name
    src = """#include <stdlib.h>
#include <mpi.h>
#include <pcu/pcu.h>
%s

int main( int argc, char* argv[] ) {
   pcu_listener_t* lsnr;

   MPI_Init( &argc, &argv );
   pcu_runner_init( argc, argv );%s

%s
   lsnr = pcu_textoutput_create();
   pcu_runner_run( lsnr );
   pcu_textoutput_destroy( lsnr );
%s
   pcu_runner_finalise();
   MPI_Finalize();
   return EXIT_SUCCESS;
}
""" % (hdrs, setup, suites, teardown)
    if not os.path.exists(os.path.dirname(target[0].abspath)):
        os.makedirs(os.path.dirname(target[0].abspath))
    f = open(target[0].abspath, "w")
    f.write(src)
    f.close()

Import('env')
env['BUILDERS']['PCUSuiteRunner']=Builder(action=build_suite_runner)
