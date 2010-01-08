import os
Import('env')

#
# Need to make a copy because SCons uses the environment.
# at it's final state, so StGermain ends up depending on
# StgDomain, etc.
#

env = env.Clone()

#
#  Build the 'pcu' sub-project.
#

SConscript('pcu/SConscript', exports='env')
env.Prepend(LIBS=['pcu'])

#
# Inside each project we will be accessing headers without the
# project name as a build_dir, so we need to let SCons know how to
# find those headers.
#

env.Append(CPPPATH=env['build_dir'] + '/include/StGermain')

#
# As a bit of a special case we need to copy these headers across
# first up. Can't use the install builder either, have to actually
# copy the bastards.
#

dir = Dir(env['build_dir'] + '/include/StGermain/Base/IO/mpirecord').abspath
if not os.path.exists(os.path.join(dir, 'mpimessaging.h')):
    if not os.path.exists(dir):
        os.makedirs(dir)
    h = Execute(Copy(File(env['build_dir'] + '/include/StGermain/Base/IO/mpirecord/mpimessaging.h'),
                     File('Base/IO/src/mpirecord/none/mpimessaging.h')))
dir = Dir(env['build_dir'] + '/include/StGermain/Base/Foundation').abspath
if not os.path.exists(os.path.join(dir, 'ClassSetup.h')):
    if not os.path.exists(dir):
        os.makedirs(dir)
    Execute(Copy(File(env['build_dir'] + '/include/StGermain/Base/Foundation/ClassSetup.h'),
                 File('Base/Foundation/src/ClassSetup.h')))
    Execute(Copy(File(env['build_dir'] + '/include/StGermain/Base/Foundation/ClassEmpty.h'),
                 File('Base/Foundation/src/ClassEmpty.h')))

#
# Keep a list of all the objects we build so we can make a library
# afterwards.
#

objs = []
suites = []
tst_exp = []
tst_input = []

#
# Process each directory uniformly.
#

dirs = Split('Base/Foundation Base/IO Base/Container Base/Automation Base/Extensibility ' \
                 'Base/Context Base/Charon/CharonSchema Base/Charon/CharonFactorySchema Base Utils libStGermain')
for d in dirs:

    # Need the module name, which is just the directory.
    mod_name = env['ESCAPE']('"' + ''.join(d.split('/')) + '"')
    cpp_defs = [('CURR_MODULE_NAME', mod_name)] + env.get('CPPDEFINES', [])

    # Setup where to look for files.
    src_dir = d + '/src'
    inc_dir = 'include/StGermain/' + d
    tst_dir = d + '/tests'
    tst_exp_dir = tst_dir + '/expected'
    tst_input_dir = tst_dir + '/input'
    tst_install_dir = 'tests/StGermain/' + d

    # Install the headers and '.def' files.
    hdrs = env.Install(inc_dir, Glob(src_dir + '/*.h'))
    defs = env.Install(inc_dir, Glob(src_dir + '/*.def'))

    # STUPID GLOBUS! Because of stupid globus we have to do this:
    cpppath = list(env['CPPPATH'])
    if 'Charon' in d:
        new_cpppath = []
        for i in range(len(cpppath)):
            if 'libxml' not in cpppath[i]:
                new_cpppath.append(cpppath[i])
        cpppath = new_cpppath

    # Build our source files.
    srcs = Glob(src_dir + '/*.c')
    srcs = [s for s in srcs if s.path.find('-meta.c') == -1]
    objs += env.SharedObject(srcs, CPPDEFINES=cpp_defs,
                             CPPPATH=cpppath)

    # Build any meta files.
    objs += env.stgSharedMeta(Glob(src_dir + '/*.meta'), CPPDEFINES=cpp_defs)

    # If we found any '.def' files make sure to register them as
    # explicit dependencies.
    if defs:
        env.Depends(hdrs + objs, defs)

    # Build any test suites we might find.
    suites += env.Object(Glob(tst_dir + '/*Suite.c'))

    # Install any test expected and input files
    tst_exp += env.Install(tst_install_dir + '/expected', Glob(tst_exp_dir + '/*'))
    tst_input += env.Install(tst_install_dir + '/input', Glob(tst_input_dir + '/*'))

# Need to install headers from libStGermain.
hdrs = env.Install('include/StGermain', Glob('libStGermain/src/*.h'))

#
# Build shared library.
#

if env['shared_libs']:
    env.SharedLibrary('lib/StGermain', objs)

#
# Build static library.
#

if env['static_libs']:

    # In order to allow other projects to register their modules,
    # we copy a template file to a standard location and direct
    # SCons to build it in to the final static library. Before the
    # build actually takes place the other projects will fill in
    # there details.
    dir = Dir(env['build_dir'] + '/StGermain').abspath
    if not os.path.exists(dir):
        os.makedirs(dir)
    f = File(env['build_dir'] + '/StGermain/stg_static_modules.c')
    Execute(Copy(f.abspath, File('libStGermain/src/stg_static_modules.c.tmpl')))

    # Now build the library.
    l = env.StaticLibrary(env['build_dir'] + '/lib/StGermain', objs)
    env.Install(env['prefix'] + '/lib', l)

#
# FlattenXML, StGermain and test runner programs.
#

libs = ['StGermain'] + env.get('LIBS', [])
env.Program('bin/FlattenXML', 'Base/FlattenXML/src/main.c', LIBS=libs)
if env['shared_libs']:
    env.Program('bin/StGermain', 'src/main.c', LIBS=libs)
env.PCUTest('tests/testStGermain', suites,
    PCU_LIBHEADERS="#include <StGermain/StGermain.h>",
    PCU_SETUP="StGermain_Init(&argc, &argv);",
    PCU_TEARDOWN="StGermain_Finalise();",
    LIBS=libs,
    PCU_EXP=tst_exp,
    PCU_INPUT=tst_input,
    PROJECT="StGermain")

#
# Copy XML validation file to correct destination.
#

xmls = env.Install('lib', 'Base/IO/src/StGermain.xsd')
