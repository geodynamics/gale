import os
Import('env')

# Need to make a copy because SCons uses the environment
# at it's final state, so StGermain ends up depending on
# StgDomain, etc.
env = env.Clone()

# Build the 'pcu' sub-project.
SConscript('pcu/SConscript', exports='env')

# Inside each project we will be accessing headers without the
# project name as a build_dir, so we need to let SCons know how to
# find those headers.
env.Append(CPPPATH=env['build_dir'] + '/include/StGermain')

# As a bit of a special case we need to copy these headers across
# first up. Can't use the install builder either, have to actually
# copy the bastards.
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

# Keep a list of all the objects we build so we can make a library
# afterwards.
objs = []
suites = []

# Process each directory uniformly.
dirs = Split('Base/Foundation Base/IO Base/Container Base/Automation Base/Extensibility ' \
                 'Base/Context Base Utils libStGermain')
for d in dirs:

    # Need the module name, which is just the directory.
    mod_name = env['ESCAPE']('"' + ''.join(d.split('/')) + '"')
    cpp_defs = [('CURR_MODULE_NAME', mod_name)] + env.get('CPPDEFINES', [])

    # Setup where to look for files.
    src_dir = d + '/src'
    inc_dir = 'include/StGermain/' + d
    tst_dir = d + '/tests'

    # Install the headers and '.def' files.
    hdrs = env.Install(inc_dir, Glob(src_dir + '/*.h'))
    defs = env.Install(inc_dir, Glob(src_dir + '/*.def'))

    # Build our source files.
    srcs = Glob(src_dir + '/*.c')
    srcs = [s for s in srcs if s.path.find('-meta.c') == -1]
    objs += env.SharedObject(srcs, CPPDEFINES=cpp_defs)

    # Build any meta files.
    objs += env.stgSharedMeta(Glob(src_dir + '/*.meta'), CPPDEFINES=cpp_defs)

    # If we found any '.def' files make sure to register them as
    # explicit dependencies.
    if defs:
        env.Depends(hdrs + objs, defs)

    # Build any test suites we might find.
    suites += env.Object(Glob(tst_dir + '/*Suite.c'))

# Need to install headers from libStGermain.
hdrs = env.Install('include/StGermain', Glob('libStGermain/src/*.h'))

# Build libraries.
if env['shared_libraries']:
    env.SharedLibrary('lib/StGermain', objs)

# FlattenXML, StGermain and test runner programs.
libs = ['StGermain'] + env.get('LIBS', [])
env.Program('bin/FlattenXML', 'Base/FlattenXML/src/main.c', LIBS=libs)
env.Program('bin/StGermain', 'src/main.c', LIBS=libs)
env.PCUTest('tests/testStGermain', suites,
            PCU_SETUP="StGermain_Init(&argc, &argv);",
            PCU_TEARDOWN="StGermain_Finalise();",
            LIBS=libs)

# Copy XML validation file to correct destination.
xmls = env.Install('lib', 'Base/IO/src/StGermain.xsd')
