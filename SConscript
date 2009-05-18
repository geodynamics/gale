import os
Import('env')

# Need to make a copy because SCons uses the environment
# at it's final state, so StGermain ends up depending on
# StgDomain, etc.
env = env.Clone()

# Inside each project we will be accessing headers without the
# project name as a prefix, so we need to let SCons know how to
# find those headers.
env.Append(CPPPATH=env['build_dir'] + '/include/glucifer')

# Keep a list of all the objects we build so we can make a library
# afterwards.
objs = []
suites = []

# Process each directory uniformly.
dirs = Split('Base Windowing RenderingEngines OutputFormats InputFormats ' \
                 'DrawingObjects WindowInteractions libglucifer')
for d in dirs:

    # Need the module name, which is just the directory.
    mod_name = env['ESCAPE']('"' + ''.join(d.split('/')) + '"')
    cpp_defs = [('CURR_MODULE_NAME', mod_name)] + env.get('CPPDEFINES', [])

    # Setup where to look for files.
    src_dir = d + '/src'
    inc_dir = 'include/glucifer/' + d
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

# Need to install headers from libglucifer.
env.Install('include/glucifer', Glob('libglucifer/src/*.h'))

# Build libraries.
if env['shared_libraries']:
    env.SharedLibrary('lib/glucifer', objs)

# Need to include the gLucifer library for binaries.
libs = ['glucifer'] + env.get('LIBS', [])

# Test runner program.
env.PCUTest('tests/testglucifer', suites,
            PCU_LIBHEADERS="#include <StGermain/StGermain.h>\n#include <StgDomain/StgDomain.h>\n"
                "#include <StgFEM/StgFEM.h>\n#include <glucifer/glucifer.h>",
            PCU_SETUP="StGermain_Init(&argc, &argv);StgDomain_Init(&argc, &argv);" \
                "StgFEM_Init(&argc, &argv);glucifer_Init(&argc, &argv);",
            PCU_TEARDOWN="glucifer_Finalise();StgFEM_Finalise();" \
                "StgDomain_Finalise();StGermain_Finalise();",
            LIBS=libs)

# Build plugins.
dirs = [('plugins/lucPlugin', 'lucPlugin')]
for d in dirs:

    if isinstance(d, tuple):
        name = d[1] + 'module'
        d = d[0]
    else:
        name = 'glucifer_' + d.split('/')[-1] + 'module'

    mod_name = env['ESCAPE']('"' + ''.join(d.split('/')) + '"')
    cpp_defs = [('CURR_MODULE_NAME', mod_name)] + env.get('CPPDEFINES', [])

    env.Install('include/glucifer/' + d.split('/')[-1],
                Glob(d + '/*.h'))

    srcs = Glob(d + '/*.c')
    srcs = [s for s in srcs if s.path.find('-meta.c') == -1]
    objs = env.SharedObject(srcs, CPPDEFINES=cpp_defs)
    objs += env.stgSharedMeta(Glob(d + '/*.meta'), CPPDEFINES=cpp_defs)

    if env['shared_libraries']:
        lib_pre = env['LIBPREFIXES']
        if not isinstance(lib_pre, list):
            lib_pre = [lib_pre]
        env.SharedLibrary('lib/' + name, objs,
                          SHLIBPREFIX='',
                          LIBPREFIXES=lib_pre + [''],
                          LIBS=libs)

# Install XML input files.
env.Install('lib/StGermain/glucifer', Glob('ModelComponents/*.xml'))
