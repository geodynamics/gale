import os
Import('env')

#
# Need to make a copy because SCons uses the environment
# at it's final state, so StGermain ends up depending on
# StgDomain, etc.
#

env = env.Clone()

#
# Inside each project we will be accessing headers without the
# project name as a prefix, so we need to let SCons know how to
# find those headers.
#

env.Append(CPPPATH=env['build_dir'] + '/include/glucifer')

#
# Need to include the gLucifer library for binaries.
#

libs = ['glucifer'] + env.get('LIBS', [])

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
    tst_exp_dir = tst_dir + '/expected'
    tst_input_dir = tst_dir + '/input'
    tst_install_dir = 'tests/glucifer/' + d

    # Install the headers and '.def' files.
    hdrs = env.Install(inc_dir, Glob(src_dir + '/*.h'))
    defs = env.Install(inc_dir, Glob(src_dir + '/*.def'))

    # Build our source files.
    srcs = Glob(src_dir + '/*.cxx')
    srcs = [s for s in srcs if s.path.find('-meta.cxx') == -1]
    objs += env.SharedObject(srcs, CPPDEFINES=cpp_defs)

    # Build any meta files.
    objs += env.stgSharedMeta(Glob(src_dir + '/*.meta'), CPPDEFINES=cpp_defs)

    # If we found any '.def' files make sure to register them as
    # explicit dependencies.
    if defs:
        env.Depends(hdrs + objs, defs)

    # Build any test suites we might find.
    suites += env.Object(Glob(tst_dir + '/*Suite.cxx'))

    # Install any test expected and input files
    tst_exp += env.Install(tst_install_dir + '/expected', Glob(tst_exp_dir + '/*'))
    tst_input += env.Install(tst_install_dir + '/input', Glob(tst_input_dir + '/*'))

# Need to install headers from libglucifer.
env.Install('include/glucifer', Glob('libglucifer/src/*.h'))

#
# Build plugins.
#

dirs = [('plugins/lucPlugin', 'lucPlugin')]
pl_objs = []
pl_regs = []
for d in dirs:

    if isinstance(d, tuple):
        name = d[1] + 'module'
        d = d[0]
    else:
        name = 'glucifer_' + d.split('/')[-1] + 'module'

    mod_name = env['ESCAPE']('"' + ''.join(d.split('/')) + '"')
    cpp_defs = [('CURR_MODULE_NAME', mod_name)] + env.get('CPPDEFINES', [])

    env.Install('include/glucifer/' + d.split('/')[-1], Glob(d + '/*.h'))

    srcs = Glob(d + '/*.cxx')
    srcs = [s for s in srcs if s.path.find('-meta.cxx') == -1]
    cur_objs = env.SharedObject(srcs, CPPDEFINES=cpp_defs)
    cur_objs += env.stgSharedMeta(Glob(d + '/*.meta'), CPPDEFINES=cpp_defs)

    # If we have shared libraries, build the dynamic plugin.
    if env['shared_libs']:
        lib_pre = env['LIBPREFIXES']
        if not isinstance(lib_pre, list):
            lib_pre = [lib_pre]
        env.SharedLibrary('lib/' + name, cur_objs,
                          SHLIBPREFIX='',
                          LIBPREFIXES=lib_pre + [''],
                          LIBS=libs)

    # If we are building static libs we need to construct a C file
    # mapping the plugin's name to its register function.
    if env['static_libs']:
        pl_regs += [name]

    # Keep track of all the plugin objects.
    pl_objs += cur_objs

#
# Build shared library.
#

if env['shared_libs']:
    env.SharedLibrary('lib/glucifer', objs)

#
# Write the static registry code.
#

if env['static_libs']:

    reg_c = '#include <StGermain/StGermain.h>\n\n'
    reg_c += 'extern int stg_num_modules;\n'
    reg_c += 'extern char **stg_module_names;\n'
    reg_c += 'extern int *stg_num_module_syms;\n'
    reg_c += 'extern char ***stg_module_syms;\n\n'
    reg_c += 'extern void ***stg_module_funcs;\n\n'
    for n in pl_regs:
        n = n[:-6]
        reg_c += 'extern void (%s_MetaAsDictionary)();\n'%n
        reg_c += 'extern void (%s_GetName)();\n'%n
        reg_c += 'extern void (%s_Register)();\n'%n
        if n.find('Toolbox') != -1:
            reg_c += 'extern void (%s_Initialise)();\n'%n
            reg_c += 'extern void (%s_Finalise)();\n'%n
    reg_c += '\n'

    reg_c += 'void glucifer_register_static_modules() {\n'
    reg_c += '   int new_num = stg_num_modules + %d;\n\n'%len(pl_regs)
    
    reg_c += '   stg_module_names = (char**)realloc(stg_module_names, new_num*sizeof(char*));\n'
    for i in range(len(pl_regs)):
        n = pl_regs[i][:-6]
        reg_c += '   stg_module_names[stg_num_modules + %d] = (char*)"%s";\n'%(i, n)
        
    reg_c += '\n   stg_num_module_syms = (int*)realloc(stg_num_module_syms, new_num*sizeof(int));\n'
    for i in range(len(pl_regs)):
        n = pl_regs[i][:-6]
        if n.find('Toolbox') != -1:
            num_syms = 5
        else:
            num_syms = 1
        reg_c += '   stg_num_module_syms[stg_num_modules + %d] = %d;\n'%(i, num_syms)

    reg_c += '\n   stg_module_syms = (char***)realloc(stg_module_syms, new_num*sizeof(char**));\n'
    for i in range(len(pl_regs)):
        n = pl_regs[i][:-6]
        if n.find('Toolbox') != -1:
            num_syms = 5
        else:
            num_syms = 1
        reg_c += '   stg_module_syms[stg_num_modules + %d] = (char**)malloc(%d*sizeof(char*));\n'%(i, num_syms)
        reg_c += '   stg_module_syms[stg_num_modules + %d][0] = (char*)"%s_Register";\n'%(i, n)
        if n.find('Toolbox') != -1:
            reg_c += '   stg_module_syms[stg_num_modules + %d][1] = (char*)"%s_MetaAsDictionary";\n'%(i, n)
            reg_c += '   stg_module_syms[stg_num_modules + %d][2] = (char*)"%s_GetName";\n'%(i, n)
            reg_c += '   stg_module_syms[stg_num_modules + %d][3] = (char*)"%s_Initialise";\n'%(i, n)
            reg_c += '   stg_module_syms[stg_num_modules + %d][4] = (char*)"%s_Finalise";\n'%(i, n)

    reg_c += '\n   stg_module_funcs = (void***)realloc(stg_module_funcs, new_num*sizeof(void**));\n'
    for i in range(len(pl_regs)):
        n = pl_regs[i][:-6]
        if n.find('Toolbox') != -1:
            num_syms = 5
        else:
            num_syms = 1
        reg_c += '   stg_module_funcs[stg_num_modules + %d] = (void**)malloc(%d*sizeof(void*));\n'%(i, num_syms)
        reg_c += '   stg_module_funcs[stg_num_modules + %d][0] = (void*)%s_Register;\n'%(i, n)
        if n.find('Toolbox') != -1:
            reg_c += '   stg_module_funcs[stg_num_modules + %d][1] = (void*)%s_MetaAsDictionary;\n'%(i, n)
            reg_c += '   stg_module_funcs[stg_num_modules + %d][2] = (void*)%s_GetName;\n'%(i, n)
            reg_c += '   stg_module_funcs[stg_num_modules + %d][3] = (void*)%s_Initialise;\n'%(i, n)
            reg_c += '   stg_module_funcs[stg_num_modules + %d][4] = (void*)%s_Finalise;\n'%(i, n)

    reg_c += '\n   stg_num_modules += %d;\n'%len(pl_regs)
    reg_c += '}\n'

    reg_filename = os.path.join(env['build_dir'], 'gLucifer', 'glucifer_static_modules.cxx')
    if not os.path.exists(os.path.dirname(reg_filename)):
        os.makedirs(os.path.dirname(reg_filename))
    reg_file = open(reg_filename, 'w')
    reg_file.write(reg_c)
    reg_file.close()
    reg_obj = env.Object(reg_filename)

    # Add our register function to the StGermain module file.
    f = open(File(env['build_dir'] + '/StGermain/stg_static_modules.cxx').abspath, 'r')
    txt = f.readlines()
    f.close()
    txt.insert(-2, '   glucifer_register_static_modules();\n')
    txt.insert(0, 'void glucifer_register_static_modules();\n')
    f = open(File(env['build_dir'] + '/StGermain/stg_static_modules.cxx').abspath, 'w')
    f.writelines(txt)
    f.close()

    # Static library.
    l = env.StaticLibrary(env['build_dir'] + '/lib/glucifer', objs + pl_objs + reg_obj)
    env.Install(env['prefix'] + '/lib', l)

#
# Test runner program.
#

env.PCUTest('tests/testglucifer', suites,
            PCU_LIBHEADERS="#include <StGermain/StGermain.h>\n#include <StgDomain/StgDomain.h>\n"
                "#include <StgFEM/StgFEM.h>\n#include <glucifer/glucifer.h>",
            PCU_SETUP="StGermain_Init(&argc, &argv);StgDomain_Init(&argc, &argv);" \
                "StgFEM_Init(&argc, &argv);glucifer_Init(&argc, &argv);",
            PCU_TEARDOWN="glucifer_Finalise();StgFEM_Finalise();" \
                "StgDomain_Finalise();StGermain_Finalise();",
            LIBS=libs,
            PCU_EXP=tst_exp,
            PCU_INPUT=tst_input,
            PROJECT="glucifer")

#
# Install XML input files.
#

env.Install('lib/StGermain/glucifer', Glob('ModelComponents/*.xml'))
