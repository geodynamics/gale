import os
from SCons.Script.SConscript import SConsEnvironment
import pprint

import stgDtd
import stgMetaXsd
import convert

#
# Helpful build utilities.
#

def check_dir_target(env, src):
    src = env.project_name + '/' + src
    if os.path.commonprefix([env['dir_target'],
                             src]) == env['dir_target']:
        return True
    return False

def build_files(env, files, dst_dir):
    for f in files:
        if not env.check_dir_target(f): continue
        env.copy_file(env.get_build_path(dst_dir) + '/' + os.path.basename(f), f)

def build_headers(env, headers, dst_dir):
    for hdr in headers:
        if not env.check_dir_target(hdr): continue
        env.Install(env.get_build_path(dst_dir), hdr)

def build_sources(env, sources, dst_dir):
    objs = []
    mod_name = ('CURR_MODULE_NAME', env['ESCAPE']('"' + ''.join(dst_dir.split('/')) + '"'))
    for src in sources:
        if not env.check_dir_target(src): continue
        tgt = env.get_build_path(dst_dir + '/' + os.path.basename(src)[:-2])
        objs += env.SharedObject(tgt, src, CPPDEFINES=[mod_name] + env.get('CPPDEFINES', []))
    return objs

def build_metas(env, metas, dst_dir):
    objs = []
    mod_name = ('CURR_MODULE_NAME', env['ESCAPE']('"' + ''.join(dst_dir.split('/')) + '"'))
    for m in metas:
        if not env.check_dir_target(m): continue
        src = env.Meta(env.get_build_path(dst_dir + '/' + os.path.basename(m)[:-5]), m)
        tgt = env.get_build_path(dst_dir + '/' + os.path.basename(src[0].path)[:-2])
        if 'tau_old_cc' in env._dict:
            objs += env.SharedObject(tgt, src,
                                     CC=env['tau_old_cc'],
                                     CPPDEFINES=[mod_name] + env.get('CPPDEFINES', []))
        else:
            objs += env.SharedObject(tgt, src, CPPDEFINES=[mod_name] + env.get('CPPDEFINES', []))
    return objs

def build_directory(env, dir, dst_dir='', with_tests=True):
    if not env.check_dir_target(dir): return
    if not dst_dir:
        dst_dir = dir
    inc_dir = 'include/' + env.project_name + '/' + dst_dir
    obj_dir = env.project_name + '/' + dst_dir
    env.build_files(env.glob(dir + '/src/*.def'), inc_dir)
    env.build_headers(env.glob(dir + '/src/*.h'), inc_dir)
    env.src_objs += env.build_sources(env.glob(dir + '/src/*.c'), obj_dir)
    env.src_objs += env.build_metas(env.glob(dir + '/src/*.meta'), obj_dir)
    env.suite_hdrs += env.glob(dir + '/tests/*Suite.h')
    env.suite_objs += env.build_sources(env.glob(dir + '/tests/*Suite.c'), obj_dir)

def build_plugin(env, dir, dst_dir='', name='', with_lib=True):
    if not env.check_dir_target(dir): return
    if not dst_dir:
        dst_dir = dir
    if not name:
        name = env.project_name + '_' + dir.split('/')[-1]
    mod_name = name + 'module'
    env.build_headers(env.glob(dir + '/*.h'), 'include/' + env.project_name + '/' + dir.split('/')[-1])
    objs = env.build_sources(env.glob(dir + '/*.c'), env.project_name + '/' + dir)
    if env['shared_libraries']:
        if with_lib:
            libs = [env.project_name] + env.get('LIBS', [])
        else:
            libs = env.get('LIBS', [])
        env.SharedLibrary(env.get_build_path('lib/' + mod_name), objs,
                          SHLIBPREFIX='',
                          LIBPREFIXES=env.make_list(env['LIBPREFIXES']) + [''],
                          LIBS=libs)
    if env['static_libraries']:
        env.src_objs += objs
        if not env['shared_libraries']:
            if not env.get('STGMODULECODE', None):
                env['STGMODULECODE'] = ''
            if not env.get('STGMODULEPROTO', None):
                env['STGMODULEPROTO'] = ''
            env['STGMODULEPROTO'] += 'int ' + env.project_name + '_' + name + '_Register( void* );\n'
            env['STGMODULECODE'] += '  ' + env.project_name + '_' + name + '_Register( mgr );\n'

def build_toolbox(env, dir, dst_dir=''):
    if not env.check_dir_target(dir): return
    if not dst_dir:
        dst_dir = env.project_name + '/' + dir
    objs = env.build_sources(env.glob(dir + '/*.c'), dst_dir)
    objs += env.build_metas(env.glob(dir + '/*.meta'), dst_dir)
    if env['shared_libraries']:
        env.SharedLibrary(env.get_target_name('lib/' + env.project_name + '_Toolboxmodule'), objs,
                          SHLIBPREFIX='',
                          LIBPREFIXES=env.make_list(env['LIBPREFIXES']) + [''],
                          LIBS=[env.project_name] + env.get('LIBS', []))
    if env['static_libraries']:
        env.src_objs += objs
        if not env['shared_libraries']:
            if not env.get('STGMODULECODE', None):
                env['STGMODULECODE'] = ''
            if not env.get('STGMODULEPROTO', None):
                env['STGMODULEPROTO'] = ''
            env['STGMODULEPROTO'] += 'int ' + env.project_name + '_Toolbox_Register( void* );\n'
            env['STGMODULEPROTO'] += 'void ' + env.project_name + '_Toolbox_Initialise( void*, int, char** );\n'
            env['STGMODULECODE'] += '  ' + env.project_name + '_Toolbox_Register( mgr );\n'
            env['STGMODULECODE'] += '  ' + env.project_name + '_Toolbox_Initialise( mgr, argc, argv );\n'

def build_module_register(env, dst, module_proto, module_code):
    src = str(module_proto)
    src += '\nvoid SingleRegister( void* mgr, int argc, char** argv ) {\n'
    src += module_code
    src += '}\n'
    if not os.path.exists(os.path.dirname(File(dst).abspath)):
        os.makedirs(os.path.dirname(File(dst).abspath))
    f = open(File(dst).abspath, "w")
    f.write(src)
    f.close()
    return File(dst)

SConsEnvironment.check_dir_target = check_dir_target
SConsEnvironment.build_files = build_files
SConsEnvironment.build_headers = build_headers
SConsEnvironment.build_sources = build_sources
SConsEnvironment.build_metas = build_metas
SConsEnvironment.build_directory = build_directory
SConsEnvironment.build_plugin = build_plugin
SConsEnvironment.build_toolbox = build_toolbox
SConsEnvironment.build_module_register = build_module_register

#
# Custom builders.
#

def loadXML( filename ):
	xml_file = file( filename )
	xml_lines = xml_file.readlines()
	xml_text = ""
	for l in xml_lines:
	        xml_text += str(l)

	try:
		dtdDict = stgDtd.readXML( xml_text )
	except:
		print 'Failed to parse as a StGermain DTD'
		raise
	try:
		return convert.dtdDict2metaXsdDict( dtdDict )
	except:
		print 'Failed to convert information from a StGermain Meta DTD to a StGermain Meta XSD'
		raise


# Builder for generating meta files (courtesy of Walter Landry).
def create_meta(target, source, env):
	output_file = file( str( target[0] ), 'wb' )
	xsdDict = loadXML( str( source[0] ) )
	# TODO: Add here things like real location (rep & path), test results & sv/hg diff status

	output_file.write( convert.metaXsdDict2stgCodeHeader() )
	output_file.write( '\n' )
	output_file.write( convert.metaXsdDict2stgStrings( xsdDict ) )
	output_file.write( '\n' )
	output_file.write( convert.metaXsdDict2stgDictionaryCode( xsdDict ) )

	output_file.close()

def gen_meta_suffix(env, sources):
    return "-meta.c"

Import('env')
env['BUILDERS']['Meta']=Builder(action=create_meta,src_suffix="meta",
                                suffix=gen_meta_suffix,single_source=True)

#
# Add object storage locations.
env.src_objs = []
env.suite_hdrs = []
env.suite_objs = []
