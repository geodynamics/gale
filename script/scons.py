import os, shutil
import glob as pyglob
from SCons.Script.SConscript import SConsEnvironment
import pprint

import stgDtd
import stgMetaXsd
import convert

#
## \file Helpful build utilities.
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

def write_pkgconfig(env, filename, name, desc='', version=0):
    """Write out a pkgconfig file."""

    # Make sure the directory structure exists.
    filename = File(filename).abspath
    dirs = os.path.dirname(filename)
    if not os.path.exists(dirs):
        os.makedirs(dirs)

    # Write the pkgconfig file.
    f = open(filename, 'w')
    build_path = env.get('build_dir', '')
    if build_path:
        f.write('prefix=%s\n' % build_path)
        f.write('exec_prefix=%s\n' % os.path.join(build_path, 'bin'))
        f.write('libdir=%s\n' % os.path.join(build_path, 'lib'))
        f.write('includedir=%s\n' % os.path.join(build_path, 'include'))
        f.write('\n')
    f.write('Name: %s\n' % name)
    f.write('Description: %s\n' % desc)
    f.write('Version: %s\n' % version)
    f.write('Requires:\n')

    # Unfortunately SCons leaves hashes in paths after calling the
    # subst command, so we'll need to expand these manually.
    old_state = {'LIBPATH': env['LIBPATH'], 'CPPPATH': env['CPPPATH']}
    env['LIBPATH'] = [Dir(p).abspath for p in env['LIBPATH']]
    env['CPPPATH'] = [Dir(p).abspath for p in env['CPPPATH']]
    f.write(env.subst('Libs: ${_LIBDIRFLAGS} ${_LIBFLAGS}') + '\n')
    f.write(env.subst('Cflags: ${_CPPINCFLAGS}') + '\n')
    env.Replace(**old_state)
    f.close()

def copy_file(env, dst, src):
    dst = File(dst).abspath
    if os.path.exists(dst):
        return
    dst_dir = os.path.dirname(dst)
    if not os.path.exists(dst_dir):
        os.makedirs(dst_dir)
    shutil.copy(src, dst)

def get_build_path(env, prefix):
    if os.path.isabs(env['build_dir']):
        bld_dir = env['build_dir']
    else:
        bld_dir = '#' + env['build_dir']
    if prefix:
        return os.path.join(bld_dir, prefix)
    else:
        return bld_dir

def get_target_name(env, source, extension=''):
    """Return the destination name for a source file with suffix 'suffix'. This
    is useful for building files into the correct build path. Returns the full
    path to the built source without extension."""
    if extension:
        src = source[:-len(extension)]
    else:
        src = source
    return env.get_build_path(src)

def glob(env, pattern):
    if not os.path.isabs(pattern):
        old = os.getcwd()
        os.chdir(Dir('.').srcnode().abspath)
        res = pyglob.glob(pattern)
        os.chdir(old)
    else:
        res = pyglob.glob(pattern)
    return res

def path_exists(env, path):
    if not os.path.isabs(path):
        old = os.getcwd()
        os.chdir(Dir('.').srcnode().abspath)
        res = os.path.exists(path)
        os.chdir(old)
    else:
        res = os.path.exists(path)
    return res

def strip_dir(env, path, subdir):
    offs = path.find(os.path.sep + subdir + os.path.sep)
    if offs != -1:
        return path[:offs] + path[offs + len(subdir) + 1:]
    offs = path.find(os.path.sep + subdir)
    if offs != -1:
        return path[:-(len(subdir) + 1)]
    return path

def make_list(self, var):
    """Convert anything into a list. Handles things that are already lists,
    tuples and strings."""

    if isinstance(var, str):
        return [var]
    elif isinstance(var, (list, tuple)):
        if not var:
            return []
        return list(var)
    elif var is None:
        return []
    else:
        return [var]

def reverse_list(self, _list):
    """Return a reversed copy of a list."""

    rev = list(_list)
    rev.reverse()
    return rev

SConsEnvironment.check_dir_target = check_dir_target
SConsEnvironment.build_files = build_files
SConsEnvironment.build_headers = build_headers
SConsEnvironment.build_sources = build_sources
SConsEnvironment.build_metas = build_metas
SConsEnvironment.build_directory = build_directory
SConsEnvironment.build_plugin = build_plugin
SConsEnvironment.build_toolbox = build_toolbox
SConsEnvironment.build_module_register = build_module_register
SConsEnvironment.write_pkgconfig = write_pkgconfig
SConsEnvironment.strip_dir = strip_dir
SConsEnvironment.copy_file = copy_file
SConsEnvironment.get_build_path = get_build_path
SConsEnvironment.get_target_name = get_target_name
SConsEnvironment.glob = glob
SConsEnvironment.path_exists = path_exists
SConsEnvironment.make_list = make_list
SConsEnvironment.reverse_list = reverse_list



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


## Builder for generating meta files (courtesy of Walter Landry).
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
