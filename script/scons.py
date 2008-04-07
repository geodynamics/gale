import os
from SCons.Script.SConscript import SConsEnvironment

#
# Helpful build utilities.
#

def build_files(env, files, dst_dir):
    for f in files:
        env.copy_file(env.get_build_path(dst_dir) + '/' + os.path.basename(f), f)

def build_headers(env, headers, dst_dir):
    for hdr in headers:
        env.Install(env.get_build_path(dst_dir), hdr)

def build_sources(env, sources, dst_dir):
    objs = []
    mod_name = ('CURR_MODULE_NAME', env['ESCAPE']('"' + ''.join(dst_dir.split('/')) + '"'))
    for src in sources:
        tgt = env.get_build_path(dst_dir + '/' + os.path.basename(src)[:-2])
        objs += env.SharedObject(tgt, src, CPPDEFINES=[mod_name] + env.get('CPPDEFINES', []))
    return objs

def build_metas(env, metas, dst_dir):
    objs = []
    mod_name = ('CURR_MODULE_NAME', env['ESCAPE']('"' + ''.join(dst_dir.split('/')) + '"'))
    for m in metas:
        src = env.Meta(env.get_build_path(dst_dir + '/' + os.path.basename(m)[:-5]), m)
        tgt = env.get_build_path(dst_dir + '/' + os.path.basename(src[0].path)[:-2])
        objs += env.SharedObject(tgt, src, CPPDEFINES=[mod_name] + env.get('CPPDEFINES', []))
    return objs

SConsEnvironment.build_files = build_files
SConsEnvironment.build_headers = build_headers
SConsEnvironment.build_sources = build_sources
SConsEnvironment.build_metas = build_metas

#
# Custom builders.
#

# Builder for generating meta files (courtesy of Walter Landry).
def create_meta(target, source, env):
    output_file = file(str(target[0]),'wb')
    output_file.write("#define XML_METADATA \"")
    xml_file = file(str(source[0]))
    xml_lines = xml_file.readlines()
    for l in xml_lines:
        output_file.write(l.replace('\"','\\\"')[:-1])
    output_file.write("\"\n#define COMPONENT_NAME ")
    for l in xml_lines:
        start=l.find('<param name="Name">')
        if start!=-1:
            end=l.find('<',start+19)
            if end==-1:
                raise RunTimeError('Malformed XML file.  The file '
                                   + str(source[0])
                                   + ' does not close off <param name="Name"> on the same line.')
            output_file.write(l[start+19:end])
            output_file.write("\n")
            break
    output_file.write('''
#define Stg_Component_Stringify( str ) #str

/* Note: Two macros are used to resolve the the extra macro level */
#define Stg_Component_Metadata_Create( name ) Stg_Component_Metadata_Create_Macro( name )
#define Stg_Component_Metadata_Create_Macro( name ) \
	const char* name ##_Meta = XML_METADATA; \
	const char* name ##_Name = #name; \
	const char* name ##_Version = VERSION; \
	const char* name ##_Type_GetMetadata() { /* hack...won't be needed when hierarchy rollout is done */\
		return name ##_Meta; \
	} \
	const char* name ##_GetMetadata() { \
		return name ##_Meta; \
	} \
	const char* name ##_GetName() { \
		return name ##_Name; \
	} \
	const char* name ##_GetVersion() { \
		return name ##_Version; \
	}

#if defined(COMPONENT_NAME) && defined(VERSION) && defined(XML_METADATA)
	Stg_Component_Metadata_Create( COMPONENT_NAME )
#endif
''')

def gen_meta_suffix(env, sources):
    return "-meta.c"

Import('env')
env['BUILDERS']['Meta']=Builder(action=create_meta,src_suffix="meta",
                                suffix=gen_meta_suffix,single_source=True)
