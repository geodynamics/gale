import os
Import('env')

env = env.Copy()
env.AppendUnique(CPPPATH=[os.path.join(env['buildPath'], # Add gLucifer include path.
                                       'include',
                                       'glucifer')])
env.project_name = 'glucifer' # Need a project name.
env.clear_all() # ... so that our structures are ready.

env.build_directory('Base')
env.build_directory('Windowing')
env.build_directory('RenderingEngines')
env.build_directory('OutputFormats')
env.build_directory('InputFormats')
env.build_directory('DrawingObjects')
env.build_directory('WindowInteractions')

env.build_headers(env.glob('libglucifer/src/*.h'), 'glucifer')
env.build_objects(env.glob('libglucifer/src/*.c'), 'libglucifer')
env.build_metadata(env.glob('libglucifer/src/*.meta'), 'libglucifer')

env.build_library(env.get_hnodes(env.SharedObject), 'glucifer')

env.build_tests(env.glob('libglucifer/tests/test*.c'),
                'glucifer', libs='glucifer')

env.build_plugin('plugins/lucPlugin', 'lucPlugin', prefix=False)

env.build_xmls(env.glob('ModelComponents/*.xml'), 'StGermain/glucifer')
