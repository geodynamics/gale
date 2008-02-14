Import('env')

# Need to copy the environment for use here.
env = env.Copy()

# Add some extra stuff.
env.proj = 'StGermain'
env.Append(CPPPATH=['#build/include/' + env.proj])

#
# Target specification.
#

# TODO: Build sources if csoap present.
env.build_headers(env.glob('Regresstor/libRegresstor/src/*.h'),
                  'StGermain/Regresstor/libRegresstor')

# Some special cases we need to copy first.
env.build_headers(['Base/Foundation/src/ClassSetup.h',
                   'Base/Foundation/src/ClassEmpty.h'],
                  'StGermain/Base/Foundation', force_copy=True)

# We have this funny 'mpirecord' business to handle now.
if env['useMpiRecord']:
    # TODO
    pass
else:
    env.build_headers('Base/IO/src/mpirecord/none/mpimessaging.h',
                      'StGermain/Base/IO/mpirecord')

env.build_directory('Base/Foundation')
env.build_directory('Base/IO')
env.build_directory('Base/Container')
env.build_directory('Base/Automation')
env.build_directory('Base/Extensibility')
env.build_directory('Base/Context')
env.build_directory('Base')
env.build_directory('Utils')

env.build_headers(env.glob('libStGermain/src/*.h'), 'StGermain')
env.build_objects(env.glob('libStGermain/src/*.c'), 'libStGermain')
env.build_metadata(env.glob('libStGermain/src/*.meta'), 'libStGermain')

env.build_library(env.get_hnodes(env.SharedObject), 'StGermain')

env.build_tests(env.glob('libStGermain/tests/test*.c'), 'StGermain',
                libs='StGermain')

env.Program('#build/bin/FlattenXML', 'Base/FlattenXML/src/main.c',
            LIBS=['StGermainBase'] + env.get('LIBS', []))
env.Program('#build/bin/StGermain', 'src/main.c',
            LIBS=['StGermain'] + env.get('LIBS', []))
