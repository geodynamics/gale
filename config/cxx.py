import checks

def apply_cxx(env, cxx):
    if cxx in ['mpicxx', 'g++']:
        pass
    elif cxx == 'xlc':
        pass

    conf = env.Configure(
        custom_tests={
            'CheckCXXFixed': checks.CheckCXXFixed
            }
        )
    if not conf.CheckCXXFixed():
        env.Exit()
    return conf.Finish()
