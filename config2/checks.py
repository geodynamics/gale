import os, SCons.SConf

def CheckLibs(context, library = None, symbol = 'main',
              header = None, language = None, autoadd = 1, extra_libs=[]):
    """
    Copied from SConf.py but modified to check all libraries at once
    instead of one at a time.
    """

    if library == []:
        library = [None]

    if not SCons.Util.is_List(library):
        library = [library]

    # ToDo: accept path for the library
    res = SCons.Conftest.CheckLib(context, [library[0]], symbol, header = header,
                                  language = language, autoadd = autoadd,
                                  extra_libs=library[1:] + extra_libs)
    context.did_show_result = 1
    return not res

def CheckLibsWithHeader(context, libs, header, language,
                       call = None, autoadd = 1, extra_libs=[]):
    # ToDo: accept path for library. Support system header files.
    """
    Copied from SConf.py, but modified to check all libraries at once
    instead of one at a time.
    """
    prog_prefix, dummy = \
        SCons.SConf.createIncludesFromHeaders(header, 0)
    if libs == []:
        libs = [None]

    if not SCons.Util.is_List(libs):
        libs = [libs]

    res = SCons.Conftest.CheckLib(context, [libs[0]], None, prog_prefix,
                                  call = call, language = language, autoadd = autoadd,
                                  extra_libs=libs[1:] + extra_libs)
    context.did_show_result = 1
    return not res

def CheckCCFixed(context):
    res = SCons.Conftest.CheckCC(context)
    context.did_show_result = 1
    return not res
