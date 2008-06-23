import os
import SConfig

def get_all_headers(cfg, hdrs):
    """For the provided configuration, return a list of all the
    headers required by it's dependencies."""

    # Process the dependencies first, to maintain the correct order.
    for dep in cfg.deps:
        get_all_headers(dep, hdrs)

    # Add the current set of headers, making sure we don't include
    # any duplicates.
    cfg_hdrs = cfg.inst.pkg.env.make_list(cfg.hdrs)
    for h in cfg_hdrs:
        if h in hdrs: continue
        if cfg.inst.fwork:
            hdrs += ['%s/%s' % (cfg.inst.fwork, h)]
        else:
            hdrs += [h]

def get_header_source(cfg):
    """From the provided configuration, return a string to use to check
    that the headers are compatible with a C compiler."""

    # Begin with standard ansi C headers.
    src = '#include<stdlib.h>\n#include<stdio.h>\n#include<string.h>\n'

    # To be safe, we have to include all the required headers from all
    # our dependencies before including any from this installation.
    hdrs = []
    get_all_headers(cfg, hdrs)

    # Now convert the list of headers into strings to be appended to the
    # source string and return.
    for h in hdrs:
        src += '#include<' + h + '>\n'
    return src

def check_headers(cfg):
    """Determine if the required headers are available by using the
    current settings in the provided configuration."""

    # Alias the installation and package for easy access.
    inst = cfg.inst
    pkg = inst.pkg

    # If there are no headers to check, automatically pass.
    if not cfg.hdrs:
        return True

    # We need to be able to locate the header files themselves before
    # trying any sanity tests. This helps make sure we're not getting
    # false positives from headers being found automatically in default
    # compiler locations.
    if not inst.fwork:
        for hdr in cfg.hdrs:
            found = False
            for hdr_dir in inst.hdr_dirs:

                # Combine the header sub-directories with the base directory
                # and header file.
                path = os.path.join(inst.base_dir, hdr_dir, hdr)

                # Check if the path exists.
                if os.path.exists(path):

                    # If so, we can break out of the directory loop and move
                    # on to the next header.
                    found = True
                    break

            # If we couldn't find this header, the configuration is invalid.
            # We can return negative straight away.
            if not found: return False

    # If we're here, we were able to locate all the header files.
    # Now we try to compile a piece of source code to ensure the headers
    # are compatible with a compiler. First we need the source code
    # to compile.
    src = get_header_source(cfg)

    # We'll need to backup the existing environment state when we enable
    # the current configuration.
    old_state = {}
    cfg.enable(pkg.env, old_state)

    # Try compiling the source and revert the environment.
    result = pkg.compile_source(src)
    pkg.env.Replace(**old_state)

    # Return the results.
    return result
