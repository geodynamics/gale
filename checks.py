import utils

def compile_test(ctx, env, com, lang, src):
    """Generate a temporary source file and try to compile it. Return
    success or failure, error code, stdout, stderr and the file if
    the compile was successful."""

    utils.log("Trying to build test source:", post_indent=1)

    src_fn = ctx.make_temp_file(src, prefix="tmp", ext=ctx.lang_exts[lang][0])
    obj_fn = utils.path.replext(src_fn, ".o")
    res, out, err = com["compile"](com, lang, env, source=src_fn, target=obj_fn)

    info = {"src_fn": src_fn, "obj_fn": obj_fn,
            "error_code": res, "stdout": out, "stderr": err}
    success = com["compile_success"](res, out, err)

    if success:
        utils.log("Success.", post_indent=-1)
    else:
        utils.log("Failure.", post_indent=-1)

    return (success, info)

def link_test(ctx, env, lnk, obj_fn):
    """Try to link the object into a program using the provided
    environment."""

    utils.log("Trying to link test object:", post_indent=1)

    prog_fn = utils.path.remext(obj_fn)
    res, out, err = lnk["link_prog"](lnk, env, source=obj_fn, target=prog_fn)

    info = {"obj_fn": obj_fn, "prog_fn": prog_fn,
            "error_code": res, "stdout": out, "stderr": err}
    success = lnk["link_success"](res, out, err)

    if success:
        utils.log("Success.", post_indent=-1)
    else:
        utils.log("Failure.", post_indent=-1)

    return (success, info)


def gen_compilers(coms, gen_coms, gen_srcs):
    coms = utils.conv.to_list(coms)
    for com in gen_coms():
        for lang in com["languages"]:
            for src_id in gen_srcs(com, lang):
                yield (com, lang, src_id)


def check_libraries(ctx, env, output, coms, lnks, apply_com_dep, apply_lnk_dep,
                    gen_deps, gen_coms, gen_lnks, gen_srcs, aux_lib_cands):
    """From the set of available compilers determine which languages we can
    use to compile the package stored in 'pkg'. We build a dictionary of the
    form
      {<compiler>:
        {<language>:
          {"source": <program test source>,
           "id":     <source id>}}}
    We keep information about every compiler/language combination so we can
    query failed combos later.
      In order to properly check compilers/languages we need to link and run
    a program, so we'll also get all the information about linking for free:
      {<linker>:
        {"aux_libs":    [<required auxilliary libraries>, ...],
         "valid":       <True/False>,
         "status_info": {...}}}"""

    # Order the dependency candidates from smallest to largest. We do this to
    # make it easier to detect supersets.
    dep_cands = list(gen_deps())
    dep_cands.sort(cmp=lambda x,y: len(x)-len(y))

    # Places to store things.
    if "valid_deps" not in output:
        output["valid_deps"] = []
    if "com" not in output:
        output["com"] = {}
    if "lnk" not in output:
        output["lnk"] = {}

    # We need a flag to indicate if we've found any valid linkers
    # for the current set of dependencies.
    found_anything = False

    # We have to try every one of these. There may be more than one
    # valid set.
    for deps in dep_cands:

        utils.log("Trialing dependencies: %s"%repr(deps), post_indent=1)

        # Don't bother testing a dependency set if we have already identified
        # a subset as valid, because adding libraries to an already valid
        # set will not give us any new information.
        okay = True
        dep_set = set(deps)
        for vd in output["valid_deps"]:
            if set(vd).issubset(dep_set):
                utils.log("Found to be a superset of an existing set.", post_indent=-1)
                okay = False
                break
        if not okay:
            continue

        # We need a flag to indicate if we've found any valid linkers
        # for the current set of dependencies.
        found_com = False

        # We use the dependencies when making a run command, so store them on the
        # environment temporarily.
        env["_deps"] = deps

        # Try each compiled set of code and compiler with each linker.
        for com, lang, src_id in gen_compilers(coms, gen_coms, gen_srcs):

            if isinstance(src_id, tuple):
                utils.log("Trialing compiler: %s, %s, %s"%(repr(com), repr(lang), repr(src_id[1])),
                          post_indent=1)
                src = src_id[0]
            else:
                utils.log("Trialing compiler: %s, %s, no source ID"%(repr(com), repr(lang)),
                          post_indent=1)
                src = src_id

            # If we've already found a valid combo for this compiler/language
            # we can skip this round.
            if com in output["com"] and lang in output["com"][com]:
                utils.log("Already done.", post_indent=-1)
                continue

            # Add all our dependency's compile environments.
            okay = True
            com_env = env.clone()
            for d in deps:
                if not apply_com_dep(d, com, lang, com_env):
                    okay = False
                    break
            if not okay:
                utils.log("Compiler %s incompatible with dependency."%repr(com), post_indent=-1)
                continue

            # Compile the source code.
            res, com_info = compile_test(ctx, com_env, com, lang, src)
            if not res:
                utils.log.unindent()
                continue

            # Add the compiler's begin and end run-time objects to the environment.
            if "rt" in com:
                rt_bak = env.backup("rt_begin", "rt_end", "rt_libs")
                env.append_unique(com["rt"]["static"])

            # Keep a flag to indicate if we found any valid linkers.
            found_lnk = False

            # Search the linkers to see if we can link and run the current
            # compiler's compiled object.
            for lnk in gen_lnks():

                utils.log("Trialing linker: %s"%repr(lnk), post_indent=1)

                # If this linker is not in the info dict, add it now.
                lnk_out = output["lnk"].get(lnk, {})

                # If we've already found a set of auxilliary libraries for which
                # this linker works use them.
                if "aux_libs" in lnk_out:
                    cur_aux_lib_cands = [lnk_out["aux_libs"]]
                else:
                    cur_aux_lib_cands = aux_lib_cands.get(lang, [])

                # Make sure we have something to try.
                if not cur_aux_lib_cands:
                    cur_aux_lib_cands = [[]]

                # Keep track of which auxilliary library sets we've already checked
                # for this round.
                done_aux_libs = []

                # Keep a flag to indicate if we found a valid set of auxilliary libraries.
                found_aux_libs = False

                # Add all our dependency's link environments.
                okay = True
                lnk_env = env.clone()
                for d in deps:
                    if not apply_lnk_dep(d, lnk, lnk_env):
                        okay = False
                        break
                if not okay:
                    utils.log("Linker %s incompatible with dependency."%repr(lnk), post_indent=-1)
                    continue

                # Try each of the auxilliary library sets.
                for aux_libs in cur_aux_lib_cands:

                    # Need two loops here to see if we need to add the current compiler's
                    # run-time libraries.
                    for iter in range(2):

                        # First loop, do things normally. Second loop add run-time environment,
                        # but only if we have something to add.
                        if iter == 1 and com["rt"]["static"]["rt_libs"]:
#                             iter_bak = env.backup("lib_dirs", "rpaths")
                            lnk_env.append_unique("lib_dirs", make_list=True, *com["rt"]["static"]["rt_lib_dirs"])
                            lnk_env.append_unique("rpaths", make_list=True, *com["rt"]["static"]["rt_rpaths"])
                            aux_libs = aux_libs + com["rt"]["static"]["rt_libs"]

                        utils.log("Trialing auxilliary libraries: %s"%repr(aux_libs), post_indent=1)

                        # Remove any auxilliary libraries that appear in this linker's static
                        # run-time library set.
                        if "rt" in lnk:
                            new_aux_libs = []
                            for l in aux_libs:
                                if l not in lnk["rt"]["static"]["rt_libs"]:
                                    new_aux_libs.append(l)
                                else:
                                    utils.log("%s is in this linker's run-time environment already."%repr(l))
                            aux_libs = new_aux_libs

                        # If we've already attempted an auxilliary library set like this then
                        # skip it.
                        if set(aux_libs) in done_aux_libs:
                            utils.log("Already tried this set of auxilliary libraries.", post_indent=-1)
#                             if iter == 1:
#                                 env.restore(iter_bak)
                            continue

                        # Add to the list of tried auxilliary libraries.
                        done_aux_libs.append(set(aux_libs))

                        # Temporarily add the auxilliray libraries to the environment.
                        lib_bak = lnk_env.backup("libs")
                        lnk_env.append_unique("libs", make_list=True, *aux_libs)

                        # Try to build the program.
                        res, lnk_info = link_test(ctx, lnk_env, lnk, com_info["obj_fn"])

                        # Restore any environment backups.
                        lnk_env.restore(lib_bak)
#                         if iter == 1:
#                             env.restore(iter_bak)

                        # Check results of the link step.
                        if not res:
                            utils.log.unindent()
                            continue

                        # Try running the program. We need to run to ensure the
                        # dependencies don't cause conflicts.
                        res, run_info = lnk_env["run_test_prog"](lnk_env, lnk_info["prog_fn"])
                        if not res:
                            utils.log.unindent()
                            continue
                        lnk_info.update(run_info)

                        # Success.
                        found_aux_libs = True
                        found_com = True
                        found_lnk = True
                        found_anything = True
                        utils.log.unindent()
                        break

                    if found_aux_libs:
                        break

                # If we found a valid set of auxilliary libraries, update the output
                # dictionary.
                utils.log.unindent()
                if found_aux_libs:

                    # If we already have a set of auxilliary libraries for this linker
                    # make sure they are consistent.
                    if "aux_libs" in lnk_out:
                        if set(lnk_out["libs"]) != set(aux_libs):
                            raise "Inconsistent auxilliary libraries."
                    else:
                        lnk_out["libs"] = aux_libs
                        if iter == 1:
                            lnk_out["lib_dirs"] = com["rt"]["static"]["rt_lib_dirs"]
                            lnk_out["rpaths"] = com["rt"]["static"]["rt_rpaths"]

                    # Store the set on the output dictionary.
                    if lnk not in output["lnk"]:
                        output["lnk"][lnk] = lnk_out

            # Restore the run-time environment.
            if "rt" in com:
                env.restore(rt_bak)

            # Update the compilers/languages.
            utils.log.unindent()
            if found_lnk:
                output["com"][com] = {}
                output["com"][com][lang] = utils.Environment()
                if isinstance(src_id, tuple):
                    output["com"][com][lang]["pp_defs"] = src_id[1]

        utils.log.unindent()
        if found_com:
            output["valid_deps"].append(deps)

    # Remove the deps from the environment.
    if "_deps" in env:
        del env["_deps"]

    return found_anything
