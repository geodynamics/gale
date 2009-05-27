import os
import config
from Module import Module
import utils

class Tool(Module):

    def __init__(self, ctx, **kw):
        Module.__init__(self, ctx, **kw)
        self.commands = []
        self.result = None

    @staticmethod
    def new_candidate(ctx, cmd, *args, **kw):
        cand = config.Config(ctx, cmd, command=cmd, *args, **kw)
        cand["options"] = utils.options.OptionSet()
        return cand

    @staticmethod
    def new_plan(mod, *args, **kw):
        return Tool.Plan(mod, *args, **kw)

    def copy_candidate(self, cand):
        cpy = Module.copy_candidate(self, cand)
        self.merge_candidates(cpy, cand)
        return cpy

    @staticmethod
    def merge_candidates(cand, trial):
        Module.merge_candidates(cand, trial)
        cand["options"].merge(trial["options"])

    @staticmethod
    def same_candidates(c0, c1):
        if not config.Module.same_candidates(c0, c1):
            return False
        return os.path.realpath(c0["command"]) == os.path.realpath(c1["command"])

    def _search(self):
        for d in self.ctx.tool_dirs:
            for c in self.commands:
                p = os.path.join(d, c)
                if not os.path.exists(p):
                    continue
                nc = self.ctx.new_candidate(self, p)
                self.add_candidate(nc)

    def setup_trial(self, cfg):
        Module.setup_trial(self, cfg)
        cfg["execute"] = self.execute
        cfg["make_command_line"] = self.make_command_line
        cfg["success"] = self.success
        cfg["string"] = self.expand_string_arg
        cfg["bool"] = self.expand_bool_arg
        cfg["strip"] = self.strip_list
        cfg["extend"] = self.extend_list

    def test_option(self, cfg, opts, err_okay=False,
                    products=[], args=[]):
        opts = utils.conv.to_list(opts)
        products = utils.conv.to_list(products)

        utils.log("Testing option: %s"%opts[0].name, post_indent=1)

        for o in opts:
            utils.log("Trying flag: %s"%str(o), post_indent=1)

            if cfg["options"].has_option(o.name):
                old_opt = cfg["options"].pop(o.name)
            else:
                old_opt = None
            cfg["options"].add_option(o)

            res, out, err = cfg["execute"](cfg, args=args)
            self.result = (res, out, err)
            if cfg["success"](res, out, err, err_okay, products):
                utils.log("Success.", post_indent=-2)
                if old_opt:
                    cfg["options"].merge_option(old_opt)
                return True

            cfg["options"].remove_option(o.name)
            if old_opt:
                cfg["options"].add_option(old_opt)

            utils.log("Failed.", post_indent=-1)
        utils.log.unindent()
        return False

    def execute(self, cfg, cmd="", env={}, args=[], **kw):
        cmd_str = self.make_command_line(cfg, cmd, env, args, **kw)

        # Run.
        utils.log("Executing: %s"%cmd_str, level=20, post_indent=1)
        ec, out, err = utils.command.run(cmd_str)

        # Log results.
        utils.log("Error code: %d"%ec, level=20)
        if out:
            utils.log("Stdout:", level=10, post_indent=1)
            utils.log(out, level=10, post_indent=-1)
        if err:
            utils.log("Stderr:", level=10, post_indent=1)
            utils.log(err, level=10, post_indent=-1)
        utils.log.unindent()

        return ec, out, err

    def make_command_line(self, cfg, cmd="", env={}, args=[], **kw):
        # If we were given a command environment base the command around that.
        if cmd:

            # Combine all the environments.
            local = cfg.clone()
            local.update(env)
            local.update(kw)

            # Produce a command string.
            cmd_str = local.subst("$%s"%cmd)

        else:
            cmd_str = ""

        # If we were given explicit arguments...
        if args:
            arg_str = cfg["options"].make_option_string(args)
        else:
            arg_str = ""

        # Combine the two strings.
        if cmd_str:
            cmd_str = "%s %s"%(cmd_str, arg_str)
        else:
            cmd_str = "%s %s"%(cfg["command"], arg_str)

        return cmd_str

    def success(self, res, out, err, err_okay=True, products=[]):
        products = utils.conv.to_list(products)
        dump = out + err
        return not (res or (not err_okay and err) or
                    self._check_output(dump) or
                    self._check_products(products))

    def expand_string_arg(self, env, opt_name, opt_val):
        if not opt_val:
            return ""
        opt_val = utils.conv.to_list(opt_val)
        cmd_args = [(opt_name, v) for v in opt_val]
        res = env["options"].make_option_string(cmd_args)
        return res

    def expand_bool_arg(self, env, opt_names):
        if not opt_names:
            return ""
        opt_names = utils.conv.to_list(opt_names)
        cmd_args = [(o, True) for o in opt_names]
        return env["options"].make_option_string(cmd_args)

    def strip_list(self, env, cmp_list, src_list):
        if src_list == "":
            src_list = []
        cmp_list = utils.conv.to_list(cmp_list)
        src_list = utils.conv.to_list(src_list)
        return [i for i in src_list if i not in cmp_list]

    def extend_list(self, env, base_list, ext_list):
        if base_list == "":
            base_list = []
        if ext_list == "":
            ext_list = []
        base_list = utils.conv.to_list(base_list)
        ext_list = utils.conv.to_list(ext_list)
        return base_list + ext_list

    def _check_output(self, output):
        return not (output.find("unrecognised") == -1 and
                    output.find("unrecognized") == -1)

    def _check_products(self, prods):
        for p in prods:
            if not os.path.exists(p): return True
        return False

    def get_summary_items(self, cfg):
        itms = []
        for k in ["command"]:
            if k in cfg:
                itms.append((k, cfg[k]))
        return itms

    def setup_options(self):
        Module.setup_options(self)
        if "location" not in self.options:
            self.options.add_option(
                utils.options.Option("location", "--" + self.opt_name, pref_sep="=",
                                     help=self.help["location"]))

    def process_options(self):
        Module.process_options(self)
        loc = self.get_option("location", None)
        if loc is not None:
            self.set_mode(3) # Need it, also don't search for others.
            self.do_search = False
            cfg = self.ctx.new_candidate(self, loc)
            self.add_candidate(cfg)

    def make_help_dict(self):
        Module.make_help_dict(self)
        self.help["location"] = "Specify the command to use for %s."%self.opt_name
