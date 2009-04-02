import re
import conv, format

class OptionVariant(object):

    def __init__(self, flgs, type="string", sep=[" ", "="], delim=",", **kw):
        flgs = conv.to_list(flgs)
        for i in range(len(flgs) - 1):
            if flgs[i] in flgs[i + 1:]:
                raise "Duplicate option string given to Option."
        sep = conv.to_list(sep)
        for i in range(len(sep) - 1):
            if sep[i] in sep[i + 1:]:
                raise "Duplicate separators in list."
        for s in sep:
            if len(s) > 1:
                raise "Separator of length > 1."
        delim = conv.to_list(delim)
        for i in range(len(delim) - 1):
            if delim[i] in delim[i + 1:]:
                raise "Duplicate delimiters in list."
        for d in delim:
            if len(d) > 1:
                raise "Delimeter of length > 1."
            if d == " ":
                raise "Whitespace delimeter."

        self.flgs = flgs
        self.type = type
        self.sep = conv.to_list(sep)
        self._pref_sep = kw.get("pref_sep", self.sep[0])
        self.sep.sort(lambda x,y: len(y)-len(x))
        self.delim = conv.to_list(delim)
        self._pref_delim = self.delim[0]
        self.delim.sort(lambda x,y: len(y)-len(x))

        if self.type == "bool":
            self.to_str = self._bool_to_str
            self.parse_val = self._parse_bool
            bf = kw.get("bool_flags", (["", "yes", "1", "t"], ["no", "0", "f"]))
            self._pref_true_flag = kw.get("pref_true_flag", conv.to_list(bf[0])[0])
            self._pref_false_flag = kw.get("pref_false_flag", conv.to_list(bf[1])[0])
            self.bool_flags = (zip(conv.to_list(bf[0]), [True for i in range(len(bf[0]))]) +
                               zip(conv.to_list(bf[1]), [False for i in range(len(bf[1]))]))
            self.bool_flags.sort(lambda x,y: len(y[0])-len(x[0]))
        elif self.type == "enum":
            self.to_str = self._enum_to_str
            self.parse_val = self._parse_enum
            self.enum = kw.get("enum", {})
            if isinstance(self.enum, list):
                self.enum = dict(zip(self.enum, self.enum))
        else:
            self.to_str = self._str_to_str
            self.parse_val = self._parse_str

    def _bool_to_str(self, val):
        if val:
            f = self._pref_true_flag
        else:
            f = self._pref_false_flag
        if f:
            s = "%s%s%s" % (self.flgs[0], self._pref_sep, f)
        else:
            s = "%s" % self.flgs[0]
        return s.strip()

    def _parse_bool(self, val_str):
        for f, v in self.bool_flags:
            if f == "":
                return (v, False)
            if f == val_str:
                return (v, True)
        return None

    def _enum_to_str(self, val):
        val = conv.to_list(val)
        s = []
        for v in val:
            e = self.enum.get(v, None)
            if e is None:
                if v not in self.enum.values():
                    raise "Invalid enumerated option."
                e = v
            s.append(e)
        return self._str_to_str(self._pref_delim.join(s))

    def _parse_enum(self, val):
        if not val:
            return None
        if val not in self.enum:
            raise "Invalid enumerated option."
        return (val, True)

    def _str_to_str(self, val):
        val = conv.to_list(val)
        return "%s%s%s" % (self.flgs[0], self._pref_sep, self._pref_delim.join(val))

    def _parse_str(self, val):
        if not val:
            return None
        return (val, True)

    def match_sep(self, val_str, b):
        for s in self.sep:
            ls = len(s)
            if ls <= len(val_str) - b and s == val_str[b:b + ls]:
                return b + ls
        return None

    def __repr__(self):
        return "/".join(self.flgs)

class Option(object):

    def __init__(self, name, *args, **kw):
        self.name = name
        vars = kw.get("variants", None)
        if vars is None:
            vars = OptionVariant(*args, **kw)
        self._vars = conv.to_list(vars)
        self.type = self._vars[0].type
        self._flg_to_var = {}
        for v in self._vars:
            if self.type != v.type:
                raise "Inconsistent OptionVariant types."
            for f in v.flgs:
                if f in self._flg_to_var:
                    raise "Duplicate option flags in Option."
                self._flg_to_var[f] = v
        self.help = kw.get("help", None)
        self.default = kw.get("default", None)

    def get_variant(self, flg):
        return self._flg_to_var[flg]

    def make_help_string(self):
        v = self._vars[0]
        if len(v.flgs) > 1:
            help = "%s" % "|".join(v.flgs)
        else:
            help = str(v.flgs[0])
        if v.type == "string":
            help += v._pref_sep + "<string>"
        elif v.type == "enum":
            help += v._pref_sep + "<%s>" % "|".join(v.enum.values())
        elif v.type == "bool":
            if v._pref_true_flag is "":
                help = "[" + help + "]"
            else:
                help += v._pref_sep + "<" + v._pref_true_flag + \
                    "|" + v._pref_false_flag + ">"
        if self.default is not None:
            help += "\n\tdefault: %s" % str(self.default)
        if self.help:
            txt = format.box(self.help, 60)
            help += "\n\t" + "\n\t".join(txt.split("\n"))
        return help

    def __repr__(self):
        s = []
        for v in self._vars:
            s.append(str(v))
        return "/".join(s)

class OptionSet(object):

    def __init__(self):
        self.options = {}
        self._order = []
        self._flg_to_opt = {}
        self._expr = None

    def __contains__(self, opt_name):
        return opt_name in self.options

    def add_option(self, opts, *args):
        opts = conv.to_list(opts)
        opts += args
        for o in opts:
            if o.name in self.options:
                raise RuntimeError, "Duplicate option name %s."%repr(o.name)
            self.options[o.name] = o
            self._order.append(o.name)
            for f in o._flg_to_var.iterkeys():
                if f in self._flg_to_opt:
                    raise "Duplicate option flag in OptionSet."
                self._flg_to_opt[f] = o
        self._expr = None

    def remove_option(self, name):
        for f in self.options[name]._flg_to_var.iterkeys():
            del self._flg_to_opt[f]
        del self.options[name]
        del self._order[self._order.index(name)]
        self._expr = None

    def pop(self, name):
        opt = self.options[name]
        self.remove_option(name)
        return opt

    def has_option(self, name):
        return name in self.options

    def merge(self, opts):
        for name in opts._order:
            if name not in self.options:
                self.add_option(opts.options[name])

    def merge_option(self, opt):
        if opt.name not in self.options:
            self.add_option(opt)
        else:
            # TODO: won't matter for now, but need to do this.
            pass

    def parse_option_string(self, opt_str):
        return self.parse_option_list(opt_str.split())

    def parse_option_list(self, words):
        found_opts = set()
        if not self._expr:
            self._make_expr()
        opt_l = []
        while len(words):
            w = words.pop(0)
            m = self._expr.match(w)
            if m is None or not m.groups():
                opt_l.append(w)
                continue
            e = m.end(m.lastindex)
            flg = w[:e]
            opt = self._flg_to_opt[flg]
            var = opt.get_variant(flg)
            if e < len(w):
                i = var.match_sep(w, e)
                if i is not None:
                    next = w[i:]
                    res = var.parse_val(next)
                    if res is None or not res[1]:
                        opt_l.append(w)
                        continue
                    opt_l.append((opt.name, res[0]))
                    found_opts.add(opt.name)
                else:
                    opt_l.append(w)
                    continue
            elif " " in var.sep:
                if not len(words):
                    next = ""
                else:
                    next = words.pop(0)
                res = var.parse_val(next)
                if res is None:
                    opt_l.append(w)
                    if next:
                        words.insert(0, next)
                    continue
                opt_l.append((opt.name, res[0]))
                found_opts.add(opt.name)
                if not res[1] and next:
                    words.insert(0, next)
            else:
                opt_l.append(w)
                continue
        opt_l += self.make_defaults(found_opts)
        return opt_l

    def make_defaults(self, found_opts):
        defs = []
        for name, opt in self.options.iteritems():
            if opt.default is None or name in found_opts:
                continue
            defs.append((opt.name, opt._vars[0].parse_val(opt.default)[0]))
        return defs

    def make_option_string(self, args):
        s = []
        for a in args:
            if isinstance(a, tuple) and len(a) == 2:
                s.append(self.options[a[0]]._vars[0].to_str(a[1]))
            elif isinstance(a, str):
                s.append(a)
            else:
                s.append(str(a))
        return " ".join(s)

    def gather(self, args):
        d = {}
        for a in args:
            if isinstance(a, tuple) and len(a) == 2:
                if a[0] in d:
                    if not isinstance(d[a[0]], list):
                        d[a[0]] = conv.to_list(d[a[0]])
                    d[a[0]].append(a[1])
                else:
                    d[a[0]] = a[1]
        return d

    def make_help_string(self):
        help = []
        for o in self._order:
            help.append(self.options[o].make_help_string())
        return "\n".join(help) + "\n"

    def _make_expr(self):
        flgs = []
        for f, o in self._flg_to_opt.iteritems():
            flgs.append((f, o))
        flgs.sort(lambda x,y: len(y[0])-len(x[0]))
        l = []
        for f in flgs:
            l.append("(%s)" % f[0])
        self._expr = re.compile("|".join(l))

    def __repr__(self):
        return str(self.options)
