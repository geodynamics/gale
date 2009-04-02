import utils

class Environment(object):

    def __init__(self):
        self._dict = {}

    def __getitem__(self, k):
        return self._dict[k]

    def __setitem__(self, k, v):
        self._dict[k] = v

    def __delitem__(self, k):
        del self._dict[k]

    def __contains__(self, k):
        return k in self._dict

    def iteritems(self):
        return self._dict.iteritems()

    def get(self, k, d):
        return self._dict.get(k, d)

    def clear(self):
        self._dict = {}

    def append(self, k, *args, **kw):
        if k in self._dict:
            if not isinstance(self._dict[k], list):
                self._dict[k] = utils.conv.to_list(self._dict[k])
            self._dict[k] += list(args)
        elif len(args) == 1:
            if kw.get("make_list", False):
                self._dict[k] = [args[0]]
            else:
                self._dict[k] = args[0]
        elif len(args) > 1:
            self._dict[k] = list(args)

    def append_unique(self, k, *args, **kw):
        self.append(k, *self._make_unique(k, args), **kw)

    def prepend(self, k, *args, **kw):
        if k in self._dict:
            if not isinstance(self._dict[k], list):
                self._dict[k] = utils.conv.to_list(self._dict[k])
            self._dict[k] = list(args) + self._dict[k]
        elif len(args) == 1:
            if kw.get("make_list", False):
                self._dict[k] = [args[0]]
            else:
                self._dict[k] = args[0]
        elif len(args) > 1:
            self._dict[k] = list(args)

    def prepend_unique(self, k, *args, **kw):
        self.prepend(k, *self._make_unique(k, args), **kw)

    def _make_unique(self, k, args):
        if k not in self._dict:
            return args
        cur = utils.conv.to_list(self._dict[k])
        new_args = []
        for v in args:
            if v not in cur:
                new_args.append(v)
        return new_args

    def subst(self, text):
        return utils.macro.subst(text, self)
