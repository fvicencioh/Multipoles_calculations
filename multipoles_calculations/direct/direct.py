# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_direct')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_direct')
    _direct = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_direct', [dirname(__file__)])
        except ImportError:
            import _direct
            return _direct
        try:
            _mod = imp.load_module('_direct', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _direct = swig_import_helper()
    del swig_import_helper
else:
    import _direct
del _swig_python_version_info

try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except __builtin__.Exception:
    class _object:
        pass
    _newclass = 0


def computeDiagonal(VL: 'double *', KL: 'double *', VY: 'double *', KY: 'double *', triangle: 'double *', centers: 'double *', kappa: 'double', K_diag: 'double', V_diag: 'double', xk: 'double *', wk: 'double *') -> "void":
    return _direct.computeDiagonal(VL, KL, VY, KY, triangle, centers, kappa, K_diag, V_diag, xk, wk)
computeDiagonal = _direct.computeDiagonal

def direct_sort(K_aux: 'double *', V_aux: 'double *', LorY: 'int', K_diag: 'double', V_diag: 'double', IorE: 'int', triangle: 'double *', tri: 'int *', k: 'int *', xi: 'double *', yi: 'double *', zi: 'double *', s_xj: 'double *', s_yj: 'double *', s_zj: 'double *', xt: 'double *', yt: 'double *', zt: 'double *', m: 'double *', mx: 'double *', my: 'double *', mz: 'double *', mKc: 'double *', mVc: 'double *', interList: 'int *', offTar: 'int *', sizeTar: 'int *', offSrc: 'int *', offTwg: 'int *', targets: 'int *', Area: 'double *', sglInt_int: 'double *', sglInt_ext: 'double *', xk: 'double *', wk: 'double *', Xsk: 'double *', Wsk: 'double *', kappa: 'double', threshold: 'double', eps: 'double', w0: 'double', aux: 'double *') -> "void":
    return _direct.direct_sort(K_aux, V_aux, LorY, K_diag, V_diag, IorE, triangle, tri, k, xi, yi, zi, s_xj, s_yj, s_zj, xt, yt, zt, m, mx, my, mz, mKc, mVc, interList, offTar, sizeTar, offSrc, offTwg, targets, Area, sglInt_int, sglInt_ext, xk, wk, Xsk, Wsk, kappa, threshold, eps, w0, aux)
direct_sort = _direct.direct_sort

def directKt_sort(Ktx_aux: 'double *', Kty_aux: 'double *', Ktz_aux: 'double *', LorY: 'int', triangle: 'double *', k: 'int *', s_xj: 'double *', s_yj: 'double *', s_zj: 'double *', xt: 'double *', yt: 'double *', zt: 'double *', m: 'double *', mKc: 'double *', interList: 'int *', offTar: 'int *', sizeTar: 'int *', offSrc: 'int *', offTwg: 'int *', Area: 'double *', Xsk: 'double *', Wsk: 'double *', kappa: 'double', threshold: 'double', eps: 'double', aux: 'double *') -> "void":
    return _direct.directKt_sort(Ktx_aux, Kty_aux, Ktz_aux, LorY, triangle, k, s_xj, s_yj, s_zj, xt, yt, zt, m, mKc, interList, offTar, sizeTar, offSrc, offTwg, Area, Xsk, Wsk, kappa, threshold, eps, aux)
directKt_sort = _direct.directKt_sort

def direct_c(K_aux: 'double *', V_aux: 'double *', LorY: 'int', K_diag: 'double', V_diag: 'double', IorE: 'int', triangle: 'double *', tri: 'int *', k: 'int *', xi: 'double *', yi: 'double *', zi: 'double *', s_xj: 'double *', s_yj: 'double *', s_zj: 'double *', xt: 'double *', yt: 'double *', zt: 'double *', m: 'double *', mx: 'double *', my: 'double *', mz: 'double *', mKc: 'double *', mVc: 'double *', targets: 'int *', Area: 'double *', sglInt_int: 'double *', sglInt_ext: 'double *', xk: 'double *', wk: 'double *', Xsk: 'double *', Wsk: 'double *', kappa: 'double', threshold: 'double', eps: 'double', w0: 'double', aux: 'double *') -> "void":
    return _direct.direct_c(K_aux, V_aux, LorY, K_diag, V_diag, IorE, triangle, tri, k, xi, yi, zi, s_xj, s_yj, s_zj, xt, yt, zt, m, mx, my, mz, mKc, mVc, targets, Area, sglInt_int, sglInt_ext, xk, wk, Xsk, Wsk, kappa, threshold, eps, w0, aux)
direct_c = _direct.direct_c

def direct_c_derivative(dKx_aux: 'double *', dKy_aux: 'double *', dKz_aux: 'double *', dVx_aux: 'double *', dVy_aux: 'double *', dVz_aux: 'double *', LorY: 'int', K_diag: 'double', V_diag: 'double', IorE: 'int', triangle: 'double *', tri: 'int *', k: 'int *', xi: 'double *', yi: 'double *', zi: 'double *', s_xj: 'double *', s_yj: 'double *', s_zj: 'double *', xt: 'double *', yt: 'double *', zt: 'double *', m: 'double *', mx: 'double *', my: 'double *', mz: 'double *', mKc: 'double *', mVc: 'double *', targets: 'int *', Area: 'double *', sglInt_int: 'double *', sglInt_ext: 'double *', xk: 'double *', wk: 'double *', Xsk: 'double *', Wsk: 'double *', kappa: 'double', threshold: 'double', eps: 'double', w0: 'double', aux: 'double *') -> "void":
    return _direct.direct_c_derivative(dKx_aux, dKy_aux, dKz_aux, dVx_aux, dVy_aux, dVz_aux, LorY, K_diag, V_diag, IorE, triangle, tri, k, xi, yi, zi, s_xj, s_yj, s_zj, xt, yt, zt, m, mx, my, mz, mKc, mVc, targets, Area, sglInt_int, sglInt_ext, xk, wk, Xsk, Wsk, kappa, threshold, eps, w0, aux)
direct_c_derivative = _direct.direct_c_derivative

def direct_c_2derivative(dKxx_aux: 'double *', dKxy_aux: 'double *', dKxz_aux: 'double *', dKyx_aux: 'double *', dKyy_aux: 'double *', dKyz_aux: 'double *', dKzx_aux: 'double *', dKzy_aux: 'double *', dKzz_aux: 'double *', dVxx_aux: 'double *', dVxy_aux: 'double *', dVxz_aux: 'double *', dVyx_aux: 'double *', dVyy_aux: 'double *', dVyz_aux: 'double *', dVzx_aux: 'double *', dVzy_aux: 'double *', dVzz_aux: 'double *', LorY: 'int', K_diag: 'double', V_diag: 'double', IorE: 'int', triangle: 'double *', tri: 'int *', k: 'int *', xi: 'double *', yi: 'double *', zi: 'double *', s_xj: 'double *', s_yj: 'double *', s_zj: 'double *', xt: 'double *', yt: 'double *', zt: 'double *', m: 'double *', mx: 'double *', my: 'double *', mz: 'double *', mKc: 'double *', mVc: 'double *', targets: 'int *', Area: 'double *', sglInt_int: 'double *', sglInt_ext: 'double *', xk: 'double *', wk: 'double *', Xsk: 'double *', Wsk: 'double *', kappa: 'double', threshold: 'double', eps: 'double', w0: 'double', aux: 'double *') -> "void":
    return _direct.direct_c_2derivative(dKxx_aux, dKxy_aux, dKxz_aux, dKyx_aux, dKyy_aux, dKyz_aux, dKzx_aux, dKzy_aux, dKzz_aux, dVxx_aux, dVxy_aux, dVxz_aux, dVyx_aux, dVyy_aux, dVyz_aux, dVzx_aux, dVzy_aux, dVzz_aux, LorY, K_diag, V_diag, IorE, triangle, tri, k, xi, yi, zi, s_xj, s_yj, s_zj, xt, yt, zt, m, mx, my, mz, mKc, mVc, targets, Area, sglInt_int, sglInt_ext, xk, wk, Xsk, Wsk, kappa, threshold, eps, w0, aux)
direct_c_2derivative = _direct.direct_c_2derivative

def coulomb_direct(xt: 'double *', yt: 'double *', zt: 'double *', m: 'double *', K_aux: 'double *') -> "void":
    return _direct.coulomb_direct(xt, yt, zt, m, K_aux)
coulomb_direct = _direct.coulomb_direct

def coulomb_energy_multipole( xt: 'double *', yt: 'double *', zt: 'double *', q: 'double *', px: 'double *', py: 'double *', pz: 'double *', px_pol: 'double *', py_pol: 'double *', pz_pol: 'double *', Qxx: 'double *', Qxy: 'double *', Qxz: 'double *', Qyx: 'double *', Qyy: 'double *', Qyz: 'double *', Qzx: 'double *', Qzy: 'double *', Qzz: 'double *', alphaxx: 'double *', thole: 'double *', polar_group: 'int *', connections_12: 'int *', pointer_connections_12: 'int *', connections_13: 'int *', pointer_connections_13: 'int *', p12scale: 'double', p13scale: 'double', K_aux: 'double *') -> "void":
    return _direct.coulomb_energy_multipole(xt, yt, zt, q, px, py, pz, px_pol, py_pol, pz_pol, Qxx, Qxy, Qxz, Qyx, Qyy, Qyz, Qzx, Qzy, Qzz, alphaxx, thole, polar_group, connections_12, pointer_connections_12, connections_13, pointer_connections_13, p12scale, p13scale, K_aux)
coulomb_energy_multipole = _direct.coulomb_energy_multipole

def compute_induced_dipole(xt: 'double *', yt: 'double *', zt: 'double *', q: 'double *', px: 'double *', py: 'double *', pz: 'double *', px_pol: 'double *', py_pol: 'double *', pz_pol: 'double *', Qxx: 'double *', Qxy: 'double *', Qxz: 'double *', Qyx: 'double *', Qyy: 'double *', Qyz: 'double *', Qzx: 'double *', Qzy: 'double *', Qzz: 'double *', alphaxx: 'double *', alphaxy: 'double *', alphaxz: 'double *', alphayx: 'double *', alphayy: 'double *', alphayz: 'double *', alphazx: 'double *', alphazy: 'double *', alphazz: 'double *', thole: 'double *', polar_group: 'int *', connections_12: 'int *', pointer_connections_12: 'int *', connections_13: 'int *', pointer_connections_13: 'int *', dphix_reac: 'double *', dphiy_reac: 'double *', dphiz_reac: 'double *', E: 'double') -> "void":
    return _direct.compute_induced_dipole(xt, yt, zt, q, px, py, pz, px_pol, py_pol, pz_pol, Qxx, Qxy, Qxz, Qyx, Qyy, Qyz, Qzx, Qzy, Qzz, alphaxx, alphaxy, alphaxz, alphayx, alphayy, alphayz, alphazx, alphazy, alphazz, thole, polar_group, connections_12, pointer_connections_12, connections_13, pointer_connections_13, dphix_reac, dphiy_reac, dphiz_reac, E)
compute_induced_dipole = _direct.compute_induced_dipole
# This file is compatible with both classic and new-style classes.


