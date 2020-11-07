import numpy as np

def myprint(mystr, val):
    print(f"{mystr}:\n")
    print(val)
    print("\n\n")

def py_gradient(f, *varargs, **kwargs):

    f = np.asanyarray(f)
    myprint("f after asanyarray", f)
    N = f.ndim  # number of dimensions


    axes = kwargs.pop('axis', None)
    myprint("Axes", axes)
    if axes is None:
        axes = tuple(range(N))
    else:
        axes = _nx.normalize_axis_tuple(axes, N)

    myprint("Axes after sorting", axes)

    len_axes = len(axes)
    n = len(varargs)
    if n == 0:
        # no spacing argument - use 1 in all axes
        dx = [1.0] * len_axes
    elif n == 1 and np.ndim(varargs[0]) == 0:
        # single scalar for all axes
        dx = varargs * len_axes
    elif n == len_axes:
        # scalar or 1d array for each axis
        dx = list(varargs)
        for i, distances in enumerate(dx):
            if np.ndim(distances) == 0:
                continue
            elif np.ndim(distances) != 1:
                raise ValueError("distances must be either scalars or 1d")
            if len(distances) != f.shape[axes[i]]:
                raise ValueError("when 1d, distances must match "
                                 "the length of the corresponding dimension")
            diffx = np.diff(distances)
            # if distances are constant reduce to the scalar case
            # since it brings a consistent speedup
            if (diffx == diffx[0]).all():
                diffx = diffx[0]
            dx[i] = diffx
    else:
        raise TypeError("invalid number of arguments")

    myprint("dx derived", dx)

    edge_order = kwargs.pop('edge_order', 1)

    myprint("edge_order", edge_order)

    if kwargs:
        raise TypeError('"{}" are not valid keyword arguments.'.format(
                                                  '", "'.join(kwargs.keys())))
    if edge_order > 2:
        raise ValueError("'edge_order' greater than 2 not supported")

    # use central differences on interior and one-sided differences on the
    # endpoints. This preserves second order-accuracy over the full domain.

    outvals = []

    # create slice objects --- initially all are [:, :, ..., :]
    slice1 = [slice(None)]*N
    slice2 = [slice(None)]*N
    slice3 = [slice(None)]*N
    slice4 = [slice(None)]*N

    myprint("slice inits", slice1)


    otype = f.dtype
    if otype.type is np.datetime64:
        # the timedelta dtype with the same unit information
        otype = np.dtype(otype.name.replace('datetime', 'timedelta'))
        # view as timedelta to allow addition
        f = f.view(otype)
    elif otype.type is np.timedelta64:
        pass
    elif np.issubdtype(otype, np.inexact):
        pass
    else:
        # all other types convert to floating point
        otype = np.double

    print("----------------------------")
    print("Loop start:\n\n\n")
    for axis, ax_dx in zip(axes, dx):
        myprint("axis", axis)
        myprint("ax_dx", ax_dx)
        if f.shape[axis] < edge_order + 1:
            raise ValueError(
                "Shape of array too small to calculate a numerical gradient, "
                "at least (edge_order + 1) elements are required.")
        # result allocation
        out = np.empty_like(f, dtype=otype)

        # spacing for the current axis
        uniform_spacing = np.ndim(ax_dx) == 0

        # Numerical differentiation: 2nd order interior
        slice1[axis] = slice(1, -1)
        slice2[axis] = slice(None, -2)
        slice3[axis] = slice(1, -1)
        slice4[axis] = slice(2, None)

        myprint("slice1 after axis assignment", slice1)
        myprint("slice2 after axis assignment", slice2)
        myprint("slice3 after axis assignment", slice3)
        myprint("slice4 after axis assignment", slice4)

        if uniform_spacing:
            out[tuple(slice1)] = (f[tuple(slice4)] - f[tuple(slice2)]) / (2. * ax_dx)
            myprint("out (uniform spacing)", out)
        else:
            dx1 = ax_dx[0:-1]
            dx2 = ax_dx[1:]
            a = -(dx2)/(dx1 * (dx1 + dx2))
            b = (dx2 - dx1) / (dx1 * dx2)
            c = dx1 / (dx2 * (dx1 + dx2))
            # fix the shape for broadcasting
            shape = np.ones(N, dtype=int)
            shape[axis] = -1
            a.shape = b.shape = c.shape = shape
            # 1D equivalent -- out[1:-1] = a * f[:-2] + b * f[1:-1] + c * f[2:]
            out[tuple(slice1)] = a * f[tuple(slice2)] + b * f[tuple(slice3)] + c * f[tuple(slice4)]


        # Numerical differentiation: 1st order edges
        if edge_order == 1:
            slice1[axis] = 0
            slice2[axis] = 1
            slice3[axis] = 0

            myprint("slice1 after edge order", slice1)
            myprint("slice2 after edge order", slice2)
            myprint("slice3 after edge order", slice3)
            myprint("slice4 after edge order", slice1clear)

            dx_0 = ax_dx if uniform_spacing else ax_dx[0]
            myprint("dx_0", dx_0)
            # 1D equivalent -- out[0] = (f[1] - f[0]) / (x[1] - x[0])
            out[tuple(slice1)] = (f[tuple(slice2)] - f[tuple(slice3)]) / dx_0

            myprint("out now", out)

            slice1[axis] = -1
            slice2[axis] = -1
            slice3[axis] = -2

            myprint("slice1 again", slice1)
            myprint("slice2 again", slice2)
            myprint("slice3 again", slice3)

            dx_n = ax_dx if uniform_spacing else ax_dx[-1]
            myprint("dx_n", dx_n)
            # 1D equivalent -- out[-1] = (f[-1] - f[-2]) / (x[-1] - x[-2])
            out[tuple(slice1)] = (f[tuple(slice2)] - f[tuple(slice3)]) / dx_n

            myprint("out now", out)

        # Numerical differentiation: 2nd order edges
        else:
            slice1[axis] = 0
            slice2[axis] = 0
            slice3[axis] = 1
            slice4[axis] = 2
            if uniform_spacing:
                a = -1.5 / ax_dx
                b = 2. / ax_dx
                c = -0.5 / ax_dx
            else:
                dx1 = ax_dx[0]
                dx2 = ax_dx[1]
                a = -(2. * dx1 + dx2)/(dx1 * (dx1 + dx2))
                b = (dx1 + dx2) / (dx1 * dx2)
                c = - dx1 / (dx2 * (dx1 + dx2))
            # 1D equivalent -- out[0] = a * f[0] + b * f[1] + c * f[2]
            out[tuple(slice1)] = a * f[tuple(slice2)] + b * f[tuple(slice3)] + c * f[tuple(slice4)]

            slice1[axis] = -1
            slice2[axis] = -3
            slice3[axis] = -2
            slice4[axis] = -1
            if uniform_spacing:
                a = 0.5 / ax_dx
                b = -2. / ax_dx
                c = 1.5 / ax_dx
            else:
                dx1 = ax_dx[-2]
                dx2 = ax_dx[-1]
                a = (dx2) / (dx1 * (dx1 + dx2))
                b = - (dx2 + dx1) / (dx1 * dx2)
                c = (2. * dx2 + dx1) / (dx2 * (dx1 + dx2))
            # 1D equivalent -- out[-1] = a * f[-3] + b * f[-2] + c * f[-1]
            out[tuple(slice1)] = a * f[tuple(slice2)] + b * f[tuple(slice3)] + c * f[tuple(slice4)]

        outvals.append(out)

        myprint("outvals", outvals)

        # reset the slice object in this dimension to ":"
        slice1[axis] = slice(None)
        slice2[axis] = slice(None)
        slice3[axis] = slice(None)
        slice4[axis] = slice(None)

        myprint("Cleaned up slice1", slice1)
        myprint("Cleaned up slice2", slice2)
        myprint("Cleaned up slice3", slice3)
        myprint("Cleaned up slice4", slice4)

    if len_axes == 1:
        return outvals[0]
    else:
        myprint("outvals returned", outvals)
        return outvals
