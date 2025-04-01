#!/usr/bin/env python3
#
# Copyright 2023 C Anthony Risinger
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy,
# modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
# is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
# IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
A15 Phase Visualizer

Overview:
    A15.py dynamically explores 3D spatial partitions. Centered on A15 (beta-tungsten), it generates honeycombs for
    numerical-stability analysis. A15.py freely shapes, layers, scales, and intermixes individual components from the
    Weaire-Phelan honeycomb, the Tetrastix prism, and the A15 phase structure itself (at the centers of both).

Usage:
    python3 A15.py [tool] [options] [params] [options] [...]

Tools:
    figure               : Create images of generated structures.
    interactive          : Interact with generated structures.

Params:
    <path>               : Read configuration from file.
    pyritohedra          : Visualize one-or-more pyritohedra.
    tetradecahedra       : Visualize one-or-more tetradecahedra.

Tool Options:
    -auto=<bool|name>    : Set or enable auto-configuration groups.
    -pop=<bool|exec>     : Open visualization in a new pop-up window.

Param Options:
    -auto=<bool|name>    : (all) Set or enable auto-configuration group.
    -h=<int|ratio>       : (pyritohedra) Adjust height. (default: "1/2")
    -planes=<str>        : (tetradecahedra) Configure planes. (default: "XYZ")
    -prescale=<ratio|int>: (all) Adjust implicit scale. (default: 20|24)

Display Options:
    -axes=<bool>         : Toggle visualization axes.
    -bars=<bool>         : Depict numeric stability as histogram bars.
    -bg=<bool>           : Toggle solid white background.
    -centers=<bool>      : Show per-type markers at shape centers.
    -edges=<bool>        : Show or hide shape edges.
    -faces=<bool>        : Show or hide shape faces.
    -lines=<bool>        : Show or hide shape lines.
    -stix=<bool>         : Show Tetrastix or Weaire-Phelan.
    -title=<bool|str>    : Show or set visualization title.
    -verts=<bool>        : Show or hide shape vertices.
    -views=<int|name>    : Set name or number of visualization views.

Scaling Options:
    -n=<int|ratio>       : Adjust the detail level of visualization.
    -rescale=<int|ratio> : Adjust the visualization scale after other operations.
    -scale=<int|ratio>   : Directly scale the visualization.

Advanced Options:
    -colormap=<name>     : Set colormap for shapes or the visualization. (default: "CMRmap")
    -savefig=<path>      : Save the visualization to a specified path. (default: "savefig.png")
    -without=<option>    : Disable previously-enabled options.

Syntax and Tips:
    - Use shorthand notations: "fig" for "figure", or "tetra" for "tetradecahedra".
    - Configuration supports various input formats: space, newline, or null-separated.
    - Add colon-separated options to tools and params: "fig:auto" or "pyrit:edges:h=G16".

Examples:
    - Random Visualization: python3 A15.py -pop -auto
    - Shape with Properties: python3 A15.py -pop,auto,edges figure -scale=0 pyritohedra -h=G06
    - Multiple Configurations: python3 A15.py -pop fig:n=10/3:edges:faces:centers:stix p++ t+:lines
    - Configuration from files: python3 A15.py - ./fig-intro.png.txt <<< :pop:title:axes:bars:views=4

Notes:
    - For deeper understanding or intricate setups, refer to the main paper or supplementary content.
    - Ensure `numpy` and `matplotlib` are properly installed.

Further Exploration:
    The properties of the A15, Weaire-Phelan, and Tetrastix structures are widely-applicable.
    Users are encouraged to experiment with various configurations and visual outputs.
"""


import sys, os, traceback
from pprint import pprint as print

import collections, functools, itertools, operator, string, numbers, random, math
import numpy, matplotlib.pyplot, matplotlib.ticker, mpl_toolkits.mplot3d, scipy.spatial


R = dict()
R["G16"] = (1+5**0.5)/2
R["G26"] = R["G16"] + 1
R["G06"] = R["G16"] - 1
R["G03"] = 2 - R["G16"]
R.update((f"-{G}", -R[G]) for G in (*R,))

V = {0: ()}
V["M"] = ((45, -90),)
V["K"] = ((45, -45),)
V["XZ"] = ((0, -90),)
V["XY"] = ((90, -90),)
V["YZ"] = ((0, -180),)
V[120] = ((30, -210),)
V["G"] = ((numpy.degrees(R["G06"]), -90),)
V["G45"] = ((numpy.degrees(R["G06"]), -45),)
V["G90"] = ((numpy.degrees(R["G06"]), -180),)
V[345] = V["XY"] + V["XZ"] + V["G45"]
V[9] = V["XY"] + V["XZ"] + V["YZ"] + V["G"] + V["G45"] + V["G90"] + V["M"] + V["K"] + V[120]
V[4] = V["XY"] + V["XZ"] + V["G45"] + V["YZ"]
V[3] = V["XY"] + V["XZ"] + V["YZ"]
V[2] = V["XY"] + V["XZ"]
V[1] = V["G45"]


def tetradecahedra(planes="XYZ", prescale=None, stix=False, at=lambda xyz: True, o=(0, 0, 0), **kwds):
    prescale = prescale if prescale is not None else 24 if stix else 20
    o1 = numpy.array((
        (  2,  0, -1),
        ( -1,  2,  0),
        (  0, -1,  2),
    ))
    o2 = numpy.roll(o1.T, 1, axis=1)
    at1 = lambda xyz, at=at: at(xyz) and (xyz[0]%4, xyz[1]%4, xyz[2]%4) == (0, 0, 0)
    at2 = lambda xyz, at=at: at(xyz) and (xyz[0]%4, xyz[1]%4, xyz[2]%4) == (2, 2, 2)
    px = pyritohedron(h=0, prescale=prescale, cube=True) if stix else numpy.array((
        ( -6,  33,  -0), (  4,  28,  20), (  4,  28, -20), (-24,  24,  -0), ( -4,  20,  28), ( -4,  20, -28),
        ( 24,  18,  15), ( 24,  18, -15), (-24,  15,  18), (-24,  15, -18), (  6,  -0,  33), ( 24,  -0,  24),
        ( 24,  -0, -24), (  6,  -0, -33), (-24, -15,  18), (-24, -15, -18), ( 24, -18,  15), ( 24, -18, -15),
        ( -4, -20,  28), ( -4, -20, -28), (-24, -24,  -0), (  4, -28,  20), (  4, -28, -20), ( -6, -33,  -0),
    ))*prescale/20
    py, pz = numpy.roll(px, 1, axis=1), numpy.roll(px, 2, axis=1)
    planes = planes.replace("X", "xxx").replace("Y", "yyy").replace("Z", "zzz")
    planes.count("x")&1 and (yield from lattice(o=o1[1], fn=lambda xyz:  px, at=at1, **kwds))
    planes.count("x")&2 and (yield from lattice(o=o2[1], fn=lambda xyz: -px, at=at2, **kwds))
    planes.count("y")&1 and (yield from lattice(o=o1[2], fn=lambda xyz:  py, at=at1, **kwds))
    planes.count("y")&2 and (yield from lattice(o=o2[2], fn=lambda xyz: -py, at=at2, **kwds))
    planes.count("z")&1 and (yield from lattice(o=o1[0], fn=lambda xyz:  pz, at=at1, **kwds))
    planes.count("z")&2 and (yield from lattice(o=o2[0], fn=lambda xyz: -pz, at=at2, **kwds))
    if kwds.get("centers") or kwds.get("lines"):
        yield from lattice(fn=lambda xyz: o1*24, at=at1, **kwds)
        yield from lattice(fn=lambda xyz: o2*24, at=at2, **kwds)


def pyritohedra(h=None, prescale=None, stix=None, cube=None, at=lambda xyz: True, **kwds):
    cube = cube if cube is not None else stix
    h = h if h is not None else 0 if stix else 1/2
    r = numpy.array(((1, 0, 0), (0, 0, -1), (0, 1, 0)))
    p = pyritohedron(h=h, prescale=prescale, cube=cube, center=kwds.get("centers"))
    fn = lambda xyz: p if xyz[0] % 4 == 0 else numpy.dot(p, r.T)
    at = lambda xyz, at=at: at(xyz) and (xyz[0]%4, xyz[1]%4, xyz[2]%4) in ((0, 0, 0), (2, 2, 2))
    return lattice(fn=fn, at=at, **kwds)


def pyritohedron(h, prescale=None, cube=None, center=None):
    prescale = prescale if prescale is not None else 20 if h else 24
    h = R[h] if h in R else operator.truediv(*h) if isinstance(h, tuple) else h
    p = numpy.unique(numpy.vstack(tuple(itertools.chain(
        () if h != 1/2 and center is None or center is False else ((0, 0, 0),),
        () if h < 0 and cube is None or cube is False else itertools.product((1, -1), repeat=3),
        () if h == 0 else (
            (0,  (1+h),  (1-h**2)), (0,  (1+h), -(1-h**2)), (0, -(1+h),  (1-h**2)), (0, -(1+h), -(1-h**2)),
            ( (1+h),  (1-h**2), 0), ( (1+h), -(1-h**2), 0), (-(1+h),  (1-h**2), 0), (-(1+h), -(1-h**2), 0),
            ( (1-h**2), 0,  (1+h)), (-(1-h**2), 0,  (1+h)), ( (1-h**2), 0, -(1+h)), (-(1-h**2), 0, -(1+h)),
        ),
    )))*prescale, axis=0)
    return p


def lattice(n=0, o=(0, 0, 0), at=lambda xyz: True, fn=lambda xyz: [(0, 0, 0)], **kwds):
    n, gen = -1, n
    if isinstance(gen, tuple):
        gen = operator.truediv(*gen) if len(gen)==2 else (gen,) if len(gen)==3 else numpy.reshape(gen, (len(gen)//3, 3))
    if isinstance(gen, numbers.Real):
        n = math.ceil(gen*2)
        if not isinstance(gen, numbers.Integral):
            r, at = gen*2, lambda xyz, at=at: at(xyz) and scipy.spatial.distance.pdist([(0, 0, 0), xyz]) <= r
        gen = itertools.product(range(-n, n+1), repeat=3)
    for xyz in filter(at, gen):
        xyz = numpy.array(xyz)
        p = numpy.array(fn(xyz))
        yield p + ((xyz + o) * 24), configuration(**kwds)


def figure(shape, *shapes, savefig:[bool|str]=False, views:tuple[tuple[int|float]]=V[1], scale:tuple[int]=(1, 96),
           title:[bool|str]=False, bars:bool=False, axes:bool=False, bg:bool=False):
    scale = tuple(map(int, scale))
    width = 96 * scale[0] / scale[1]
    base1 = collections.defaultdict(set)
    base2 = collections.defaultdict(set)
    rgb2d = collections.defaultdict(list)
    supheaders = collections.defaultdict(int)
    supheader = collections.namedtuple("header", "volume vertices count edges config scale")
    mathtex = collections.namedtuple("mathtex", "approx cdot mathbf epsilon frac dfrac delta Delta")
    normalize = matplotlib.colors.Normalize(vmin=0)
    colormap = matplotlib.colormaps["CMRmap"]
    figsize = (
        (10/4, 10) if not views else
        (10, 10-bars) if (len(views)**0.5).is_integer() else
        (10/len(views), 10) if bars else
        (10, 10/len(views))
    )
    figrats = (
        operator.truediv(*figsize),
        operator.truediv(*figsize[::-1]),
    )
    gridspec_kw = dict(
        left=0, right=1, bottom=0, top=1, hspace=0, wspace=0,
    )
    per_subplot_kw = dict(filter(None, (
        ((*map(str, range(0, len(views))),), dict(projection="3d")),
        (("B",), dict(projection="rectilinear")) if bars else None,
    )))
    mosaic = numpy.array(
        [[*"012"], [*"345"], [*"678"]] if len(views) == 9 else
        [[*"01"], [*"23"]] if len(views) == 4 else
        [*([str(i)] for i in range(len(views)))] if bars else
        [[str(i) for i in range(len(views))]]
    )
    mosaic_kwds = dict(dpi=300, mosaic=mosaic, figsize=figsize, gridspec_kw=gridspec_kw, per_subplot_kw=per_subplot_kw)
    if bars:
        gridspec_kw.update(width_ratios=[figrats[1]*2**-3, *[1/mosaic.shape[1]]*mosaic.shape[1]] if views else [1])
        mosaic = mosaic_kwds["mosaic"] = numpy.insert(mosaic, 0, ["B"]*mosaic.shape[0], axis=1) if views else numpy.array([["B"]])
    fig, axs = (None, dict()) if 0 in mosaic.shape else matplotlib.pyplot.subplot_mosaic(**mosaic_kwds)
    fs = (figsize[0]**2 + figsize[1]**2)**0.5 + min(figsize)/(len(views)+1)
    mtx = mathtex(*(fr"\{field}" for field in mathtex._fields))
    alpha = (3*math.log2(2+len(shapes)))**-1
    bbox_extra_artists = list()
    for ident, c in (shape, *shapes):
        shape = ident * c.rescale[0] * scale[0] / c.rescale[1] / scale[1]
        for vn in numpy.nditer(ident):
            vn = float(vn).as_integer_ratio()
            base1[int(math.log2(vn[1]))].add(abs(vn[0]))
        if c.rescale[1] == 1 and (isinstance(c.rescale[0], numbers.Integral) or c.rescale[0].is_integer()):
            for vn in numpy.nditer(shape):
                vn = float(vn).as_integer_ratio()
                base2[int(math.log2(vn[1]))].add(abs(vn[0]))
        res = operator.truediv(*c.rescale)
        if len(shape) < 4:
            fsz = fs * math.log(1 + res)
            for ax in axs:
                if not ax.isdigit():
                    continue
                if len(shape) != 3:
                    if c.verts:
                        xs, ys, zs = shape[:, 0], shape[:, 1], shape[:, 2]
                        axs[ax].scatter(xs=xs, ys=ys, zs=zs, s=fsz, c="#1a1a1a", alpha=1/4, marker="o")
                elif c.lines:
                    rgb2d[(fsz, "B" if c.centers else "b", shape[0][0], shape[0][1])].append(tuple(shape[0]))
                    rgb2d[(fsz, "R" if c.centers else "r", shape[1][1], shape[1][2])].append(tuple(shape[1]))
                    rgb2d[(fsz, "G" if c.centers else "g", shape[2][0], shape[2][2])].append(tuple(shape[2]))
                elif c.centers:
                    axkwds = dict(ms=fsz/6, linestyle="none", alpha=1/4)
                    axs[ax].plot(xs=shape[0][0], ys=shape[0][1], zs=shape[0][2], c="b", marker="^", **axkwds)
                    axs[ax].plot(xs=shape[1][0], ys=shape[1][1], zs=shape[1][2], c="r", marker=">", **axkwds)
                    axs[ax].plot(xs=shape[2][0], ys=shape[2][1], zs=shape[2][2], c="g", marker="<", **axkwds)
            continue
        hull = scipy.spatial.ConvexHull(shape)
        color = matplotlib.colormaps[c.colormap]((normalize(hull.volume)+0.25)/1.50)
        normals = collections.defaultdict(functools.partial(collections.defaultdict, int))
        for xyz in hull.simplices:
            if c.faces:
                for ax in axs:
                    if not ax.isdigit():
                        continue
                    axs[ax].add_collection3d(mpl_toolkits.mplot3d.art3d.Poly3DCollection(
                        [shape[xyz]], color=color, alpha=alpha, linewidths=0,
                    ))
            curr = shape[xyz[1]] - shape[xyz[2]]
            for l, r in itertools.combinations([0, 1, 2], 2):
                last, curr = curr, shape[xyz[l]] - shape[xyz[r]]
                norm = (numpy.cross(last, curr)/numpy.linalg.norm(numpy.cross(last, curr))).round(4)
                normals[frozenset([tuple(norm), tuple(-norm)])][(xyz[l], xyz[r], round(sum(curr**2)**0.5, 4))] += 1
        dists = set()
        for k in normals:
            es = normals[k]
            for k2 in es:
                if es[k2] != 1:
                    continue
                start, end, dist = k2
                dists.add(dist)
                if not c.edges:
                    continue
                color, zorder, lw = "#040404", 1, 1
                xs, ys, zs = shape[[start, end], 0], shape[[start, end], 1], shape[[start, end], 2]
                for ax in axs:
                    if not ax.isdigit():
                        continue
                    axs[ax].plot(xs=xs, ys=ys, zs=zs, zorder=zorder, color=color, lw=lw, alpha=alpha*3/4)
        if c.centers or c.verts or (c.faces and not c.edges):
            ls = len(shape)
            vertices = shape[hull.vertices]
            xs, ys, zs = vertices[:, 0], vertices[:, 1], vertices[:, 2]
            centered = c.centers and ls % 2 == 1 and ls != len(vertices)
            for ax in axs:
                if not ax.isdigit():
                    continue
                if c.verts:
                    axs[ax].scatter(xs=xs, ys=ys, zs=zs, s=fs*math.log(1+res)/3, marker=".")
                elif c.faces:
                    axs[ax].scatter(xs=xs, ys=ys, zs=zs, s=0, marker=".")
                if centered:
                    axs[ax].scatter(
                        xs=shape[ls//2, 0:1], ys=shape[ls//2, 1:2], zs=shape[ls//2, 2:3], s=fs*math.log(1+res),
                        c="#1a1a1a", alpha=1/4, marker="h" if shape[ls//2, 0:1] % 2 else "H",
                    )
        hscale = (c.rescale[0]*scale[0], c.rescale[1]*scale[1])
        if isinstance(c.rescale[0], numbers.Integral):
            hscale = (hscale[0]//math.gcd(*hscale), hscale[1]//math.gcd(*hscale))
        supheaders[supheader(
            count=-1,
            scale=hscale,
            edges=tuple(sorted(dists)),
            vertices=len(hull.vertices),
            volume=round(hull.volume, 4),
            config=c,
        )] += 1
    for i, ax in axs.items():
        ax.grid(False)
        ax.tick_params(pad=fs/8, labelsize=fs/2)
        ax.margins(*(0,)*len(ax._axis_names))
        try:
            i, view = int(i), views[int(i)]
        except (ValueError, IndexError):
            ax.set_xticklabels([])
            ax.set_xlabel("")
            ax.set_xticks([])
            ax.set_ylabel("")
            ax.set_yticks([])
            continue
        formatter = lambda value, xy: f"${value/width:.0f}$"
        locator = matplotlib.ticker.MultipleLocator(base=width)
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)
        ax.zaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_locator(locator)
        ax.yaxis.set_major_locator(locator)
        ax.zaxis.set_major_locator(locator)
        ax.xaxis.set_pane_color((0, 0, 0, 0))
        ax.yaxis.set_pane_color((0, 0, 0, 0))
        ax.zaxis.set_pane_color((0, 0, 0, 1/30))
        ax.xaxis.line.set_color("r")
        ax.yaxis.line.set_color("g")
        ax.zaxis.line.set_color("b")
        ax.set_proj_type("ortho")
        ax.set_zorder(len(views)-i)
        ax.set_box_aspect(aspect=(1, 1, 1))
        ax.set_xlabel("$N_{X}$", fontweight="bold", labelpad=fs/2, fontsize=fs*0.75, color="r", clip_on=False)
        ax.set_ylabel("$N_{Y}$", fontweight="bold", labelpad=fs/2, fontsize=fs*0.75, color="g", clip_on=False)
        ax.set_zlabel("$N_{Z}$", fontweight="bold", labelpad=fs/2, fontsize=fs*0.75, color="b", clip_on=False)
        bbox_extra_artists.extend((ax.xaxis.label, ax.yaxis.label, ax.zaxis.label))
        if view[0] in (90, 270):
            ax.set_box_aspect(aspect=(1, 1, 1), zoom=1+2**-2)
            ax.yaxis._axinfo['juggled'] = (2, 1, 0)
            ax.set_zticklabels([])
            ax.set_zlabel("")
            ax.set_zticks([])
        elif view[0] in (0, 180) and view[1] in (90, 270):
            ax.set_box_aspect(aspect=(1, 1, 1), zoom=1+2**-2)
            ax.zaxis._axinfo['juggled'] = (1, 2, 0)
            ax.set_yticklabels([])
            ax.set_ylabel("")
            ax.set_yticks([])
        elif view[0] in (0, 180) and view[1] in (0, 180):
            ax.set_box_aspect(aspect=(1, 1, 1), zoom=1+2**-2)
            ax.set_xticklabels([])
            ax.set_xlabel("")
            ax.set_xticks([])
        for axinfo in (ax.xaxis._axinfo, ax.yaxis._axinfo, ax.zaxis._axinfo):
            axinfo["tick"]["inward_factor"], axinfo["tick"]["outward_factor"] = 0.0, 0.0
        if rgb2d:
            for (sz, rgb, _, _), line in rgb2d.items():
                line = numpy.unique(line, axis=0)
                marker = "^" if rgb == "B" else ">" if rgb == "R" else "<" if rgb == "G" else "none"
                ax.plot(xs=line[:, 0], ys=line[:, 1], zs=line[:, 2], linewidth=sz/12, alpha=1/4,
                        marker=marker, ms=sz/6, c=rgb.lower())
        if axes:
            golden = numpy.degrees(R["G06"])
            bbox_extra_artists.append(ax.text2D(
                x=1+2**-8, y=1-2**-5, transform=ax.transAxes, va="top",
                ha="right", fontsize=fs*0.75, color="#1a1a1a",
                s=r"$%s_{azim}\ %s_{elev}\ \mathbf{%s}$" % (
                    r"\Phi-1" if view[1] == golden else fr"{view[1]}^\circ",
                    r"\Phi-1" if view[0] == golden else fr"{view[0]}^\circ",
                    ax.name if ax.name != "3d" else "ortho" if ax._focal_length == numpy.inf else "persp",
                ),
            ))
        else:
            ax.zaxis.set_pane_color((0, 0, 0, 0))
            ax.xaxis.line.set_color("none")
            ax.yaxis.line.set_color("none")
            ax.zaxis.line.set_color("none")
            ax.set_xticklabels([])
            ax.set_xlabel("")
            ax.set_xticks([])
            ax.set_yticklabels([])
            ax.set_ylabel("")
            ax.set_yticks([])
            ax.set_zticklabels([])
            ax.set_zlabel("")
            ax.set_zticks([])
        xyzmm = max(map(abs, itertools.chain(ax.get_xlim(), ax.get_ylim(), ax.get_zlim())))
        ax.set_xlim(-xyzmm, xyzmm)
        ax.set_ylim(-xyzmm, xyzmm)
        ax.set_zlim(-xyzmm, xyzmm)
        ax.view_init(*view)
    if len(base2) != 1 and 0 in base2 and len(base2[0]) == 1 and 0.0 in base2[0]:
        origin = base2.pop(0).pop()
        base2[min(base2)].add(origin)
    base1max, base2max = max((0, *base1)), max((0, *base2))
    base2rat = operator.truediv(*scale).as_integer_ratio()
    base2min = max(int(math.log2(base2rat[1])), base1max)
    base2set = {2**(base2max-b)*n for b in base2 for n in base2[b]}
    base2gcd = 1 if len(base2set) == 1 else math.gcd(*base2set)
    base2mm0 = base2rat[0] - base2gcd*2**(base2min-base2max)
    if bars and base2:
        ax = axs["B"]
        ax.set_zorder(10)
        ax.set_xlim(-1-2**-3, 1+2**-3)
        for spine in ax.spines.values():
            spine.set_visible(False)
        if axes:
            ax.set_title(
                x=0, y=-3*2**-7, loc="left", ma="left", va="top", fontsize=fs*3*2**-2, pad=1,
                label="${}$".format("$\n$".join(list(filter(None, [
                    fr"\mathbf{{N_{{1}}}} = {width:.55}_{{mm}}",
                ] + [
                    r"N_{%s} %s %s_{m} %s %s" % (
                        f"2^{{{math.log2(n):.0f}}}" if math.log2(n).is_integer() else n,
                        "=" if (n*10/254).is_integer() else mtx.approx,
                        mm2m(n*width),
                        "=" if (n*10/254).is_integer() else mtx.approx,
                        fr"{ft}_{{ft}}\,{ins:.0f}_{{in}}" if ins != 0 else
                            fr"{ft}_{{ft}}",
                    )
                    for mm, n, ft, ins in (
                        (mm/1, int(round(mm/width, 0)), int(mm*10//3048), (mm*10-(mm*10//3048*3048))/254)
                        for mm in sorted((2**16*width, 1524, 1524*2, 1524*30, 1524*40))
                    )
                ] + [(" " if base2mm0==0 else "$\n$\\ldots ").join(filter(None, [
                    fr"\mathbf{{\epsilon_{{\Delta}}}} = \epsilon_{{\delta}} - \epsilon_{{N}}",
                    fr"= \frac{{{base2rat[0]} - {base2gcd} \cdot 2^{{{base2min}\!-\!{base2max}}}}}{{2^{{{base2min}}}}}",
                    "= 0" if base2mm0==0 else fr"= \frac{{{base2rat[0]} - {base2gcd * 2**(base2min-base2max)}}}{{2^{{{base2min}}}}}",
                    None if base2mm0==0 else f"{mtx.approx} {(base2rat[0] - base2gcd * 2**(base2min-base2max))/2**base2min:.11}"
                ]))])))),
            )
        stepc, counts, bins, gaps, hist = 4, [], [], [], numpy.histogram(
            a=[k for k in sorted(base2) for _ in base2[k]],
            bins=sorted({*base2, max(base2min, base2max)+1}),
        )
        for i, (c, b) in enumerate(itertools.zip_longest(*hist, fillvalue=0)):
            cgap = 0 if i == 0 or bins[len(bins)-1] == b - 1 else b - 1 - bins[len(bins)-1]
            if cgap != 0:
                counts.append(0), bins.append(b-1), gaps.append(cgap)
            if b != max(base2min, base2max)+1:
                counts.append(int(c)), bins.append(b), gaps.append(1)
        cnil, cnix = sum(count==0 for count in counts), max((1, *gaps))
        csum, cmax = sum(counts), max(counts)
        ax.set_title(
            x=0.5, y=1+2**-6, loc="center", ma="left", va="bottom", fontsize=fs*5*2**-3, fontweight="bold",
            label="\n".join(filter(None, [
                f"${csum}$ $binary_{{64}}$ float{'s'*(csum>1)}",
                f"${len(base2)}$ $rational$ epsilon{'s'*(len(base2)>1)}",
            ])),
        )
        cmap = colormap.resampled(len(set(counts)) * stepc + 2*stepc)
        cmap.set_extremes(under=cmap(0), over=cmap(cmap.N-1))
        for i, (count, power, gap) in enumerate(zip(counts, bins, gaps)):
            y = ax.get_ylim()[1]
            w = max(1/8, count/cmax) if count else max(1/16, gap/cnix)
            bc = cmap(stepc + (count % (cmap.N - 2*stepc))) if count else cmap.get_over()
            lc = cmap(stepc + (count + ((cmap.N - 2*stepc) // 2) % (cmap.N - 2*stepc))) if count else cmap.get_under()
            ax.hlines(xmin=-1, xmax=1, y=y+2, colors="#e8e8e8", clip_on=False, linewidth=0.5, capstyle="butt")
            bax = ax.barh(
                left=-1, width=2*w, y=y+1, height=2, color=bc,
                hatch=None if count else "///" if w==1 and cnil>1 else "//",
            )
            count and ax.bar_label(bax, color=lc, label_type="center", fontsize=fs*3*2**-2, fmt=f"${count:.0f}$")
            bbox_extra_artists.append(ax.annotate(
                textcoords="offset points", xytext=(-fs/4, 0), xy=(-1, y+1.625), clip_on=False,
                va="center", ha="right", fontsize=fs*3*2**-2, color="#1a1a1a",
                text=r"$2^{%s}$" % -power,
            ))
    if fig and title:
        nscale = (base2gcd*2**(base2min-base2max)).as_integer_ratio()
        nscale = (nscale[0]//math.gcd(*nscale), (nscale[1]*2**base2min)//math.gcd(*nscale))
        bbox_extra_artists.append(fig.suptitle(
            x=(0.7 if figrats[1]>=3 else 0.55) if bars and views else 0.5,
            y=1+2**-5 if axes and views else 1+3*2**-5 if axes or bars else 1,
            ma="right", va="baseline" if not axes and mosaic.shape[0]==1 else "bottom",
            fontsize=fs*1.125, fontfamily="monospace",
            t="\n".join(list(filter(None, [
                hasattr(title, "title") and title.title(),
                r"$\mathbf{%s}_{\left(\epsilon_{\delta} = %s\right)}$" % (
                    "integer" if base2rat[1] == 1 else
                        "binary" if scale[0] == 1 and scale == base2rat else
                        "stable" if scale == base2rat else
                        "unstable",
                    scale[0] if scale[1] == 1 else
                        fr"2^{{{-base2min}}} = \epsilon_{{N}}" if scale[0] == 1 and scale == base2rat else
                        fr"\frac{{{scale[0]}}}{{2^{{{base2min}}}}} = \epsilon_{{N}}" if scale == base2rat else
                        fr"\frac{{{scale[0]}}}{{{scale[1]}}} {mtx.approx} \frac{{{base2rat[0]}}}{{2^{{{base2min}}}}}"
                            fr"\right) > \left(\frac{{{nscale[0]}}}{{2^{{{math.log2(nscale[1]):.0f}}}}} = \epsilon_{{N}}",
                ),
            ] + [pretty(h) for h in sorted((h._replace(count=c) for h, c in supheaders.items()))] + [
            ]))),
        ))
    if fig and savefig:
        metadata = {"Title": hasattr(title, "title") and title.title() or None,
                    "Software": "https://github.com/infimalabs/space/",
                    "Artist": "https://infima.space/"}
        fig.savefig(bbox_extra_artists=bbox_extra_artists, bbox_inches="tight", transparent=not bg, pad_inches=0,
                    fname=savefig, metadata=metadata)


def mm2m(mm:float):
    # We can't divide by 1000 (binary float) and this needs way more than the
    # default decimal.Decimal precision (28); visually move the decimal point
    # and :heart: metric instead.
    l, _, r = f"{float(abs(mm)):.111}".partition(".")
    l = l if len(l) >= 4 else ("0"*(4-len(l)) + l)
    return f"{l[:-3]}.{l[-3:]}{r}".rstrip("0")


def pretty(h):
    desc = None
    power = math.log2(h.scale[1])
    denom = fr"2^{{{power:.0f}}}" if power.is_integer() else h.scale[1]
    epsilon = (
        h.scale[0] if h.scale[1]==1 else
        fr"2^{{{-power:.0f}}}" if h.scale[0]==1 and power.is_integer() else
        fr"\frac{{{h.scale[0]}}}{{{denom}}}"
    )
    if h.vertices == 24 and len(h.edges) == 4:
        desc = fr"tetradecahedra_{{\left(\epsilon={epsilon}\right)}}"
    if h.vertices == 20:
        desc = "pyritohedra"
        if len(h.edges) == 1:
            desc = fr"dodecahedra_{{\left(\epsilon={epsilon}\right)}}"
        elif round(h.edges[0]/h.edges[1], 2) == 0.76:
            desc += fr"_{{\left(\epsilon={epsilon}\right)}}"
        else:
            desc += fr"_{{\left(edges={h.edges}, \epsilon={epsilon}\right)}}"
    if h.vertices == 12:
        desc = "pseudoicosahedra"
        if len(h.edges) == 1:
            desc = fr"icosahedra_{{\left(\epsilon={epsilon}\right)}}"
        elif round(h.edges[0]/h.edges[1], 2) == 0.65:
            desc += fr"_{{\left(h=\frac{{7}}{{5}}, \epsilon={epsilon}\right)}}"
        else:
            desc += fr"_{{\left(edges={h.edges}, \epsilon={epsilon}\right)}}"
    if h.vertices == 8:
        desc = "hexahedra"
        if len(h.edges) == 1:
            desc = fr"cube_{{\left(\epsilon={epsilon}\right)}}"
        else:
            desc += fr"_{{\left(edges={h.edges}, \epsilon={epsilon}\right)}}"
    if desc is None:
        desc = fr"polyhedra_{{\left(edges={h.edges}, vertices={h.vertices}, \epsilon={epsilon}\right)}}"
    return fr"${h.count} \cdot {desc}$"


def configuration(*args, **kwds):
    make, args = (args[0], args[1:]) if args and callable(args[0]) else (configuration._make, args)
    c = (auto, axes, bars, bg, centers, colormap, edges, faces, lines, n,
     pop, rescale, savefig, scale, stix, title, verts, views) = make()(**kwds)
    colormap = "CMRmap" if colormap is None else colormap
    auto = False if auto is None else auto
    factor = random.randint(1, 4) if auto else 1
    scale = random.choice((random.randint(-7, 22), (factor, 100), (factor, 96))) if scale is None else scale
    if isinstance(scale, numbers.Integral) and scale >= 0:
        scale = (factor*(2**scale//96+1)/2**scale).as_integer_ratio()
    elif isinstance(scale, tuple) and scale[1] < 0:
        scale = (scale[0]*2**scale[1]).as_integer_ratio()
    elif isinstance(scale, numbers.Integral):
        scale = (factor*2**scale).as_integer_ratio()
    elif isinstance(scale, numbers.Real):
        scale = scale.as_integer_ratio()
    scale = (scale[0]//math.gcd(*scale), scale[1]//math.gcd(*scale))
    astype = float
    if rescale is None:
        astype = int
        rescale = (1, 1)
    elif isinstance(rescale, numbers.Real):
        astype = int if isinstance(rescale, numbers.Integral) else float
        rescale = (1, 1) if rescale in {-1, 1, -1.0, 1.0} else rescale.as_integer_ratio()
    elif isinstance(rescale, tuple) and isinstance(rescale[0], numbers.Integral):
        astype = int if operator.truediv(*rescale).as_integer_ratio()[1] == 1 else float
    rescale = (
        (astype(1), 1) if 0 in rescale else
        (astype(abs(rescale[0])), abs(rescale[1])) if operator.truediv(*rescale) > 0 else
        (astype(abs(rescale[1])), abs(rescale[0]))
    )
    shuffle = auto and views is None
    views = [
        (float(elev % 360), float(azim % 360))
        for elev, azim in (
            V[random.choice((0, *(random.randint(1, 4),)*9))] if views is None else
            V[views] if views.__hash__ and views in V else
            views
        )
    ]
    shuffle and random.shuffle(views)
    views = (*views,)
    n = (
        auto and factor*10/9 or 1 if n is None else
        n if isinstance(n, numbers.Real) else
        n[0] if n[1] == 1 else
        operator.truediv(*n)
    )
    savefig = (
        "/".join(p or "" for p in savefig) if isinstance(savefig, tuple) else
        os.environ.get("SAVEFIG", "savefig.png") if savefig is None else
        savefig
    )
    pop = (
        () if pop is None else
        (*pop, savefig) if isinstance(pop, tuple) else
        ("open" if pop is True else pop, savefig)
    )
    axes = auto and random.random() < 0.9 if axes is None else axes
    bars = auto and (True if not views else random.random()<0.9) if bars is None else bars
    title = auto and (axes or bars or len(views) > 1) if title is None else title
    centers = auto and random.choice((True, getattr(n, "real", 0)>2)) if centers is None else centers
    edges = auto and random.choice((True, False, centers)) if edges is None else edges
    lines = auto and random.choice((True, False, not edges)) if lines is None else lines
    faces = auto and random.choice((True, False, edges, lines)) if faces is None else faces
    verts = auto and random.choice((True, False, not edges, not lines, not faces)) if verts is None else verts
    stix = auto and random.choice((True, False, edges, not lines, not faces, verts)) if stix is None else stix
    bg = auto and True if bg is None else bg
    ns, configs = locals().copy(), list()
    ks = {k for k in kwds if ns[k] is not None}
    for k in args:
        if isinstance(k, tuple):
            configs.append(make(k)(*map(ns.pop, k)))
            ks.difference_update(k)
            continue
        configs.append(ns.pop(k))
        ks.discard(k)
    key = (*sorted(ks),)
    konfig = make(key)(*map(ns.pop, key))
    return (*configs, konfig) if configs else konfig

configuration._field_defaults = dict.fromkeys((*sorted((
    "auto", "axes", "bars", "bg", "centers", "colormap", "edges", "faces", "lines", "n",
    "pop", "rescale", "savefig", "scale", "stix", "title", "verts", "views",
)),))

configuration._fields = (*sorted(configuration._field_defaults),)

configuration._make = (
    lambda t=configuration._fields, *, name="configuration", defaults=configuration._field_defaults, cache=dict(): (
        cache[t] if t in cache else cache.setdefault(t, collections.namedtuple(name, t, defaults=map(defaults.get, t)))
    )
)


def flags(name, *args, **kwds):
    if name is not None:
        name, *more = name.split(":")
        rescale = kwds.get("rescale", (1, 1))
        rescale = rescale if isinstance(rescale, tuple) else (rescale, 1)
        rescale = rescale[0]*2**(name.count("+")-name.count("-")), rescale[1]
        name = name.lower().translate(str.maketrans(dict.fromkeys("+-_", "")))
        n = name.translate(str.maketrans(dict.fromkeys(string.ascii_lowercase, "")))
        name = name.translate(str.maketrans(dict.fromkeys(string.digits+".,/", "")))
        args = ("-rescale=%s/%s" % rescale, *args) if rescale != (1, 1) else args
        args = (f"-n={n}", *args) if n else args
        args = (*(f"-{m}" for m in more if m), *args)
    yays = {"1", "true", "yes", "enable", "with"}
    nays = {"0", "false", "no", "disable", "without"}
    while args:
        arg = "/dev/stdin" if args[0] == "-" else args[0]
        if not arg.strip("-:./"):
            args = () if arg else args
            break
        if arg[0] == ":":
            more = arg.replace(":", ":-")[1:].split(":")
            args = (*more, *args[1:])
            continue
        if arg[0] in (".", "/"):
            txt = sys.stdin.read() if arg=="/dev/stdin" else open(arg).read()
            sep = "\0" if "\0" in txt else "\n" if txt.count("\n") > 1 else " "
            more = txt.strip(string.whitespace if sep == " " else sep).split(sep)
            args = (*more, *args[1:])
            continue
        if arg[0] != "-":
            break
        arg, args = args[0], args[1:]
        keys, eq, vals = arg.partition("=")
        vals = vals.split("/" if "/" in vals else ",")
        *nay, keys = keys.lower().lstrip("-").split("-", 1)
        nay = bool(nays & set(nay))
        for i, v in enumerate(vals):
            if not v:
                vals[i] = None if eq else not nay
                continue
            try:
                vals[i] = int(v) if v.removeprefix("-").isdecimal() else float(v)
            except ValueError:
                v = v.lower()
            else:
                continue
            if v in yays:
                vals[i] = not nay
            elif v in nays:
                vals[i] = nay
        for key in keys.split("/" if "/" in keys else ","):
            if len(key) > 3 and key in yays:
                kwds.update(dict.fromkeys(vals, True))
            elif len(key) > 3 and key in nays:
                kwds.update(dict.fromkeys(vals, False))
            else:
                kwds[key] = vals[0] if len(vals)==1 else (*vals,)
    kwds = {k: v for k, v in kwds.items() if v is not None}
    return (
        (name, args, kwds) if name is not None or not args else
        flags(args[0], *args[1:], **kwds)
    )


def interactive(*args, **kwds):
    figure(*args, **kwds)
    with matplotlib.pyplot.ion():
        matplotlib.pyplot.show(block=True)


def help(*args, **kwds):
    sys.stdout.write(f"---{__doc__}---\n")
    sys.exception() and traceback.print_exception(sys.exception())


def main(*args, **kwds):
    tool, args, kwds = (
        flags("figure", "-auto", "-pop", **kwds) if not args else
        flags(None, *args, **kwds) if args[0] and args[0][0] in "-:./" else
        flags(*args)
    )
    toolc, shapec, pop, more = configuration(
        ("scale", "views", "bars", "title", "axes", "savefig", "bg"),
        ("auto", "rescale", "n"),
        "pop", **kwds,
    )
    args = args or (*filter(None, (
        # Largest first.
        "pyritohedra+" if random.random() < 3/4 else None,
        "pyritohedra",
        "tetradecahedra+" if random.random() < 3/4 else None,
        "tetradecahedra",
    )),)
    shapers = {
        "pyritohedra": pyritohedra,
        "tetradecahedra": tetradecahedra,
    }
    tooling = {
        None: figure,
        "help": help,
        "figure": figure,
        "interactive": interactive,
    }
    shapes = list()
    configs = dict()
    while args:
        shape, args, kwds = flags(args[0], *args[1:], **shapec._asdict(), **more._asdict())
        if shape not in shapers:
            shape = sorted(
                (len(os.path.commonprefix((shape, s))), s.startswith(shape) and s or shape) for s in filter(None, shapers))[-1][1]
        for ident, config in shapers[shape](**kwds):
            k = shape, config._fields, config
            c = configs[k] = configs[k] if k in configs else configuration(
                **dict(dict.fromkeys(("colormap", "edges", "faces", "verts", "centers", "lines")), **config._asdict()))
            shapes.append((ident, c))
    if tool not in tooling:
        tool = sorted(
            (len(os.path.commonprefix((tool, t))), t.startswith(tool) and t or tool) for t in filter(None, tooling))[-1][1]
    tooling[tool](*shapes, **toolc._asdict())
    pop and os.spawnvp(os.P_WAIT, pop[0], pop)


if __name__ == "__main__":
    try: status = main(*sys.argv[1:])
    except: status = help() or 1
    sys.exit(status)
