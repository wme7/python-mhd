#!/usr/bin/env python


if __name__ == "__main__":

    from pylab import cm, matplotlib
    from matplotlib.colors import ListedColormap

    palette_file = open('/Users/jzrake/Work/codes/classics/athena/src/palette.h')

    cmap_name = None
    cmap_adds = { }

    for line in palette_file:

        if line.startswith('static float'):
            cmap_name = line.split()[2].split('[')[0]
            cmap_list = [ ]

        elif cmap_name is not None:

            if line.startswith('}'):
                cmap_adds[cmap_name] = ListedColormap(cmap_list, name=cmap_name)
                cmap_name = None
            else:
                cmap_list.append([float(x) for x in line.split(',')[0:3]])

    from numpy import linspace, newaxis
    from pylab import imshow, colorbar, show, subplot, title, setp

    X = linspace(-1,1,256)
    Y = linspace(-1,1,256)
    test_data = X[:,newaxis] + Y[newaxis,:]

    for n,k in enumerate(cmap_adds):

        ax = subplot(4,2,n+1)
        imshow(test_data, cmap=cmap_adds[k])
        colorbar(shrink=0.8)
        title(k)
        setp(ax.get_xticklabels(), visible=False)
        setp(ax.get_yticklabels(), visible=False)

    show()

    from pickle import dump
    dump(cmap_adds, open('cmaps.pk', 'w'))
