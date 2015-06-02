#!/usr/bin/env python


if __name__ == "__main__":
    from matplotlib import pyplot as plt
    import numpy as np
    import json
    
    fpass = open("c2p_pass.json", "r")
    dpass = np.array(json.load(fpass))

    ftime = open("c2p_time.json", "r")
    dtime = np.array(json.load(ftime))

    G = dtime[:,:,0]
    B = dtime[:,:,1]
    D = {
         'noble2dzt' : dtime[:,:,2],
         'duffell3d' : dtime[:,:,3],
         'noble1dw'  : dtime[:,:,4],
         'anton2dzw' : dtime[:,:,5] }

    kwargs = { 'interpolation': 'bilinear',
               'origin': 'bottom',
               }

    for k,v in D.items():
        N = v.shape[0]/5
        plt.figure()
        plt.imshow(v.T, **kwargs)
        plt.colorbar()
        plt.xlabel(r"$\gamma$")
        plt.ylabel(r"$\log_{10} B$")
        plt.title(k)
        plt.xticks(range(len(G[:,0]))[N::N], [r"$10^{%1.1f}$"%np.log10(x) for x in G[N::N,0]])
        plt.yticks(range(len(B[0,:]))[N::N], [r"%2.1f"%y for y in B[0,N::N]])

    plt.show()

