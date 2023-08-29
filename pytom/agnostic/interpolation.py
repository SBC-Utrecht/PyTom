from numba import njit


@njit
def fill_values_real_spline(src, dst, mtx):
    # dont need to pass the shape (dims_src, dims) numba works with numpy
    dims = dst.shape

    for dx in range(dims[0]):
        for dy in range(dims[1]):
            for dz in range(dims[2]):

                x = float(mtx[0, 0] * dx + mtx[0, 1] * dy + mtx[0, 2] * dz + mtx[0, 3])
                y = float(mtx[1, 0] * dx + mtx[1, 1] * dy + mtx[1, 2] * dz + mtx[1, 3])
                z = float(mtx[2, 0] * dx + mtx[2, 1] * dy + mtx[2, 2] * dz + mtx[2, 3])

                dst[dx, dy, dz] = splineInterpolation(src, x, y, z)


@njit
def CUB_INT(x,y,xj,x2,x3):
    return y*(x-x2)/(xj-x2)*(x-x3)/(xj-x3)


@njit
def CSPL_INT(f1,f2,f3,f4,a,b,c,d):
    return 3*(f1*(b-a)+f2*(c-a)+f3*(d-b)+f4*(d-c))


@njit
def CSPL_CALC(c2,c3,D1,D2,off):
    return c2+D1*off+((3)*(c3-c2)-(2)*D1-D2)*off*off+((2)*(c2-c3)+D1+D2)*off*off*off


@njit
def linearInterpolation(data,x,y,z):

    sizex, sizey, sizez = data.shape
    if x < 0. or y < 0. or z < 0. or x > (sizex-1) or y > (sizey-1) or z > (sizez-1):
        return 0.

    floorx = int(x)  # integer position along x
    floory = int(y)  # integer position along y
    floorz = int(z)  # integer position along z
    xoffseth = x - floorx  # floating point offset along x
    yoffseth = y - floory  # floating point offset along y
    zoffseth = z - floorz  # floating point offset along z

    case = (xoffseth > 0.000001)*1 + (yoffseth > 0.000001)*2 + (zoffseth>0.000001)*4
    # print(xoffseth, yoffseth, zoffseth, case)
    if case == 7:
        x00x = data[floorx][floory][floorz] + (data[floorx+1][floory][floorz]-data[floorx][floory][floorz])*xoffseth
        #src += this->stridey
        floory += 1
        x01x = data[floorx][floory][floorz] + (data[floorx+1][floory][floorz]-data[floorx][floory][floorz])*xoffseth
        #src += this->stridez;
        floorz += 1
        x11x = data[floorx][floory][floorz] + (data[floorx+1][floory][floorz]-data[floorx][floory][floorz])*xoffseth
        #src -= this->stridey;
        floory += 1
        x10x = data[floorx][floory][floorz] + (data[floorx+1][floory][floorz]-data[floorx][floory][floorz])*xoffseth #src[0] + (src[1]-src[0])*xoffseth;
        x0yx = x00x + (x01x - x00x) * yoffseth

        a =  x0yx + ((x10x + (x11x - x10x) * yoffseth) - x0yx) * zoffseth
        # print(a)
        return

    elif case == 6:
        x0y0 = data[floorx][floory][floorz] + (data[floorx][floory+1][floorz]-data[floorx][floory][floorz])*yoffseth
        #src += this->stridez;
        floorz += 1
        return x0y0 + (data[floorx][floory][floorz] + (data[floorx][floory+1][floorz]-data[floorx][floory][floorz])*yoffseth - x0y0)*zoffseth

    elif case == 5:
        x00x = data[floorx][floory][floorz] + (data[floorx+1][floory][floorz]-data[floorx][floory][floorz])*xoffseth
        #src += this->stridez;
        floorz += 1
        return x00x + (data[floorx][floory][floorz] + (data[floorx+1][floory][floorz]-data[floorx][floory][floorz])*xoffseth - x00x)*zoffseth
    elif case == 3:
        x00x = data[floorx][floory][floorz] + (data[floorx+1][floory][floorz]-data[floorx][floory][floorz])*xoffseth
        #src += this->stridey;
        floory += 1
        return x00x + (data[floorx][floory][floorz] + (data[floorx+1][floory][floorz]-data[floorx][floory][floorz])*xoffseth - x00x)* yoffseth

    elif case == 4:
        return data[floorx][floory][floorz] + (data[floorx][floory][floorz+1] - data[floorx][floory][floorz])*zoffseth
    elif case == 2:
        return data[floorx][floory][floorz] + (data[floorx][floory+1][floorz] - data[floorx][floory][floorz])*yoffseth
    elif case == 1:
        return data[floorx][floory][floorz] + (data[floorx+1][floory][floorz] - data[floorx][floory][floorz])*xoffseth
    else:
        return data[floorx, floory, floorz]


@njit
def cubicInterpolation(data, x, y, z):
    is3D = data.shape[2] > 1
    sizex, sizey, sizez = data.shape

    # is the current position in the data or outside. return default value if outside
    if is3D and (x < 1. or y < 1. or z < 1. or x > (sizex-2) or y > (sizey-2) or z > (sizez-2)):
        return linearInterpolation(data, x, y, z)

    if not is3D and (x < 1. or y < 1. or x > (sizex-2) or y > (sizey-2)):
        return linearInterpolation(data, x, y, z)

    floorx = int(x)  # integer position along x
    floory = int(y)  # integer position along y
    floorz = int(z)  # integer position along z
    xoffseth = x - floorx  # floating point offset along x
    yoffseth = y - floory  # floating point offset along y
    zoffseth = z - floorz  # floating point offset along z

    if xoffseth < 0.0000001 and yoffseth < 0.00000001 and zoffseth < 0.0000001:
        return data[floorx, floory, floorz]

    px_pl_1 = floorx + 1  # all voxels plus 1 from x -> p3,p6,p9
    px_mi_1 = floorx - 1  # all voxels minus 1 from x -> p1,p4,p7

    py_pl_1 = floory + 1  # all voxels plus 1 from y -> p1,p2,p3
    py_mi_1 = floory - 1  # all voxels minus 1 from x -> p7,p8,p9

    lowerLayerBound = -1
    upperLayerBound = 1
    layerOffset = 1

    layerValues = [0., 0., 0.]

    if not is3D:
        lowerLayerBound = 0
        upperLayerBound = 0
        layerOffset = 0

    #//interpolation values for each layer (z)
    for zIteration in range(lowerLayerBound, upperLayerBound + 1):
        #//current position in memory plus current z layer offset in voxels (of type T)
        #//first will be negative (-1), second 0 (same layer), third is 1, next layer

        #//load the pixel values
        v1 = data[floorx-1][floory-1][floorz+zIteration] #*(src-1-this->stridey + zIteration*this->stridez); //one line up in y direction, one position back in x
        v2 = data[floorx  ][floory-1][floorz+zIteration] #*(src  -this->stridey + zIteration*this->stridez); //one line up in y direction
        v3 = data[floorx+1][floory-1][floorz+zIteration] #*(src+1-this->stridey + zIteration*this->stridez); //one line up in y direction, one position forward in x
        v4 = data[floorx-1][floory  ][floorz+zIteration] #*(src-1 + zIteration*this->stridez); //same line in y
        v5 = data[floorx  ][floory  ][floorz+zIteration] #*(src + zIteration*this->stridez); //...
        v6 = data[floorx+1][floory  ][floorz+zIteration] #*(src+1 + zIteration*this->stridez);
        v7 = data[floorx-1][floory+1][floorz+zIteration] #*(src-1+this->stridey + zIteration*this->stridez);
        v8 = data[floorx  ][floory+1][floorz+zIteration] #*(src  +this->stridey + zIteration*this->stridez);
        v9 = data[floorx+1][floory+1][floorz+zIteration] #*(src+1+this->stridey + zIteration*this->stridez);

        #print("Value %f %f %f %f %f %f %f %f %f \n".format(v1,v2,v3,v4,v5,v6,v7,v8,v9))

        #//interpolate first row 1 2 3
        line1 = CUB_INT(x,v1,px_mi_1,floorx,px_pl_1) #; //px1,px2,px3
        line2 = CUB_INT(x,v2,floorx,px_pl_1,px_mi_1) #; //px2,px3,px1
        line3 = CUB_INT(x,v3,px_pl_1,px_mi_1,floorx) #; //px3,px1,px2
        #//store values into v1
        #//printf("Line 1 %f %f %f\n",line1,line2,line3);
        v1 = line1+line2+line3


        #//same for the next rows
        line1 = CUB_INT(x,v4,px_mi_1,floorx,px_pl_1)
        line2 = CUB_INT(x,v5,floorx,px_pl_1,px_mi_1)
        line3 = CUB_INT(x,v6,px_pl_1,px_mi_1,floorx)
        #//printf("Line 2 %f %f %f\n",line1,line2,line3);
        v2 = line1+line2+line3

        line1 = CUB_INT(x,v7,px_mi_1,floorx,px_pl_1)
        line2 = CUB_INT(x,v8,floorx,px_pl_1,px_mi_1)
        line3 = CUB_INT(x,v9,px_pl_1,px_mi_1,floorx)
        v3 = line1+line2+line3

        #//interpolate col 2 5 8 in y direction
        line1 = CUB_INT(y,v1,py_mi_1,floory,py_pl_1)
        line2 = CUB_INT(y,v2,floory,py_pl_1,py_mi_1)
        line3 = CUB_INT(y,v3,py_pl_1,py_mi_1,floory)
        #//printf("Row 1%f %f %f\n",line1,line2,line3);

        layerValues[zIteration + layerOffset] = line1 + line2 + line3

    #//printf("Layer Values %f %f %f \n",layerValues[0],layerValues[1],layerValues[2]);
    #//printf("FloorZ %d %d %d \n",floorz-1,floorz,floorz +1);

    if (is3D):
        line1 = CUB_INT(z,layerValues[0],(floorz-1),floorz,(floorz+1))
        line2 = CUB_INT(z,layerValues[1],floorz,(floorz+1),(floorz-1))
        line3 = CUB_INT(z,layerValues[2],(floorz+1),(floorz-1),floorz)

    #//printf("Layer 1 %f %f %f\n",line1,line2,line3);

    return line1 + line2 + line3


@njit
def splineInterpolation(data, x, y, z):
    sizex, sizey = data.shape[:2]

    is3D = len(data.shape) == 3

    if is3D:
        sizez = data.shape[2]

    # is the current position in the data or outside. return default value if outside
    if is3D and (x < 2. or y < 2. or z < 2. or x > sizex-3 or y > sizey-3 or z > sizez-3):
        return cubicInterpolation(data, x, y, z)

    if not is3D and (x < 2. or y < 2. or x > (sizex-3) or y > (sizey-3)):
        return cubicInterpolation(data, x, y, z)

    floorx = int(x)  # integer position along x
    floory = int(y)  # integer position along y
    floorz = int(z)  # integer position along z
    xoffseth = x - floorx  # floating point offset along x
    yoffseth = y - floory  # floating point offset along y
    zoffseth = z - floorz  # floating point offset along z

    if xoffseth < 0.0000001 and yoffseth < 0.00000001 and zoffseth < 0.0000001:
        return data[floorx, floory, floorz]

    f1 = -0.1556
    f2 = 0.3111
    f3 = -0.0889
    f4 = 0.0444

    lowerLayerBound = -1
    upperLayerBound = 2
    layerOffset = 1
    layerValues = [0., 0., 0., 0.]

    if not is3D:
        lowerLayerBound = 0
        upperLayerBound = 0
        layerOffset = 0

    # interpolation values for each layer (z)
    for zIteration in range(lowerLayerBound, upperLayerBound +1):
        # load the pixel values
        v1 = data[floorx-1][floory-1][floorz+zIteration] #one line up in y direction, one position back in x
        v2 = data[floorx  ][floory-1][floorz+zIteration] #*(src  -this->stridey + zIteration*this->stridez); #one line up in y direction
        v3 = data[floorx+1][floory-1][floorz+zIteration] #*(src+1-this->stridey + zIteration*this->stridez); #one line up in y direction, one position forward in x
        v4 = data[floorx+2][floory-1][floorz+zIteration] #*(src+2-this->stridey + zIteration*this->stridez);

        v5 = data[floorx-1][floory][floorz+zIteration] #*(src-1 + zIteration*this->stridez); //same line in y
        v6 = data[floorx  ][floory][floorz+zIteration] #*(src + zIteration*this->stridez); //...
        v7 = data[floorx+1][floory][floorz+zIteration] #*(src+1 + zIteration*this->stridez);
        v8 = data[floorx+2][floory][floorz+zIteration] #*(src+2 + zIteration*this->stridez);

        v9  = data[floorx-1][floory+1][floorz+zIteration] #*(src-1+this->stridey + zIteration*this->stridez);
        v10 = data[floorx  ][floory+1][floorz+zIteration] #*(src  +this->stridey + zIteration*this->stridez);
        v11 = data[floorx+1][floory+1][floorz+zIteration] #*(src+1+this->stridey + zIteration*this->stridez);
        v12 = data[floorx+2][floory+1][floorz+zIteration] #*(src+2+this->stridey + zIteration*this->stridez);

        v13 = data[floorx-1][floory+2][floorz+zIteration] #*(src-1+2*this->stridey + zIteration*this->stridez);
        v14 = data[floorx  ][floory+2][floorz+zIteration] #*(src  +2*this->stridey + zIteration*this->stridez);
        v15 = data[floorx+1][floory+2][floorz+zIteration] #*(src+1+2*this->stridey + zIteration*this->stridez);
        v16 = data[floorx+2][floory+2][floorz+zIteration] #*(src+2+2*this->stridey + zIteration*this->stridez);

        #print(f'{v1:.5f} {v2:.5f} {v3:.5f} {v4:.5f} {v5:.2f} {v6:.2f} {v7:.2f} {v8:.2f} {v9:.2f} {v10:.2f} {v11:.2f} {v12:.2f} {v13:.2f} {v14:.2f} {v15:.2f} {v16:.2f}')

        #calculate spline value for line 1 2 3 4 above pixel of interest */
        D1 = CSPL_INT(f1,f2,f3,f4,v1,v2,v3,v4)
        D2 = CSPL_INT(f4,f3,f2,f1,v1,v2,v3,v4)
        line1 = CSPL_CALC(v2,v3,D1,D2,xoffseth)

        #print(D1, D2, line1, xoffseth)

        #calculate spline value for line 5 6 7 8 above pixel of interest */
        D1 = CSPL_INT(f1,f2,f3,f4,v5,v6,v7,v8)
        D2 = CSPL_INT(f4,f3,f2,f1,v5,v6,v7,v8)
        line2 = CSPL_CALC(v6,v7,D1,D2,xoffseth)

        #*calculate spline value for line 9 10 11 12 above pixel of interest */
        D1 = CSPL_INT(f1,f2,f3,f4,v9,v10,v11,v12)
        D2 = CSPL_INT(f4,f3,f2,f1,v9,v10,v11,v12)
        line3 = CSPL_CALC(v10,v11,D1,D2,xoffseth)

        #/*calculate spline value for line 13 14 15 16 above pixel of interest */
        D1 = CSPL_INT(f1,f2,f3,f4,v13,v14,v15,v16)
        D2 = CSPL_INT(f4,f3,f2,f1,v13,v14,v15,v16)
        line4 = CSPL_CALC(v14,v15,D1,D2,xoffseth)

        #/*finaly, calculate spline into y direction and save into value[z]*/
        D1 = CSPL_INT(f1,f2,f3,f4,line1,line2,line3,line4)
        D2 = CSPL_INT(f4,f3,f2,f1,line1,line2,line3,line4)
        layerValues[zIteration + layerOffset] = CSPL_CALC(line2,line3,D1,D2,yoffseth)

    if (is3D):
        #/*calculate spline value for z direction*/
        D1 = CSPL_INT(f1,f2,f3,f4,layerValues[0],layerValues[1],layerValues[2],layerValues[3])
        D2 = CSPL_INT(f4,f3,f2,f1,layerValues[0],layerValues[1],layerValues[2],layerValues[3])
        D1 = CSPL_CALC(layerValues[1],layerValues[2],D1,D2,zoffseth)
        return D1
    else:
        return layerValues[0]
