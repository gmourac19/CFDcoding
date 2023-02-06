########################################################################################################################
########################################################################################################################
########################################################################################################################
#***********************************************************************************************************************
#
# WARNING: The file named '2Dmesh.vtk' has to be located in the same folder as this solver ('2Dheatdiffusion_nonOOP.vtk')
# is for proper functioning.
#
# REQUIRED LIBRARIES: matplotlib, meshio and numpy
#
# COMMENTS: this is an attempt of solving a 2D temperature distribuition problem using the FVM discretization. A coarse
# mesh was generated for faster computing.
#
# PROBLEM FORMULATION:
#                         250 °C
#                 .....................
#                 | q_gen = 1000 W/m³ |
#                 |                   |
#       100 °C    |                   |    200 °C
#                 |  k = 100 W/(mK)   |
#                 |                   |
#                 ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
#                         150 °C
#
#***********************************************************************************************************************
########################################################################################################################
########################################################################################################################
########################################################################################################################

import matplotlib.pyplot as plt
import meshio
import matplotlib.pyplot as mpl
import numpy as np

def readmesh(filename):
    return meshio.read(filename)


def nodeinfo(mesh):
    mesh_pointcount = len(mesh.points)  # talvez tenha de somar uma unidade
    mesh_points = mesh.points
    x_pointloc = []
    y_pointloc = []
    nodecount = 0
    for nodecount in range(0, mesh_pointcount):
        x_pointloc.append(mesh_points[nodecount, 0])
        y_pointloc.append(mesh_points[nodecount, 1])
    return nodecount, x_pointloc, y_pointloc


#       print("node {}: (x,y) = ({},{})".format(nodeCount,mesh_points[nodeCount,0],mesh_points[nodeCount,1]))
#       print("node {}: (x,y) = ({},{})".format(nodeCount, x_pointloc[nodeCount], y_pointloc[nodeCount]))

def cellinfo(mesh, x_pointloc, y_pointloc):
    # x_pointloc = x_pointloc
    # y_pointloc = y_pointloc
    mesh_cellcount = len(mesh.cells_dict['quad'])  # talvez tenha de somar uma unidade
    meshcells = mesh.cells_dict['quad']
    cellnode1 = []
    cellnode2 = []
    cellnode3 = []
    cellnode4 = []
    cellcount = 0
    for cellcount in range(0, mesh_cellcount):
        cellnode1.append(meshcells[cellcount][0])
        cellnode2.append(meshcells[cellcount][1])
        cellnode3.append(meshcells[cellcount][2])
        cellnode4.append(meshcells[cellcount][3])

    #        print("cell {} has nodes {} {} {} {}".format(cellCount, cellnodes1[cellCount], cellnodes2[cellCount],
    #                                                     cellnodes3[cellCount], cellnodes4[cellCount]))
    #        #return cellcount, cellnode1, cellnode2, cellnode3, cellnode4
    node1_x = []
    node1_y = []
    node2_x = []
    node2_y = []
    node3_x = []
    node3_y = []
    node4_x = []
    node4_y = []
    for ind_aux in range(0, cellcount + 1):
        node1_index = cellnode1[ind_aux]
        node1_x.append(x_pointloc[node1_index]), node1_y.append(y_pointloc[node1_index])
        node2_index = cellnode2[ind_aux]
        node2_x.append(x_pointloc[node2_index]), node2_y.append(y_pointloc[node2_index])
        node3_index = cellnode3[ind_aux]
        node3_x.append(x_pointloc[node3_index]), node3_y.append(y_pointloc[node3_index])
        node4_index = cellnode4[ind_aux]
        node4_x.append(x_pointloc[node4_index]), node4_y.append(y_pointloc[node4_index])
    #        #return node1_x, node1_y, node2_x, node2_y, node3_x, node3_y, node4_x, node4_y
    centroid_x = []
    centroid_y = []
    node1_xy = []  # node (x,y): x1, y1
    node2_xy = []  # node (x,y): x2, y2
    node3_xy = []  # node (x,y): x3, y3
    node4_xy = []  # node (x,y): x4, y4

    face12 = [0, 0, 0, 0]

    for cell_index in range(0, cellcount + 1):
        centroid_x.append(
            (node1_x[cell_index] + node2_x[cell_index] + node3_x[cell_index] + node4_x[cell_index]) / 4)
        centroid_y.append(
            (node1_y[cell_index] + node2_y[cell_index] + node3_y[cell_index] + node4_y[cell_index]) / 4)
        # print("cell {} has centroid located in (x,y): {} {}".format(cell_index, centroid_x[cell_index],
        # centroid_y[cell_index]))
        face12[0] = node1_x[cell_index]
        face12[1] = node1_y[cell_index]
        face12[2] = node2_x[cell_index]
        face12[3] = node2_y[cell_index]

    return centroid_x, centroid_y, node1_xy, node1_x, node1_y, node2_x, node2_y, node3_x, node3_y, node4_x, node4_y

def connectdots(n1x, n1y, n2x, n2y, n3x, n3y, n4x, n4y, nelem):
    posx = np.zeros([nelem, 5])
    posy = np.zeros([nelem, 5])
    for index in range(0, nelem):
        posx[index, 0] = n1x[index]
        posx[index, 1] = n2x[index]
        posx[index, 2] = n3x[index]
        posx[index, 3] = n4x[index]
        posx[index, 4] = n1x[index]

        posy[index, 0] = n1y[index]
        posy[index, 1] = n2y[index]
        posy[index, 2] = n3y[index]
        posy[index, 3] = n4y[index]
        posy[index, 4] = n1y[index]
    return posx, posy


def definingfaces(posx, posy, nel, lx, ly, cx, cy):
    face = np.zeros([1, 6])  ### usar o vstack
    mymat = np.zeros([1, 6])
    xmat = []
    ymat = []
    nxw = np.zeros([nel, 1])
    nxe = np.zeros([nel, 1])
    nyn = np.zeros([nel, 1])
    nys = np.zeros([nel, 1])
    nodecomponents = np.zeros([nel, 4])
    maxval = []
    minval = []
    Lx = lx
    Ly = ly
    cx = cx
    cy = cy
    neighbour_cells = np.zeros([len(cx), 4])
    for auxind in range(0, nelem):
        for alfa in range(0, 4):
            x1 = posx[auxind][alfa]
            y1 = posy[auxind][alfa]
            x2 = posx[auxind][alfa + 1]
            y2 = posy[auxind][alfa + 1]
            xc = (x1 + x2) * 0.5
            yc = (y1 + y2) * 0.5
            mymat = [x1, y1, x2, y2, xc, yc]
            face = np.vstack((face, mymat))
    face = np.delete(face, obj=0, axis=0)
    facedesc = ["" for x in range(len(face))]

    for i in range(0, len(face)):
        if face[i, 4] == 0:
            facedesc[i] = 'wall_1'

        elif face[i, 5] == Ly:
            facedesc[i] = 'wall_2'

        elif face[i, 4] == Lx:
            facedesc[i] = 'wall_3'

        elif face[i, 5] == 0:
            facedesc[i] = 'wall_4'

        elif face[i, 4] > 0 and face[i, 5] > 0:
            facedesc[i] = 'internal'

        elif face[i, 4] < Lx or face[i, 5] < Ly:
            facedesc[i] = 'internal'

    for t in range(0, nel):
        for p in range(0, nel):
            if cy[t] == cy[p]:
                xx = abs(cx[t] - cx[p])
                if xx != 0:
                    xmat.append(xx)
                else:
                    continue

            if cx[t] == cx[p]:
                yy = abs(cy[t] - cy[p])
                if yy != 0:
                    ymat.append(abs(cy[t] - cy[p]))
                else:
                    continue

        cdx = min(xmat)
        cdy = min(ymat)
        nxw[t] = cx[t] - cdx
        nxe[t] = cx[t] + cdx
        if nxw[t] < 0:
            nxw[t] = -11.000
        if nxe[t] > lx:
            nxe[t] = -33.000

        nyn[t] = cy[t] + cdy
        nys[t] = cy[t] - cdy
        if nyn[t] > ly:
            nyn[t] = -22.000
        if nys[t] < 0:
            nys[t] = -44.000

    for yy in range(0, nel):
        nodecomponents[yy, 0] = nyn[yy]
        nodecomponents[yy, 1] = nys[yy]
        nodecomponents[yy, 2] = nxe[yy]
        nodecomponents[yy, 3] = nxw[yy]

    elensew = np.zeros([nel, 4])
    cx = np.around(cx, 3)
    cy = np.around(cy, 3)
    nodecomponents = np.around(nodecomponents, 3)
    for cc in range(0, nel):
        if nodecomponents[cc, 3] != -11:
            stocx1 = nodecomponents[cc, 3]
            for aa in range(0, nel):
                if cx[aa] == stocx1 and cy[aa] == cy[cc]:
                    elensew[cc, 3] = aa
        if nodecomponents[cc, 3] == -11:
            elensew[cc, 3] = -11

        if nodecomponents[cc, 2] != -33:
            stocx3 = nodecomponents[cc, 2]
            for aa in range(0, nel):
                if cx[aa] == stocx3 and cy[aa] == cy[cc]:
                    elensew[cc, 2] = aa
        if nodecomponents[cc, 2] == -33:
            elensew[cc, 2] = -33

        if nodecomponents[cc, 1] != -44:
            stocy4 = nodecomponents[cc, 1]
            for aa in range(0, nel):
                if cy[aa] == stocy4 and cx[aa] == cx[cc]:
                    elensew[cc, 1] = aa
        if nodecomponents[cc, 1] == -44:
            elensew[cc, 1] = -44

        if nodecomponents[cc, 0] != -22:
            stocy1 = nodecomponents[cc, 0]
            for aa in range(0, nel):
                if cy[aa] == stocy1 and cx[aa] == cx[cc]:
                    elensew[cc, 0] = aa
        if nodecomponents[cc, 0] == -22:
            elensew[cc, 0] = -22
    print(cdx)
    return face, facedesc, elensew, cdx


def solver2D(nel, elensew, dx):
    # Ax = b
    area = 0.1
    k = 100
    DA = k * area / dx
    S_vol = 1000
    SV = S_vol * area * dx
    Twall1 = 100
    Twall2 = 250
    Twall3 = 200
    Twall4 = 150
    A = np.zeros([nel, nel])
    b = np.zeros([nel, 1])
    wall1 = Twall1 * (2 * DA) + SV
    wall2 = Twall2 * (2 * DA) + SV
    wall3 = Twall3 * (2 * DA) + SV
    wall4 = Twall4 * (2 * DA) + SV
    spwall1 = -2 * DA
    spwall2 = -2 * DA
    spwall3 = -2 * DA
    spwall4 = -2 * DA
    spinternal = 0
    internal = SV
    elcoef = np.zeros([nel, 5])  ## an, as, ae, aw, sp
    for cc in range(0, nel):
        if elensew[cc, 3] == -11 and elensew[cc, 0] == -22:
            elcoef[cc, 0] = 0
            elcoef[cc, 1] = -(k * area / dx)
            elcoef[cc, 2] = -(k * area / dx)
            elcoef[cc, 3] = 0
            elcoef[cc, 4] = spwall1 + spwall2
            rr1 = elensew[cc][1]
            rr2 = elensew[cc][2]
            A[cc, int(rr1)] = elcoef[cc, 1]
            A[cc, int(rr2)] = elcoef[cc, 2]
            A[cc, cc] = - elcoef[cc, 0] - elcoef[cc, 1] - elcoef[cc, 2] - elcoef[cc, 3] - elcoef[cc, 4]
            b[cc] = wall1 + wall2

        if elensew[cc, 2] == -33 and elensew[cc, 0] == -22:
            elcoef[cc, 0] = 0
            elcoef[cc, 1] = -(k * area / dx)
            elcoef[cc, 2] = 0
            elcoef[cc, 3] = -(k * area / dx)
            elcoef[cc, 4] = spwall3 + spwall2
            rr1 = elensew[cc][1]
            rr3 = elensew[cc][3]
            A[cc, int(rr1)] = elcoef[cc, 1]
            A[cc, int(rr3)] = elcoef[cc, 3]
            A[cc, cc] = - elcoef[cc, 0] - elcoef[cc, 1] - elcoef[cc, 2] - elcoef[cc, 3] - elcoef[cc, 4]
            b[cc] = wall2 + wall3

        if elensew[cc, 1] == -44 and elensew[cc, 2] == -33:
            elcoef[cc, 0] = -(k * area / dx)
            elcoef[cc, 1] = 0
            elcoef[cc, 2] = 0
            elcoef[cc, 3] = -(k * area / dx)
            elcoef[cc, 4] = spwall4 + spwall3
            rr0 = elensew[cc][0]
            rr3 = elensew[cc][3]
            A[cc, int(rr0)] = elcoef[cc, 0]
            A[cc, int(rr3)] = elcoef[cc, 3]
            A[cc, cc] = - elcoef[cc, 0] - elcoef[cc, 1] - elcoef[cc, 2] - elcoef[cc, 3] - elcoef[cc, 4]
            b[cc] = wall3 + wall4

        if elensew[cc, 1] == -44 and elensew[cc, 3] == -11:
            elcoef[cc, 0] = -(k * area / dx)
            elcoef[cc, 1] = 0
            elcoef[cc, 2] = -(k * area / dx)
            elcoef[cc, 3] = 0
            elcoef[cc, 4] = spwall4 + spwall1
            rr0 = elensew[cc][0]
            rr2 = elensew[cc][2]
            A[cc, int(rr0)] = elcoef[cc, 0]
            A[cc, int(rr2)] = elcoef[cc, 2]
            A[cc, cc] = - elcoef[cc, 0] - elcoef[cc, 1] - elcoef[cc, 2] - elcoef[cc, 3] - elcoef[cc, 4]
            b[cc] = wall1 + wall4

        # elementos dos quatro cantos

        if elensew[cc, 0] == -22 and elensew[cc, 2] != -33:
            if elensew[cc, 0] == -22 and elensew[cc, 3] == -11:
                continue
            elcoef[cc, 0] = 0
            elcoef[cc, 1] = -(k * area / dx)
            elcoef[cc, 2] = -(k * area / dx)
            elcoef[cc, 3] = -(k * area / dx)
            elcoef[cc, 4] = spwall2
            rr1 = elensew[cc][1]
            rr2 = elensew[cc][2]
            rr3 = elensew[cc][3]
            A[cc, int(rr1)] = elcoef[cc, 1]
            A[cc, int(rr2)] = elcoef[cc, 2]
            A[cc, int(rr3)] = elcoef[cc, 3]
            A[cc, cc] = - elcoef[cc, 0] - elcoef[cc, 1] - elcoef[cc, 2] - elcoef[cc, 3] - elcoef[cc, 4]
            b[cc] = wall2

        if elensew[cc, 1] == -44 and elensew[cc, 3] != -11:
            if elensew[cc, 1] == -44 and elensew[cc, 2] == -33:
                continue

            elcoef[cc, 0] = -(k * area / dx)
            elcoef[cc, 1] = 0
            elcoef[cc, 2] = -(k * area / dx)
            elcoef[cc, 3] = -(k * area / dx)
            elcoef[cc, 4] = spwall4
            rr0 = elensew[cc][0]
            rr2 = elensew[cc][2]
            rr3 = elensew[cc][3]
            A[cc, int(rr0)] = elcoef[cc, 0]
            A[cc, int(rr2)] = elcoef[cc, 2]
            A[cc, int(rr3)] = elcoef[cc, 3]
            A[cc, cc] = - elcoef[cc, 0] - elcoef[cc, 1] - elcoef[cc, 2] - elcoef[cc, 3] - elcoef[cc, 4]
            b[cc] = wall4

        if elensew[cc, 2] == -33 and elensew[cc, 0] != -22:
            if elensew[cc, 2] == -33 and elensew[cc, 1] == -44:
                continue

            elcoef[cc, 0] = -(k * area / dx)
            elcoef[cc, 1] = -(k * area / dx)
            elcoef[cc, 2] = 0
            elcoef[cc, 3] = -(k * area / dx)
            elcoef[cc, 4] = spwall3
            rr0 = elensew[cc][0]
            rr1 = elensew[cc][1]
            rr3 = elensew[cc][3]
            A[cc, int(rr0)] = elcoef[cc, 0]
            A[cc, int(rr1)] = elcoef[cc, 1]
            A[cc, int(rr3)] = elcoef[cc, 3]
            A[cc, cc] = - elcoef[cc, 0] - elcoef[cc, 1] - elcoef[cc, 2] - elcoef[cc, 3] - elcoef[cc, 4]
            b[cc] = wall3

        if elensew[cc, 3] == -11 and elensew[cc, 1] != -44:
            if elensew[cc, 3] == -11 and elensew[cc, 0] == -22:
                continue

            elcoef[cc, 0] = -(k * area / dx)
            elcoef[cc, 1] = -(k * area / dx)
            elcoef[cc, 2] = -(k * area / dx)
            elcoef[cc, 3] = 0
            elcoef[cc, 4] = spwall1
            rr0 = elensew[cc][0]
            rr1 = elensew[cc][1]
            rr2 = elensew[cc][2]
            A[cc, int(rr0)] = elcoef[cc, 0]
            A[cc, int(rr1)] = elcoef[cc, 1]
            A[cc, int(rr2)] = elcoef[cc, 2]
            A[cc, cc] = - elcoef[cc, 0] - elcoef[cc, 1] - elcoef[cc, 2] - elcoef[cc, 3] - elcoef[cc, 4]
            b[cc] = wall1

        if elensew[cc, 0] != -22 and elensew[cc, 1] != -44 and elensew[cc, 2] != -33 and elensew[cc, 3] != -11:
            elcoef[cc, 0] = -(k * area / dx)
            elcoef[cc, 1] = -(k * area / dx)
            elcoef[cc, 2] = -(k * area / dx)
            elcoef[cc, 3] = -(k * area / dx)
            elcoef[cc, 4] = spinternal

            rr0 = elensew[cc][0]
            rr1 = elensew[cc][1]
            rr2 = elensew[cc][2]
            rr3 = elensew[cc][3]
            A[cc, int(rr0)] = elcoef[cc, 0]
            A[cc, int(rr1)] = elcoef[cc, 1]
            A[cc, int(rr2)] = elcoef[cc, 2]
            A[cc, int(rr3)] = elcoef[cc, 3]
            A[cc, cc] = - elcoef[cc, 0] - elcoef[cc, 1] - elcoef[cc, 2] - elcoef[cc, 3] - elcoef[cc, 4]
            b[cc] = internal

            ## an, as, ae, aw, sp 2 5 11 12

    # Gauss-Seidel method as presented in Malalasekera & Versteeg
    T1 = np.zeros([nel, nel])
    T2 = np.zeros([nel, nel])
    c = np.zeros([nel, 1])
    # T = np.random.uniform(low=99, high=253, size=[nel, 1])
    T = np.zeros([nel, 1])
    res = 0
    # T = np.zeros([nel, 1])
    for iteration in range(0, 10000):
        T_ant = T
        for i in range(0, nel):
            for j in range(0, nel):
                if i > j:
                    T1[i, j] = -A[i, j] / A[i, i]
                    T2[i, j] = 0
                elif i == j:
                    T1[i, j] = 0
                    T2[i, j] = 0
                elif i < j:
                    T1[i, j] = 0
                    T2[i, j] = -A[i, j] / A[i, i]
                c[i] = b[i] / A[i, i]
        T = np.dot(np.array(T1), np.array(T)) + np.dot(np.array(T2), np.array(T_ant)) + np.array(c)
        res = abs(np.dot(np.array(np.transpose(T)), np.array(T)) - np.dot(np.array(np.transpose(T_ant)), np.array(T_ant)))
        print(res)
        if res < 10 ** -1:
            break

    Tmat = np.zeros([int(nel ** 0.5), int(nel ** 0.5)])
    ee = 0
    for iii in range(0, int(nel ** 0.5)):
        for jjj in range(0, int(nel ** 0.5)):
            Tmat[iii, jjj] = T[ee]
            ee = ee + 1

    return A, b, Tmat, res


# np.set_printoptions(edgeitems=30, linewidth=100000, formatter=dict(float=lambda x: "%.3g" % x))
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)

Lx = 4.0
Ly = 4.0
# print("number of processors: ", mp.cpu_count())

mymesh = readmesh('2Dmesh.vtk')
totalnodes, x_pointpos, y_pointpos = nodeinfo(mymesh)
centr_x, centr_y, cell, n1x, n1y, n2x, n2y, n3x, n3y, n4x, n4y = cellinfo(mymesh, x_pointpos, y_pointpos)
px, py = connectdots(n1x, n1y, n2x, n2y, n3x, n3y, n4x, n4y, len(centr_x))
nlines = len(mymesh.cells_dict['line'])
nelem = len(centr_x)
print(nlines)
faces, facedesc, nsew, deltax = definingfaces(px, py, nelem, Lx, Ly, centr_x, centr_y)
print(nsew)
print('dx = {}'.format(deltax))
A_mat, b_mat, T, res = solver2D(nelem, nsew, deltax)
# print(A_mat)
# print(b_mat)
print(T)


mpl.figure()
mpl.scatter(n1x, n1y, marker='.', color='b')
mpl.scatter(n2x, n2y, marker='.', color='b')
mpl.scatter(n3x, n3y, marker='.', color='b')
mpl.scatter(n4x, n4y, marker='.', color='b')
mpl.scatter(centr_x, centr_y, marker='.', color='r')
daplot = mpl.subplot()
for i, txt in enumerate(centr_x):
    daplot.annotate(i, (centr_x[i], centr_y[i]))  # i+1

for a in range(0, len(centr_x)):
    x12 = [px[a][0], px[a][1]]
    y12 = [py[a][0], py[a][1]]

    x23 = [px[a][1], px[a][2]]
    y23 = [py[a][1], py[a][2]]

    x34 = [px[a][2], px[a][3]]
    y34 = [py[a][2], py[a][3]]

    x41 = [px[a][3], px[a][0]]
    y41 = [py[a][3], py[a][0]]

    mpl.plot(x12, y12, 'b')
    mpl.plot(x23, y23, 'b')
    mpl.plot(x34, y34, 'b')
    mpl.plot(x41, y41, 'b')

mpl.figure()
nodes = int((nelem ** 0.5) + 0)
print(nodes)

# ---------------------------------------------------------------------

teste = np.zeros([nodes + 2, nodes + 2])
for aa in range(0, nodes + 2):
    for bb in range(0, nodes + 2):
        if aa == 0:
            teste[aa, bb] = 250

        if bb == 0:
            teste[aa, bb] = 100

        if bb == nodes + 1:
            teste[aa, bb] = 200

        if aa == nodes + 1:
            teste[aa, bb] = 150

for xx in range(0, nodes):
    for yy in range(0, nodes):
        teste[xx + 1, yy + 1] = T[xx, yy]

# ---------------------------------------------------------------------

nodes = nodes + 2
X = (np.linspace(0, Lx, nodes))
print(X)
Y = np.flip(np.linspace(0, Ly, nodes))
print(Y)

Xmesh, Ymesh = np.meshgrid(X, Y)


mpl.contourf(Xmesh, Ymesh, teste)
mpl.colorbar()
mpl.show()

print(A_mat)