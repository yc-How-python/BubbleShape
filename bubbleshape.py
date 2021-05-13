import ctypes
from time import time

# -*- coding: utf-8 -*-
import multiprocessing as mp
import os
import time
from math import fabs
import numpy as np
import matplotlib.pyplot as plt

from goto import with_goto


def checkBoundary(x, y, z):
    (x1, x2) = x.split(maxsplit=2)

    (y1, y2) = y.split(maxsplit=2)

    (z1, z2) = z.split(maxsplit=2)
    # read boundary
    LX = float(x2) - float(x1)
    LY = float(y2) - float(y1)
    LZ = (float(z2) - float(z1))
    return LX, LY, LZ, x2, x1, y2, y1, z2, z1


def checkPBC(matrix, box):
    [LX, LY, LZ, x2, x1, y2, y1, z2, z1] = box
    listcopy = []
    dO_O = 6
    for atomsarray in matrix:

        dx = fabs(atomsarray[1] - float(x2))
        dy = fabs(atomsarray[2] - float(y1))
        dz = fabs(atomsarray[3] - float(z1))


        if (dx <= dO_O or dy <= dO_O or dz <= dO_O):
            if dx <= dO_O:
                movex = -LX

                k = [atomsarray[0], atomsarray[1] + movex, atomsarray[2], atomsarray[3]]



                listcopy.append( k)
            if dy <= dO_O:
                movey = LY
                k = [atomsarray[0], atomsarray[1], atomsarray[2] + movey, atomsarray[3]]





                listcopy.append( k)
            if dz <= dO_O:
                movez = LZ
                k = [atomsarray[0], atomsarray[1], atomsarray[2], atomsarray[3] + movez]





                listcopy.append( k)
            if dx <= dO_O and dy <= dO_O:
                movex = -LX
                movey = LY
                movez = 0
                k = [atomsarray[0], atomsarray[1] + movex, atomsarray[2] + movey, atomsarray[3] + movez]





                listcopy.append( k)
            if dx <= dO_O and dz <= dO_O:
                movex = -LX
                movez = LZ
                movey = 0
                k = [atomsarray[0], atomsarray[1] + movex, atomsarray[2] + movey, atomsarray[3] + movez]




                listcopy.append( k)
            if dy <= dO_O and dz <= dO_O:
                movex = 0
                movey = LY
                movez = LZ
                k = [atomsarray[0], atomsarray[1] + movex, atomsarray[2] + movey, atomsarray[3] + movez]





                listcopy.append( k)
            if dx <= dO_O and dy <= dO_O and dz <= dO_O:
                movex = -LX
                movey = LY
                movez = LZ
                k = [atomsarray[0], atomsarray[1] + movex, atomsarray[2] + movey, atomsarray[3] + movez]





                listcopy.append( k)


        else:
            pass

    return np.array(listcopy)




def data_extractor(line):
    (aid, atype, ax, ay, az) = line.split(maxsplit=4)
    return [float(aid), float(ax), float(ay), float(az)]


def Loading_atom(targetO,targetC, datafile):
    listW = []
    listC=[]
    while 1:
        line = datafile.readline()
        if line == 'ITEM: NUMBER OF ATOMS\n':
            Noatom_str = datafile.readline()
            atomnumber = int(Noatom_str)

        elif line == 'ITEM: TIMESTEP\n':
            title = datafile.readline()
        elif line == 'ITEM: BOX BOUNDS pp pp pp\n':
            (LX, LY, LZ, x2, x1, y2, y1, z2, z1) = checkBoundary(datafile.readline(), datafile.readline(),
                                                                 datafile.readline())
        elif line == 'ITEM: ATOMS id type xu yu zu\n':  # atomdata
            roundnum = 1
            while roundnum <= atomnumber:
                atomstr = datafile.readline()
                if atomstr == '':
                    break
                else:

                    (aid, atype, ax, ay, az) = atomstr.split(maxsplit=4)
                if atype == targetO:
                    O = [float(aid), float(ax), float(ay), float(az)]

                    H1 = data_extractor(datafile.readline())
                    H2 = data_extractor(datafile.readline())

                    listW.append(O)
                    roundnum = roundnum + 3
                    continue
                elif atype == targetC:
                    listC.append([float(aid), float(ax), float(ay), float(az)])
                    roundnum += 1
                    continue
                else:
                    roundnum +=1
                    continue
            # print('Total: ', len(listW), ' Water',len(listC), 'Methane')

            return (np.array(listW),np.array(listC), title.split()[0], [LX, LY, LZ, x2, x1, y2, y1, z2, z1], atomnumber,listC)

            break  # while1

def Radius_cal(metheneGroup,upperBound,lowerBound):
    # print(metheneGroup)
    Radiusmin=99
    for each in range(int(np.abs(float(upperBound) - float(lowerBound)))+1):
        halfPoint =  each + float(lowerBound)
    # print(halfPoint)
        upperGroup= metheneGroup[metheneGroup >= halfPoint]
        lowerGroup = metheneGroup[metheneGroup < halfPoint]
        # print('upperGroup',upperGroup)
        # print('lowerGroup',lowerGroup)
        if len(upperGroup)>0 :
            upperRadius =upperGroup.max()-upperGroup.min()
        else:
            upperRadius =0
        if len(lowerGroup) > 0:
            lowerRadius = lowerGroup.max() - lowerGroup.min()
        else:
            lowerRadius = 0
        Radius=upperRadius+lowerRadius
        if Radiusmin >=round(Radius, 2):
            Radiusmin = round(Radius, 2)
        else:
            continue
    return  Radiusmin
def Bubble_process(partial_data):#atomarray, title, box

    listresult = [[],[],[],[]]
    liststep = []
    # print('Initiated F4_proc. ')
    for data in partial_data:
        listcopyC = checkPBC(data[1], data[3])
        listcopyW = checkPBC(data[0], data[3])
        listC_O = search_neighbor(data[0], data[1],listcopyC,listcopyW)
        dissovled_methane=[]
        bubble=[]
        for order in range(len(listC_O)):
            if len(listC_O[order]) >=11:
                dissovled_methane.append(listC_O[order])
            else:

                bubble.append(order)
        [LX, LY, LZ, x2, x1, y2, y1, z2, z1]=data[3]

        Rz=Radius_cal(data[1][bubble][:,3],z2, z1)
        Ry=Radius_cal(data[1][bubble][:,2], y2,y1)
        Rx=Radius_cal(data[1][bubble][:,1], x2, x1)
        listresult[0].append(Rx)
        listresult[1].append(Ry)
        listresult[2].append(Rz)
        listresult[3].append(Rx*Ry*Rz*4/3*3.1415926)
        liststep.append(int(data[2]))
    # print('Finished F4_proc. ')
    return listresult,liststep


def search_neighbor(listW,listC, listcopyC,listcopyW):
    listj = []
    listneighboring = []
    arr, col= listC.shape
    dO_O = 6  # The distance between 2 oxygen atoms which are connected by  a hydrogen bond
    dO_O2 = 6 ** 2

    anglecut = np.radians(30)
    for i in range(arr):  # each Carbon atom's...
        idi = listC[i][0]
        dist = np.square(listC[i][1:] - listW[:, 1:]).sum(axis=1)# ...distance to Oxygen atoms

        for j in np.where(dist <= dO_O2)[0]:
            da = dist[j]

            if j not in listj:
                listj.append(j)

            else:
                continue

        whereCopy = np.where(listcopyC[:,  0] == idi) # ...'s copes'...
        for copyorder in whereCopy[0]:

            dist = np.square(listcopyC[copyorder][1:] - listW[:,  1:]).sum(axis=1)# ... distance to Oxygen atoms

            for order in np.where(dist <= dO_O2)[0]:
                distic = dist[order]



                if order not in listj:
                    listj.append(order)

                else:
                    continue

        for copyorder in whereCopy[0]:

            dist = np.square(listcopyC[copyorder][1:] - listcopyW[:, 1:]).sum(axis=1) #... distance to Oxygen atoms' copes

            for N in np.where(dist <= dO_O2)[0]:
                distic = dist[N]





                wherem = np.where(listW[:, 0] == listcopyW[N][0])
                for order in wherem[0]:
                    # print('add copy\'s neighbor\'s original atom')
                    if order not in listj:
                        listj.append(order)



        dist = np.square(listC[i][1:] - listcopyW[:, 1:]).sum(axis=1)# distance to Oxygen atoms' copes

        for orderC2 in np.where(dist <= dO_O2)[0]:
            distci = dist[orderC2]

            if distci >= 3 * dO_O2:
                print('wrong!/n')


            wherem = np.where(listcopyW[orderC2][0] == listW[:,  0])
            for order in wherem[0]:
                # print('add_neighbor\'s original atom')
                if order not in listj:
                    listj.append(order)




            else:
                continue

        listneighboring.append(listj.copy())
        listj = []

    return listneighboring





def is_suffix_lmp(suffix: str):
    if suffix == 'lammpstrj':
        return True
    return True


#



def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


import argparse
from matplotlib.pyplot import MultipleLocator
from matplotlib import cm
from scipy.interpolate import interp1d

if __name__ == '__main__':

    parser = argparse.ArgumentParser(

        description='Notice:\n' + '\n 1.The code and the folder containing the trajectories to be analyzed should be in the same directory.\n' + ' 2.trajectories must be in lammpstrj format and contain only water molecules.')
    parser.add_argument('-i', type=str, default=r'F:\TrjData\TPE_290K_run2', help="Path of folder containing the trjs to be analysed")
    parser.add_argument('-t', type=str, default='1', help="Symbol of oxygen atom")
    parser.add_argument('-n', type=int, default=11, help="Core")
    args = parser.parse_args()
    # mp.cpu_count()
    try:
        os.mkdir('./Bubble_results_'+os.path.split(args.i)[-1])

    except FileExistsError:
        print('已有该文件夹存在')
    else:
        pass

    # Constant
    foldername = args.i
    t1 = time.time()

    PROJECT_DIR_PATH = os.path.dirname(os.path.abspath(os.path.abspath(__file__)))
    # print('Current Directory', PROJECT_DIR_PATH)
    DIR_PATH = os.path.join(PROJECT_DIR_PATH, foldername)
    files = os.listdir(DIR_PATH)
    # print(files)
    pool = mp.Pool(args.n)
    V0=[]

    # plot
    figtotal=plt.figure(figsize=(3,3),dpi=300)
    axtotal=figtotal.gca()

    plt.rc('font', family='Times New Roman', )

    font = {'family': 'Times New Roman',
            'color': 'black',
            'weight': 'normal',
            'size': 12,
            }

    plt.xlabel(r'θ$_i$ (degree)', fontsize=12, fontdict=font)
    plt.ylabel('Kinetic/Potential Energy (J)', fontsize=12)
    b=0.75
    # plt.title('T = '+filename+' K')
    axtotal.tick_params(labelsize=8, direction='in', width=0.5)
    axtotal.spines['bottom'].set_linewidth(b)
    axtotal.spines['left'].set_linewidth(b)
    axtotal.spines['top'].set_linewidth(b)
    axtotal.spines['right'].set_linewidth(b)
    plt.gca().yaxis.set_ticks_position('left')


    listResults=[]
    listName=[]
    
    for filename in files:
        if float(os.path.splitext(os.path.split(filename)[-1])[0].replace('_','.'))>=0.9:
            continue

        filename = os.path.join(DIR_PATH, filename)


        foldername = args.i

        name, suffix = os.path.splitext(filename)
        name = name.split('/')[-1]

        CORE = args.n
        listdata = []

        listResult = [[],[],[],[]]

        listStep = []

        if is_suffix_lmp(suffix):
            os.system('clear')

            datafile = open(filename, 'r')
            ending = datafile.seek(0, 2)
            datafile.seek(0, 0)
            print('Loading...')
            while 1:
                if datafile.tell() == ending:
                    break
                else:
                    pass
                listdata.append(Loading_atom('1', '4', datafile))
            print('Done!\n\n\n')

            multi_process = []
            totalTasks = len(listdata)
            averTasks = totalTasks // CORE
            print("Number of data segments:", totalTasks)
            print("AverTasks:", averTasks)
            # print('Identifying cages...')

            for i in range(CORE):
                if i == CORE - 1:
                    remain = totalTasks % CORE + 1
                else:
                    remain = 0
                multi_process.append(
                    pool.apply_async(Bubble_process, (listdata[i * averTasks:(i + 1) * averTasks + remain],)))

            for each_process in multi_process:
                # print('Finished')
                if each_process.get() == []:

                    r = []

                    s = []

                else:
                    r, s = each_process.get()

                listResult[0].extend(r[0]) #x
                listResult[1].extend(r[1]) #y
                listResult[2].extend(r[2]) #z
                listResult[3].extend(r[3]) #v
                listStep.extend(s)

            name = os.path.split(name)[-1]
            # fig=plt.figure(figsize=(3,3),dpi=300)
            # ax=plt.gca()
            # listResult[0][2:101]=moving_average(listResult[0])
            # listResult[1][2:101] = moving_average(listResult[1])
            # listResult[2][2:101] = moving_average(listResult[2])
            # listResult[3][2:101] = moving_average(listResult[3])
            # ax.plot(listStep, listResult[0], label='x')
            # ax.plot(listStep, listResult[1], label='y')
            # ax.plot(listStep, listResult[2], label='z')
            listResults.append(listResult)
            listName.append(float(os.path.splitext(os.path.split(filename)[-1])[0].replace('_','.')))
            # axtotal.plot(listStep,listResult[3],lw=1,label=float(os.path.splitext(os.path.split(filename)[-1])[0].replace('_','.')))
            # plt.title(float(os.path.splitext(os.path.split(filename)[-1])[0].replace('_','.')))
            # ax.legend()
            # plt.show()

    listName=np.array(listName)
    nameorder = listName.argsort()
    listName = listName[nameorder]
    listResults = np.array(listResults)[nameorder]


    def draw(listName, listResults):
        fig = plt.figure(figsize=(2, 3), dpi=600)
        ax = fig.gca()
        C1 = cm.get_cmap("Greens", 50)
        C2 = cm.get_cmap("Blues", 50)
        C3 = cm.get_cmap("Purples", 50)
        C4 = cm.get_cmap("Oranges", 50)
        colorbar = np.vstack((C1(np.linspace(0.6, 0.9, 6)), C2(np.linspace(0.6, 0.9, 5)), C3(np.linspace(0.6, 0.9, 5)),
                              C4(np.linspace(0.6, 0.9, 5))))
        plt.rc('font', family='Times New Roman', )
        ax.tick_params(labelsize=8, direction='in', width=0.5)
        ax.spines['bottom'].set_linewidth(b)
        ax.spines['left'].set_linewidth(b)
        ax.spines['top'].set_linewidth(b)
        ax.spines['right'].set_linewidth(b)
        plt.gca().yaxis.set_ticks_position('left')
        plt.xlabel('Time (ns)', fontsize=12)
        plt.ylabel('Ellipticity', fontsize=12)
        markerlist = ['*', '^', '>', 'v', '<', '^', '>', 'v', '<', '^', '>', 'v', '<', '^', '>', 'v', '<', '^', '>',
                      'v',
                      '<']
        for eorder in range(len(listName)):
            x = range(len(listResults[eorder][0]))
            ys = listResults[eorder]
            # ys[0][2:101]=moving_average(ys[0])
            # ys[1][2:101] = moving_average(ys[1])
            # ys[2][2:101] = moving_average(ys[2])
            # ys[3][2:101] = moving_average(ys[3])
            y = ys[2]/np.min([ys[0], ys[1]], axis=0)
            f = interp1d(x, y, kind='cubic')
            if eorder == 0:
                ax.plot(np.arange(0, len(listResults[0][0]) - 1, 0.1), f(np.arange(0, len(listResults[0][0]) - 1, 0.1)),
                        ls='--', label='E = %s V/nm' % str(name), lw=0.7, c='black')
            else:
                ax.plot(np.arange(0, len(listResults[0][0]) - 1, 0.1), f(np.arange(0, len(listResults[0][0]) - 1, 0.1)),
                        c=colorbar[eorder], lw=0.7, markersize=4, marker=markerlist[eorder], markevery=100,
                        label='E = %s V/nm' % str(name))
        plt.tight_layout(pad=0)
        plt.savefig('bubble.tiff',dpi=600)
    draw(listName, listResults)

    import xvgreader.xvgreader as xre
    titleList = []
    X = []
    Y = []
    Z = []
    #
    PROJECT_DIR_PATH = os.path.abspath(os.path.abspath(r'F:\Python\xvgreader\TPE_290K_run2\Mtot'))

    print('当前路径', PROJECT_DIR_PATH)

    print(DIR_PATH)
    all_files = os.listdir(PROJECT_DIR_PATH)
    titleList = []
    X = []
    Y = []
    Z = []
    for eachfile in all_files:
        name, suffix = os.path.splitext(eachfile)
        if suffix != '.xvg':
            continue
        else:
            pass
        print(DIR_PATH + '\\' + eachfile)
        file, legends = xre.read_xvg(os.path.join(PROJECT_DIR_PATH,eachfile), unpack=True)

        groupnumber = len(file)
        X.append(file[0])
        if len(Y) == 0:
            for order in range(groupnumber):
                Y.append([file[order]])
        else:
            for order in range(groupnumber):
                Y[order].append(file[order])
        Z.append(float(name.split('K')[0].replace('_', '.')))
        # Z.append(float(name.split('V')[0].replace('_','.')))#X,sI,sII,Total,nameList
        print(legends)
    Ylist = []
    for each in Y:
        Ylist.append(np.array(each))
    X = np.array(X)
    Z = np.array(Z)
    sortorder = Z.argsort()
    Y=Ylist[3][sortorder]
    def drawsub(xlist, ylist, namelist,Y):
        b = 0
        fig, ax = plt.subplots(5, 1, facecolor='w', figsize=(1, 2), dpi=600, sharex=True, )
        fig2, ax2 = plt.subplots(4, 1, facecolor='w', figsize=(1, 2), dpi=600, sharex=True)
        plt.rc('font', family='Times New Roman', )
        
        C1 = cm.get_cmap("Greens", 50)
        C2 = cm.get_cmap("Blues", 50)
        C3 = cm.get_cmap("Purples", 50)
        C4 = cm.get_cmap("Oranges", 50)
        markerlist = ['*', '^', '>', 'v', '<', '^', '>', 'v', '<', '^', '>', 'v', '<', '^', '>', 'v', '<', '^', '>',
                      'v',
                      '<']
        cmap = np.vstack((C1(np.linspace(0.6, 0.9, 6)), C2(np.linspace(0.6, 0.9, 5)), C3(np.linspace(0.6, 0.9, 5)),
                          C4(np.linspace(0.6, 0.9, 5))))


        # ax.xaxis.set_major_locator(x_major_locator)
        # ax2.xaxis.set_major_locator(x_major_locator)

        for order in list(range(len(ylist))):
            if order == 0:
                # s = fig2.add_subplot(510 + order+1,sharex=True,)
                x = range(len(listResults[order][0]))
                ys = listResults[order]
                y = ys[2] / np.min([ys[0], ys[1]], axis=0)
                axtwin = ax2[order].twinx()
                axtwin.plot(np.array(range(len(Y[order])))/10, Y[order],c='black', lw=0.5, label='E = ' + str(namelist[order]) + 'V/nm')
                axtwin.tick_params(labelsize=4, direction='in', width=0.5)
                axtwin.set_ylim([-600,8000])
                ax2[order].plot(x, y,c='magenta', lw=0.7, label='E = ' + str(namelist[order]) + 'V/nm')
                ax2[order].tick_params(labelsize=4, direction='in', width=0.5)
                ax2[order].set_ylim([0,4])

            elif order <= 8:



                if order < 4:
                    x = range(len(listResults[order][0]))
                    ys = listResults[order]
                    y = ys[2] / np.min([ys[0], ys[1]], axis=0)
                    # s = fig2.add_subplot(510 + order+1,sharex=True)
                    axtwin = ax2[order].twinx()
                    axtwin.plot(np.array(range(len(Y[order])))/10, Y[order], c=cmap[order], lw=0.5, label='E = ' + str(namelist[order]) + 'V/nm')
                    axtwin.tick_params(labelsize=4, direction='in', width=0.5)
                    axtwin.set_ylim([-600,8000])
                    ax2[order].plot(x, y, lw=0.7, label='E = ' + str(namelist[order - 1]) + 'V/nm',
                                    c='magenta',
                                    )
                    ax2[order].tick_params(labelsize=4, direction='in', width=0.5)
                    ax2[order].set_ylim([0, 4])


                else:
                    # s = fig.add_subplot(510 + order-3,sharex=True)
                    # ax2.plot(range(len(y)),y*100,label=str(namelist[order]) + 'V/nm')
                    x = range(len(listResults[order][0]))
                    ys = listResults[order]
                    y = ys[2] / np.min([ys[0], ys[1]], axis=0)
                    axtwin = ax[order - 4].twinx()
                    axtwin.plot(np.array(range(len(Y[order]))) / 10, Y[order], c=cmap[order], lw=0.5,
                                label='E = ' + str(namelist[order]) + 'V/nm')
                    axtwin.tick_params(labelsize=4, direction='in', width=0.5)
                    axtwin.set_ylim([-600, 8000])
                    ax[order - 4].plot(x, y, lw=0.7, label='E = ' + str(namelist[order - 1]) + 'V/nm',
                                      c='magenta')
                    ax[order - 4].tick_params(labelsize=4, direction='in', width=0.5)
                    ax[order - 4].set_ylim([0, 4])


                    # ax[order - 4].annotate(round(radlimit, 1), xy=(radlimit, 0.5), xytext=(radlimit, 0.5),
                    #                        arrowprops=dict(facecolor='magenta', shrink=1, width=0.1), fontsize=4)
        # ax.legend(fontsize=6, )
        # ax2.legend(fontsize=6,)
        # ax.tight_layout(pad=0)
        # ax2.tight_layout(pad=0)

        fig.tight_layout(pad=0)
        fig2.tight_layout(pad=0)
        fig.subplots_adjust(wspace=0, hspace=0)
        fig2.subplots_adjust(wspace=0, hspace=0)

        fig.tight_layout(pad=0)
        fig2.tight_layout(pad=0)

        fig.savefig('alpha1.tiff',dpi=600)
        fig2.savefig('alpha2.tiff', dpi=600)
    drawsub([], listResults, listName,Y)
