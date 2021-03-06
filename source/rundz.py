#!/usr/bin/env python

# @author: Bo Xin
# @      Large Synoptic Survey Telescope

# main function

import os
import sys
# import glob
import argparse
import subprocess
import multiprocessing

import numpy as np
from astropy.io import fits
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from cwfsAlgo import cwfsAlgo
from cwfsInstru import cwfsInstru
from cwfsImage import cwfsImage


def main():
    parser = argparse.ArgumentParser(
        description='-----Wavefront sensor separation analysis------')

    parser.add_argument('obsID', help='observation ID number')
    parser.add_argument('-nsnap', dest='nsnap', default=100, type=int,
                        help='Number of Image pairs (snaps)')
    parser.add_argument('-mag', dest='mag', default=-1, type=int,
                        help='Star magnitude')
    parser.add_argument('-stampSize', dest='stampSize', default=-1, type=int,
                        help='Size of Image Stamp')
    parser.add_argument('-phosimoff', help='w/o running Phosim',
                        action='store_true')
    parser.add_argument('-e', dest='eimage', default=0, type=int,
                        help='eimage only? 0: eimage (default), \
                        1, amplifier image')
    parser.add_argument('-cwfsoff', help='w/o running cwfs',
                        action='store_true')
    parser.add_argument('-p', dest='numproc', default=1, type=int,
                        help='Number of Processors Phosim uses')
    parser.add_argument('-d', dest='debugLevel', type=int,
                        default=0, choices=(-1, 0, 1),
                        help='debug level, -1=quiet, 0=normal, \
                        1=verbose, default=0')
    args = parser.parse_args()
    if args.debugLevel >= 1:
        print(args)

    dz, instFile, cmdFile, filter, field = createPertFiles(
        args.obsID, args.nsnap, args.mag, args.debugLevel)

    if args.stampSize == -1:
        # stampSize = 2**np.ceil(np.log(dz/1.2335/10e-3+50)/np.log(2))
        if dz<3.0:
            stampSize = np.ceil((dz/1.2335/10e-3+100)/10)*10
        else:
            stampSize = np.ceil((2.0/1.2335/10e-3+100)/10)*10
    else:
        stampSize = args.stampSize
        
    if not args.phosimoff:
        runPhosim(args.obsID, dz, instFile, cmdFile, args.nsnap, filter, field,
                  args.eimage, stampSize, args.numproc, args.debugLevel)

    if dz < 3.0:
        instruFile = 'lsst%02d' % (dz * 10)
    else:
        instruFile = 'lsst20' 
    algoFile = 'exp'
    if not args.cwfsoff:
        parallelCwfs(args.obsID, args.eimage, instruFile, algoFile,
                     stampSize, args.nsnap, field, args.numproc,
                     args.debugLevel)

    plotMeanSTD(args.obsID, args.nsnap, args.debugLevel)


def plotMeanSTD(obsID, nsnap, debugLevel):
    znmax = 22
    zer = np.zeros((znmax - 3, nsnap))

    x = range(4, znmax + 1)
    dz, r0seeing500, vKseeing500, seed, teleState, filter, field, exptime, \
      ccdMode = \
        parseObsID(obsID, -1)

    # get the truth
    if filter == 6:
        wavelength = 500  # in nm
        intrinsic35 = np.loadtxt(
            '../../simulation/activeoptics/data/intrinsic_zn.txt')
    elif filter == 7:
        wavelength = 770  # in nm
        if teleState == 'design':
            intrinsic35 = np.loadtxt(
                '../../simulation/activeoptics/data/intrinsic_zn_770nm.txt')
        elif teleState.startswith('M2xp'):
            intrinsic35 = np.loadtxt(
                'data/M2_r2_%4.2f_zn_770nm.txt' % (int(teleState[4:6])/10))
        elif teleState == 'M2rxp001deg':
            intrinsic35 = np.loadtxt(
                'data/M2_r4_0.010_zn_770nm.txt')
    elif filter == 2: #Chuck: r-band should be good enough
        wavelength = 622 #use Leff
        if teleState == 'design':
            intrinsic35 = np.loadtxt(
                '../../simulation/activeoptics/data/intrinsic_zn_r.txt')        
        elif teleState.startswith('M2xp'):
            intrinsic35 = np.loadtxt(
                'data/M2_r2_%4.2f_zn_r.txt' % (int(teleState[4:6])/10))
        elif teleState == 'M2rxp001deg':
            intrinsic35 = np.loadtxt(
                'data/M2_r4_0.010_zn_r.txt')
    else:
        wavelength = 500  # in nm
        intrinsic35 = np.loadtxt(
            '../../simulation/activeoptics/data/intrinsic_zn.txt')
        
    intrinsic35 = intrinsic35 * wavelength
    if field == 'center':
        ztrue = intrinsic35[0, 3:znmax]
    elif field == 'UR':
        ztrue = intrinsic35[31, 3:znmax] 
    elif field == 'LL':
        ztrue = intrinsic35[33, 3:znmax] 

    ax = plt.subplot(2, 1, 1)
    plt.plot(x, ztrue, label='Truth (Optics only)',
             marker='o', color='b', markersize=5)
    goodIdx = np.ones(nsnap)==1
    for isnap in range(nsnap):
        zFile = 'output/%s/wfs_%s_%03d.txt' % (obsID, obsID, isnap)
        zer[:, isnap] = np.loadtxt(zFile)

        if np.std(zer[:,isnap])>1000: #larger than 1um
            print('cwfs has a problem with snap# %d (idx starts from 0)\n'%(
                isnap))
            print(zer[:, isnap])
            goodIdx[isnap] = False
        else:
            if isnap == 0:
                plt.plot(x, zer[:, isnap], label='CWFS results (%d pairs)' % nsnap,
                         marker='.', color='r', markersize=10, linestyle='--')
            else:
                plt.plot(x, zer[:, isnap],  # label = '',
                         marker='.', color='r', markersize=10, linestyle='--')

    plt.plot(x, ztrue, marker='o', color='b', markersize=5)
    ax.set_xlim(3.5, znmax + 0.5)
    # plt.xlabel('Zernike index')
    plt.ylabel('Coefficients (nm)')
    plt.legend(loc="best", shadow=True, fancybox=True)
    plt.grid()

    if dz <3.0:
        plt.title('%3.1fmm, %4.2f arcsec, %s, field: %s, %6.1f s' % (
            dz, vKseeing500, teleState, field, exptime))
    else:
        plt.title('mid-%3.1fum, %4.2f arcsec, %s, field: %s, %6.1f s' % (
            (dz-2.0)*10, vKseeing500, teleState, field, exptime))        
    ax = plt.subplot(2, 1, 2)
    plt.plot(x, ztrue, label='Truth (Optics only)',
             marker='o', color='b', markersize=5)
    plt.errorbar(x, np.mean(zer[:, goodIdx], axis=1),
                 yerr=np.std(zer[:, goodIdx], axis=1),
                 linestyle='--', marker='.', color='r', markersize=10,
                 linewidth=2, label='CWFS Mean and STD')
    plt.plot(x, ztrue, marker='o', color='b', markersize=5)
    ax.set_xlim(3.5, znmax + 0.5)
    plt.legend(loc="best",
               shadow=True, fancybox=True)
    plt.xlabel('Zernike index')
    plt.ylabel('Mean and STD (nm)')
    plt.grid()

    plt.savefig('output/%s/wfs_%s_%d.png' % (obsID, obsID, nsnap))
    # plt.show()

    outFile = 'output/%s/wfs_%s_%d_sum.txt' % (obsID, obsID, nsnap)
    np.savetxt(outFile, np.vstack((
        ztrue, np.mean(zer, axis=1), np.std(zer, axis=1),
        np.sqrt(np.sum((zer-np.tile(ztrue.reshape(-1,1),
                                    (1,nsnap)))**2,axis=1)/nsnap) )))

# def preproc():


def parallelCwfs(obsID, eimage, instruFile, algoFile, stampSize, nsnap,
                 field, numproc, debugLevel):

    inst = cwfsInstru(instruFile, stampSize)
    algo = cwfsAlgo(algoFile, inst, debugLevel)
    if field == 'center':
        I1Field = [0, 0]
        I2Field = [0, 0]
        model = 'onAxis'
    elif field == 'UR':
        I1Field = [1.166, 1.166]
        I2Field = [1.186, 1.166]
        model = 'offAxis'
    elif field == 'LL':
        I1Field = [-1.166, -1.166]
        I2Field = [-1.186, -1.166]
        model = 'offAxis'

    jobs = []
    counter = 0
    for isnap in range(nsnap):
        # runcwfs(obsID, eimage, isnap, I1Field, I2Field, inst, algo, model)
        p = multiprocessing.Process(
            target=runcwfs, name='cwfs%d' % isnap, args=(
                obsID, eimage, isnap, I1Field, I2Field, inst, algo, model))
        jobs.append(p)
        p.start()
        counter += 1
        if (counter == numproc) or (isnap == nsnap - 1):
            for p in jobs:
                p.join()
            counter = 0
            jobs = []


def runcwfs(obsID, eimage, isnap, I1Field, I2Field, inst, algo, model):
    if eimage == 0:
        I1File = 'image/%s/wfe_%s_%03d_1.fits' % (obsID, obsID, isnap)
        I2File = 'image/%s/wfe_%s_%03d_0.fits' % (obsID, obsID, isnap)
    else:
        I1File = 'image/%s/wfs_%s_%03d_1.fits' % (obsID, obsID, isnap)
        I2File = 'image/%s/wfs_%s_%03d_0.fits' % (obsID, obsID, isnap)

    I1 = cwfsImage(I1File, I1Field, 'intra')
    I2 = cwfsImage(I2File, I2Field, 'extra')
    algo.reset(I1, I2)
    algo.runIt(inst, I1, I2, model)

    outputDir = 'output/%s' % obsID
    try:
        os.stat(outputDir)
    except:
        os.makedirs(outputDir)    
    zFile = '%s/wfs_%s_%03d.txt' % (outputDir, obsID, isnap)
    np.savetxt(zFile, algo.zer4UpNm)
    return


def runPhosim(obsID, dz, instFile, cmdFile, nsnap, filter, field, eimage,
              stampSize, numproc, debugLevel):

    if obsID[0] == '0':  #Phosim (cpp code, not phosim.py) ignores the leading '0'
        obsIDPhosim = '9'+obsID[1:]
    else:
        obsIDPhosim = obsID
       
    phosimDir = '../../simulation/phosimSE/'
    logDir = 'image/%s'% obsID
    try:
        os.stat(logDir)
    except:
        os.makedirs(logDir)
            
    phosimLog = '%s/wfs_%s.log' % (logDir, obsID)

    isc = 'lsst%02d' % (dz * 10)
    iscDir = '%s/data/%s/' % (phosimDir, isc)
    try:
        os.stat(iscDir)
    except:
        os.makedirs(iscDir)
        runProgram('ln %s/data/lsst/* %s/' % (phosimDir, iscDir))
        runProgram('rm %s/focalplanelayout.txt' % (iscDir))
        runProgram('rm %s/segmentation.txt' % (iscDir))
        runProgram('cp data/focalplanelayout_%02d.txt \
        %s/focalplanelayout.txt' % (
            dz * 10, iscDir))
        runProgram('cp data/segmentation_splitR22_S11.txt \
        %s/segmentation.txt' % iscDir)

    myargs = '%s -c %s -i %s -e %d -p %d > %s 2>&1' % (
        instFile, cmdFile, isc, eimage, numproc, phosimLog)
    if debugLevel >= 1:
        print('********Runnnig PHOSIM with following parameters********')
        print('Check the log file below for progress')
        print('%s' % myargs)

    try:
        runProgram('python %s/phosim.py' % phosimDir, argstring=myargs)
    except RuntimeError:
        print('Phosim RuntimeError')
        sys.exit()
        
    if filter == 6:
        filterStr = '1' #500, g-band
    elif filter == 7:
        filterStr = '3' #770, i-band
    else:
        filterStr = '%d' % filter

    for itra in range(2):
        if field == 'center':
            chipStr = 'R22_S11_C%d' % itra  # 0 is extra, 1 is intra
            ampStr = 'R22_S11_C%d4' % itra
        elif field == 'UR':
            chipStr = 'R44_S00_C%d' % itra  # 0 is extra, 1 is intra
            ampStr = 'R44_S00_C%d4' % itra
        elif field == 'LL':
            chipStr = 'R00_S22_C%d' % itra  # 0 is extra, 1 is intra
            ampStr = 'R00_S22_C%d4' % itra
        if itra == 0:
            ampCenter = [366, 1862]
            eCenter = [1820, 2180]
        else:
            ampCenter = [366, 1862]
            eCenter = [180, 2180]
            
        for isnap in range(nsnap):

            # get stamp from amplifer image
            if eimage == 1:
                src = '%s/output/%s_a_%s_f%s_%s_E%03d.fits.gz' % (
                    phosimDir, isc, obsIDPhosim, filterStr, ampStr, isnap)
                runProgram('gunzip -f %s' % src)
                src = src.replace('.gz', '')
                IHDU = fits.open(src)
                amp = IHDU[0].data
                IHDU.close()
                stamp = amp[
                    ampCenter[1] - stampSize / 2:ampCenter[1] + stampSize / 2,
                    ampCenter[0] - stampSize / 2:ampCenter[0] + stampSize / 2]
                stampFile = 'image/wfs_%s_%03d_%d.fits' % (obsID, isnap, itra)
                if os.path.isfile(stampFile):
                    os.remove(stampFile)
                hdu = fits.PrimaryHDU(stamp)
                hdu.writeto(stampFile)

            # get stamp from e-image
            src = '%s/output/%s_e_%s_f%s_%s_E%03d.fits.gz' % (
                phosimDir, isc, obsIDPhosim, filterStr, chipStr, isnap)
            runProgram('gunzip -f %s' % src)
            src = src.replace('.gz', '')
            IHDU = fits.open(src)
            amp = IHDU[0].data
            IHDU.close()

            nPreCut = 2
            stamp0 = amp[max(0, eCenter[1] - nPreCut* stampSize):eCenter[1] + nPreCut* stampSize,
                        max(0, eCenter[0] - nPreCut*stampSize):eCenter[0] + nPreCut*stampSize]
            centroid = ndimage.measurements.center_of_mass(stamp0)
            offsety = centroid[0] - nPreCut * stampSize + 1
            offsetx = centroid[1] - nPreCut * stampSize + 1
            if (eCenter[1] - nPreCut*stampSize<0):
                offsety -= eCenter[1] - nPreCut*stampSize
            if (eCenter[0] - nPreCut*stampSize<0):
                offsetx -= eCenter[0] - nPreCut*stampSize
            stamp = amp[
                eCenter[1] - stampSize / 2 + offsety:
                eCenter[1] + stampSize / 2 + offsety,
                eCenter[0] - stampSize / 2 + offsetx:
                eCenter[0] + stampSize / 2 + offsetx]

            stampFile = 'image/%s/wfe_%s_%03d_%d.fits' % (obsID, obsID, isnap, itra)
            if os.path.isfile(stampFile):
                os.remove(stampFile)
            stamp = np.rot90(stamp,2)
            if (stamp.shape[1] < stamp.shape[0]):
                stamp = np.hstack((stamp,np.zeros((
                    stamp.shape[0],stamp.shape[0]-stamp.shape[1]))))
            hdu = fits.PrimaryHDU(stamp)
            hdu.writeto(stampFile)

    # for f in glob.glob('%s/output/*'%phosimDir):
        # os.remove(f)


def createPertFiles(obsID, nsnap, mag, debugLevel):

    dz, r0seeing500, vKseeing500, seed, teleState, filter, field, exptime, \
      ccdMode = \
        parseObsID(obsID, debugLevel)
    if mag == -1: #5mag = 100x; we use 14 mag for dz=1.0mm
        if exptime == 15 or exptime == 10:
            mag = 14 #it takes too long if we increase intensity for 2.0mm, etc.
            # mag = 14 - np.log((dz/1.0)**2)/np.log(100**0.2)
        elif exptime == 150:
            mag = 16
        elif exptime == 1:
            mag = 12
            
    if obsID[0] == '0':  #Phosim ignores the leading '0'
        obsIDPhosim = '9'+obsID[1:]
    else:
        obsIDPhosim = obsID
                
    # create inst file
    if field == 'center':
        source = 'data/fieldCenter.inst'
    else:
        source = 'data/field%s.inst' % field
        
    instFile = 'pert/wfs_%s.inst' % obsID
    fidr = open(source, 'r')
    fidw = open(instFile, 'w')
    if teleState.startswith('M2xp'):
        fidw.write('move 6 %d\n'% (int(teleState[4:6])*100)) #in micron
    elif teleState == 'M2rxp001deg':
        fidw.write('move 8 36\n')  #in arcsec
    for line in fidr:            
        if line.startswith('Opsim_obshistid'):
            line = 'Opsim_obshistid %s\n' % obsIDPhosim
        elif line.startswith('Opsim_filter'):
            if (filter == 6):
                line = 'Opsim_filter %d\n' % 1 #500, g-band
            elif (filter == 7):
                line = 'Opsim_filter %d\n' % 3 #770, i-band
            else:
                line = 'Opsim_filter %d\n' % filter
        elif line.startswith('Opsim_rawseeing'):
            line = 'Opsim_rawseeing %6.4f\n' % r0seeing500
        elif line.startswith('SIM_SEED'):
            line = 'SIM_SEED %d\n' % seed
        elif line.startswith('SIM_VISTIME'):
            line = 'SIM_VISTIME %d\n' % (nsnap * exptime + (nsnap - 1) * 3)
        elif line.startswith('SIM_NSNAP'):
            line = 'SIM_NSNAP %d\n' % nsnap
        elif line.startswith('object'):
            line = line.replace('17.000000', '%9.6f' % mag)
            if filter < 6:
                line = line.replace('sed_500.txt', 'sed_flat.txt')
            elif filter==7:
                line = line.replace('sed_500.txt', 'sed_770.txt')
                
        fidw.write(line)
    fidr.close()
    fidw.close()

    source = 'data/designOptics.cmd'
    cmdFile = 'pert/wfs_%s.cmd' % (obsID)
    fidr = open(source, 'r')
    fidw = open(cmdFile, 'w')
    for line in fidr:
        if teleState != 'design' and line.startswith('perturbationmode'):
            fidw.write('perturbationmode 1\n')
        else:
            fidw.write(line)
    if ccdMode == 0:
        fidw.write('detectormode 0\n')
    elif ccdMode == 1:
        fidw.write('cleardefects\n')
    if r0seeing500 < 1e-5:
        fidw.write('clearturbulence\n')
        fidw.write('clearopacity\n')
        fidw.write('atmosphericdispersion 0\n')
        if ccdMode == 0:
            fidw.write('opticsonlymode 1\n')
    fidr.close()
    fidw.close()

    return dz, instFile, cmdFile, filter, field


def parseObsID(obsID, debugLevel):
    dz = int(obsID[0:2]) / 10
    r0seeing500 = int(obsID[2])*0.2
    if abs(r0seeing500 -  0.4)<1e-5:
        r0 = 0.2002
    elif abs(r0seeing500 -  0.6)<1e-5:
        r0 = 0.1382
    elif abs(r0seeing500 -  0.8)<1e-5:
        r0 = 0.1058
    elif abs(r0seeing500 -  1.0)<1e-5:
        r0 = 0.0859
    elif abs(r0seeing500 -  1.2)<1e-5:
        r0 = 0.0724
    elif abs(r0seeing500 -  1.4)<1e-5:
        r0 = 0.0626
    elif abs(r0seeing500 -  0.0)<1e-5:
        r0 = 0

    #below, make seed independent of dz.
    seed = int(obsID[3]) * 1000 + 7 + sum(int(i) for i in obsID[2:])
    if obsID[4] == '0':
        teleState = 'design'
    elif obsID[4] == '1':
        teleState = 'M2xp05mm'
    elif obsID[4] == '2':
        teleState = 'M2rxp001deg'
    elif obsID[4] == '3':
        teleState =  'M2xp02mm'
    elif obsID[4] == '4':
        teleState =  'M2xp04mm'
    elif obsID[4] == '5':
        teleState =  'M2xp06mm'
    elif obsID[4] == '6':
        teleState =  'M2xp08mm'
    elif obsID[4] == '7':
        teleState =  'M2xp10mm'
               
    # bands = 'ugrizy'
    # band = bands[int(obsID[5])]
    filter = int(obsID[5])
    if int(obsID[6]) == 0:
        field = 'center'
    elif int(obsID[6]) == 1:
        field = 'UR'
    elif int(obsID[6]) == 3:
        field = 'LL'

    exptimeList = [0.1, 1, 10, 15, 150, 1500]
    exptime = exptimeList[int(obsID[7])]
    ccdMode = int(obsID[8])
    
    L0 = 30
    if r0 > 0:
        vKseeing500 = 0.976 * 500e-9 / (r0 / 3600 / 180 * np.pi) * \
            np.sqrt(1 - 2.183 * (r0 / L0)**0.356)
        r0seeing500 = 0.976 * 500e-9 / (r0 / 3600 / 180 * np.pi)
    else:
        vKseeing500 = 0
        r0seeing500 = 0

    if debugLevel >= 0:
        print('--------------------------')
        if dz<3.0:
            print('dz=%3.1fmm' % dz)
        else:
            print('dz=2.0mm, midpoint = %dum' % ((dz-2.0)*10))
            
        print('vKseeing500=%4.2f arcsec; r0seeing500=%6.4f arcsec' % (
            vKseeing500, r0seeing500))
        print('seed=%d' % seed)
        print('teleState = %s' % teleState)
        if (filter == 6):
            print('wavelength = 500nm\n')
        if (filter == 7):
            print('wavelength = 770nm\n')
        else:
            print('filter = %s' % filter)
        print('field = %s' % field)
        print('--------------------------')

    return dz, r0seeing500, vKseeing500, seed, teleState, filter, field, exptime, ccdMode


def runProgram(command, binDir=None, argstring=None):
    myCommand = command
    if binDir is not None:
        myCommand = os.path.join(binDir, command)
    if argstring is not None:
        myCommand += (' ' + argstring)
    if subprocess.call(myCommand, shell=True) != 0:
        raise RuntimeError("Error running %s" % myCommand)

if __name__ == "__main__":
    main()
