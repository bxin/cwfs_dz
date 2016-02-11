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
        stampSize = 2**np.ceil(np.log(dz/1.2335/10e-3+50)/np.log(2))
    else:
        stampSize = args.stampSize
        
    if not args.phosimoff:
        runPhosim(args.obsID, dz, instFile, cmdFile, args.nsnap, filter, field,
                  args.eimage, stampSize, args.numproc, args.debugLevel)

    instruFile = 'lsst%02d' % (dz * 10)
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
        parseObsID(obsID, debugLevel)

    # get the truth
    if filter == 6:
        wavelength = 500  # in nm
        intrinsic35 = np.loadtxt(
            '../../simulation/activeoptics/data/intrinsic_zn.txt')
    elif filter == 7:
        wavelength = 770  # in nm
        intrinsic35 = np.loadtxt(
            '../../simulation/activeoptics/data/intrinsic_zn_770nm.txt')
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
    for isnap in range(nsnap):
        zFile = 'output/wfs_%s_%03d.txt' % (obsID, isnap)
        zer[:, isnap] = np.loadtxt(zFile)

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

    plt.title('%3.1fmm, %4.2f arcsec, %s, field: %s, %6.1f s' % (
        dz, vKseeing500, teleState, field, exptime))

    ax = plt.subplot(2, 1, 2)
    plt.plot(x, ztrue, label='Truth (Optics only)',
             marker='o', color='b', markersize=5)
    plt.errorbar(x, np.mean(zer, axis=1), yerr=np.std(zer, axis=1),
                 linestyle='--', marker='.', color='r', markersize=10,
                 linewidth=2, label='CWFS Mean and STD')
    plt.plot(x, ztrue, marker='o', color='b', markersize=5)
    ax.set_xlim(3.5, znmax + 0.5)
    plt.legend(loc="best",
               shadow=True, fancybox=True)
    plt.xlabel('Zernike index')
    plt.ylabel('Mean and STD (nm)')
    plt.grid()

    plt.savefig('output/wfs_%s_%d.png' % (obsID, nsnap))
    # plt.show()

    outFile = 'output/wfs_%s_%d_sum.txt' % (obsID, nsnap)
    np.savetxt(outFile, np.vstack((
        ztrue, np.mean(zer, axis=1), np.std(zer, axis=1))))

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
        I1Field = [1.168, 1.168]
        I2Field = [1.184, 1.168]
        model = 'offAxis'
    elif field == 'LL':
        I1Field = [-1.168, -1.168]
        I2Field = [-1.184, -1.168]
        model = 'offAxis'

    jobs = []
    counter = 0
    for isnap in range(nsnap):
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
        I1File = 'image/wfe_%s_%03d_1.fits' % (obsID, isnap)
        I2File = 'image/wfe_%s_%03d_0.fits' % (obsID, isnap)
    else:
        I1File = 'image/wfs_%s_%03d_1.fits' % (obsID, isnap)
        I2File = 'image/wfs_%s_%03d_0.fits' % (obsID, isnap)

    I1 = cwfsImage(I1File, I1Field, 'intra')
    I2 = cwfsImage(I2File, I2Field, 'extra')
    algo.reset(I1, I2)
    algo.runIt(inst, I1, I2, model)

    zFile = 'output/wfs_%s_%03d.txt' % (obsID, isnap)
    np.savetxt(zFile, algo.zer4UpNm)
    return


def runPhosim(obsID, dz, instFile, cmdFile, nsnap, filter, field, eimage,
              stampSize, numproc, debugLevel):

    if obsID[0] == '0':  #Phosim ignores the leading '0'
        obsIDPhosim = '9'+obsID[1:]
    else:
        obsIDPhosim = obsID
       
    phosimDir = '../../simulation/phosimSE/'
    phosimLog = 'image/log/wfs_%s.log' % (obsID)

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

    myargs = '%s -c %s -i %s -e %d -p %d > %s' % (
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
            if itra == 0:
                ampCenter = [366, 1862]
                eCenter = [1855, 2180]
            else:
                ampCenter = [366, 1862]
                eCenter = [147, 2180]
        elif field == 'UR':
            chipStr = 'R44_S00_C%d' % itra  # 0 is extra, 1 is intra
            ampStr = 'R44_S00_C%d4' % itra
            if itra == 0:
                ampCenter = [366, 1862]
                eCenter = [1133, 1880]
            else:
                ampCenter = [366, 1862]
                eCenter = [863, 1897]
        elif field == 'LL':
            chipStr = 'R00_S22_C%d' % itra  # 0 is extra, 1 is intra
            ampStr = 'R00_S22_C%d4' % itra
            if itra == 0:
                ampCenter = [366, 1862]
                eCenter = [865, 2180]
            else:
                ampCenter = [366, 1862]
                eCenter = [1147, 2177]
            
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

            stamp = amp[eCenter[1] - stampSize / 2:eCenter[1] + stampSize / 2,
                        eCenter[0] - stampSize / 2:eCenter[0] + stampSize / 2]
            stampFile = 'image/wfe_%s_%03d_%d.fits' % (obsID, isnap, itra)
            if os.path.isfile(stampFile):
                os.remove(stampFile)
            hdu = fits.PrimaryHDU(stamp)
            hdu.writeto(stampFile)

    # for f in glob.glob('%s/output/*'%phosimDir):
        # os.remove(f)


def createPertFiles(obsID, nsnap, mag, debugLevel):

    dz, r0seeing500, vKseeing500, seed, teleState, filter, field, exptime, \
      ccdMode = \
        parseObsID(obsID, debugLevel)
    if mag == -1: #5mag = 100x; we use 17 mag for dz=1.0mm
        mag = 17 - np.log((dz/1.0)**2)/np.log(100**0.2)

    if obsID[0] == '0':  #Phosim ignores the leading '0'
        obsIDPhosim = '9'+obsID[1:]
    else:
        obsIDPhosim = obsID
                
    # create inst file
    if int(obsID[7]) == 0:
        source = 'data/fieldCenter.inst'
    else:
        source = 'data/field%s.inst' % field
        
    instFile = 'pert/wfs_%s.inst' % obsID
    fidr = open(source, 'r')
    fidw = open(instFile, 'w')
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

    if obsID[5] == '0':
        source = 'data/designOptics.cmd'
    cmdFile = 'pert/wfs_%s.cmd' % (obsID)
    fidr = open(source, 'r')
    fidw = open(cmdFile, 'w')
    for line in fidr:
        fidw.write(line)
    if ccdMode == 0:
        fidw.write('detectormode 0\n')
    elif ccdMode == 1:
        fidw.write('cleardefects\n')
    if r0seeing500 == 0:
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
    r0seeing500 = obsID[2:4]
    if r0seeing500 == '04':
        r0 = 0.2002
    elif r0seeing500 == '06':
        r0 = 0.1382
    elif r0seeing500 == '08':
        r0 = 0.1058
    elif r0seeing500 == '10':
        r0 = 0.0859
    elif r0seeing500 == '12':
        r0 = 0.0724
    elif r0seeing500 == '14':
        r0 = 0.0626

    seed = int(obsID[4]) * 1000 + 7 + sum(int(i) for i in obsID)
    if obsID[5] == '0':
        teleState = 'design'
    elif obsID[5] == '1':
        teleState = 'perturbed'
    # bands = 'ugrizy'
    # band = bands[int(obsID[5])]
    filter = int(obsID[6])
    if int(obsID[7]) == 0:
        field = 'center'
    elif int(obsID[7]) == 1:
        field = 'UR'
    elif int(obsID[7]) == 3:
        field = 'LL'

    exptimeList = [0.1, 1, 10, 15, 150, 1500]
    exptime = exptimeList[int(obsID[8])]
    ccdMode = int(obsID[9])
    
    L0 = 30
    vKseeing500 = 0.976 * 500e-9 / (r0 / 3600 / 180 * np.pi) * \
        np.sqrt(1 - 2.183 * (r0 / L0)**0.356)
    r0seeing500 = 0.976 * 500e-9 / (r0 / 3600 / 180 * np.pi)

    if debugLevel >= 0:
        print('--------------------------')
        print('dz=%3.1fmm' % dz)
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
