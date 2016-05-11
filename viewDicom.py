#!/ccnc_bin/venv/bin/python

from __future__ import division
import os
import sys
import matplotlib.pyplot as plt
import dicom
import numpy as np
import re
import argparse
import textwrap
from progressbar import AnimatedMarker,ProgressBar,Percentage,Bar


def filter_modality(modality, dicomList):
    newDicomList = []
    for dicomName in dicomList:
        if re.search(modality, dicomName):
            newDicomList.append(dicomName)
    return newDicomList


def getName(dicomDir):
    try:
        searched = re.search('([A-Z]{3}\d{2,3}_[A-Z]{2,3})/(baseline|followup)', dicomDir, re.IGNORECASE)
        name  = '_'.join([searched.group(1), searched.group(2)])
    except:
        name = re.sub('\/','_',dicomDir)
    return name

def prac():
    #dataLocation = '/Volumes/promise/CCNC_MRI_3T/NOR'
    dataLocation = os.getcwd()
    dirList = returnDir(dataLocation)
    dtiList = filter_modality('DTI', dirList)
    for dicomDir in dtiList:
        viewDicom(dicomDir, 0, '/Volumes/promise/cerebellum')


def noMixDicom(dicomList):
    instanceNumList = []
    for i in dicomList:
        ds = dicom.read_file(i)
        instanceNum = ds.InstanceNumber
        if instanceNum in instanceNumList:
            print 'Not all dicoms are from one modality !'
            print 'Please check'
            return False
            break
        else:
            instanceNumList.append(instanceNum)
    print 'All from one modality'
    return True


def returnDir(rootDir):
    dirList = []
    for root, dirs, files in os.walk(rootDir):
        if 'DTI_FA' in root:
            pass
        elif 'DTI_EXP' in root:
            pass
        elif 'DTI_COLFA' in root:
            pass
        elif 'DKI_FA' in root:
            pass
        elif 'DKI_EXP' in root:
            pass
        elif 'DKI_COLFA' in root:
            pass
        elif 'scout' in root.lower():
            pass
        else:
            for oneFile in files:
                if oneFile.endswith('dcm') or \
                        oneFile.endswith('IMA') or \
                        oneFile.endswith('DCM'):
                    dirList.append(root)
                    print root
                    break
    return dirList 

def viewDicom(dicomDir, sliceNum, outdir, verbose=False):
    name = getName(dicomDir)

    dicomFiles = [os.path.join(dicomDir,x) for x in os.listdir(dicomDir) if x.endswith('dcm') or \
                                                     x.endswith('IMA') or \
                                                     x.endswith('DCM')]

    #if noMixDicom(dicomFiles):
        #pass
    #else:
        #os.error('There are mixture of dicoms')


    # Get ref file
    try:
        # RefDs = dicom.read_file(dicomFiles[0])
        RefDs = dicom.read_file(dicomFiles[sliceNum])
        dataArray = RefDs.pixel_array
        # Load dimensions based on the number of rows, columns, and slices (along the Z axis)
        ConstPixelDims = (int(RefDs.Rows), 
                          int(RefDs.Columns), 
                          len(dicomFiles))

        # Load spacing values (in mm)
        ConstPixelSpacing = (float(RefDs.PixelSpacing[0]), 
                             float(RefDs.PixelSpacing[1]), 
                             float(RefDs.SliceThickness))

        # The array is sized based on 'ConstPixelDims'
        ArrayDicom = np.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)

    except:
        dataArray = forcedRead(dicomFiles[sliceNum])
        ArrayDicom = np.zeros((dataArray.shape[-2],
                               dataArray.shape[-1],
                               len(dicomFiles)), dtype=float)


    fig = plt.figure('Dicom viewer', dpi=300)
    plt.axes().set_aspect('equal', 'datalim')
    plt.set_cmap(plt.gray())

    if verbose:
        pass
    else:
        plt.pcolormesh(np.flipud(dataArray))
        plt.show()
        sys.exit('Finished')


    # loop through all the DICOM files
    num=0
    pbar=ProgressBar().start()
    for filenameDCM in dicomFiles:
        pbar.update(num)
        num+=100 / len(dicomFiles)
        # read the file
        try:
            ds = dicom.read_file(filenameDCM)
            # store the raw image data
            ArrayDicom[:, :, dicomFiles.index(filenameDCM)] = ds.pixel_array  
        except:
            ArrayDicom[:, :, dicomFiles.index(filenameDCM)] = forcedRead(filenameDCM)
    pbar.finish()



    #plt.pcolormesh(RefDs.pixel_array)
    try:
        if 'MOSAIC' in ds.ImageType:
            # DTI
            mosaic_data = np.flipud(ArrayDicom[:, :, sliceNum])
            #slices = blockshaped(mosaic_data, 128, 128)#$[:,:,:]
            #plt.pcolormesh(stacked)
            plt.pcolormesh(mosaic_data)

        else:
            if sliceNum < 0:
                sliceNums = [sliceNum-x for x in range(5)]
            else:
                sliceNums = [sliceNum+x for x in range(5)]

            stacked = np.hstack([np.flipud(ArrayDicom[x,:,:]) for x in sliceNums])
            plt.pcolormesh(stacked)

    except:
        plt.pcolormesh(np.flipud(ArrayDicom[:,:,sliceNum]))

    #cfm = plt.get_current_fig_manager('Dicom viewer')
    #cfm.window.raise_Tue 10 May 2016 10:27:35()

    #print os.path.join(outdir, name+'.png')
    #plt.savefig(os.path.join(outdir, name+'.png'))
    plt.show()


def blockshaped(arr, nrows, ncols):
    """
    Return an array of shape (n, nrows, ncols) where
    n * nrows * ncols = arr.size

    If arr is a 2D array, the returned array should look like n subblocks with
    each subblock preserving the "physical" layout of arr.
    """
    h, w = arr.shape
    return (arr.reshape(h//nrows, nrows, -1, ncols)
               .swapaxes(1,2)
               .reshape(-1, nrows, ncols))

def forcedRead(dicomFile):
    '''
    When normal pydicom does not work,
    "force" option is given True.

    Returns numpy array
    '''
    ds = dicom.read_file(dicomFile, force=True)

    rawInfo = ds.values()[0]

    if 'repval' in dir(rawInfo):
        rawData = rawInfo.repval
    elif 'value' in dir(rawInfo):
        rawData = rawInfo.value
    else:
        rawData = RefDs.values()[0][3]

    matrixLenDict = {'dti':(2, 896,896), 
                    'dki':(2, 896,896),
                    'T2':(320,288),
                    'rest':(128,128),
                    'T1':(256,256)}
    for modality, matShape in matrixLenDict.iteritems():
        if modality in dicomFile:
            matrixLen = reduce(lambda x, y: x*y, matShape) 
            data = map(ord,rawData[-matrixLen:])
            data_array = np.array(data).reshape(matShape)

    return data_array


    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
            {codeName} : Visualize the dicoms using matplotlib
            ========================================
            '''.format(codeName=os.path.basename(__file__))))

    parser.add_argument(
        '-i', '--dicomDir',
        help='Data directory location, default=pwd',
        default=os.getcwd())
    parser.add_argument(
        '-s', '--sliceNumber',
        help='Number of horizontal slice',
        )
    parser.add_argument(
        '-o', '--out',
        help='png output location',
        default=os.getcwd(),
        )
    parser.add_argument(
        '-v', '--verbose',
        help='Calculate whole image array',
        action='store_true',
        default=False
        )

    args = parser.parse_args()

    viewDicom(args.dicomDir, int(args.sliceNumber), args.out, args.verbose)



