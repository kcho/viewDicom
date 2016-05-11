from __future__ import division
import os
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

def viewDicom(dicomDir, sliceNum, outdir):
    #name = re.sub(dicomDir, '\/', '_')
    name = getName(dicomDir)
    print 'haha'


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
        RefDs = dicom.read_file(dicomFiles[0], force=True)
        #rawData = RefDs.values()[0][3]
        try:
            rawData = RefDs.values()[0].repval
        except:
            rawData = RefDs.values()[0][3][-896*896*2:]

        #num = rawData.index('A1/PFP/FS') + 18
        imageMatrixList = map(ord, rawData)
        print len(imageMatrixList)
        data1 = np.array(imageMatrixList[0::2]).reshape((896,896))
        data2 = np.array(imageMatrixList[1::2]).reshape((896,896))
        ArrayDicom = np.zeros((896,896,65), dtype=float)


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
            ds = dicom.read_file(filenameDCM, force=True)
            try:
                rawData = ds.values()[0].repval
            except:
                rawData = RefDs.values()[0][3][-896*896*2:]
            imageMatrixList = map(ord, rawData)
            data1 = np.array(imageMatrixList[0::2]).reshape((896,896))
            data2 = np.array(imageMatrixList[1::2]).reshape((896,896))
            ArrayDicom[:, :, dicomFiles.index(filenameDCM)] = data1  
    pbar.finish()


    fig = plt.figure('Dicom viewer', dpi=300)
    plt.axes().set_aspect('equal', 'datalim')
    plt.set_cmap(plt.gray())

    #plt.pcolormesh(RefDs.pixel_array)
    try:
        if 'MOSAIC' in ds.ImageType:
            # DTI
            mosaic_data = np.flipud(ArrayDicom[:, :, 10])
            slices = blockshaped(mosaic_data, 128, 128)#$[:,:,:]
            stacked = np.hstack([
                np.flipud(slices[38,:,:]), 
                np.flipud(slices[39,:,:]),
                np.flipud(slices[40,:,:])
                ])
            #plt.pcolormesh(stacked)
            plt.pcolormesh(mosaic_data)

        else:
            plt.pcolormesh(np.flipud(ArrayDicom[:, :, sliceNum]))
            #plt.pcolormesh(np.flipud(ArrayDicom[:, -sliceNum, :]))
        #plt.imshow(ds.pixel_array)
    except:
        plt.pcolormesh(np.flipud(ArrayDicom[:,:,sliceNum]))

    #cfm = plt.get_current_fig_manager('Dicom viewer')
    #cfm.window.raise_Tue 10 May 2016 10:27:35()

    print os.path.join(outdir, name+'.png')
    plt.savefig(os.path.join(outdir, name+'.png'))



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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
            {codeName} : Visualize the dicoms
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

    args = parser.parse_args()

    #viewDicom(args.dicomDir, int(args.sliceNumber), args.out)
    prac()



