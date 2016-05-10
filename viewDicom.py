from __future__ import division
import os
from matplotlib import pyplot, cm
import dicom
import numpy
import re
import argparse
import textwrap
from progressbar import AnimatedMarker,ProgressBar,Percentage,Bar


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


def filter_modality(modality, dicomList):
    newDicomList = []
    for dicomName in dicomList:
        if re.search(modality, dicomName):
            newDicomList.append(dicomName)
    return newDicomList


def getName(dicomName):
    try:
        name = re.search('([A-Z]{3}\d{2,3}_[A-Z]{2,3}', dicomName).group(1)
    except:
        name = re.sub('\/','_',dicomName)
    return name

def prac():
    dataLocation = '/Volumes/promise/CCNC_MRI_3T/NOR'
    dirList = returnDir(dataLocation)
    dtiList = filter_modality('DTI', dirList)
    for dicomDir in dtiList:
        try:
            viewDicom(dicomDir, 0, '/Volumes/promise/cerebellum')
        except:
            print dicomDir, 'does not work'

def viewDicom(dataLocation, sliceNum, outdir):
    print dataLocation
    #name = re.sub(dataLocation, '\/', '_')
    name = getName(dataLocation)

    dicomFiles = [os.path.join(dataLocation,x) for x in os.listdir(dataLocation) if x.endswith('dcm') or \
                                                     x.endswith('IMA') or \
                                                     x.endswith('DCM')]

    # Get ref file
    #RefDs = dicom.read_file(dicomFiles[0])

    RefDs = dicom.read_file(dicomFiles[sliceNum])

    # Load dimensions based on the number of rows, columns, and slices (along the Z axis)
    ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(dicomFiles))

    # Load spacing values (in mm)
    ConstPixelSpacing = (float(RefDs.PixelSpacing[0]), float(RefDs.PixelSpacing[1]), float(RefDs.SliceThickness))


    # The array is sized based on 'ConstPixelDims'
    ArrayDicom = numpy.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)

    # loop through all the DICOM files

    instance_number = {}
    num=0
    pbar=ProgressBar().start()
    for filenameDCM in dicomFiles:
        pbar.update(num)
        num+=100 / len(dicomFiles)
        # read the file
        ds = dicom.read_file(filenameDCM)
        instance_number[filenameDCM] = ds.InstanceNumber - 1
        # store the raw image data
        #ArrayDicom[:, :, dicomFiles.index(filenameDCM)] = ds.pixel_array  
        ArrayDicom[:, :, instance_number[filenameDCM]] = ds.pixel_array  
    pbar.finish()


    fig = pyplot.figure('Dicom viewer', dpi=300)
    pyplot.axes().set_aspect('equal', 'datalim')
    pyplot.set_cmap(pyplot.gray())

    #pyplot.pcolormesh(RefDs.pixel_array)
    if 'MOSAIC' in ds.ImageType:
        # DTI
        mosaic_data = ArrayDicom[:, :, 10]
        slices = blockshaped(mosaic_data, 128, 128)#$[:,:,:]
        print slices.shape, sliceNum
        stacked = numpy.hstack([
            numpy.flipud(slices[38,:,:]), 
            numpy.flipud(slices[39,:,:]),
            numpy.flipud(slices[40,:,:])
            ])
        #pyplot.pcolormesh(stacked)
        pyplot.pcolormesh(mosaic_data)

    else:
        pyplot.pcolormesh(numpy.flipud(ArrayDicom[:, :, sliceNum]))
        #pyplot.pcolormesh(numpy.flipud(ArrayDicom[:, -sliceNum, :]))
    #pyplot.imshow(ds.pixel_array)

    #cfm = pyplot.get_current_fig_manager('Dicom viewer')
    #cfm.window.raise_Tue 10 May 2016 10:27:35()

    print 'ha'
    print os.path.join(outdir, name+'.png')
    pyplot.savefig(os.path.join(outdir, name+'.png'))



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



