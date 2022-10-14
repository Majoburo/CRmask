import numpy as np
from astropy import stats
from astropy.io import fits
import glob

#Value for resistant mean
def cr_rej(im, lo_rej=5, hi_rej=5, rot=1):
    sigma = 3
    #if rot == 1:
    #    im = rotate(im)
    nim,nx,ny = im.shape

    mask = np.full((nim, nx, ny),True)
    clean =  np.full((nx, ny),True)
    good =  np.full((nx, ny),True)

    #Rescale spectra
    #Could subtract floor, but assuming roughly equal exposures
    if nim > 1:
        scale0 = np.median(im[0,:,:])
        for i in range(1,nim):
            scale = np.median(im[i,:,:])
            im[i,:,:] = im[i,:,:]*scale0/scale

    #Ideally, we would look row-by-row first, doing pixel-by-pixel

    mean, median, meansig = stats.sigma_clipped_stats(im, sigma=sigma, maxiters=5,axis=0)
    mask = ((im - mean) > -1*lo_rej*meansig)*((im - mean) < hi_rej*meansig)
    good = np.sum(mask,axis=0)
    clean = np.sum(im*mask/good,axis=0)


    return clean

def main():
    files=["r1090.fits",
          "r1091.fits",
          "r1092.fits",
          "r1093.fits"]
    im=[]
    for file in files:
        hdu = fits.open(file)
        im.append(hdu[0].data)
    im = np.array(im)

    clean = cr_rej(im, lo_rej=5, hi_rej=5, rot=1)
    hdu[0].data = clean
    hdu.writeto('test.fits',overwrite=True)
    hdu.close()

if __name__ == "__main__":
    main()


