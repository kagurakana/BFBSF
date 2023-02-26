# Progressive Nonuniformity Correction for Aero-optical Thermal Radiation Images via BFBSF

This code is based on the paper *Progressive Nonuniformity Correction for Aero-optical Thermal Radiation Images via Bilateral Filtering and Bézier Surface Fitting*.

The `main.m` is in the `Code_BFBSF/raw/src` file, to run the raw code, Matlab should be installed, the code was only tested on Matlab R2021. 

## RAW set
The BFBSF Algorithm is in the `src/IR_correction.m`, the BFBSF+ code is in the `src/IR_correction_complex.m`. Data of TABLE I in our paper was test on `raw` set.

## MSP set
The Muti-stage Parallel (MSP) version is in the MSP folder. The first run will open the Matlab parallel pool automatically, which may takes a few time.


Before testing the PSNR and SSIM, the clear image and the correction image should be normlized to the `[0, 1]` interval.

If need any help with these codes, please send email to peko*stu.wit.edu.cn or 1278820830*qq.com * -> @
