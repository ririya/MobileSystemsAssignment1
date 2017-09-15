#include <jni.h>
#include <string>
#include <arm_neon.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <android/log.h>

#define PI 3.14159

//
//void readChannels(int height, int width, int pixels[], int Gpixels[], int Bpixels[], int Rpixels[]);
//void buildRGBimage(int height, int width, int pixels[], int Bpixels[], int Gpixels[], int Rpixels[]);

struct GaussianKernel {

    int r;
    double weightSum;
//    double kernel[15][15];
    double **kernel;
};

void readChannels(jint height, jint width, jint *pixels, jint *Gpixels, jint *Bpixels, jint *Rpixels)
{
    jint x,y,p, BB, GG, RR;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            // From Google Developer: int color = (A & 0xff) << 24 | (R & 0xff) << 16 | (G & 0xff) << 16 | (B & 0xff);
            p = pixels[y * width + x];
            BB = p & 0xff;
            GG = (p >> 8) & 0xff;
            RR = (p >> 16) & 0xff;
//                int AA = (p>>24)& 0xff;
            Gpixels[y * width + x] = GG;
            Rpixels[y * width + x] = RR;
            Bpixels[y * width + x] = BB;
        }
    }
}

void buildRGBimage(jint height, jint width, jint pixels[], jint Bpixels[], jint Gpixels[], jint Rpixels[])
{

    jint x,y, BB,GG,RR,AA;

    AA = 0xff;

    for (y=0; y<height; y++)
    {
        for (x = 0; x < width; x++)
        {
            BB = Bpixels[y*width+x];
            GG = Gpixels[y*width+x];
            RR = Rpixels[y*width+x];
//                int AA = 0xff;
            int color = (AA & 0xff) << 24 | (RR & 0xff) << 16 | (GG & 0xff) << 8 | (BB & 0xff);
            pixels[y*width+x] = color;
        }
    }
}


GaussianKernel initializeKernel2D(jdouble sigma)
{
    int r = (int) ceil(3 * sigma);
    int kernelLen = 2 * r + 1;
//    double G[kernelLen][kernelLen];
    double weightSum = 0;
    double temp;

    struct GaussianKernel GKernel;
    GKernel.kernel = (double**) malloc(kernelLen * sizeof(double *));

    for (int y = -r; y <= r; y++)
    {
        GKernel.kernel[y+r] = (double*)  malloc(kernelLen * sizeof(double));
        for (int x = -r; x <= r; x++)
        {
            temp =  exp(-(pow(x, 2) + pow(y, 2)) / (2 * pow(sigma, 2))) / (2 * PI * pow(sigma, 2));
            GKernel.kernel[y + r][x + r] = temp;
            weightSum = weightSum + temp;
        }
    }

    GKernel.r = r;
    GKernel.weightSum = weightSum;

    return GKernel;
}


void gaussianBlur2D(jint x, jint y, jint width, jint height, GaussianKernel GKernel, jint *p, jint *blurredP)
{
    double pBlur = 0;
    long pixelIndex;
    long relPixelIndex;

    pixelIndex = y * width + x;

    int r = GKernel.r;
    double weightSum = GKernel.weightSum;

    for(int ky=-r; ky<=r; ky++)
    {
        for (int kx = -r; kx <= r; kx++)
        {
            int relY = y + ky;
            int relX = x + kx;

            if (relY < 0 || relY >= height || relX < 0 || relX >= width)
            {
                continue;
            }
            else
            {
                relPixelIndex = relY * width + relX;
                pBlur = pBlur + GKernel.kernel[ky+r][kx+r] * p[relPixelIndex];
            }
        }
    }

    blurredP[pixelIndex] = (int) round(pBlur/weightSum);

}



jint* applyGaussianBlurToAllPixels(jint *pixels, long length, jint width, jint height, jint a0, jint a1, jint a2, jint a3, jfloat s_far, jfloat s_near)
{
    jint* blurredPixels = (jint *) malloc(sizeof(jint) * length);
    memcpy(blurredPixels,pixels,length);

    GaussianKernel Gfar = initializeKernel2D(s_far);
    GaussianKernel Gnear = initializeKernel2D(s_near);
    double sigma10,sigma32;

    jint y,x;

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "0 to a0");

    for (y=0; y<= a0; y++)
    {
        for (x = 0; x<width; x++)
        {
            gaussianBlur2D( x,  y,  width,  height, Gfar, pixels,blurredPixels);
        }
    }

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "a0 to a1");



    for (y=a0+1; y<a1; y++)
    {
        sigma10 = s_far*(double)(a1-y)/(double)(a1-a0);

        if (sigma10<0.7)
        {
            continue;
        }

        GaussianKernel G10 = initializeKernel2D(sigma10);
        for (x = 0; x<width; x++)
        {
            gaussianBlur2D( x,  y,  width,  height,  G10, pixels,blurredPixels);
        }

        free(G10.kernel);
    }

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "a2 to a3");

    for (y=a2+1; y<a3; y++)
    {
        sigma32 = s_near*(double)(y-a2)/(double)(a3-a2);
        if (sigma32<0.7)
        {
            continue;
        }
        GaussianKernel G32 = initializeKernel2D(sigma32);
        for (x = 0; x<width; x++)
        {
            gaussianBlur2D( x,  y,  width,  height, G32, pixels,blurredPixels);
        }

        free(G32.kernel);
    }

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "a3 to height");

    for (y=a3; y< height; y++)
    {
        for (x = 0; x<width; x++)
        {
            gaussianBlur2D( x,  y,  width,  height, Gnear, pixels,blurredPixels);
        }
    }

    free(Gfar.kernel);
    free(Gnear.kernel);

    return blurredPixels;
}




extern "C" {
JNIEXPORT jintArray JNICALL
Java_meteor_asu_edu_speedytiltshift_SpeedyTiltShift_nativeTiltShift(JNIEnv *env,
                                                                    jobject This,
                                                                    jintArray pixels_,
                                                                    jint width, jint height,
                                                                    jint a0, jint a1, jint a2, jint a3, jfloat s_far, jfloat s_near
                                                                    ) {
    int32x4_t sum_vec = vdupq_n_s32(0);
    jint *pixels = env->GetIntArrayElements(pixels_, NULL);
    long length = env->GetArrayLength(pixels_);
    jintArray pixelsOut = env->NewIntArray(length);

    jint *Gpixels =  (jint *) malloc(sizeof(jint) * length);
    jint *Bpixels = (jint *) malloc(sizeof(jint) * length);
    jint *Rpixels = (jint *) malloc(sizeof(jint) * length);

    jint *GpixelsBlur;
    jint *BpixelsBlur;
    jint *RpixelsBlur;

    //obtain separate channels
    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "Reading Channels");
    readChannels(height, width, pixels, Gpixels, Bpixels, Rpixels);

    //blur each channel

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "Applying R blur");
    RpixelsBlur = applyGaussianBlurToAllPixels(Rpixels, length, width, height, a0, a1, a2, a3, s_far, s_near);
    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "Applying G blur");
    GpixelsBlur = applyGaussianBlurToAllPixels(Gpixels, length, width, height, a0, a1, a2, a3, s_far, s_near);
    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "Applying B blur");
    BpixelsBlur = applyGaussianBlurToAllPixels(Bpixels, length, width, height, a0, a1, a2, a3, s_far, s_near);


//
//        //Build blurred RGB image
    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "Building final RGB Image");
    buildRGBimage(height, width, pixels, BpixelsBlur, GpixelsBlur, RpixelsBlur);

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "Freeing Memory");
    free(Gpixels);
    free(Bpixels);
    free(Rpixels);
    free(GpixelsBlur);
    free(BpixelsBlur);
    free(RpixelsBlur);

    env->SetIntArrayRegion(pixelsOut, 0, length, pixels);
    env->ReleaseIntArrayElements(pixels_, pixels, 0);
    return pixelsOut;
}


JNIEXPORT jintArray JNICALL
Java_meteor_asu_edu_speedytiltshift_SpeedyTiltShift_nativeTiltShiftNeon(JNIEnv *env,
                                                                    jobject This,
                                                                    jintArray pixels_,
                                                                    jint width, jint height,
                                                                    jint a0, jint a1, jint a2, jint a3, jfloat s_far, jfloat s_near
) {
    int32x4_t sum_vec = vdupq_n_s32(0);
    jint *pixels = env->GetIntArrayElements(pixels_, NULL);
    long length = env->GetArrayLength(pixels_);
    jintArray pixelsOut = env->NewIntArray(length);

    env->SetIntArrayRegion(pixelsOut, 0, length, pixels);
    env->ReleaseIntArrayElements(pixels_, pixels, 0);
    return pixelsOut;
}
}