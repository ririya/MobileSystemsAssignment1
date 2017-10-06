#include <jni.h>
#include <string>
#include <arm_neon.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <android/log.h>

#define PI 3.14159
#define Vertical 0
#define Horizontal 1

#define neighborElSize 8
#define kernelElSize 4

#define kernelSize = 31;

struct GaussianKernel2D {

    int r;
    double weightSum;

    double **kernel;
};

struct GaussianKernelNeon
{
    float32_t* G;
    int kernelLen;
    int padLen;
    int r;
    int offset;
};

void readChannels(jint height, jint width, jint *pixels, int *Gpixels, int *Bpixels, int *Rpixels)
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

void buildRGBimage(jint height, jint width, jint pixels[], int Bpixels[], int Gpixels[], int Rpixels[])
{

    jint x,y, BB,GG,RR,AA;
    jint color;
    AA = 0xff;

    for (y=0; y<height; y++)
    {
        for (x = 0; x < width; x++)
        {
            BB = Bpixels[y*width+x];
            GG = Gpixels[y*width+x];
            RR = Rpixels[y*width+x];
//                int AA = 0xff;
            color = (AA & 0xff) << 24 | (RR & 0xff) << 16 | (GG & 0xff) << 8 | (BB & 0xff);
            pixels[y*width+x] = color;
        }
    }
}

double** initializeKernel2D(jdouble sigma)
{
    int r = (int) ceil(3 * sigma);
    int kernelLen = 2 * r + 1;
    double weightSum = 0;
    double temp;

//    struct GaussianKernel2D GKernel;
    double** G = (double**) malloc(kernelLen * sizeof(double *));

    int y, x;

    for ( y = -r; y <= r; y++)
    {
        G[y+r] = (double*)  malloc(kernelLen * sizeof(double));
        for (x = -r; x <= r; x++)
        {
            temp =  exp(-(pow(x, 2) + pow(y, 2)) / (2 * pow(sigma, 2))) / (2 * PI * pow(sigma, 2));
            G[y + r][x + r] = temp;
            weightSum = weightSum + temp;
        }
    }

//    GKernel.r = r;
//    GKernel.weightSum = weightSum;

//    return GKernel;
    return G;
}


void gaussianBlur2D(jdouble sigma, jint x, jint y, jint width, jint height, double **G, int *p, int *blurredP)
{
    double pBlur = 0;
    long pixelIndex;
    long relPixelIndex;
    int r = (int) ceil(3 * sigma);

    pixelIndex = y * width + x;

//    int r = (int) ceil(3 * sigma);
//    int r = GKernel.r;
//    double weightSum = GKernel.weightSum;

    int ky, kx, relY,relX;

    for(ky=-r; ky<=r; ky++)
    {
        for (kx = -r; kx <= r; kx++)
        {
            relY = y + ky;
            if (relY < 0){
                continue;
            }
            if(relY >= height)
            {
                continue;
            }
            relX = x + kx;
            if(relX < 0)
            {
                continue;
            }
            if(relX >= width)
            {
                continue;
            }

            relPixelIndex = relY * width + relX;
            pBlur = pBlur + G[ky+r][kx+r] * p[relPixelIndex];

        }
    }

    blurredP[pixelIndex] = (int) round(pBlur);

}



int* applyGaussianBlurToAllPixels(int *pixels, long length, jint width, jint height, jint a0, jint a1, jint a2, jint a3, jfloat s_far, jfloat s_near)
{
    int* blurredPixels = (int *) malloc(sizeof(int) * length);

    memcpy(blurredPixels,pixels,sizeof(int) * length);

    double **Gfar;
    double **G10;
    double **G32;
    double **Gnear;

    double sigma10,sigma32;
    jint y,x;

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "0 to a0");
    Gfar = initializeKernel2D(s_far);

    for (y=0; y<= a0; y++)
    {
        for (x = 0; x<width; x++)
        {
            gaussianBlur2D(s_far, x,  y,  width,  height, Gfar, pixels,blurredPixels);
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

        else
        {
            G10 = initializeKernel2D(sigma10);
            for (x = 0; x < width; x++) {
                gaussianBlur2D(sigma10, x, y, width, height, G10, pixels, blurredPixels);
            }

            free(G10);
        }
    }

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "a2 to a3");

    for (y=a2+1; y<a3; y++)
    {
        sigma32 = s_near*(double)(y-a2)/(double)(a3-a2);
        if (sigma32<0.7)
        {
            continue;
        }
        else
        {

            G32 = initializeKernel2D(sigma32);
            for (x = 0; x < width; x++) {
                gaussianBlur2D(sigma32, x, y, width, height, G32, pixels, blurredPixels);
            }

            free(G32);

        }

    }

    Gnear = initializeKernel2D(s_near);

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "a3 to height");

    for (y=a3; y< height; y++)
    {
        for (x = 0; x<width; x++)
        {
            gaussianBlur2D(s_near, x,  y,  width,  height, Gnear, pixels,blurredPixels);
        }
    }

    free(Gnear);
    free(Gfar);

    return blurredPixels;
}


extern "C" {
JNIEXPORT jintArray JNICALL
Java_meteor_asu_edu_speedytiltshift_SpeedyTiltShift_nativeTiltShift(JNIEnv *env,
                                                                    jobject This,
                                                                    jintArray pixels_,
                                                                    jint width, jint height,
                                                                    jint a0, jint a1, jint a2,
                                                                    jint a3, jfloat s_far,
                                                                    jfloat s_near
) {
    jint *pixels = env->GetIntArrayElements(pixels_, NULL);
    const long length = env->GetArrayLength(pixels_);
    jintArray pixelsOut = env->NewIntArray(length);

    int *Gpixels = (int *) malloc(sizeof(int) * length);
    int *Bpixels = (int *) malloc(sizeof(int) * length);
    int *Rpixels = (int *) malloc(sizeof(int) * length);

    int *GpixelsBlur = (int *) malloc(sizeof(int) * length);
    int *BpixelsBlur = (int *) malloc(sizeof(int) * length);
    int *RpixelsBlur = (int *) malloc(sizeof(int) * length);

    //obtain separate channels
    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "Reading Channels");
    readChannels(height, width, pixels, Gpixels, Bpixels, Rpixels);

    //blur each channel

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "Applying R blur");
    RpixelsBlur = applyGaussianBlurToAllPixels(Rpixels, length, width, height, a0, a1, a2, a3,
                                               s_far, s_near);
    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "Applying G blur");
    GpixelsBlur = applyGaussianBlurToAllPixels(Gpixels, length, width, height, a0, a1, a2, a3,
                                               s_far, s_near);
    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "Applying B blur");
    BpixelsBlur = applyGaussianBlurToAllPixels(Bpixels, length, width, height, a0, a1, a2, a3,
                                               s_far, s_near);

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
}

GaussianKernelNeon initializeKernel1DNeon(double sigma)
{
    int r = (int) ceil(3 * sigma);
    int kernelLen = 2 * r + 1;

    struct GaussianKernelNeon GKernel;

    int numNeonVecs = ceil((float) kernelLen/kernelElSize);
    int padLen =  numNeonVecs*kernelElSize;

    float32_t * G = (float32_t*) calloc(padLen, sizeof(float32_t*));   //calloc initializes all elements with 0's

    int k;

    for(k=-r; k<=r; k++)
    {
        G[k+r] =  exp(-pow(k/sigma,2)/2)/sqrt(2*PI*pow(sigma,2));
    }

    GKernel.r = r;
    GKernel.G = G;
    GKernel.kernelLen = kernelLen;
    GKernel.padLen = padLen;
    GKernel.offset = 0;

    return GKernel;
}

float32_t gaussianBlurChannelNeon(uint8x8_t vec_8, GaussianKernelNeon GKernel, float32_t blurredPixel)
{
    float32x4_t sum_vec = vdupq_n_f32(0);
    float32x4_t curr_kernel;

    uint16x8_t vec_u16 = vmovl_u8(vec_8);

    int tempOffset = GKernel.offset;   // saves temporary offset ( will be used for other channels)

    uint16x4_t vec_16_low = vget_low_u16(vec_u16);
    uint32x4_t vec_u32_low = vmovl_u16(vec_16_low);
    float32x4_t vec_f32_low = vcvtq_f32_u32(vec_u32_low);

    curr_kernel = vld1q_f32(GKernel.G + GKernel.offset);

    sum_vec = vmlaq_f32(sum_vec,curr_kernel,vec_f32_low);

    GKernel.offset += kernelElSize;

    if (GKernel.offset >= GKernel.padLen) {    // pixel block length might be greater than kernel length because pixels are read in blocks of 8 and kernel in blocks of 4

        blurredPixel += vgetq_lane_f32(sum_vec,0);
        blurredPixel += vgetq_lane_f32(sum_vec,1);
        blurredPixel += vgetq_lane_f32(sum_vec,2);
        blurredPixel += vgetq_lane_f32(sum_vec,3);

        GKernel.offset = tempOffset;          //resets temporary offset to be used for other channels

        return blurredPixel;
    }

    uint16x4_t vec_16_high = vget_high_u16(vec_u16);
    uint32x4_t vec_u32_high = vmovl_u16(vec_16_high);
    float32x4_t vec_f32_high = vcvtq_f32_u32(vec_u32_high);

    curr_kernel = vld1q_f32(GKernel.G + GKernel.offset);

    sum_vec = vmlaq_f32(sum_vec,curr_kernel,vec_f32_high);

    blurredPixel += vgetq_lane_f32(sum_vec,0);
    blurredPixel += vgetq_lane_f32(sum_vec,1);
    blurredPixel += vgetq_lane_f32(sum_vec,2);
    blurredPixel += vgetq_lane_f32(sum_vec,3);

    GKernel.offset = tempOffset;       //resets temporary offset to be used for other channels

    return blurredPixel;

}

jint gaussianBlurNeon(int x,  int y,  int width,  int height, int direction, GaussianKernelNeon GKernel, jint* pixels)
{
    int i, pixels_lowerLimit, pixels_upperLimit, temp_lowerOffset, pixelRange;

    int numNeighborVecs = ceil((float)GKernel.padLen/neighborElSize);

    int neighborPadLen = numNeighborVecs*neighborElSize;

    jint* neighborPixels = (jint*) calloc(neighborPadLen, sizeof(jint*));   // calloc initializes all elements with 0's

    if (direction == Vertical)
    {
        int temp = y;
        y = x;
        x = temp;

        temp = height;
        height = width;
        width = temp;
    }

    pixels_upperLimit = width < x+GKernel.r ? width : x+GKernel.r;
    pixels_lowerLimit =  0 > x-GKernel.r? 0 : x-GKernel.r;
    pixelRange = pixels_upperLimit - pixels_lowerLimit + 1;

    temp_lowerOffset = pixels_lowerLimit - (x-GKernel.r);

    memcpy(neighborPixels+temp_lowerOffset, pixels+y*width + pixels_lowerLimit, sizeof(jint) * pixelRange);

    uint8_t* neighborByteArray = (uint8_t*) (neighborPixels);

    int offset = 0;

    float32_t  R_blurredPixel = 0;
    float32_t G_blurredPixel = 0;
    float32_t B_blurredPixel = 0;

        for (int i = 0; i < numNeighborVecs; i++)
        {
            uint8x8x4_t vecs = vld4_u8(neighborByteArray+4*neighborElSize*i);
            uint8x8_t B = vecs.val[0];
            uint8x8_t G = vecs.val[1];
            uint8x8_t R = vecs.val[2];

            R_blurredPixel = gaussianBlurChannelNeon(R, GKernel, R_blurredPixel);
            G_blurredPixel = gaussianBlurChannelNeon(G, GKernel, G_blurredPixel);
            B_blurredPixel = gaussianBlurChannelNeon(B, GKernel, B_blurredPixel);

            GKernel.offset += neighborElSize;

            if (GKernel.offset >= GKernel.padLen){
                break;
            }
        }


    int RR = roundf(R_blurredPixel);
    int GG = roundf(G_blurredPixel);
    int BB = roundf(B_blurredPixel);

    jint colorBlurredPixel = (0xff) << 24 | (RR & 0xff) << 16 | (GG & 0xff) << 8 |  BB & 0xff;

    free(neighborPixels);

    return colorBlurredPixel;

}


jint* transposeImage(jint* pixels, long length, jint width, jint height)
{
    int x,y;

    jint* transPixels = (jint *) malloc(sizeof(jint) * length);

    for (y=0;y<height;y++)
    {
        for (x = 0; x < width; x++)
        {
            transPixels[x*height + y] = pixels[y*width + x];
        }
    }

    return transPixels;
}


jint* gaussianBlurSingleDirection(jint* pixels, int direction, long length, jint width, jint height, jint a0, jint a1, jint a2, jint a3, jfloat s_far, jfloat s_near)
{
    int x,y;

    jint* blurredPixels =  (jint *) malloc(sizeof(jint) * length);
    memcpy(blurredPixels,pixels,sizeof(jint) * length);

    jint* pixelsOrg;
    jint* transPixels;

    if (direction == Vertical)   //transpose image for easier reading of the columns
    {
        transPixels = transposeImage(pixels, length, width, height);
        pixelsOrg = transPixels;
    }
    else
    {
        pixelsOrg = pixels;
    }

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "0 to a0");

    double sigma10;
    double sigma32;
    struct GaussianKernelNeon Gfar = initializeKernel1DNeon(s_far);


    for (y=0; y<= a0; y++)
    {
        for (x = 0; x<width; x++)
        {
            blurredPixels[y * width + x] = gaussianBlurNeon(x,  y,  width,  height, direction, Gfar, pixelsOrg);
        }

    }

    free(Gfar.G);

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "a3 to height");

    struct GaussianKernelNeon Gnear = initializeKernel1DNeon(s_near);

    for (y=a3; y< height; y++)
    {
        for (x = 0; x<width; x++)
        {
            blurredPixels[y * width + x] = gaussianBlurNeon(x,  y,  width,  height, direction, Gnear, pixelsOrg);
        }
    }

    free(Gnear.G);

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "a0 to a1");

    for (y=a0+1; y<a1; y++)
    {
        sigma10 = s_far*(double)(a1-y)/(double)(a1-a0);

        if (sigma10<0.7)
        {
            continue;
        }

        else
        {
            struct GaussianKernelNeon G10 = initializeKernel1DNeon(sigma10);
            for (x = 0; x < width; x++) {
                blurredPixels[y * width + x] = gaussianBlurNeon(x,  y,  width,  height, direction, G10, pixelsOrg);
            }

            free(G10.G);
        }
    }

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "a2 to a3");

    for (y=a2+1; y<a3; y++)
    {
        sigma32 = s_near*(double)(y-a2)/(double)(a3-a2);
        if (sigma32<0.7)
        {
            continue;
        }
        else
        {
            struct GaussianKernelNeon G32 = initializeKernel1DNeon(sigma32);
            for (x = 0; x < width; x++)
            {
                blurredPixels[y * width + x] = gaussianBlurNeon(x,  y,  width,  height, direction, G32, pixelsOrg);
            }

            free(G32.G);
        }
    }

    if (direction == Vertical)
    {
        free(transPixels);
    }

    return blurredPixels;

}


extern "C" {
JNIEXPORT jintArray JNICALL
Java_meteor_asu_edu_speedytiltshift_SpeedyTiltShift_nativeTiltShiftNeon(JNIEnv *env,
                                                                        jobject This,
                                                                        jintArray pixels_,
                                                                        jint width, jint height,
                                                                        jint a0, jint a1, jint a2,
                                                                        jint a3, jfloat s_far,
                                                                        jfloat s_near
) {
    int32x4_t sum_vec = vdupq_n_s32(0);
    jint *pixels = env->GetIntArrayElements(pixels_, NULL);
    long length = env->GetArrayLength(pixels_);
    jintArray pixelsOut = env->NewIntArray(length);

    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "Vertical blur");
    jint * blurredVertical = gaussianBlurSingleDirection(pixels, Vertical, length, width, height, a0, a1, a2, a3, s_far, s_near);
    __android_log_print(ANDROID_LOG_DEBUG, "MyTag", "Horizontal blur");
    jint * blurredHorizontal = gaussianBlurSingleDirection(blurredVertical,Horizontal, length, width, height, a0, a1, a2, a3, s_far, s_near);

    env->SetIntArrayRegion(pixelsOut, 0, length, blurredHorizontal);
    env->ReleaseIntArrayElements(pixels_, pixels, 0);

    free(blurredVertical);
    free(blurredHorizontal);

    return pixelsOut;

}
}
