#include <jni.h>
#include <string>
#include <arm_neon.h>

//
//void readChannels(int height, int width, int pixels[], int Gpixels[], int Bpixels[], int Rpixels[]);
//void buildRGBimage(int height, int width, int pixels[], int Bpixels[], int Gpixels[], int Rpixels[]);

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

    jint Gpixels[length];
    jint Bpixels[length];
    jint Rpixels[length];

    //obtain separate channels
    readChannels(height, width, pixels, Gpixels, Bpixels, Rpixels);

    //blur each channel

//    Rpixels = applyGaussianBlurToAllPixels(Rpixels, width, height, a0, a1, a2, a3, s_far, s_near);
//    Gpixels = applyGaussianBlurToAllPixels(Gpixels, width, height, a0, a1, a2, a3, s_far, s_near);
//    Bpixels = applyGaussianBlurToAllPixels(Bpixels, width, height, a0, a1, a2, a3, s_far, s_near);
//
//        //Build blurred RGB image


    buildRGBimage(height, width, pixels, Bpixels, Gpixels, Rpixels);
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