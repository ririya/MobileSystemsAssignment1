package meteor.asu.edu.speedytiltshift;

import android.graphics.Bitmap;
import android.graphics.BitmapFactory;
import android.util.Log;

import java.io.Console;
import java.lang.*;

/**
 * Created by roblkw on 7/26/17.
 */

class GaussianKernel{

    private double [][] kernel;
    private double weightSum;
    private int r;

    public GaussianKernel(double [][] kernel,  double weightSum, int r)
    {
        this.kernel = kernel;
        this.weightSum = weightSum;
        this.r = r;
    }

    public double[][] GetKernel()
    {
        return this.kernel;
    }

    public double GetWeightSum()
    {
        return this.weightSum;
    }

    public int Get_r()
    {
        return this.r;
    }

}


public class SpeedyTiltShift {

    static {
        System.loadLibrary("nativetiltshift-lib");
    }

    public static Bitmap tiltshift_java(Bitmap in, int a0, int a1, int a2, int a3, float s_far, float s_near){
        Bitmap out;
        out=in.copy(in.getConfig(),true);

        int width=in.getWidth();
        int height=in.getHeight();

        Log.d("TILTSHIFT_JAVA","hey:"+width+","+height);

        int pixelsLength = width*height;

        int[] pixels = new int[pixelsLength];
        int[] Gpixels = new int[pixelsLength];
        int[] Bpixels = new int[pixelsLength];
        int[] Rpixels = new int[pixelsLength];


        int offset=0;
        int stride = width;
        in.getPixels(pixels,offset,stride,0,0,width,height);
        //obtain separate channels

        Log.d("myTag", "Reading Channels");
        readChannels(height, width, pixels, Gpixels, Bpixels, Rpixels);

        //blur each channel


        Log.d("myTag", "Blur R Channel");
        Rpixels = applyGaussianBlurToAllPixels(Rpixels, width, height, a0, a1, a2, a3, s_far, s_near);
        Log.d("myTag", "Blur G Channel");
        Gpixels = applyGaussianBlurToAllPixels(Gpixels, width, height, a0, a1, a2, a3, s_far, s_near);
        Log.d("myTag", "Blur B Channel");
        Bpixels = applyGaussianBlurToAllPixels(Bpixels, width, height, a0, a1, a2, a3, s_far, s_near);
//
//        //Build blurred RGB image

        Log.d("myTag", "Buliding final RGB image");
        buildRGBimage(height, width, pixels, Bpixels, Gpixels, Rpixels);

        out.setPixels(pixels,offset,stride,0,0,width,height);

        Log.d("TILTSHIFT_JAVA","hey2");
        return out;
    }

    private static void buildRGBimage(int height, int width, int[]pixels, int[]Bpixels, int[] Gpixels, int[] Rpixels)
    {

        int BB,GG,RR,AA;

        AA = 0xff;

        for (int y=0; y<height; y++)
        {
            for (int x = 0; x < width; x++)
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

    private static void readChannels(int height, int width, int[] pixels, int[] Gpixels, int[] Bpixels, int[] Rpixels)
    {
        int p, BB,GG,RR;

        for (int y=0; y<height; y++){
            for (int x = 0; x<width; x++){
                // From Google Developer: int color = (A & 0xff) << 24 | (R & 0xff) << 16 | (G & 0xff) << 16 | (B & 0xff);
                p = pixels[y*width+x];
                BB = p & 0xff;
                GG = (p>>8)& 0xff;
                RR = (p>>16)& 0xff;
//                int AA = (p>>24)& 0xff;
                Gpixels[y*width+x] = GG;
                Rpixels[y*width+x] = RR;
                Bpixels[y*width+x] = BB;
            }
        }
    }

    private static int[] applyGaussianBlurToAllPixels(int[] pixels, int width, int height, int a0, int a1, int a2, int a3, float s_far, float s_near)
    {
        int[] blurredPixels = pixels.clone();

        GaussianKernel Gfar = initializeKernel2D(s_far);
        GaussianKernel Gnear = initializeKernel2D(s_near);
        double sigma10, sigma32;

        Log.d("myTag", "0-a0");
        for (int y=0; y<= a0; y++)
        {
            for (int x = 0; x<width; x++)
            {
                gaussianBlur2D( x,  y,  width,  height, Gfar, pixels,blurredPixels);
            }
        }

        Log.d("myTag", "a0-a1");
        for (int y=a0+1; y<a1; y++)
        {
            sigma10 = s_far*(double)(a1-y)/(double)(a1-a0);
            if (sigma10 < 0.7)
            {
                continue;
            }
            GaussianKernel G10 = initializeKernel2D(sigma10);
            for (int x = 0; x<width; x++)
            {
                gaussianBlur2D( x,  y,  width,  height,  G10, pixels,blurredPixels);
            }
        }

        Log.d("myTag", "a2-a3");
        for (int y=a2+1; y<a3; y++)
        {
            sigma32 = s_near*(double)(y-a2)/(double)(a3-a2);
            if (sigma32 < 0.7)
            {
                continue;
            }
            GaussianKernel G32 = initializeKernel2D(sigma32);
            for (int x = 0; x<width; x++)
            {
                gaussianBlur2D( x,  y,  width,  height, G32, pixels,blurredPixels);
            }
        }
        Log.d("myTag", "a3-height");
        for (int y=a3; y< height; y++)
        {
            for (int x = 0; x<width; x++)
            {
                gaussianBlur2D( x,  y,  width,  height, Gnear, pixels,blurredPixels);
            }
        }

        return blurredPixels;
    }

    private static double[] initializeKernel(double sigma)
    {
        int r = (int) Math.ceil(3*sigma);
        int kernelLen = 2*r+1;
        double[] G = new double[kernelLen];

        for(int k=-r; k<=r; k++)
        {

            G[k+r] =  Math.exp(-Math.pow(k/sigma,2)/2)/Math.sqrt(2*Math.PI*Math.pow(sigma,2));

        }

        return G;
    }

    private static GaussianKernel initializeKernel2D(double sigma)
    {
        int r = (int) Math.ceil(3 * sigma);
        int kernelLen = 2 * r + 1;
        double[][] G = new double[kernelLen][kernelLen];
        double weightSum = 0;
        double temp;

        for (int y = -r; y <= r; y++)
        {
            for (int x = -r; x <= r; x++)
            {
                temp =  Math.exp(-(Math.pow(x, 2) + Math.pow(y, 2)) / (2 * Math.pow(sigma, 2))) / (2 * Math.PI * Math.pow(sigma, 2));
                G[y + r][x + r] = temp;
                weightSum = weightSum + temp;
            }
        }

        GaussianKernel GKernel = new GaussianKernel(G,weightSum, r);

        return GKernel;
    }

    private static void gaussianBlur2D(int x, int y, int width, int height, GaussianKernel GKernel, int[] p, int[] blurredP)
    {
        double pBlur = 0;
        int pixelIndex;
        int relPixelIndex;

        pixelIndex = y * width + x;

        int r = GKernel.Get_r();
        double[][] G = GKernel.GetKernel();
        double weightSum = GKernel.GetWeightSum();

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
                    pBlur = pBlur + G[ky+r][kx+r] * p[relPixelIndex];
//                    pBlur = 0xff;
                }
            }
        }

        blurredP[pixelIndex] = (int)Math.round(pBlur/weightSum);

    }


    public static Bitmap tiltshift_cpp(Bitmap in, int a0, int a1, int a2, int a3, float s_far, float s_near){

        Bitmap out;
        out=in.copy(in.getConfig(),true);

        int width=in.getWidth();
        int height=in.getHeight();

        Log.d("TILTSHIFT_JAVA","hey:"+width+","+height);

        int[] pixels = new int[width*height];

        int offset=0;
        int stride = width;
        in.getPixels(pixels,offset,stride,0,0,width,height);

        nativeTiltShift(pixels, width, height, a0, a1, a2, a3, s_far, s_near);

        return out;
    }
    public static Bitmap tiltshift_neon(Bitmap in, int a0, int a1, int a2, int a3, float s_far, float s_near){
        return in;
    }
    private static native int[] nativeTiltShift(int[] pixels, int imgW, int imgH, int a0, int a1, int a2, int a3, float s_far, float s_near);
    private static native int[] nativeTiltShiftNeon(int[] pixels, int imgW, int imgH, int a0, int a1, int a2, int a3, float s_far, float s_near);

}
