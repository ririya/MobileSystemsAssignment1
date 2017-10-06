package meteor.asu.edu.speedytiltshift;

import android.content.res.ColorStateList;
import android.graphics.Bitmap;
import android.graphics.BitmapFactory;
import android.graphics.Color;
import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import android.view.View;
import android.widget.ImageView;
import android.widget.TextView;
import android.util.Log;

public class MainActivity extends AppCompatActivity {


    private Bitmap bmpIn;
    private Bitmap bmpOut;
    private ImageView imageView;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);



        BitmapFactory.Options opts = new BitmapFactory.Options();
        opts.inPreferredConfig = Bitmap.Config.ARGB_8888; // Each pixel is 4 bytes: Alpha, Red, Green, Blue
        bmpIn = BitmapFactory.decodeResource(getResources(), R.drawable.input, opts);
        imageView = (ImageView) findViewById(R.id.imageView);
    }
    public void tiltshiftjava(View view){
        TextView javaElapsedView = (TextView) findViewById(R.id.javaElapsed);
        javaElapsedView.setText("Blurring");
        javaElapsedView.setTextColor(Color.BLUE);
//        int a = 1000/0;
        long startTime = System.currentTimeMillis();
        Log.d("myTag", "Java Blur");
//        bmpOut = SpeedyTiltShift.tiltshift_java(bmpIn, 100, 200, 300, 400, 0.5f, 2.1f);
        bmpOut = SpeedyTiltShift.tiltshift_java(bmpIn, 300, 600, 1000, 1300, 5.0f, 5.0f);
        imageView.setImageBitmap(bmpOut);

        double javaTime = (double) (System.currentTimeMillis() - startTime)/1000;
        String javaTimeString = "Java Elapsed Time: " + Double.toString(javaTime) + " seconds";
        javaElapsedView.setText(javaTimeString);
        javaElapsedView.setTextColor(Color.WHITE);

    }
    public void tiltshiftcpp(View view){
        TextView cppElapsedView = (TextView) findViewById(R.id.cppElapsed);
        cppElapsedView.setText("Blurring");
        cppElapsedView.setTextColor(Color.BLUE);

        long startTime = System.currentTimeMillis();
        Log.d("myTag", "C Blur");
//        bmpOut = SpeedyTiltShift.tiltshift_cpp(bmpIn, 100, 200, 300, 400, 0.5f, 2.1f);
        bmpOut = SpeedyTiltShift.tiltshift_cpp(bmpIn, 300, 600, 1000, 1300, 5.0f, 5.0f);
        imageView.setImageBitmap(bmpOut);


        double cppTime = (double) (System.currentTimeMillis() - startTime)/1000;
        String cppTimeString = "Cpp Elapsed Time: " + Double.toString(cppTime) + " seconds";
        cppElapsedView.setText(cppTimeString);
        cppElapsedView.setTextColor(Color.WHITE);
    }
    public void tiltshiftneon(View view){
        TextView neonElapsedView = (TextView) findViewById(R.id.neonElapsed);
        neonElapsedView.setText("Blurring");
        neonElapsedView.setTextColor(Color.BLUE);

        long startTime = System.currentTimeMillis();
        Log.d("myTag", "Neon Blur");
//        bmpOut = SpeedyTiltShift.tiltshift_neon(bmpIn, 100, 200, 300, 400, 0.5f, 2.1f);
        bmpOut = SpeedyTiltShift.tiltshift_neon(bmpIn, 300, 600, 1000, 1300, 5.0f, 5.0f);
        imageView.setImageBitmap(bmpOut);

        double neonTime = (double) (System.currentTimeMillis() - startTime)/1000;
        String neonTimeString = "Neon Elapsed Time: " + Double.toString(neonTime) + " seconds";
        neonElapsedView.setText(neonTimeString);
        neonElapsedView.setTextColor(Color.WHITE);
    }
}

