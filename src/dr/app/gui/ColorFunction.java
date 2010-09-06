package dr.app.gui;

import java.awt.*;

/**
 * @author Alexei Drummond
 */
public class ColorFunction {

    Color[] colors;
    float[] points;

    float[] rgba1 = new float[4];
    float[] rgba2 = new float[4];


    public ColorFunction(Color[] colors, float[] points) {
        this.colors = colors;
        this.points = points;

        if (points.length != colors.length) throw new IllegalArgumentException();
        //if (points[0] != 0.0) throw new IllegalArgumentException();
        //if (points[points.length - 1] != 1.0) throw new IllegalArgumentException();
        for (int i = 0; i < points.length - 1; i++) {
            if (points[i + 1] < points[i]) {
                throw new IllegalArgumentException();
            }
        }
    }

    public Color getColor(float I) {

        for (int i = 0; i < points.length - 1; i++) {
            if (I >= points[i] && I <= points[i + 1]) {
                return interpolate(colors[i], colors[i + 1], I - points[i], points[i + 1] - points[i]);
            }
        }
        return Color.BLACK;
    }

    private Color interpolate(Color x, Color y, float s, float t) {

        rgba1 = x.getRGBComponents(rgba1);
        rgba2 = y.getRGBComponents(rgba2);
        float[] rgba3 = new float[4];

        for (int i = 0; i < rgba1.length; i++) {
            rgba3[i] = (s * rgba2[i] + (t - s) * rgba1[i]) / t;
        }

        return new Color(rgba3[0], rgba3[1], rgba3[2], rgba3[3]);
    }

    public void setAlpha(float alpha) {

        for (int i = 0; i < colors.length; i++) {
            float[] rgba = colors[i].getRGBComponents(rgba1);
            colors[i] = new Color(rgba[0], rgba[1], rgba[2], alpha);
        }
    }
}
