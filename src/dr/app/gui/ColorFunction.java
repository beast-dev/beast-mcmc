/*
 * ColorFunction.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

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
