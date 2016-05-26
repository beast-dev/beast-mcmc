/*
 * ChannelColorScheme.java
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

package dr.geo.color;

import java.awt.*;
import java.util.ArrayList;

/**
 * Created by msuchard on 4/18/14.
 */

public interface ChannelColorScheme {
    Color getColor(java.util.List<Double> input, java.util.List<Double> min, java.util.List<Double> max);


    static class MultipleChannelColorScheme implements ChannelColorScheme {
        private final ColorScheme[] schemes;

        MultipleChannelColorScheme(ColorScheme[] schemes) {
            this.schemes = schemes;
        }

        public Color getColor(java.util.List<Double> input, java.util.List<Double> min, java.util.List<Double> max) {
            java.util.List<Color> colors = new ArrayList<Color>();

            final int channels = schemes.length; // assumes the same length as input, min, max
            for (int i = 0; i < channels; ++i) {
                colors.add(schemes[i].getColor(input.get(i), min.get(i), max.get(i)));
            }

            return blend(colors);
        }

        private static Color blend(java.util.List<Color> colors) {
            double totalAlpha = 0.0;
            for (Color color : colors) {
                totalAlpha += color.getAlpha();
            }

            double r = 0.0;
            double g = 0.0;
            double b = 0.0;
            double a = 0.0;
            for (Color color : colors) {
                double weight = color.getAlpha();
                r += weight * color.getRed();
                g += weight * color.getGreen();
                b += weight * color.getBlue();
                a = Math.max(a, weight);
            }
            r /= totalAlpha;
            g /= totalAlpha;
            b /= totalAlpha;

            return new Color((int) r, (int) g, (int) b, (int) a);
        }
    }

    public static final ChannelColorScheme CHANNEL_RED = new MultipleChannelColorScheme(
            new ColorScheme[]{ColorScheme.WHITE_RED}
    );

    public static final ChannelColorScheme CHANNEL_RED_BLUE = new MultipleChannelColorScheme(
            new ColorScheme[]{ColorScheme.WHITE_RED, ColorScheme.WHITE_BLUE}
    );

}