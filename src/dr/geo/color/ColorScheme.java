/*
 * ColorScheme.java
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

/**
 * Created by msuchard on 4/18/14.
 */
public interface ColorScheme {

    Color getColor(double input, double min, double max);

    public static final ColorScheme HEATMAP = new ColorScheme() {
        double getRampValue(double input, double min, double max) {
            double end = 1.0 / 6.0;
            double start = 0.0;
            return (input - min) / (max - min) * (end - start);
        }

        public Color getColor(double input, double min, double max) {
            float hue = (float) getRampValue(input, min, max);
            float saturation = 0.85f;
            float alpha = 1.0f;
            return Color.getHSBColor(hue, saturation, alpha);
        }
    };

    interface Ramp {
        int mix(double ratio, int c1, int c2);
    }

    public static final Ramp LINEAR = new Ramp() {
        public int mix(double ratio, int c1, int c2) {
            return (int) (c1 * ratio + c2 * (1.0 - ratio));
        }
    };

    abstract class RampedColorScheme implements ColorScheme {

        protected abstract Color getMaxColor();

        protected abstract Color getMinColor();

        private final Ramp ramp;

        private RampedColorScheme(Ramp ramp) {
            this.ramp = ramp;
        }

        public Color getColor(double input, double min, double max) {
            double value = (input - min) / (max - min);
            int red = ramp.mix(value, getMaxColor().getRed(), getMinColor().getRed());
            int green = ramp.mix(value, getMaxColor().getGreen(), getMinColor().getGreen());
            int blue = ramp.mix(value, getMaxColor().getBlue(), getMinColor().getBlue());
            return new Color(red, green, blue);
        }
    }

    abstract class TransparentColorScheme implements ColorScheme {

        TransparentColorScheme(ColorScheme scheme) {
            this.scheme = scheme;
        }

        protected abstract double getTransparentValue();

        protected Color getBaseTransparentColor() {
            return Color.WHITE;
        }

        public Color getColor(double input, double min, double max) {
            if (input == getTransparentValue()) {
                int rgb = getBaseTransparentColor().getRGB();
                rgb = 0x00FFFFFF & rgb; // make complete transparent
                return new Color(rgb, true);
            } else {
                return scheme.getColor(input, min, max);
            }
        }

        private final ColorScheme scheme;
    }

    public static final ColorScheme WHITE_RED = new RampedColorScheme(LINEAR) {
        protected Color getMaxColor() {
            return Color.RED;
        }

        protected Color getMinColor() {
            return Color.WHITE;
        }
    };

    public static final ColorScheme WHITE_BLUE = new RampedColorScheme(LINEAR) {
        protected Color getMaxColor() {
            return Color.BLUE;
        }

        protected Color getMinColor() {
            return Color.WHITE;
        }
    };

    public static final ColorScheme TRANPARENT_WHITE_RED = new TransparentColorScheme(WHITE_RED) {

        protected double getTransparentValue() {
            return Double.NaN;
        }
    };

    public static final ColorScheme TRANPARENT_HEATMAP = new TransparentColorScheme(HEATMAP) {

        protected double getTransparentValue() {
            return Double.NaN;
        }
    };

    public static final ColorScheme TRANPARENT0_HEATMAP = new TransparentColorScheme(HEATMAP) {

        protected double getTransparentValue() {
            return 0.0;
        }
    };

    public static final ColorScheme TRANPARENT_HEATMAP2 = new TransparentColorScheme(
            new RampedColorScheme(LINEAR) {

                protected Color getMaxColor() {
                    return Color.YELLOW;
                }

                protected Color getMinColor() {
                    return Color.RED;
                }
            }
    ) {

        protected double getTransparentValue() {
            return Double.NaN;
        }
    };
}