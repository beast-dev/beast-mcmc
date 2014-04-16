/*
 * TIFFWriter.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.util;

/**
 * @author Marc A. Suchard
 * @author Liya Thomas -- most code taken from Liya's November 2001 free source code
 */


//import java.awt.*;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

//import java.awt.image.IndexColorModel;

public class TIFFWriter {

    public static final short MAXROWS = 6000;    // maximum # of rows
    public static final short MAXCOLUMNS = 3000; // maximum # of columns

    public static void writeDoubleArray(String fileName, double[][] inputImageInt) {
        writeDoubleArray(fileName, inputImageInt, "png", HEATMAP);
    }

    public static void writeDoubleArrayMultiChannel(String fileName, List<double[][]> inputImageIntList, String format,
                                                    MultipleChannelColorScheme scheme) {
        // Get size, assumes the same for all matrix in list
        int dim1 = inputImageIntList.get(0).length;
        int dim2 = (inputImageIntList.get(0))[0].length;
        BufferedImage image = new BufferedImage(dim1, dim2, BufferedImage.TYPE_INT_ARGB);

        List<Double> max = new ArrayList<Double>();
        List<Double> min = new ArrayList<Double>();

        final int channels = inputImageIntList.size();

        for (int c = 0; c < channels; ++c) {
            double[][] inputImageInt = inputImageIntList.get(c);
            double tmax = Double.NEGATIVE_INFINITY;
            double tmin = Double.POSITIVE_INFINITY;
            for (int i = 0; i < dim1; ++i) {
                for (int j = 0; j < dim2; ++j) {
                    double value = inputImageInt[i][j];
                    if (value > tmax) tmax = value;
                    else if (value < tmin) tmin = value;
                }
            }
            max.add(tmax);
            min.add(tmin);
        }

        for (int i = 0; i < dim1; ++i) {
            for (int j = 0; j < dim2; ++j) {
                List<Double> input = new ArrayList<Double>();
                for (int c = 0; c < channels; ++c) {
                    double value = (inputImageIntList.get(c))[i][j];
                    input.add(value);
                }
                image.setRGB(i, j, scheme.getColor(input, min, max).getRGB());
            }
        }

        try {
            javax.imageio.ImageIO.write(image, format, new File(fileName));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private interface ColorScheme {
        Color getColor(double input, double min, double max);
    }

    private static Color blend(List<Color> colors) {
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

    private interface Ramp {
        int mix(double ratio, int c1, int c2);
    }

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

    public static final Ramp LINEAR = new Ramp() {
        public int mix(double ratio, int c1, int c2) {
            return (int) (c1 * ratio + c2 * (1.0 - ratio));
        }
    };

    private static abstract class RampedColorScheme implements ColorScheme {
        protected abstract Color getMaxColor();

        protected abstract Color getMinColor();

        private final Ramp ramp;

        private RampedColorScheme(Ramp ramp) {
            this.ramp = ramp;
        }

        private int mix(double ratio, int c1, int c2) {
            return (int) (c1 * ratio + c2 * (1.0 - ratio));
        }

        public Color getColor(double input, double min, double max) {
            double value = (input - min) / (max - min);
            int red = ramp.mix(value, getMaxColor().getRed(), getMinColor().getRed());
            int green = ramp.mix(value, getMaxColor().getGreen(), getMinColor().getGreen());
            int blue = ramp.mix(value, getMaxColor().getBlue(), getMinColor().getBlue());
            return new Color(red, green, blue);
        }
    }

    public static final RampedColorScheme WHITE_RED = new RampedColorScheme(LINEAR) {
        protected Color getMaxColor() {
            return Color.RED;
        }

        protected Color getMinColor() {
            return Color.WHITE;
        }
    };

    public static final RampedColorScheme WHITE_BLUE = new RampedColorScheme(LINEAR) {
        protected Color getMaxColor() {
            return Color.BLUE;
        }

        protected Color getMinColor() {
            return Color.WHITE;
        }
    };

    private interface ChannelColorScheme {
        Color getColor(List<Double> input, List<Double> min, List<Double> max);
    }

    public static class MultipleChannelColorScheme implements ChannelColorScheme {
        private final ColorScheme[] schemes;

        MultipleChannelColorScheme(ColorScheme[] schemes) {
            this.schemes = schemes;
        }

        public Color getColor(List<Double> input, List<Double> min, List<Double> max) {
            List<Color> colors = new ArrayList<Color>();

            final int channels = schemes.length; // assumes the same length as input, min, max
            for (int i = 0; i < channels; ++i) {
                colors.add(schemes[i].getColor(input.get(i), min.get(i), max.get(i)));
            }

            return blend(colors);
        }
    }


    public static final MultipleChannelColorScheme CHANNEL_RED = new MultipleChannelColorScheme(
            new ColorScheme[]{WHITE_RED}
    );

    public static final MultipleChannelColorScheme CHANNEL_RED_BLUE = new MultipleChannelColorScheme(
            new ColorScheme[]{WHITE_RED, WHITE_BLUE}
    );


//     public enum ColorSchemeList implements ColorScheme {
//         HEATMAP {
//             double getRampValue(double input, double min, double max) {
//                 double end = 1.0 / 6.0;
//                 double start = 0.0;
//                 return (input - min) / (max - min) * (end - start);
//             }
//
//             public Color getColor(double input, double min, double max) {
//                float hue = (float) getRampValue(input, min, max);
//                float saturation = 0.85f;
//                float alpha = 1.0f;
//                return Color.getHSBColor(hue, saturation, alpha);
//             }
//         },
//
//         WHITE_RED {
//            Color COLOR_MAX = Color.RED;
//            Color COLOR_MIN = Color.WHITE;
//
//             private int mix(double ratio, int c1, int c2) {
//                 return (int) (c1 * ratio + c2 * (1.0 - ratio));
//             }
//
//             public Color getColor(double input, double min, double max) {
//                double ratio = (input - min) / (max - min);
//                int red = mix(ratio, COLOR_MAX.getRed(), COLOR_MIN.getRed());
//                int green = mix(ratio, COLOR_MAX.getGreen(), COLOR_MIN.getGreen());
//                 int blue = mix(ratio, COLOR_MAX.getBlue(), COLOR_MIN.getBlue());
//                 return new Color(red, green, blue);
//             }
//         },
//
//
//
//         HEATMAP_TRANSPARENT {
//             double getRampValue(double input, double min, double max) {
//                 double end = 1.0 / 6.0;
//                 double start = 0.0;
//                 return (input - min) / (max - min) * (end - start);
//             }
//
//             public Color getColor(double input, double min, double max) {
//                float hue = (float) getRampValue(input, min, max);
////                float saturation = (input == 0.0) ? 0.0f : 0.85f;
//                float saturation = 0.85f;
//                float alpha = (input == 0.0) ? 0.0f : 1.0f;
////                float alpha = 1.0f;
//                return Color.getHSBColor(hue, saturation, alpha);
//             }
//         }
//     };

    public static void writeDoubleArray(String fileName, double[][] inputImageInt, String format, ColorScheme scheme) {
        BufferedImage image = new BufferedImage(inputImageInt.length, inputImageInt[0].length, BufferedImage.TYPE_INT_ARGB);

        double max = Double.NEGATIVE_INFINITY;
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < inputImageInt.length; ++i) {
            for (int j = 0; j < inputImageInt[i].length; ++j) {
                double value = inputImageInt[i][j];
                if (value > max) max = value;
                else if (value < min) min = value;
            }
        }

        for (int i = 0; i < inputImageInt.length; ++i) {
            for (int j = 0; j < inputImageInt[i].length; ++j) {
                double value = inputImageInt[i][j];
                image.setRGB(i, j, scheme.getColor(value, min, max).getRGB());
            }
        }

        try {
            javax.imageio.ImageIO.write(image, format, new File(fileName));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

//    public static double getRampValue(double input, double min, double max) {
//        double end = 1.0 / 6.0;
//        double start = 0.0;
//        return (input - min) / (max - min) * (end - start);
//    }

//    public static Color getColor(double input, double min, double max) {
//        float hue = (float) getRampValue(input, min, max);
//        float saturation = 0.85f;
//        float alpha = 1.0f;
//        return Color.getHSBColor(hue, saturation, alpha);
//    }

    // Create TIFF image of integer array
    public static void writeDoubleArray(
            DataOutputStream dataOut,
            double[][] inputImageInt) {

        final int rows = inputImageInt.length;
        final int columns = inputImageInt[0].length;

        if (rows < 0 || rows > MAXROWS || columns < 0 || columns > MAXCOLUMNS)
            throw new RuntimeException("Invalid # rows and # columns");

        // offset to the end of data (gray values) in file
        int pos = 8 + rows * columns;

        try {

			/*
                               *  Write the header
			 */
            short i, j;
            i = (short) 'I';
            j = (short) (i * 256 + i);
            fputword(dataOut, j);
            fputword(dataOut, (short) 42);
            fputlong(dataOut, pos);

			/*
                               * Write the bitmap
			 */
            for (i = 0; i < rows; i++)
                for (j = 0; j < columns; j++) {
                    int datum = (int) inputImageInt[i][j];
                    dataOut.writeByte((byte) datum);
                }

			/*
			 * Write the tags
			 */

            fputword(dataOut, (short) 8);                                                                                                    // # of tags
            writetiftag(dataOut, SubFileType, TIFFshort, 1, 1);
            writetiftag(dataOut, ImageWidth, TIFFshort, 1, columns);
            writetiftag(dataOut, ImageLength, TIFFshort, 1, rows);
            writetiftag(dataOut, BitsPerSample, TIFFshort, 1, 8);
            writetiftag(dataOut, Compression, TIFFshort, 1, 1);
            writetiftag(dataOut, PhotoMetricInterp, TIFFshort, 1, 1);                      // for gray values only
            writetiftag(dataOut, StripOffsets, TIFFlong, 1, 8);                                        // beginning of image data
            writetiftag(dataOut, PlanarConfiguration, TIFFshort, 1, 1);

            fputlong(dataOut, 0);

        } catch (java.io.IOException read) {
            System.out.println("Error occured while writing output file.");
        }

    }

    /*
     * write one TIFF tag to the IFD
     */
    static void writetiftag(DataOutputStream dataOut, short tag, short type, int length, int offset) {
        fputword(dataOut, tag);
        fputword(dataOut, type);
        fputlong(dataOut, length);
        fputlong(dataOut, offset);
    } /* writetiftag */

    /*
     * function: fputword
     */
    static void fputword(DataOutputStream dataOut, short n) {
        try {
            dataOut.writeByte((byte) n);
            dataOut.writeByte((byte) (n >> 8));
        } catch (java.io.IOException read) {
            System.out.println("Error occured while writing output file.");
        }

    } /* fputword */

    /*
     * function: fputlong
     */
    static void fputlong(DataOutputStream dataOut, int n) {
        try {
            dataOut.writeByte((byte) n);
            dataOut.writeByte((byte) (n >> 8));
            dataOut.writeByte((byte) (n >> 16));
            dataOut.writeByte((byte) (n >> 24));
        } catch (java.io.IOException read) {
            System.out.println("Error occured while writing output file.");
        }

    } /* fputlong */


    public static final short GOOD_WRITE = 0;
    public static final short BAD_WRITE = 1;
    public static final short BAD_READ = 2;
    public static final short MEMORY_ERROR = 3;
    public static final short WRONG_BITS = 4;

    public static final short RGB_RED = 0;
    public static final short RGB_GREEN = 1;
    public static final short RGB_BLUE = 2;
    public static final short RGB_SIZE = 3;

    /*
     * TIFF object sizes
     */
    public static final short TIFFbyte = 1;
    public static final short TIFFascii = 2;
    public static final short TIFFshort = 3;
    public static final short TIFFlong = 4;
    public static final short TIFFrational = 5;

    /*
     * TIFF tag names
     */
    public static final short NewSubFile = 254;
    public static final short SubFileType = 255;
    public static final short ImageWidth = 256;
    public static final short ImageLength = 257;
    public static final short RowsPerStrip = 278;
    public static final short StripOffsets = 273;
    public static final short StripByteCounts = 279;
    public static final short SamplesPerPixel = 277;
    public static final short BitsPerSample = 258;
    public static final short Compression = 259;
    public static final short PlanarConfiguration = 284;
    public static final short Group3Options = 292;
    public static final short Group4Options = 293;
    public static final short FillOrder = 266;
    public static final short Threshholding = 263;
    public static final short CellWidth = 264;
    public static final short CellLength = 265;
    public static final short MinSampleValue = 280;
    public static final short MaxSampleValue = 281;
    public static final short PhotoMetricInterp = 262;
    public static final short GrayResponseUnit = 290;
    public static final short GrayResponseCurve = 291;
    public static final short ColorResponseUnit = 300;
    public static final short ColorResponseCurves = 301;
    public static final short XResolution = 282;
    public static final short YResolution = 283;
    public static final short ResolutionUnit = 296;
    public static final short Orientation = 274;
    public static final short DocumentName = 269;
    public static final short PageName = 285;
    public static final short XPosition = 286;
    public static final short YPosition = 287;
    public static final short PageNumber = 297;
    public static final short ImageDescription = 270;
    public static final short Make = 271;
    public static final short Model = 272;
    public static final short FreeOffsets = 288;
    public static final short FreeByteCounts = 289;
    public static final short ColorMap = 320;
    public static final short Artist = 315;
    public static final short DateTime = 306;
    public static final short HostComputer = 316;
    public static final short Software = 305;

}




