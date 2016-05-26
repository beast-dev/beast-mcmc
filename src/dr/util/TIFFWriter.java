/*
 * TIFFWriter.java
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

package dr.util;

/**
 * @author Marc A. Suchard
 * @author Liya Thomas -- most code taken from Liya's November 2001 free source code
 */


//import java.awt.*;

import dr.geo.color.ChannelColorScheme;
import dr.geo.color.ColorScheme;

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
        writeDoubleArray(fileName, inputImageInt, "png", ColorScheme.HEATMAP);
    }

    public static void writeDoubleArrayMultiChannel(String fileName, List<double[][]> inputImageIntList, String format,
                                                    ChannelColorScheme scheme) {
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




