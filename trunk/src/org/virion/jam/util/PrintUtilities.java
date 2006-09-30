/*
 * PrintUtilities.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package org.virion.jam.util;

import javax.swing.*;
import java.awt.*;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;
import java.lang.reflect.InvocationTargetException;

/**
 * A simple utility class that lets you very simply print
 * an arbitrary component. Just pass the component to the
 * PrintUtilities.printComponent. The component you want to
 * print doesn't need a print method and doesn't have to
 * implement any interface or do anything special at all.
 * <P>
 * If you are going to be printing many times, it is marginally more
 * efficient to first do the following:
 * <PRE>
 * PrintUtilities printHelper = new PrintUtilities(theComponent);
 * </PRE>
 * then later do printHelper.print(). But this is a very tiny
 * difference, so in most cases just do the simpler
 * PrintUtilities.printComponent(componentToBePrinted).
 * <p/>
 * 7/99 Marty Hall, http://www.apl.jhu.edu/~hall/java/
 * May be freely used or adapted.
 */

public class PrintUtilities implements Printable {
    private Component componentToBePrinted;
    private boolean scaled;

    public static void printComponent(Component c) {
        new PrintUtilities(c).print();
    }

    public static void printComponentScaled(Component c) {
        new PrintUtilities(c, true).print();
    }

    public PrintUtilities(Component componentToBePrinted) {
        this.componentToBePrinted = componentToBePrinted;
    }
    public PrintUtilities(Component componentToBePrinted, boolean scaled) {
        this(componentToBePrinted);
        this.scaled = scaled;
    }

    public void print() {
        PrinterJob printJob = PrinterJob.getPrinterJob();
        printJob.setPrintable(this);
        if (printJob.printDialog()) {
//            RepaintManager.currentManager(componentToBePrinted).paintDirtyRegions();

            try {
                printJob.print();
            } catch (PrinterException pe) {
                System.out.println("Error printing: " + pe);
            }
    }
    }

    public int print(Graphics g, PageFormat pageFormat, int pageIndex) {
        if(scaled) return printScaled(g, pageFormat, pageIndex);
        if (pageIndex > 0) {
            return (NO_SUCH_PAGE);
        } else {
            Graphics2D g2d = (Graphics2D) g;
            g2d.translate(pageFormat.getImageableX(), pageFormat.getImageableY());
            disableDoubleBuffering(componentToBePrinted);
            componentToBePrinted.paint(g2d);
            enableDoubleBuffering(componentToBePrinted);
            return (PAGE_EXISTS);
        }
    }

    public int printScaled(Graphics g, PageFormat pageFormat, int pageIndex) {
        return printScaled(componentToBePrinted,g, pageFormat, pageIndex);
    }

    public static int printScaled(Component componentToBePrinted, Graphics g, PageFormat pageFormat, int pageIndex) {
        if (pageIndex > 0) {
            return (NO_SUCH_PAGE);
        } else {
            Graphics2D g2d = (Graphics2D) g;

            double x0 = pageFormat.getImageableX();
            double y0 = pageFormat.getImageableY();

            double w0 = pageFormat.getImageableWidth();
            double h0 = pageFormat.getImageableHeight();

            double w1 = componentToBePrinted.getWidth();
            double h1 = componentToBePrinted.getHeight();

            double scale;

            if (w0 / w1 < h0 / h1) {
                scale = w0 / w1;
            } else {
                scale = h0 / h1;
            }

            g2d.translate(x0, y0);
            g2d.scale(scale, scale);

            disableDoubleBuffering(componentToBePrinted);
            componentToBePrinted.paint(g2d);
            enableDoubleBuffering(componentToBePrinted);

            return (PAGE_EXISTS);
        }
    }


    /**
     * The speed and quality of printing suffers dramatically if
     * any of the containers have double buffering turned on.
     * So this turns if off globally.
     *
     * @see enableDoubleBuffering
     */
    public static void disableDoubleBuffering(Component c) {
        RepaintManager currentManager = RepaintManager.currentManager(c);
        currentManager.setDoubleBufferingEnabled(false);
    }

    /**
     * Re-enables double buffering globally.
     */

    public static void enableDoubleBuffering(Component c) {
        RepaintManager currentManager = RepaintManager.currentManager(c);
        currentManager.setDoubleBufferingEnabled(true);
    }
}
