/*
 * TabDelimitedFormatter.java
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

package dr.inference.loggers;

import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;

/**
 * A class that writes a log in tab delimited format.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TabDelimitedFormatter.java,v 1.5 2005/05/24 20:25:59 rambaut Exp $
 */
public class TabDelimitedFormatter implements LogFormatter {

    protected PrintWriter printWriter;

    public TabDelimitedFormatter(PrintWriter printWriter) {
        this.printWriter = printWriter;
    }

    public TabDelimitedFormatter(OutputStream stream) {
        this(new PrintWriter(new OutputStreamWriter(stream)));
    }

    boolean outputLabels = true;

    public TabDelimitedFormatter(PrintWriter printWriter, boolean labels) {

        this(printWriter);
        outputLabels = labels;
    }

    public void startLogging(String title) {
        // DO NOTHING    
    }

    public void logHeading(String heading) {
        if (heading != null) {
            String[] lines = heading.split("[\r\n]");
            for (String line : lines) {
                printWriter.println("# " + line);
            }
        }
        printWriter.flush();
    }

    public void logLine(String line) {
        printWriter.println(line);
        printWriter.flush();
    }

    public void logLabels(String[] labels) {
        if (outputLabels) {
            if (labels.length > 0) {
                printWriter.print(labels[0]);
            }

            for (int i = 1; i < labels.length; i++) {
                printWriter.print('\t');
                printWriter.print(labels[i]);
            }

            printWriter.println();
            printWriter.flush();
        }
    }

    public void logValues(String[] values) {

        if (values.length > 0) {
            printWriter.print(values[0]);
        }

        for (int i = 1; i < values.length; i++) {
            printWriter.print('\t');
            printWriter.print(values[i]);
        }

        printWriter.println();
        printWriter.flush();
    }

    public void stopLogging() {
        // Nothing to do...
    }

}
