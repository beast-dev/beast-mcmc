/*
 * TreeLogger.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.tree;

import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;

import java.io.*;

/**
 * @author Guy Baele
 */
public class TreeDebugLogger extends MCLogger {

    private DebugTabDelimitedFormatter formatter;
    private OutputStream fileOut;
    private PrintStream out;
    private String fileName;
    private File file;
    private TreeLogger treeLogger;

    public TreeDebugLogger(String fileName, TreeLogger treeLogger) {
        super(true, 0);

        this.fileName = fileName;
        this.file = new File(fileName);
        try {
            this.fileOut = new FileOutputStream(this.file);
            this.out = new PrintStream(this.fileOut);
        } catch (IOException ioe) {
            throw new RuntimeException("File '" + fileName +
                    "' can not be opened to perform tree logging (TreeDebugLogger).");
        }
        //boolean argument ensures that a new file is written every time
        this.formatter = new DebugTabDelimitedFormatter(this.out);

        //make a TreeLogger object based on arguments
        //TODO Perform a more proper check on the final arguments, could/should copy them from TreeLogger
        this.treeLogger = new TreeLogger(treeLogger.getTree(), treeLogger.getBranchRates(), treeLogger.getTreeAttributeProviders(),
                treeLogger.getTreeTraitProviders(), this.formatter, 1, true, true, true, null, null);
    }

    public void writeToFile(long state) {

        try {
            //make sure to reopen output file
            this.fileOut = new FileOutputStream(this.file);
            this.out = new PrintStream(this.fileOut);
            this.formatter.setPrintStream(this.out);

            //make sure start of NEXUS file is written
            this.treeLogger.startLogging();

            //perform custom log
            this.treeLogger.log(state);

            //make sure NEXUS file has a proper end
            this.treeLogger.stopLogging();

        } catch (IOException ioe) {
            throw new RuntimeException("File '" + fileName +
                    "' can not be opened to perform tree logging (TreeDebugLogger).");
        }

    }

    @Override
    public void log(long state) {
        //do not log to file as a regular logger in the XML file
        //only write to file at dump_every x, by calling writeToFile()
        //hence do nothing here
    }

    @Override
    public void startLogging() {
        //do nothing
    }

    @Override
    public void stopLogging() {
        //do nothing
    }

    public String getFileName() {
        return this.fileName;
    }

    private class DebugTabDelimitedFormatter implements LogFormatter {
        private PrintStream out;

        public DebugTabDelimitedFormatter(PrintStream out) {
            this.out = out;
        }

        public void setPrintStream(PrintStream set) {
            this.out = set;
        }

        @Override
        public void startLogging(String title) {
            //do nothing
        }

        @Override
        public void logHeading(String heading) {
            throw new RuntimeException("TreeDebugLogger: logHeading(String heading) not implemented.");
        }

        @Override
        public void logLine(String line) {
            this.out.println(line);
        }

        @Override
        public void logLabels(String[] labels) {
            throw new RuntimeException("TreeDebugLogger: logLabels(String[] labels) not implemented.");
        }

        @Override
        public void logValues(String[] values) {
            throw new RuntimeException("TreeDebugLogger: logValues(String[] values) not implemented.");
        }

        @Override
        public void stopLogging() {
            try {
                out.close();
                fileOut.close();
            } catch (IOException ioe) {
                throw new RuntimeException(ioe);
            }
        }

    }

}
