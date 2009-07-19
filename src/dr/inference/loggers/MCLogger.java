/*
 * MCLogger.java
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

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * A class for a general purpose logger.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: MCLogger.java,v 1.18 2005/05/24 20:25:59 rambaut Exp $
 */
public class MCLogger implements Logger {
    
	// the file that this logger is logging to or 
	// null if the LogFormatter is logging to a non-file print stream
	String fileName;
	
	/**
     * Output performance stats in this log
     */
     private boolean performanceReport;

    /**
     * Constructor. Will log every logEvery.
     *
     * @param formatter the formatter of this logger
     * @param logEvery  logging frequency
     */
    public MCLogger(String fileName, LogFormatter formatter, int logEvery, boolean performanceReport) {

        addFormatter(formatter);
        this.logEvery = logEvery;
        this.performanceReport = performanceReport;
        this.fileName = fileName;
    }
    
    /**
     * Constructor. Will log every logEvery.
     *
     * @param formatter the formatter of this logger
     * @param logEvery  logging frequency
     */
    public MCLogger(LogFormatter formatter, int logEvery, boolean performanceReport) {
        this(null, formatter, logEvery, performanceReport);
    }

    public final void setTitle(String title) {
        this.title = title;
    }

    public final String getTitle() {
        return title;
    }

    public int getLogEvery() {
        return logEvery;
    }

    public void setLogEvery(int logEvery) {
        this.logEvery = logEvery;
    }

    public final void addFormatter(LogFormatter formatter) {

        formatters.add(formatter);
    }

    public final void add(Loggable loggable) {

        LogColumn[] columns = loggable.getColumns();

        for (LogColumn column : columns) {
            addColumn(column);
        }
    }

    public final void addColumn(LogColumn column) {

        columns.add(column);
    }

    public final void addColumns(LogColumn[] columns) {

        for (LogColumn column : columns) {
            addColumn(column);
        }
    }

    public final int getColumnCount() {
        return columns.size();
    }

    public final LogColumn getColumn(int index) {

        return columns.get(index);
    }

    public final String getColumnLabel(int index) {

        return columns.get(index).getLabel();
    }

    public final String getColumnFormatted(int index) {

        return columns.get(index).getFormatted();
    }
    
    /**
     * @return file name or null if this logger is logging to a non-file print stream
     */
    public final String getFileName() {
    	return fileName;
    }

    protected void logHeading(String heading) {
        for (LogFormatter formatter : formatters) {
            formatter.logHeading(heading);
        }
    }

    protected void logLine(String line) {
        for (LogFormatter formatter : formatters) {
            formatter.logLine(line);
        }
    }

    protected void logLabels(String[] labels) {
        for (LogFormatter formatter : formatters) {
            formatter.logLabels(labels);
        }
    }

    protected void logValues(String[] values) {
        for (LogFormatter formatter : formatters) {
            formatter.logValues(values);
        }
    }

    public void startLogging() {

        for (LogFormatter formatter : formatters) {
            formatter.startLogging(title);
        }

        if (title != null) {
            logHeading(title);
        }

        if (logEvery > 0) {
            final int columnCount = getColumnCount();
            String[] labels = new String[columnCount + 1];

            labels[0] = "state";

            for (int i = 0; i < columnCount; i++) {
                labels[i + 1] = getColumnLabel(i);
            }

            logLabels(labels);
        }
    }

    public void log(int state) {

        if (state == 0) {
            startTime = System.currentTimeMillis();
            formatter.setMaximumFractionDigits(2);
        }

        if (logEvery > 0 && (state % logEvery == 0)) {

            final int columnCount = getColumnCount();

            String[] values = new String[columnCount + (performanceReport ? 2 : 1)];

            values[0] = Integer.toString(state);

            for (int i = 0; i < columnCount; i++) {
                values[i + 1] = getColumnFormatted(i);
            }

            if( performanceReport ) {
                if (state > 0) {

                    long timeTaken = System.currentTimeMillis() - startTime;

                    double hoursPerMillionStates = (double) timeTaken / (3.6 * (double) state);

                    values[columnCount + 1] = formatter.format(hoursPerMillionStates) + " hours/million states";
                } else {
                    values[columnCount + 1] = "-";
                }
            }

            logValues(values);
        }
    }

    public void stopLogging() {

        for (LogFormatter formatter : formatters) {
            formatter.stopLogging();
        }
    }

    private String title = null;

    private ArrayList<LogColumn> columns = new ArrayList<LogColumn>();

    protected int logEvery = 0;

    public List<LogFormatter> getFormatters() {
        return formatters;
    }

    public void setFormatters(List<LogFormatter> formatters) {
        this.formatters = formatters;
    }

    protected List<LogFormatter> formatters = new ArrayList<LogFormatter>();

    private long startTime;
    private NumberFormat formatter = NumberFormat.getNumberInstance();

}
