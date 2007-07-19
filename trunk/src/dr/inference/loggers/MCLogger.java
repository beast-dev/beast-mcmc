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

import dr.util.Identifiable;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * A class for a general purpose logger.
 *
 * @version $Id: MCLogger.java,v 1.18 2005/05/24 20:25:59 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class MCLogger implements Logger {

	public static final String LOG = "log";
	public static final String ECHO = "echo";
	public static final String ECHO_EVERY = "echoEvery";	
	public static final String TITLE = "title";
	public static final String FILE_NAME = "fileName";
    public static final String FORMAT = "format";
	public static final String TAB = "tab";	
	public static final String HTML = "html";	
	public static final String PRETTY = "pretty";	
	public static final String LOG_EVERY = "logEvery";	

	public static final String COLUMNS = "columns";
	public static final String COLUMN = "column";
	public static final String LABEL = "label";
	public static final String SIGNIFICANT_FIGURES = "sf";
	public static final String DECIMAL_PLACES = "dp";
	public static final String WIDTH = "width";
	
	/**
	 * Constructor. Will log every logEvery.
     * @param formatter
     * @param logEvery  logging frequency
     */
	public MCLogger(LogFormatter formatter, int logEvery) {
		
		addFormatter(formatter);
		this.logEvery = logEvery;
	}

	public final void setTitle(String title) {
	
		this.title = title;
	}

	public final String getTitle() {
	
		return title;
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

	public final int getColumnCount() { return columns.size(); }
	
	public final LogColumn getColumn(int index) { 
	
		return columns.get(index);
	}
	
	public final String getColumnLabel(int index) { 
	
		return columns.get(index).getLabel();
	}
	
	public final String getColumnFormatted(int index) { 
	
		return columns.get(index).getFormatted();
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
				labels[i+1] = getColumnLabel(i);
			}
			
			logLabels(labels);
		}
	}
		
	public void log(int state) {
	
		if (logEvery > 0 && (state % logEvery == 0)) {

            final int columnCount = getColumnCount();
            String[] values = new String[columnCount + 1];
			
			values[0] = Integer.toString(state);
			
			for (int i = 0; i < columnCount; i++) {
				values[i+1] = getColumnFormatted(i);
			}
			
			logValues(values);
		}
	}
	
	public void stopLogging() {

        for (LogFormatter formatter : formatters) {
            formatter.stopLogging();
        }
    }

    // directory where beast xml file resides
    private static File masterBeastDirectory = null;

    public static void setMasterDir(File fileName) {
        masterBeastDirectory = fileName;
    }

    /**
     * Allow a file relative to beast xml file with a prefix of ./
     * @param xo log element
     * @param parserName for error messages
     * @return logger object
     * @throws XMLParseException if file can't be created for some reason
     */
    protected static PrintWriter getLogFile(XMLObject xo, String parserName) throws XMLParseException {
        if (xo.hasAttribute(FILE_NAME)) {

            String fileName = xo.getStringAttribute(FILE_NAME);
            
            final boolean localFile = fileName.startsWith("./");
            final boolean relative = masterBeastDirectory != null && localFile;
            if( localFile ) {
                fileName = fileName.substring(2);
            }

            final File file = new File(fileName);
            final String name = file.getName();
            String parent = file.getParent();

            if (!file.isAbsolute()) {
                if( relative ) {
                    parent = masterBeastDirectory.getAbsolutePath();
                } else {
                    parent = System.getProperty("user.dir");
                }
            }
            final File logFile = new File(parent, name);

//	         System.out.println("Writing log file to "+parent+System.getProperty("path.separator")+name);
            try {
                return new PrintWriter(new FileOutputStream(logFile));
            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + logFile.getAbsolutePath() + "' can not be opened for " + parserName + " element.");
            }

        }
        return new PrintWriter(System.out);
    }

    public static LoggerParser PARSER = new LoggerParser();
	
	public static class LoggerParser extends AbstractXMLObjectParser {


        public String getParserName() { return LOG; }

        /**
         * @return an object based on the XML element it was passed.
         */
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
	
			// You must say how often you want to log
			int logEvery = xo.getIntegerAttribute(LOG_EVERY);
			
			PrintWriter pw = getLogFile(xo, getParserName());

            LogFormatter formatter = new TabDelimitedFormatter(pw);
			
			MCLogger logger = new MCLogger(formatter, logEvery);
								
			if (xo.hasAttribute(TITLE)) {
				logger.setTitle(xo.getStringAttribute(TITLE)); 			
			}
			
			for (int i =0; i < xo.getChildCount(); i++) {
			
				Object child = xo.getChild(i);
					
				if (child instanceof Columns) {
				
					logger.addColumns(((Columns)child).getColumns());
					
				} else if (child instanceof Loggable) {
				
					logger.add((Loggable)child);
					
				} else if (child instanceof Identifiable) {

					logger.addColumn(new LogColumn.Default(((Identifiable)child).getId(), child));
					
				} else {

					logger.addColumn(new LogColumn.Default(child.getClass().toString(), child));
				}
			}

			return logger;
		}


        //************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
	
		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			AttributeRule.newIntegerRule(LOG_EVERY),
			new StringAttributeRule(FILE_NAME, 
				"The name of the file to send log output to. " + 
				"If no file name is specified then log is sent to standard output", true),
			new StringAttributeRule(TITLE,
				 "The title of the log", true),
			new OrRule( 
				new XMLSyntaxRule[] {
					new ElementRule(Columns.class, 1, Integer.MAX_VALUE),
					new ElementRule(Loggable.class, 1, Integer.MAX_VALUE),
					new ElementRule(Object.class, 1, Integer.MAX_VALUE)	
				}
			)
		};
		
		public String getParserDescription() { 
			return "Logs one or more items at a given frequency to the screen or to a file";
		}
	
		public Class getReturnType() { return MLLogger.class; }
		
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
}
