/*
 * ErrorLogHandler.java
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

import java.util.logging.Formatter;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.StreamHandler;

public class ErrorLogHandler extends StreamHandler {

    public ErrorLogHandler(int maxErrorCount) {
        setOutputStream(System.err);
        setFormatter(new MessageLogFormatter());

        this.maxErrorCount = maxErrorCount;
    }


    public void publish(LogRecord record) {
        super.publish(record);
        flush();

        if (record.getLevel() == Level.SEVERE) {
            errorCount++;

            if (errorCount > maxErrorCount) {
                if (errorCount > 1) {
                    throw new RuntimeException("ErrorLog: Maximum number of errors (" + (maxErrorCount + 1) + ") reached. Terminating BEAST");
                } else {
                    throw new RuntimeException("An error was encountered. Terminating BEAST");
                }
            }
        }
    }

    public void close() {
        flush();
    }

    public int getErrorCount() {
        return errorCount;
    }

    private class MessageLogFormatter extends Formatter {

        // Line separator string.  This is the value of the line.separator
        // property at the moment that the SimpleFormatter was created.
        private final String lineSeparator = System.getProperty("line.separator");

        // AR - is there a reason why this was used? It causes warnings at compile
//        private final String lineSeparator = (String) java.security.AccessController.doPrivileged(
//                new sun.security.action.GetPropertyAction("line.separator"));

        /**
         * Format the given LogRecord.
         *
         * @param record the log record to be formatted.
         * @return a formatted log record
         */
        public synchronized String format(LogRecord record) {
            StringBuffer sb = new StringBuffer();
            String message = formatMessage(record);
            sb.append(message);
            sb.append(lineSeparator);
            return sb.toString();
        }
    }


    private final int maxErrorCount;
    private int errorCount = 0;
}