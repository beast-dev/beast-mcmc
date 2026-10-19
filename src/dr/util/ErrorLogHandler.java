/*
 * ErrorLogHandler.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.util;

import java.util.concurrent.atomic.AtomicInteger;
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
            // each chain of a multi-chain run counts its own errors
            AtomicInteger counter = threadErrorCount.get();
            int errorCount = (counter != null ? counter : this.errorCount).incrementAndGet();

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
        return errorCount.get();
    }

    /**
     * Counts the errors logged from the current thread with the given counter instead of the
     * handler's own. Pass null to go back to the handler's counter.
     */
    public static void setThreadErrorCounter(AtomicInteger counter) {
        if (counter == null) {
            threadErrorCount.remove();
        } else {
            threadErrorCount.set(counter);
        }
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
    private final AtomicInteger errorCount = new AtomicInteger(0);

    private static final ThreadLocal<AtomicInteger> threadErrorCount = new ThreadLocal<AtomicInteger>();
}