/*
 * CitationLogHandler.java
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

import java.io.OutputStream;
import java.util.logging.Formatter;
import java.util.logging.LogRecord;
import java.util.logging.StreamHandler;

public class CitationLogHandler extends StreamHandler {
    private static final CitationLogHandler INSTANCE = new CitationLogHandler();
    private final MessageLogFormatter citationFormatter = new MessageLogFormatter();

    public static CitationLogHandler getHandler(OutputStream stream) {
        INSTANCE.setFormatter(INSTANCE.citationFormatter);
        INSTANCE.setOutputStream(stream);
        return INSTANCE;
    }

    public static void closeHandler() {
        INSTANCE.flush();
        INSTANCE.close();
    }

    public void publish(LogRecord record) {
        super.publish(record);
        INSTANCE.flush();
        //System.err.println("published: " + record.getMessage());
    }

    private class MessageLogFormatter extends Formatter {

        // Line separator string. This is the value of the line.separator
        // property at the moment that the SimpleFormatter was created.
        private final String lineSeparator = System.getProperty("line.separator");

        // AR - is there a reason why this was used? It causes warnings at compile
//        private final String lineSeparator = (String) java.security.AccessController.doPrivileged(
//                new sun.security.action.GetPropertyAction("line.separator"));

        /**
         * Format the given LogRecord.
         * @param record the log record to be formatted.
         * @return a formatted log record
         */
        public synchronized String format(LogRecord record) {
            final StringBuffer sb = new StringBuffer();
            sb.append(formatMessage(record));
            sb.append(lineSeparator);
            return sb.toString();
        }
    }

}
