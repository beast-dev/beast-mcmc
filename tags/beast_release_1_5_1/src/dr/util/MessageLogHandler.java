/*
 * MessageLogHandler.java
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

package dr.util;

import java.util.logging.*;

public class MessageLogHandler extends StreamHandler {

	public MessageLogHandler() {
		setOutputStream(System.out);
		setFormatter(new MessageLogFormatter());
	}


	public void publish(LogRecord record) {
		super.publish(record);
		flush();
	}

	public void close() {
		flush();
	}

	private class MessageLogFormatter extends Formatter {
				
		// Line separator string.  This is the value of the line.separator
		// property at the moment that the SimpleFormatter was created.
		private String lineSeparator = (String) java.security.AccessController.doPrivileged(
		        new sun.security.action.GetPropertyAction("line.separator"));

		/**
		 * Format the given LogRecord.
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
}

