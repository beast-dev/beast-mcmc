/*
 * ResumedLogWriter.java
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

package dr.inference.loggers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Writer;

/**
 * Writes to the end of an existing log file of a chain that is resumed from a saved state.
 * What the file already has must not be written again: its header, and the state the chain is
 * resumed from. While the thread that is logging has asked for it (see setSuppressed), what is
 * written is discarded. Writers of other files are not affected.
 */
public final class ResumedLogWriter {

    private static final ThreadLocal<Boolean> SUPPRESSED = new ThreadLocal<Boolean>();

    private ResumedLogWriter() {
    }

    /**
     * @return a writer that appends to the file
     */
    public static PrintWriter open(File file) throws IOException {
        final Writer out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file, true)));
        return new PrintWriter(new Writer() {
            @Override
            public void write(char[] buffer, int offset, int length) throws IOException {
                if (!isSuppressed()) {
                    out.write(buffer, offset, length);
                }
            }

            @Override
            public void flush() throws IOException {
                out.flush();
            }

            @Override
            public void close() throws IOException {
                out.close();
            }
        });
    }

    /**
     * Turns the discarding of what this thread writes to resumed log files on or off. On while
     * the loggers write their headers and while they log the state the chain is resumed from.
     */
    public static void setSuppressed(boolean suppressed) {
        if (suppressed) {
            SUPPRESSED.set(Boolean.TRUE);
        } else {
            SUPPRESSED.remove();
        }
    }

    private static boolean isSuppressed() {
        return SUPPRESSED.get() != null;
    }
}
