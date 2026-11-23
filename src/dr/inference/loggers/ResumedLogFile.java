/*
 * ResumedLogFile.java
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

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;

/**
 * Helpers for the log files of a chain that is resumed from a saved state: naming a new file
 * for the states that follow, and checking and preparing an existing file to be appended to.
 */
public final class ResumedLogFile {

    private static final String COMPRESSED_EXTENSION = ".zst";
    private static final int BLOCK_SIZE = 1 << 16;
    private static final int PREFIX_LENGTH = 64;

    private ResumedLogFile() {
    }

    /**
     * @return the file for the states logged after resuming from the given state: the name
     * with ".from&lt;state&gt;" before its extension (run.log gives run.from3000.log)
     */
    public static File renamed(File file, long state) {
        String name = file.getName();
        String compressed = "";
        if (name.toLowerCase().endsWith(COMPRESSED_EXTENSION)) {
            compressed = name.substring(name.length() - COMPRESSED_EXTENSION.length());
            name = name.substring(0, name.length() - COMPRESSED_EXTENSION.length());
        }
        final String addition = ".from" + state;
        final int dot = name.lastIndexOf('.');
        if (dot <= 0) {
            name = name + addition;
        } else {
            name = name.substring(0, dot) + addition + name.substring(dot);
        }
        return new File(file.getParentFile(), name + compressed);
    }

    /**
     * @return the last state that a log written every logEvery states holds when the chain
     * has reached the given state
     */
    public static long expectedLastState(long state, long logEvery) {
        return (logEvery > 0 ? (state / logEvery) * logEvery : state);
    }

    /**
     * Checks that an existing log or trees file ends at the state that a saved state continues
     * from, and leaves it ready to be appended to: the closing "End;" of a trees file is
     * removed, and so are blank lines at the end. If the file goes past the saved state it is
     * only cut back to it when force is true.
     *
     * @param file            the log or trees file
     * @param checkpointState the state the chain is resumed from
     * @param logEvery        how often the file is logged to
     * @param force           remove states logged after the saved state rather than fail
     * @return the number of states removed from the end of the file
     * @throws IOException           if the file cannot be read or changed
     * @throws IllegalStateException if the file cannot be appended to, with the reason
     */
    public static int prepareForAppend(File file, long checkpointState, long logEvery, boolean force)
            throws IOException {

        final long expected = expectedLastState(checkpointState, logEvery);

        try (RandomAccessFile raf = new RandomAccessFile(file, "rw")) {
            final long length = raf.length();

            long end = length;      // the end (exclusive) of the lines not yet looked at
            int surplus = 0;        // states found after the expected one
            long lastState = -1;    // the last state in the file
            boolean partial = false;

            if (length > 0) {
                raf.seek(length - 1);
                if (raf.read() != '\n') {
                    // the run was stopped while a line was being written
                    partial = true;
                }
            }

            while (end > 0) {
                // a complete line ends with a line feed, which is not part of its content
                final long contentEnd = (partial && end == length) ? end : end - 1;
                final long start = lineStart(raf, contentEnd);
                final String prefix = readPrefix(raf, start, contentEnd);
                final long state = parseState(prefix);

                if (partial && end == length) {
                    if (!force) {
                        throw new IllegalStateException("The file " + file.getPath() +
                                " ends with an incomplete line, so it cannot safely be appended to.");
                    }
                    if (state >= 0) {
                        surplus++;
                        lastState = state;
                    }
                } else if (state >= 0) {
                    if (lastState < 0) {
                        lastState = state;
                    }
                    if (state == expected) {
                        if (surplus > 0 && !force) {
                            throw ahead(file, lastState, checkpointState, expected);
                        }
                        if (end != length) {
                            raf.setLength(end);
                        }
                        return surplus;
                    }
                    if (state < expected) {
                        if (surplus > 0) {
                            throw new IllegalStateException("The file " + file.getPath() +
                                    " does not have state " + expected + ", which the saved state at " +
                                    checkpointState + " continues from, so it cannot be appended to.");
                        }
                        throw new IllegalStateException("The file " + file.getPath() + " ends at state " +
                                state + " but the saved state is at " + checkpointState + ": the states up to " +
                                expected + " are missing, so it cannot be appended to.");
                    }
                    surplus++;
                } else if (!isIgnorable(prefix)) {
                    // the header: there are no more states above it
                    break;
                }

                end = start;
            }

            if (lastState > expected && !force) {
                throw ahead(file, lastState, checkpointState, expected);
            }
            throw new IllegalStateException("The file " + file.getPath() + " does not have state " + expected +
                    ", which the saved state at " + checkpointState + " continues from, so it cannot be appended to.");
        }
    }

    private static IllegalStateException ahead(File file, long lastState, long checkpointState, long expected) {
        return new IllegalStateException("The file " + file.getPath() + " goes on to state " + lastState +
                " but the saved state is at " + checkpointState + ": appending would repeat the states after " +
                expected + "." +
                "\nUse --force_append to remove those states from the file and append, or resume without" +
                "\n--append to write the new states to new files.");
    }

    /**
     * @return the state of a line of a log file ("3000&lt;tab&gt;...") or of a trees file
     * ("tree STATE_3000 ..."), or -1 if it is not such a line
     */
    static long parseState(String line) {
        int i = 0;
        if (line.startsWith("tree STATE_")) {
            i = "tree STATE_".length();
        }
        final int first = i;
        while (i < line.length() && Character.isDigit(line.charAt(i))) {
            i++;
        }
        if (i == first || i - first > 18) {
            return -1;
        }
        if (first == 0) {
            // a line of a log file: the state is followed by a tab
            if (i >= line.length() || line.charAt(i) != '\t') {
                return -1;
            }
        } else if (i < line.length() && !Character.isWhitespace(line.charAt(i)) && line.charAt(i) != '=' &&
                line.charAt(i) != '[') {
            return -1;
        }
        return Long.parseLong(line.substring(first, i));
    }

    /**
     * @return whether a line at the end of the file can be removed: blank, or the "End;" that
     * closes the trees block
     */
    private static boolean isIgnorable(String prefix) {
        final String trimmed = prefix.trim();
        return trimmed.isEmpty() || trimmed.equalsIgnoreCase("End;");
    }

    /**
     * @return where the line that ends (exclusive) at end starts
     */
    private static long lineStart(RandomAccessFile raf, long end) throws IOException {
        final byte[] block = new byte[BLOCK_SIZE];
        long position = end;
        while (position > 0) {
            final int size = (int) Math.min(BLOCK_SIZE, position);
            position -= size;
            raf.seek(position);
            raf.readFully(block, 0, size);
            for (int i = size - 1; i >= 0; i--) {
                if (block[i] == '\n') {
                    return position + i + 1;
                }
            }
        }
        return 0;
    }

    private static String readPrefix(RandomAccessFile raf, long start, long end) throws IOException {
        final int size = (int) Math.min(PREFIX_LENGTH, Math.max(0, end - start));
        final byte[] bytes = new byte[size];
        raf.seek(start);
        raf.readFully(bytes, 0, size);
        String prefix = new String(bytes, StandardCharsets.ISO_8859_1);
        if (prefix.endsWith("\r")) {
            prefix = prefix.substring(0, prefix.length() - 1);
        }
        return prefix;
    }
}
