/*
 * FileHelpers.java
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

import java.io.*;

/**
 * @author Joseph Heled
 *         Date: 15/04/2008
 */
public class FileHelpers {

    public static final String FILE_NAME = "fileName";

    /**
     * @param file the file to read the line numbers from
     * @return Number of lines in file
     * @throws IOException low level file error
     */
    public static int numberOfLines(File file) throws IOException {
        RandomAccessFile randFile = new RandomAccessFile(file, "r");
        long lastRec = randFile.length();
        randFile.close();
        FileReader fileRead = new FileReader(file);
        LineNumberReader lineRead = new LineNumberReader(fileRead);
        lineRead.skip(lastRec);
        int count = lineRead.getLineNumber() - 1;
        fileRead.close();
        lineRead.close();
        return count;
    }

    /**
     * Resolve file from name.
     * <p/>
     * Keep A fully qualified (i.e. absolute path) as is. A name starting with a "./" is
     * relative to the master directory (set by FileHelpers.setMasterDir).
     * Any other name is stripped of any directory
     * component and placed in the "user.dir" directory.
     *
     * @param fileName an absolute or relative file name
     * @return a File object resolved from provided file name
     */
    public static File getFile(String fileName, String prefix) {
        final boolean localFile = fileName.startsWith("./");
        final boolean relative = masterDirectory != null && localFile;
        if (localFile) {
            fileName = fileName.substring(2);
        }

        if (prefix != null) {
            fileName = prefix + fileName;
        }

        final File file = new File(fileName);
        final String name = file.getName();
        String parent = file.getParent();

        if (!file.isAbsolute()) {
            String p;
            if (relative) {
                p = masterDirectory.getAbsolutePath();
            } else {
                p = System.getProperty("user.dir");
            }
            if (parent != null && parent.length() > 0) {
                parent = p + '/' + parent;
            } else {
                parent = p;
            }
        }
        return new File(parent, name);
    }

    public static File getFile(String fileName) {
        return getFile(fileName, null);
    }

    // directory where beast xml file resides
    private static File masterDirectory = null;

    public static void setMasterDir(File fileName) {
        masterDirectory = fileName;
    }

    // Zstandard compression: output files are compressed when their name ends in ".zst";
    // input files are decompressed when they start with the Zstandard magic number.

    public static final String ZSTD_EXTENSION = ".zst";

    private static final String ZSTD_LEVEL_PROPERTY = "beast.zstd.level";
    private static final int DEFAULT_ZSTD_LEVEL = 3;
    private static final String ZSTD_FLUSH_PROPERTY = "beast.zstd.flush.seconds";
    private static final long DEFAULT_ZSTD_FLUSH_SECONDS = 10;
    private static final byte[] ZSTD_MAGIC = {(byte) 0x28, (byte) 0xB5, (byte) 0x2F, (byte) 0xFD};

    public static boolean isZstdFileName(File file) {
        return file.getName().toLowerCase().endsWith(ZSTD_EXTENSION);
    }

    /**
     * Opens a file for writing, compressing with Zstandard if its name ends in ".zst".
     */
    public static OutputStream getOutputStream(File file) throws IOException {
        OutputStream out = new FileOutputStream(file);
        if (isZstdFileName(file)) {
            int level = Integer.getInteger(ZSTD_LEVEL_PROPERTY, DEFAULT_ZSTD_LEVEL);
            try {
                out = ZstdStreams.compress(out, level);
            } catch (IOException | RuntimeException | LinkageError e) {
                out.close();
                throw new IOException("Unable to open Zstandard compressed file '" + file.getName() + "': " + e.getMessage(), e);
            }
        }
        return out;
    }

    /**
     * Opens a file for writing as a PrintWriter, compressing with Zstandard if its name ends
     * in ".zst". A compressed writer is closed at JVM shutdown if it has not been already, so
     * that the file ends with a complete frame even when a run is interrupted.
     */
    public static PrintWriter getPrintWriter(File file) throws IOException {
        if (!isZstdFileName(file)) {
            return new PrintWriter(new FileOutputStream(file));
        }
        final CompressedPrintWriter writer = new CompressedPrintWriter(
                new OutputStreamWriter(getOutputStream(file)));
        try {
            Runtime.getRuntime().addShutdownHook(new Thread(writer::close));
        } catch (IllegalStateException ise) {
            // already shutting down
        }
        return writer;
    }

    /**
     * Opens a file for reading, transparently decompressing it if it is Zstandard compressed.
     */
    public static InputStream getInputStream(File file) throws IOException {
        BufferedInputStream in = new BufferedInputStream(new FileInputStream(file));
        if (!startsWith(in, ZSTD_MAGIC)) {
            return in;
        }
        try {
            return ZstdStreams.decompress(in);
        } catch (IOException | RuntimeException | LinkageError e) {
            in.close();
            throw new IOException("Unable to read Zstandard compressed file '" + file.getName() + "': " + e.getMessage(), e);
        }
    }

    public static BufferedReader getReader(File file) throws IOException {
        return new BufferedReader(new InputStreamReader(getInputStream(file)));
    }

    public static BufferedReader getReader(String fileName) throws IOException {
        return getReader(new File(fileName));
    }

    private static boolean startsWith(BufferedInputStream in, byte[] magic) throws IOException {
        in.mark(magic.length);
        try {
            for (byte b : magic) {
                if (in.read() != (b & 0xFF)) {
                    return false;
                }
            }
            return true;
        } finally {
            in.reset();
        }
    }

    /**
     * A writer whose stream must be closed (not just flushed) to finish the file. Loggers
     * flush after every line, and every flush ends a compressed block, which costs a lot of
     * ratio when the lines are short; so flushes are passed on at most once per interval.
     */
    public static class CompressedPrintWriter extends PrintWriter {
        private final long flushIntervalNanos;
        private long lastFlush = System.nanoTime();

        CompressedPrintWriter(Writer out) {
            super(out);
            flushIntervalNanos = Long.getLong(ZSTD_FLUSH_PROPERTY, DEFAULT_ZSTD_FLUSH_SECONDS) * 1000000000L;
        }

        @Override
        public void flush() {
            long now = System.nanoTime();
            if (now - lastFlush >= flushIntervalNanos) {
                lastFlush = now;
                super.flush();
            }
        }
    }
}
