/*
 * OutputFileHelper.java
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

import java.io.File;
import java.io.IOException;

/**
 * Resolves where output files (logs, trees, checkpoints, citations) are written.
 * <p>
 * Two directories are tracked:
 * <ul>
 * <li>An <em>explicit</em> output directory, requested by the user with the
 * {@code -output} command line option (stored in the {@value #OUTPUT_DIR_PROPERTY}
 * system property) or by {@link #setOutputDirectory(File)}. Every relative output
 * file name is placed in this directory.</li>
 * <li>A <em>default</em> output directory, inferred from the first log file the XML
 * parser resolves. It is only used by outputs that must follow the log files (the
 * checkpoint file) when no explicit directory was given.</li>
 * </ul>
 * The explicit directory always wins over the inferred one.
 *
 * @author BEAST Development Team
 */
public class OutputFileHelper {

    public static final String OUTPUT_DIR_PROPERTY = "beast.output.dir";

    private static File outputDirectory = null;
    private static File defaultOutputDirectory = null;

    /**
     * Set the explicit output directory. Relative output file names will be placed here.
     *
     * @param directory the directory where output files should be written, or null to unset
     */
    public static synchronized void setOutputDirectory(File directory) {
        outputDirectory = directory;
    }

    /**
     * Record the directory of a resolved log file as the default output directory.
     * Only the first call has an effect so that all outputs follow the first log file.
     *
     * @param directory the directory the log file was resolved to
     */
    public static synchronized void setDefaultOutputDirectory(File directory) {
        if (defaultOutputDirectory == null) {
            defaultOutputDirectory = directory;
        }
    }

    /**
     * @return the explicitly requested output directory (set programmatically or via the
     * {@value #OUTPUT_DIR_PROPERTY} system property), or null if none was requested
     */
    public static synchronized File getExplicitOutputDirectory() {
        if (outputDirectory != null) {
            return outputDirectory;
        }

        String outputDirProperty = System.getProperty(OUTPUT_DIR_PROPERTY);
        if (outputDirProperty != null && !outputDirProperty.trim().isEmpty()) {
            return new File(outputDirProperty);
        }

        return null;
    }

    /**
     * @return the directory output files should be written to: the explicit directory if one
     * was requested, otherwise the directory inferred from the log files, otherwise null
     */
    public static synchronized File getOutputDirectory() {
        File explicit = getExplicitOutputDirectory();
        if (explicit != null) {
            return explicit;
        }
        return defaultOutputDirectory;
    }

    /**
     * Clear both the explicit and the inferred output directory.
     */
    public static synchronized void clearOutputDirectory() {
        outputDirectory = null;
        defaultOutputDirectory = null;
    }

    /**
     * Resolve an output file name against an output directory.
     * <p>
     * An absolute file name is returned as is. If {@code directory} is null the name is
     * resolved exactly as {@link FileHelpers#getFile(String, String)} would. Otherwise a
     * leading {@code ./} is stripped, the optional prefix is applied to the file name and
     * the result (including any relative sub-directory) is placed inside {@code directory}.
     *
     * @param fileName  an absolute or relative file name
     * @param prefix    an optional prefix for the file name (may be null)
     * @param directory the output directory to resolve relative names against (may be null)
     * @return the resolved file
     */
    public static File resolveOutputFile(String fileName, String prefix, File directory) {
        if (directory == null || new File(fileName).isAbsolute()) {
            return FileHelpers.getFile(fileName, prefix);
        }

        if (fileName.startsWith("./")) {
            fileName = fileName.substring(2);
        }

        File file = new File(fileName);
        String name = (prefix != null ? prefix : "") + file.getName();
        String parent = file.getParent();

        File parentDir = (parent != null && !parent.isEmpty()) ? new File(directory, parent) : directory;
        return new File(parentDir, name);
    }

    /**
     * Check that a directory exists, is a directory and can have files created in it.
     * The writability check creates and deletes a temporary file because
     * {@link File#canWrite()} only inspects permission bits and misses read-only
     * volumes and platform ACLs.
     *
     * @param directory the directory to check
     * @throws IOException with a user-readable message if the directory is unusable
     */
    public static void validateOutputDirectory(File directory) throws IOException {
        if (!directory.exists()) {
            throw new IOException("Output directory does not exist: " + directory.getAbsolutePath());
        }
        if (!directory.isDirectory()) {
            throw new IOException("Output path is not a directory: " + directory.getAbsolutePath());
        }

        File probe = null;
        try {
            probe = File.createTempFile("beast_output_check", ".tmp", directory);
        } catch (IOException | SecurityException e) {
            throw new IOException("Output directory is not writable: " + directory.getAbsolutePath());
        } finally {
            if (probe != null && !probe.delete()) {
                probe.deleteOnExit();
            }
        }
    }
}
