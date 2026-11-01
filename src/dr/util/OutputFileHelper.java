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

/**
 * Helper class to store and retrieve the output file directory
 * This allows output files (logs, checkpoints, etc.) to be written to the same directory
 *
 * @author BEAST Development Team
 */
public class OutputFileHelper {

    private static File outputDirectory = null;
    public static final String OUTPUT_DIR_PROPERTY = "beast.output.dir";

    /**
     * Set the output directory
     * @param directory the directory where output files should be written
     */
    public static synchronized void setOutputDirectory(File directory) {
        outputDirectory = directory;
    }

    /**
     * Get the output directory
     * Checks first for an explicitly set directory, then for the system property
     * @return the output directory, or null if not set
     */
    public static synchronized File getOutputDirectory() {
        if (outputDirectory != null) {
            return outputDirectory;
        }
        
        // Check for system property (set via command line --output option)
        String outputDirProperty = System.getProperty(OUTPUT_DIR_PROPERTY);
        if (outputDirProperty != null && !outputDirProperty.trim().isEmpty()) {
            return new File(outputDirProperty);
        }
        
        return null;
    }

    /**
     * Clear the output directory
     */
    public static synchronized void clearOutputDirectory() {
        outputDirectory = null;
    }
}
