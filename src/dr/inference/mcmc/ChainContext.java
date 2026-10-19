/*
 * ChainContext.java
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

package dr.inference.mcmc;

import dr.inference.model.ModelRegistry;
import dr.math.MathUtils;
import dr.util.ErrorLogHandler;

import java.io.*;
import java.nio.file.Files;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * What belongs to one chain of a multi-chain run: its random number stream, the registry of
 * its models and parameters (used, e.g., for checkpointing), its count of numerical errors and
 * how its output files are named. A context is bound to a thread while the chain's model is being
 * built and on the thread that runs the chain. With nothing bound (a plain single-chain
 * run) everything behaves as it always has.
 *
 * @author Andrew Rambaut
 */
public final class ChainContext {

    /**
     * @param name        a label for messages, e.g. "chain 2"
     * @param fileSuffix  inserted before the extension of every output file name ("" for none)
     * @param writesFiles false for the heated members of a tempering set, whose output is
     *                    discarded (the cold member's loggers are moved to whichever chain is
     *                    cold)
     * @param writesScreen false if the chain's screen log should be discarded (only one of
     *                     several independent chains shows its screen log)
     * @param seed        the seed of this chain's random number stream
     */
    public ChainContext(String name, String fileSuffix, boolean writesFiles, boolean writesScreen, long seed) {
        this.name = name;
        this.fileSuffix = fileSuffix;
        this.writesFiles = writesFiles;
        this.writesScreen = writesScreen;
        this.seed = seed;
        this.random = new MathUtils.RandomStream(seed);
    }

    /**
     * Binds this context (its random number stream, registry and error count) to the
     * current thread.
     */
    public void bind() {
        CURRENT.set(this);
        MathUtils.setThreadRandom(random);
        ModelRegistry.setCurrent(registry);
        ErrorLogHandler.setThreadErrorCounter(errorCount);
    }

    /**
     * Removes any context from the current thread.
     */
    public static void unbind() {
        CURRENT.remove();
        MathUtils.setThreadRandom(null);
        ModelRegistry.setCurrent(null);
        ErrorLogHandler.setThreadErrorCounter(null);
    }

    public ModelRegistry getRegistry() {
        return registry;
    }

    /**
     * @return the context bound to the current thread or null if there is none
     */
    public static ChainContext current() {
        return CURRENT.get();
    }

    public String getName() {
        return name;
    }

    public long getSeed() {
        return seed;
    }

    public boolean writesFiles() {
        return writesFiles;
    }

    public boolean writesScreen() {
        return writesScreen;
    }

    public String getFileSuffix() {
        return fileSuffix;
    }

    /**
     * @return the file with this chain's suffix inserted before the extension
     */
    public File addSuffix(File file) {
        if (fileSuffix.isEmpty()) {
            return file;
        }
        return new File(file.getParentFile(), addSuffix(file.getName()));
    }

    /**
     * @return the file name with this chain's suffix inserted before the extension
     */
    public String addSuffix(String fileName) {
        if (fileSuffix.isEmpty()) {
            return fileName;
        }
        int dot = fileName.lastIndexOf('.');
        if (dot <= 0) {
            return fileName + fileSuffix;
        }
        return fileName.substring(0, dot) + fileSuffix + fileName.substring(dot);
    }

    /**
     * @return a file for output that is not wanted. It is in a temporary directory that is
     * removed when the program exits.
     */
    public static synchronized File getDiscardFile(String fileName) {
        try {
            if (discardDirectory == null) {
                discardDirectory = Files.createTempDirectory("beast-discarded").toFile();
                discardDirectory.deleteOnExit();
            }
            // files registered after the directory are deleted before it
            File file = File.createTempFile("discard", "-" + new File(fileName).getName(), discardDirectory);
            file.deleteOnExit();
            return file;
        } catch (IOException e) {
            throw new RuntimeException("Unable to create a temporary file: " + e.getMessage());
        }
    }

    /**
     * @return a writer that discards everything
     */
    public static PrintWriter getDiscardWriter() {
        return new PrintWriter(new Writer() {
            @Override
            public void write(char[] buffer, int offset, int length) {
            }

            @Override
            public void flush() {
            }

            @Override
            public void close() {
            }
        });
    }

    private static final ThreadLocal<ChainContext> CURRENT = new ThreadLocal<ChainContext>();
    private static File discardDirectory = null;

    private final String name;
    private final String fileSuffix;
    private final boolean writesFiles;
    private final boolean writesScreen;
    private final long seed;
    private final MathUtils.RandomStream random;
    private final ModelRegistry registry = new ModelRegistry();
    private final AtomicInteger errorCount = new AtomicInteger(0);
}
