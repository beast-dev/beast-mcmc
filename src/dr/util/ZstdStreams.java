/*
 * ZstdStreams.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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

import com.github.luben.zstd.ZstdInputStream;
import com.github.luben.zstd.ZstdOutputStream;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * The only class that touches the zstd-jni library, so that code paths which never
 * see a Zstandard file do not need the library on the classpath. Use the methods in
 * FileHelpers rather than calling this directly.
 *
 * @author Philippe Lemey
 */
final class ZstdStreams {

    private ZstdStreams() {
    }

    static OutputStream compress(OutputStream out, int level) throws IOException {
        // A flush ends the current block but keeps the compression window, so the
        // per-line flushes of the loggers leave a readable file without costing much ratio.
        return new ZstdOutputStream(out, level);
    }

    static InputStream decompress(InputStream in) throws IOException {
        ZstdInputStream zin = new ZstdInputStream(in);
        // read files whose last frame was never closed (e.g., a run that was killed)
        zin.setContinuous(true);
        return zin;
    }
}
