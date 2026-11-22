/*
 * NexusExporter.java
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

package dr.app.tools;

import java.io.IOException;
import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.List;
import java.util.Map;

import dr.evolution.alignment.Alignment;
import dr.evolution.tree.Tree;

/**
 * Deprecated: Use {@link dr.evolution.io.NexusExporter} instead.
 * This class is kept for backward compatibility only.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @deprecated Use {@link dr.evolution.io.NexusExporter}
 */
@Deprecated
public class NexusExporter extends dr.evolution.io.NexusExporter {

    public NexusExporter(PrintStream out) {
        super(out);
    }

    public NexusExporter(PrintStream out, boolean useNumericalTaxa) {
        super(out, useNumericalTaxa);
    }

    public NexusExporter(PrintStream out, AttributeType writeAttributesAs) {
        super(out, writeAttributesAs);
    }

    public NexusExporter(PrintStream out, AttributeType writeAttributesAs, boolean useNumericalTaxa) {
        super(out, writeAttributesAs, useNumericalTaxa);
    }

}// END: class
