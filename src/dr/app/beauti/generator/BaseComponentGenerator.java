/*
 * BaseComponentGenerator.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beauti.generator;

import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.util.XMLWriter;

/**
 * This is an abstract base class for component generators. Currently it simply
 * wraps each insertion with a pair of comments.
 * @author Andrew Rambaut
 * @version $Id$
 */
public abstract class BaseComponentGenerator extends Generator implements ComponentGenerator {
    private Generator callingGenerator;

    protected BaseComponentGenerator(final BeautiOptions options) {
        super(options, null);
    }

    @Override
    public void checkOptions() throws GeneratorException {
        // default is to do nothing
    }

    public void generateAtInsertionPoint(Generator generator, final InsertionPoint point, final Object item, final XMLWriter writer) {
        callingGenerator = generator;
        writer.writeComment("START " + getCommentLabel());
        generate(point, item, "", writer);
        writer.writeComment("END " + getCommentLabel());
        writer.writeBlankLine();
    }

    public void generateAtInsertionPoint(Generator generator, final InsertionPoint point, final Object item, final String prefix, final XMLWriter writer) {
        callingGenerator = generator;
        writer.writeComment("START " + getCommentLabel());
        generate(point, item, prefix, writer);
        writer.writeComment("END " + getCommentLabel());
        writer.writeBlankLine();
    }

    public Generator getCallingGenerator() {
        return callingGenerator;
    }

    protected abstract void generate(final InsertionPoint point, final Object item, final String prefix, final XMLWriter writer);

    protected abstract String getCommentLabel();
}
