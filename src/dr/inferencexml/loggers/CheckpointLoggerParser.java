/*
 * LoggerParser.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.loggers;

import dr.app.checkpoint.BeastCheckpointer;
import dr.util.FileHelpers;
import dr.xml.*;

/**
 * @author Guy Baele
 */

public class CheckpointLoggerParser extends AbstractXMLObjectParser {

    public static final String LOG_CHECKPOINT = "logCheckpoint";
    public static final String CHECKPOINT_EVERY = "checkpointEvery";
    public static final String CHECKPOINT_FINAL = "checkpointFinal";
    public static final String FILE_NAME = FileHelpers.FILE_NAME;
    public static final String ALLOW_OVERWRITE_LOG = "overwrite";

    /**
     * @return an object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final String fileName = xo.getStringAttribute(FILE_NAME);

        final int checkpointEvery = xo.getIntegerAttribute(CHECKPOINT_EVERY);

        final int checkpointFinal = xo.getIntegerAttribute(CHECKPOINT_FINAL);

        final boolean overwrite = xo.getBooleanAttribute(ALLOW_OVERWRITE_LOG);

        return BeastCheckpointer.getInstance(fileName, checkpointEvery, checkpointFinal, overwrite);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserName() {
        return LOG_CHECKPOINT;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newLongIntegerRule(CHECKPOINT_EVERY),
            AttributeRule.newLongIntegerRule(CHECKPOINT_FINAL),
            AttributeRule.newBooleanRule(ALLOW_OVERWRITE_LOG, true),
            new StringAttributeRule(FILE_NAME, "The name of the file to checkpoint to.", false)
    };

    public String getParserDescription() {
        return "Logs the state of the Markov chain at a given frequency to a checkpoint file";
    }

    public Class getReturnType() {
        return BeastCheckpointer.class;
    }

}
