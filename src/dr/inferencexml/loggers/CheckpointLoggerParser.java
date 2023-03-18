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
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Guy Baele
 */

public class CheckpointLoggerParser extends AbstractXMLObjectParser {

    public static final String LOG_CHECKPOINT = "logCheckpoint";
    public static final String CHECKPOINT_EVERY = "checkpointEvery";
    public static final String FILE_NAME = FileHelpers.FILE_NAME;

    /**
     * @return an object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        //TODO actually write some code
        return null;
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

    };

    public String getParserDescription() {
        return "Logs the state of the Markov chain at a given frequency to a checkpoint file";
    }

    public Class getReturnType() {
        //TODO is this right?
        return BeastCheckpointer.class;
    }

}
