/*
 * PosteriorPredictiveSequenceLogger.java
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

package dr.evomodel.treedatalikelihood.pps;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;

/**
 * Logs posterior predictive sequence datasets
 */
public class PosteriorPredictiveSequenceLogger extends MCLogger {

    public static final String TAXA_ID = "taxa";

    private final PredictiveDataGenerator generator;
    private final PosteriorPredictiveXMLWriter xmlWriter;
    private int replicateCount = 0;

    public PosteriorPredictiveSequenceLogger(PredictiveDataGenerator generator, LogFormatter formatter, long logEvery) {
        super(formatter, logEvery, false);
        this.generator = generator;
        this.xmlWriter = new PosteriorPredictiveXMLWriter();
    }

    @Override
    public void startLogging() {
        writeBlock(xmlWriter.writeRootOpen());
        writeBlock(xmlWriter.writeTaxaBlock(generator.getTaxa(), TAXA_ID));

        DataType dataType = generator.getDataType();
        if (dataType instanceof GeneralDataType) {
            GeneralDataType generalDataType = (GeneralDataType) dataType;
            writeBlock(xmlWriter.writeGeneralDataTypeBlock(generalDataType, generalDataType.getId()));
        }
    }

    @Override
    public void log(long state) {
        SimpleAlignment alignment = generator.simulate();

        replicateCount++;
        String id = generator.getParentId() + "_" + replicateCount;

        writeBlock(xmlWriter.writeAlignmentBlock(alignment, id, state, generator.getParentId()));
    }

    @Override
    public void stopLogging() {
        writeBlock(xmlWriter.writeRootClose());
        super.stopLogging();
    }

    private void writeBlock(String xmlFragment) {
        String lineSeparator = System.lineSeparator();
        if (xmlFragment.endsWith(lineSeparator)) {
            xmlFragment = xmlFragment.substring(0, xmlFragment.length() - lineSeparator.length());
        }
        logLine(xmlFragment);
    }
}
