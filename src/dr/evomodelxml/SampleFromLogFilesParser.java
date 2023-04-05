/*
 * SampleFromLogFilesParser.java
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

package dr.evomodelxml;

import dr.evomodel.SampleFromLogFiles;
import dr.evomodel.tree.EmpiricalTreeDistributionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeTraceAnalysis;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.Logger;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.xml.*;

import java.io.*;
import java.util.*;

import static dr.inferencexml.loggers.LoggerParser.getLogFile;

/**
 * @author Marc Suchard
 */
public class SampleFromLogFilesParser extends AbstractXMLObjectParser {

    private final static String PARSER_NAME = "sampleFromLogFiles";
    private final static String SAMPLE_BLOCK = "sample";
    private final static String EXECUTE = "execute";
    private final static String FILE_NAME = "fileName";
    private final static String COLUMN_NAME = "columnName";
    private final static String FIRST_SAMPLE = "firstSample";
    private final static String LAST_SAMPLE = "lastSample";
    private final static String NUMBER_SAMPLES = "numberOfSamples";
    private final static String CHECK_BLOCK = "check";
    private final static String TOLERANCE = "tolerance";

    public String getParserName() {
        return PARSER_NAME;
    }

    private final Map<String, LogFileTraces> logCache = new HashMap<>();

    private LogFileTraces getCachedLogFile(String fileName) throws XMLParseException {

        if (logCache.containsKey(fileName)) {
            return logCache.get(fileName);
        } else {
            File file = new File(fileName);
            LogFileTraces traces = new LogFileTraces(fileName, file);
            traces.setBurnIn(0);
            try {
                traces.loadTraces();
            } catch (TraceException | IOException e) {
                throw new XMLParseException(e.getMessage());
            }

            logCache.put(fileName, traces);
            return traces;
        }
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxoExecute = xo.getChild(EXECUTE);
        List<Loggable> loggable = new ArrayList<>();
        ArrayList<Logger> loggers = new ArrayList<>();

        for (int i = 0; i < cxoExecute.getChildCount(); ++i) {
            Object obj = cxoExecute.getChild(i);
            if (obj instanceof Loggable) {
                loggable.add((Loggable) obj);
            }

            if (obj instanceof Logger) {
                loggers.add((Logger) obj);
            }
        }

        final PrintWriter pw = getLogFile(xo, getParserName());

        boolean printStatus = xo.hasAttribute(FILE_NAME);

        SampleFromLogFiles action = new SampleFromLogFiles(loggable, loggers, pw, printStatus);
        
        for (XMLObject cxo : xo.getAllChildren(SAMPLE_BLOCK)) {

            EmpiricalTreeDistributionModel treeModel = (EmpiricalTreeDistributionModel) cxo.getChild(TreeModel.class);
            Parameter parameter = (Parameter) cxo.getChild(Parameter.class);

            if (treeModel != null) {

                action.addTreeBinding(new SampleFromLogFiles.TreeBinding(treeModel));

            } else if (parameter != null) {

                LogIndex logIndex = parseBlock(cxo);
                action.addParameterBinding(new SampleFromLogFiles.ParameterBinding(
                        parameter, logIndex.traces, logIndex.traceIndexParameter));
            }
        }

        for (XMLObject cxo : xo.getAllChildren(CHECK_BLOCK)) {

            Likelihood likelihood = (Likelihood) cxo.getChild(Likelihood.class);
            LogIndex logIndex = parseBlock(cxo);
            double tolerance = cxo.getAttribute(TOLERANCE, 1E-6);

            action.addCheckBinding(new SampleFromLogFiles.CheckBinding(likelihood, logIndex.traces,
                    logIndex.traceIndexParameter, tolerance));
        }

        long firstSample = xo.getAttribute(FIRST_SAMPLE, -1L);
        long lastSample = xo.getAttribute(LAST_SAMPLE, -1L);
        int numberSamples = xo.getAttribute(NUMBER_SAMPLES, 100);

        action.run(firstSample, lastSample, numberSamples);

        return null;
    }

    class LogIndex {
        LogFileTraces traces;
        int traceIndexParameter;

        LogIndex(LogFileTraces traces, int traceIndexParameter) {
            this.traces = traces;
            this.traceIndexParameter = traceIndexParameter;
        }
    }

    private LogIndex parseBlock(XMLObject cxo) throws XMLParseException {

        String fileName = cxo.getStringAttribute(FILE_NAME);
        String columnName = cxo.getStringAttribute(COLUMN_NAME);
        java.util.logging.Logger.getLogger("dr.evomodelxml").info("Reading " + columnName
                + " from " + fileName);

        LogFileTraces traces = getCachedLogFile(fileName);

        int traceIndexParameter = -1;
        for (int j = 0; j < traces.getTraceCount(); j++) {
            String traceName = traces.getTraceName(j);
            if (traceName.trim().equals(columnName)) {
                traceIndexParameter = j;
            }
        }

        if (traceIndexParameter == -1) {
            throw new XMLParseException("Column '" + columnName + "' can not be found in " + fileName);
        }

        return new LogIndex(traces, traceIndexParameter);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Analyses and reports on a trace consisting of trees.";
    }

    public Class getReturnType() {
        return TreeTraceAnalysis.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newLongIntegerRule(FIRST_SAMPLE, true),
            AttributeRule.newLongIntegerRule(LAST_SAMPLE, true),
            AttributeRule.newIntegerRule(NUMBER_SAMPLES, true),
            AttributeRule.newStringRule(FILE_NAME, true),
            new ElementRule(SAMPLE_BLOCK,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(EmpiricalTreeDistributionModel.class),
                                    new ElementRule(Parameter.class)),
                            new StringAttributeRule(FILE_NAME, "File name", true),
                            new StringAttributeRule(COLUMN_NAME, "Column name", true),
                    }, 1, Integer.MAX_VALUE),
            new ElementRule(EXECUTE, new XMLSyntaxRule[] {
                    new ElementRule(Loggable.class, 1, Integer.MAX_VALUE),
            }),
            new ElementRule(CHECK_BLOCK,
                    new XMLSyntaxRule[]{
                            new ElementRule(Likelihood.class),
                            new StringAttributeRule(FILE_NAME, "File name"),
                            new StringAttributeRule(COLUMN_NAME, "Likelihood column name"),
                            AttributeRule.newDoubleRule(TOLERANCE, true),
                    }, 0, Integer.MAX_VALUE),
    };
}