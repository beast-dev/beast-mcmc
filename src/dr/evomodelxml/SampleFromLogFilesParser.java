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

import dr.app.beast.BeastVersion;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.EmpiricalTreeDistributionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeTraceAnalysis;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.model.Parameter;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.math.MathUtils;
import dr.xml.*;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;

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

    public String getParserName() {
        return PARSER_NAME;
    }

    class TreeBinding {

        Tree[] trees;
        EmpiricalTreeDistributionModel treeModel;

        public TreeBinding(EmpiricalTreeDistributionModel treeModel) {
            this.treeModel = treeModel;
            this.trees = treeModel.getTrees();
        }
    }

    class ParameterBinding {
        Parameter parameter;
        LogFileTraces traces;
        int index;

        public ParameterBinding(Parameter parameter, LogFileTraces traces, int index) {
            this.parameter = parameter;
            this.traces = traces;
            this.index = index;
        }
    }

    class Action {

        final private List<TreeBinding> treeBindings = new ArrayList<>();
        final private List<ParameterBinding> parameterBindings = new ArrayList<>();

        private final LogColumn[] log;
        private final TabDelimitedFormatter formatter;
        private final boolean printToScreen;

        public Action(List<Loggable> loggable, PrintWriter printWriter,
                      boolean printStatus) {

            int dim = 0;
            for (Loggable log : loggable) {
                dim += log.getColumns().length;
            }

            log = new LogColumn[dim];

            int index = 0;
            for (Loggable log : loggable) {
                LogColumn[] columns = log.getColumns();
                for (LogColumn column : columns) {
                    this.log[index] = column;
                    ++index;
                }
            }

            formatter = new TabDelimitedFormatter(printWriter);
            this.printToScreen = printStatus;
//            printTitle();
        }

        private void printTitle() {
            final BeastVersion version = new BeastVersion();

            String title = "BEAST " + version.getVersionString() + "\n" +
                    "Generated " + (new Date()) + " [seed=" + MathUtils.getSeed() + "]\n" +
                    System.getProperty("command_line", "");

            formatter.logHeading(title);
        }

        private void setParameter(ParameterBinding binding, int sample) {
            Parameter parameter = binding.parameter;
            LogFileTraces traces = binding.traces;
            int index = binding.index;

            for (int i = 0; i < parameter.getDimension(); ++i) {
                parameter.setParameterValueQuietly(i,
                        traces.getTrace(index + i).getValue(sample));
            }
            parameter.fireParameterChangedEvent();
        }

        private void setTree(TreeBinding binding, int sample) {
            EmpiricalTreeDistributionModel treeModel = binding.treeModel;
            treeModel.setTree(sample);
        }

        private class StateInfo {
            long firstState;
            long maxState;
            long stepSize;
            int stateCount;

            public StateInfo(long firstState, long maxState, long stepSize) {
                this.firstState = firstState;
                this.maxState = maxState;
                this.stepSize = stepSize;
                this.stateCount = (int) ((maxState - firstState) / stepSize) + 1;
            }

            public long getState(int sample) {
                return sample * stepSize + firstState;
            }
        }

        private StateInfo ensureCompatibleLogs() {

            StateInfo info = null;

            if (parameterBindings.size() > 0) {
                ParameterBinding pb = parameterBindings.get(0);
                long maxState = pb.traces.getMaxState();
                long stepSize = pb.traces.getStepSize();
                long firstState = pb.traces.getFirstState();

                info = new StateInfo(firstState, maxState, stepSize);

                for (int i = 1; i < parameterBindings.size(); ++i) {
                    ParameterBinding x = parameterBindings.get(i);
                    if (x.traces.getMaxState() != info.maxState ||
                            x.traces.getStepSize() != info.stepSize ||
                            x.traces.getFirstState() != info.firstState) {
                        throw new RuntimeException("Mis-balanced log file");
                    }
                }
            }

            if (treeBindings.size() > 0) {
                for (TreeBinding tb : treeBindings) {

                    int stateCount = tb.trees.length;

                    String id0 = tb.trees[0].getId();
                    String id1 = tb.trees[1].getId();
                    String idN = tb.trees[stateCount - 1].getId();

                    long state0 = Long.parseLong(id0.replace("STATE_", ""));
                    long state1 = Long.parseLong(id1.replace("STATE_", ""));
                    long stateN = Long.parseLong(idN.replace("STATE_", ""));

                    if (info == null) {
                        info = new StateInfo(state0, stateN, state1 - state0);
                    }

                    if (stateN != info.maxState ||
                            state0 != info.firstState ||
                            (state1 - state0) != info.stepSize ||
                            stateCount != info.stateCount) {
                        throw new RuntimeException("Mis-balanced tree log file");
                    }
                }
            }

            return info;
        }

        public void run(long firstSample,
                        long lastSample,
                        int numberSamples) {

            StateInfo info = ensureCompatibleLogs();

            if (firstSample == -1L) {
                firstSample = info.firstState;
            }

            if (lastSample == -1L) {
                lastSample = info.maxState;
            }

            int firstIndex = (int)((firstSample - info.firstState) / info.stepSize);
            int lastIndex = (int)((lastSample - info.firstState) / info.stepSize);
            int range = lastIndex - firstIndex + 1;

            printTitle();

            List<String> labels = new ArrayList<>();
            labels.add("sample");
            labels.add("state");

            for (LogColumn column : log) {
                labels.add(column.getLabel());
            }
            formatter.logLabels(labels.toArray(new String[0]));

            for (int i = 0; i < numberSamples; ++i) {

                int sample = firstIndex + MathUtils.nextInt(range);

                for (ParameterBinding binding : parameterBindings) {
                    setParameter(binding, sample);
                }

                for (TreeBinding binding : treeBindings) {
                    setTree(binding, sample);
                }

                List<String> values = new ArrayList<>();
                values.add(Integer.toString(i));
                values.add(Long.toString(info.getState(sample)));

                for (LogColumn column : log) {
                    values.add(column.getFormatted());
                }
                formatter.logValues(values.toArray(new String[0]));

                if (printToScreen) {
                    Logger.getLogger("dr.evomodelxml").info("Iteration " + i + " completed.");
                }
            }

            formatter.stopLogging();
        }

        public void addTreeBinding(TreeBinding treeBinding) {
            treeBindings.add(treeBinding);
        }

        public void addParameterBinding(ParameterBinding parameterBinding) {
            parameterBindings.add(parameterBinding);
        }
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
            } catch (TraceException e) {
                throw new XMLParseException(e.getMessage());
            } catch (IOException e) {
                throw new XMLParseException(e.getMessage());
            }

            logCache.put(fileName, traces);
            return traces;
        }
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxoExecute = xo.getChild(EXECUTE);
        List<Loggable> loggable = new ArrayList<>();

        for (int i = 0; i < cxoExecute.getChildCount(); ++i) {
            Object obj = cxoExecute.getChild(i);
            if (obj instanceof Loggable) {
                loggable.add((Loggable) obj);
            }
        }

        final PrintWriter pw = getLogFile(xo, getParserName());

        boolean printStatus = xo.hasAttribute(FILE_NAME);

        Action action = new Action(loggable, pw, printStatus);
        
        for (XMLObject cxo : xo.getAllChildren(SAMPLE_BLOCK)) {

            EmpiricalTreeDistributionModel treeModel = (EmpiricalTreeDistributionModel) cxo.getChild(TreeModel.class);
            Parameter parameter = (Parameter) cxo.getChild(Parameter.class);

            if (treeModel != null) {

                action.addTreeBinding(new TreeBinding(treeModel));

            } else if (parameter != null) {
                
                String fileName = cxo.getStringAttribute(FILE_NAME);
                String columnName = cxo.getStringAttribute(COLUMN_NAME);
                Logger.getLogger("dr.evomodelxml").info("Reading " + columnName + " from " + fileName);

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

                action.addParameterBinding(new ParameterBinding(parameter, traces, traceIndexParameter));
            }
        }

        long firstSample = xo.getAttribute(FIRST_SAMPLE, -1L);
        long lastSample = xo.getAttribute(LAST_SAMPLE, -1L);
        int numberSamples = xo.getAttribute(NUMBER_SAMPLES, 100);

        action.run(firstSample, lastSample, numberSamples);

        return null;
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
    };
}