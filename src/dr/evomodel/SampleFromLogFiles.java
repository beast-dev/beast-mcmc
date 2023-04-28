/*
 * SampleFromLogFiles.java
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

package dr.evomodel;

import dr.app.beast.BeastVersion;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.EmpiricalTreeDistributionModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.Logger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.trace.LogFileTraces;
import dr.math.MathUtils;
import dr.util.Timer;
import dr.util.Transform;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

/**
 * @author Marc Suchard
 */

public class SampleFromLogFiles {

    private final List<TreeBinding> treeBindings = new ArrayList<>();
    private final List<ParameterBinding> parameterBindings = new ArrayList<>();
    private final List<CheckBinding> checkBindings = new ArrayList<>();

    private final List<Logger> loggers;
    private final LogColumn[] log;
    private final TabDelimitedFormatter formatter;
    private final boolean printToScreen;

    private final Timer timer = new Timer();

    public SampleFromLogFiles(List<Loggable> loggable, List<Logger> loggers, PrintWriter printWriter,
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

        this.loggers = loggers;
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
        Transform transform = binding.transform;

        for (int i = 0; i < parameter.getDimension(); ++i) {
            double value = traces.getTrace(index + i).getValue(sample);
            if (transform != null) {
                value = transform.transform(value);
            }
            parameter.setParameterValueQuietly(i, value);
        }
        parameter.fireParameterChangedEvent();
    }

    private void setTree(TreeBinding binding, int sample) {
        EmpiricalTreeDistributionModel treeModel = binding.treeModel;
        treeModel.setTree(sample);
    }

    private void check(CheckBinding binding, int sample) {
        Likelihood likelihood = binding.likelihood;
        LogFileTraces traces = binding.traces;
        int index = binding.index;
        double tolerance = binding.tolerance;

        double currentValue = likelihood.getLogLikelihood();
        double logValue = traces.getTrace(index).getValue(sample);

        if (Math.abs(currentValue - logValue) > tolerance) {
            throw new RuntimeException("Check failed: restored value (" + currentValue +
                    ") != log value (" + logValue + ")");
        }
    }

    private static class StateInfo {
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

        for (Logger logger : loggers) {
            logger.startLogging();
        }

        timer.start();

        for (int i = 0; i < numberSamples; ++i) {

            int sample = firstIndex + MathUtils.nextInt(range);

            for (ParameterBinding binding : parameterBindings) {
                setParameter(binding, sample);
            }

            for (TreeBinding binding : treeBindings) {
                setTree(binding, sample);
            }

            for (CheckBinding binding : checkBindings) {
                check(binding, sample);
            }

            List<String> values = new ArrayList<>();
            values.add(Integer.toString(i));
            values.add(Long.toString(info.getState(sample)));

            for (LogColumn column : log) {
                values.add(column.getFormatted());
            }
            formatter.logValues(values.toArray(new String[0]));

            for (Logger logger : loggers) {
                logger.log(info.getState(sample));
            }

            if (printToScreen) {
                java.util.logging.Logger.getLogger("dr.evomodelxml").info("Iteration " + i + " completed.");
            }
        }

        timer.stop();

        formatter.stopLogging();

        for (Logger logger : loggers) {
            logger.stopLogging();
        }
    }
           //dr.util.Timer
    public Timer getTimer() {
        return timer;
    }

    public void addTreeBinding(TreeBinding treeBinding) {
        treeBindings.add(treeBinding);
    }

    public void addParameterBinding(ParameterBinding parameterBinding) {
        parameterBindings.add(parameterBinding);
    }

    public void addCheckBinding(CheckBinding checkBinding) {
        checkBindings.add(checkBinding);
    }

    public static class TreeBinding {

        Tree[] trees;
        EmpiricalTreeDistributionModel treeModel;

        public TreeBinding(EmpiricalTreeDistributionModel treeModel) {
            this.treeModel = treeModel;
            this.trees = treeModel.getTrees();
        }
    }

    public static class ParameterBinding {
        Parameter parameter;
        LogFileTraces traces;
        int index;
        Transform transform;

        public ParameterBinding(Parameter parameter, LogFileTraces traces, int index, Transform transform) {
            this.parameter = parameter;
            this.traces = traces;
            this.index = index;
            this.transform = transform;
        }
    }

    public static class CheckBinding {
        Likelihood likelihood;
        LogFileTraces traces;
        int index;
        double tolerance;

        public CheckBinding(Likelihood likelihood, LogFileTraces traces, int index, double tolerance) {
            this.likelihood = likelihood;
            this.traces = traces;
            this.index = index;
            this.tolerance = tolerance;
        }
    }
}
