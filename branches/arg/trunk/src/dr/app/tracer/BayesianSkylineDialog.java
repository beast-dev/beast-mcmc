/*
 * BayesianSkylineDialog.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.app.tracer;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.io.Importer;
import dr.evolution.tree.Tree;
import dr.util.Variate;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.util.LongTask;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashSet;
import java.util.Iterator;

public class BayesianSkylineDialog {

	private JFrame frame;

	private String[][] argumentGuesses = {
			{ "populationsize", "population", "popsize" },
			{ "groupsize", "groups" } };

	private String[] argumentNames = new String[] {
		"Population Size", "Group Size"
	};

    private final JButton button = new JButton("Choose File...");
    private ActionListener buttonListener;

    private final JTextField fileNameText = new JTextField("not selected", 16);
    private File treeFile = null;
    private WholeNumberField binCountField;

	private String[] argumentTraces = new String[argumentNames.length];
	private JComboBox[] argumentCombos = new JComboBox[argumentNames.length];
    private JComboBox maxHeightCombo = new JComboBox(new String[] {"Lower 95% HPD", "Median", "Mean", "Upper 95% HPD"});
	private JComboBox rootHeightCombo;
    private JComboBox stepwiseLinearCombo = new JComboBox(new String[] {"Stepwise (Constant)", "Linear Change"});
	private String rootHeightTrace = "None selected";
	private OptionsPanel optionPanel;

	public BayesianSkylineDialog(JFrame frame) {
		this.frame = frame;

		for (int i = 0; i < argumentNames.length; i++) {
			argumentCombos[i] = new JComboBox();
			argumentTraces[i] = "None selected";
		}

		rootHeightCombo = new JComboBox();

        binCountField = new WholeNumberField(2, 2000);
        binCountField.setValue(100);
        binCountField.setColumns(4);

	  	optionPanel = new OptionsPanel();
	}

	private int findArgument(JComboBox comboBox, String argument) {
		for (int i = 0; i < comboBox.getItemCount(); i++) {
			String item = ((String)comboBox.getItemAt(i)).toLowerCase();
			if (item.indexOf(argument) != -1) return i;
		}
		return -1;
	}

    private String getNumericalSuffix(String argument) {
        int i = argument.length() - 1;

        if (i < 0) return "";

        char ch = argument.charAt(i);

        if (!Character.isDigit(ch)) return "";

        while (i > 0 && Character.isDigit(ch)) {
            i -= 1;
            ch = argument.charAt(i);
        }

        String suffix = argument.substring(i+1, argument.length());
        return suffix;
    }

    private int getTraceRange(TraceList traceList, int first) {
        int i = 1;
        int k = first;

        String name = traceList.getTraceName(first);
        String root = name.substring(0, name.length() - 1);
        while (k < traceList.getTraceCount() && traceList.getTraceName(k).equals(root + i)) {
            i++;
            k++;
        }

        return i - 1;
    }

 	public int showDialog(TraceList traceList) {

         HashSet roots = new HashSet();
         for (int j = 0; j < traceList.getTraceCount(); j++) {
             String statistic = traceList.getTraceName(j);
             String suffix = getNumericalSuffix(statistic);
             if (suffix.equals("1")) {
                 roots.add(statistic.substring(0, statistic.length() - 1));
             }
         }

         if (roots.size() == 0) {
             throw new IllegalArgumentException("no traces found with a range of numerical suffixes (1-n).");
         }

		 for (int i = 0; i < argumentCombos.length; i++) {
			argumentCombos[i].removeAllItems();

             Iterator iter = roots.iterator();
             while (iter.hasNext()) {
                argumentCombos[i].addItem((String)iter.next());
             }

             int index = findArgument(argumentCombos[i], argumentTraces[i]);

             for (int j = 0; j < argumentGuesses[i].length; j++) {
                 if (index != -1) break;

                 index = findArgument(argumentCombos[i], argumentGuesses[i][j]);
             }
             if (index == -1) index = 0;

             argumentCombos[i].setSelectedIndex(index);
		}

		setArguments();

		for (int j = 0; j < traceList.getTraceCount(); j++) {
			String statistic = traceList.getTraceName(j);
			rootHeightCombo.addItem(statistic);
		}
		int index = findArgument(rootHeightCombo, rootHeightTrace);
		if (index == -1) index = findArgument(rootHeightCombo, "root");
		if (index == -1) index = findArgument(rootHeightCombo, "height");
		if (index == -1) index = 0;
		rootHeightCombo.setSelectedIndex(index);

        final JOptionPane optionPane = new JOptionPane(optionPanel,
                                    JOptionPane.QUESTION_MESSAGE,
                                    JOptionPane.OK_CANCEL_OPTION,
                                    null,
                                    null,
                                    null);
		optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        button.removeActionListener(buttonListener);
        buttonListener = new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                FileDialog dialog = new FileDialog(frame,
                        "Open Trees Log File...",
                        FileDialog.LOAD);

                dialog.setVisible(true);
                if (dialog.getFile() == null) {
                    // the dialog was cancelled...
                    return;
                }

                treeFile = new File(dialog.getDirectory(), dialog.getFile());
                fileNameText.setText(treeFile.getName());
            }};

        button.addActionListener(buttonListener);

		final JDialog dialog = optionPane.createDialog(frame, "Bayesian Skyline Analysis");
   		dialog.pack();

        int result = JOptionPane.CANCEL_OPTION;

        boolean done = true;
        do {
            done = true;
		dialog.show();

		Integer value = (Integer)optionPane.getValue();
		if (value != null && value.intValue() != -1) {
			result = value.intValue();
		}

		if (result == JOptionPane.OK_OPTION) {
                if (treeFile == null) {
                    JOptionPane.showMessageDialog(frame, "A tree file was not selected",
                            "Error parsing file",
                            JOptionPane.ERROR_MESSAGE);
                    done = false;
                } else {
			for (int i = 0; i < argumentCombos.length; i++) {
				argumentTraces[i] = (String)argumentCombos[i].getSelectedItem() + "1";
			}
			rootHeightTrace = (String)rootHeightCombo.getSelectedItem();
		}
            }
        } while (!done);

		return result;
	}

	private void setArguments() {
 		optionPanel.removeAll();

		JLabel label = new JLabel("<html>Warning! This analysis should only be run on traces where<br>" +
 									"the Bayesian Skyline plot was specified as the demographic in BEAST.<br>" +
 									"<em>Any other model will produce meaningless results.</em></html>");
		label.setFont(label.getFont().deriveFont(((float)label.getFont().getSize() - 2)));
 		optionPanel.addSpanningComponent(label);
		optionPanel.addSeparator();

        if (treeFile != null) {
            fileNameText.setText(treeFile.getName());
        }
        fileNameText.setEditable(false);

        JPanel panel = new JPanel(new BorderLayout(0,0));
        panel.add(fileNameText, BorderLayout.CENTER);
        panel.add(button, BorderLayout.EAST);
        optionPanel.addComponentWithLabel("Trees Log File: ", panel);

        optionPanel.addSeparator();

        optionPanel.addComponentWithLabel("Bayesian skyline variant: ", stepwiseLinearCombo);

        optionPanel.addSeparator();

        optionPanel.addLabel("Select the traces to use for the arguments:");

 		for (int i = 0; i < argumentNames.length; i++) {
			optionPanel.addComponentWithLabel(argumentNames[i] + ":",
													argumentCombos[i]);
		}

		optionPanel.addSeparator();

        optionPanel.addComponentWithLabel("Maximum time obtained from root height's:", maxHeightCombo);

        optionPanel.addComponentWithLabel("Select the trace of the root height:", rootHeightCombo);

        optionPanel.addSeparator();
        optionPanel.addComponentWithLabel("Number of bins:", binCountField);
	}

	javax.swing.Timer timer = null;

    public void createBayesianSkylineFrame(TraceList traceList, TracerFrame parent) {

		DemographicFrame frame = new DemographicFrame(parent);
		frame.initialize();

		int firstPopSize = traceList.getTraceIndex(argumentTraces[0]);
        int popSizeCount = getTraceRange(traceList, firstPopSize);
        int firstGroupSize = traceList.getTraceIndex(argumentTraces[1]);
        int groupSizeCount = getTraceRange(traceList, firstGroupSize);
        int binCount = binCountField.getValue().intValue();
        boolean isLinear = stepwiseLinearCombo.getSelectedIndex() > 0;

        if (isLinear) {
            if (groupSizeCount != popSizeCount - 1) {
                JOptionPane.showMessageDialog(frame,
                        "For the linear change Bayesian skyline model there should\n" +
                                "one fewer group size than population size parameters. Either\n" +
                                "this is a stepwise (constant) model or the wrong parameters\n" +
                                "were specified. Please try again and check.",
                        "Error creating Bayesian skyline",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }
        } else {
            if (groupSizeCount != popSizeCount) {
                JOptionPane.showMessageDialog(frame,
                        "For the stepwise (constant) Bayesian skyline model there should\n" +
                        "be the same number of group size as population size parameters. \n" +
                        "Either this is a linear change model or the wrong parameters\n" +
                        "were specified. Please try again and check.",
                        "Error creating Bayesian skyline",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }
        }
		final AnalyseBayesianSkylineTask analyseTask = new AnalyseBayesianSkylineTask(traceList,
                                                                                        treeFile,
                                                                                        firstPopSize,
                popSizeCount,
                                                                                        firstGroupSize,
                groupSizeCount,
                                                                                        binCount,
                                                                                        isLinear,
                                                                                        frame);

		final ProgressMonitor progressMonitor = new ProgressMonitor(parent,
												"Analysing Bayesian Skyline",
												"", 0, analyseTask.getLengthOfTask());
		progressMonitor.setProgress(0);
		progressMonitor.setMillisToDecideToPopup(0);
		progressMonitor.setMillisToPopup(0);

		timer = new javax.swing.Timer(1000, new ActionListener() {
				public void actionPerformed(ActionEvent evt) {
					progressMonitor.setProgress(analyseTask.getCurrent());
					if (progressMonitor.isCanceled() || analyseTask.done()) {
					    progressMonitor.close();
					    analyseTask.stop();
					    timer.stop();
					}
				}
			});

		analyseTask.go();
		timer.start();
	}

	class AnalyseBayesianSkylineTask extends LongTask {

		TraceList traceList;
		DemographicFrame frame;
        File treeFile;
        int firstPopSize;
        int firstGroupSize;
        int popSizeCount;
        int groupSizeCount;
        int binCount;

        int stateCount;

        double[][] popSizes;
        double[][] groupSizes;

		private int lengthOfTask = 0;
		private int current = 0;
		private String message;
        private boolean isLinear;

        public AnalyseBayesianSkylineTask(TraceList traceList, File treeFile, int firstPopSize, int popSizeCount,
                                          int firstGroupSize, int groupSizeCount, int binCount,
                                          boolean isLinear, DemographicFrame frame) {
			this.traceList = traceList;
			this.frame = frame;
            this.treeFile = treeFile;
			this.firstPopSize = firstPopSize;
            this.firstGroupSize = firstGroupSize;
            this.popSizeCount = popSizeCount;
            this.groupSizeCount = groupSizeCount;
            this.binCount = binCount;
            this.isLinear = isLinear;

            lengthOfTask = traceList.getStateCount() + binCount;

            stateCount = traceList.getStateCount();

            popSizes = new double[popSizeCount][stateCount];
            for (int i = 0; i < popSizeCount; i++) {
                traceList.getValues(firstPopSize + i, popSizes[i]);
            }
            groupSizes = new double[groupSizeCount][stateCount];
            for (int i = 0; i < groupSizeCount; i++) {
                traceList.getValues(firstGroupSize + i, groupSizes[i]);
            }

		}

		public int getCurrent() { return current; }
		public int getLengthOfTask() { return lengthOfTask; }

		public String getDescription() { return "Calculating Bayesian skyline..."; }
		public String getMessage() { return message; }

		public Object doWork() {

            double[] heights = new double[stateCount];
            traceList.getValues(traceList.getTraceIndex(rootHeightTrace), heights);

            TraceDistribution distribution = new TraceDistribution(heights, traceList.getStepSize());

            double timeMean = distribution.getMean();
            double timeMedian = distribution.getMedian();
            double timeUpper = distribution.getUpperHPD();
            double timeLower = distribution.getLowerHPD();
            
            double maxHeight = 0.0;
            switch (maxHeightCombo.getSelectedIndex()) {
                // setting a timeXXXX to -1 means that it won't be displayed...
                case 0: maxHeight = timeLower; timeLower = -1; timeMean = -1; timeMedian = -1; timeUpper = -1; break;
                case 1: maxHeight = timeMedian; timeMean = -1; timeMedian = -1; timeUpper = -1; break;
                case 2: maxHeight = timeMean; timeMean = -1; timeMedian = -1; timeUpper = -1; break;
                case 3: maxHeight = timeUpper; timeUpper = -1; break;
            }
            double delta = maxHeight / (binCount - 1);

            try {
                BufferedReader reader = new BufferedReader(new FileReader(treeFile));

                String line = reader.readLine();

                TreeImporter importer = null;
                if (line.toUpperCase().startsWith("#NEXUS")) {
                    importer = new NexusImporter(reader);
                } else {
                    importer = new NewickImporter(reader);
                }

                int burnin = traceList.getBurnIn();
                int skip = burnin / traceList.getStepSize();
                int state = 0;

                while (importer.hasTree() && state < skip) {
                    importer.importNextTree();
                    state += 1;
                }

                current = 0;

                // the age of the end of this group
                double[][] groupTimes = new double[stateCount][];
                int tips = 0;
                state = 0;

                while (importer.hasTree()) {
                    Tree tree = importer.importNextTree();
                    TreeIntervals intervals = new TreeIntervals(tree);
                    int intervalCount = intervals.getIntervalCount();
                    tips = tree.getExternalNodeCount();

                    // get the coalescent intervales only
                    groupTimes[state] = new double[tips-1];
                    double totalTime = 0.0;
                    int groupSize = 1;
                    int groupIndex = 0;
                    int subIndex = 0;
                    if (firstGroupSize > 0) {
                        double g = groupSizes[groupIndex][state];
                        if (g != Math.round(g)) {
                            throw new RuntimeException("Group size " + groupIndex + " should be integer but found:" + g);
                        } else groupSize = (int)Math.round(g);
                    }

                    for (int j = 0; j < intervalCount; j++) {

                        totalTime += intervals.getInterval(j);

                        if (intervals.getIntervalType(j) == IntervalType.COALESCENT) {
                            subIndex += 1;
                            if (subIndex == groupSize) {
                                groupTimes[state][groupIndex] = totalTime;
                                subIndex = 0;
                                groupIndex += 1;
                                if (groupIndex < groupSizeCount) {
                                    double g = groupSizes[groupIndex][state];
                                    if (g != Math.round(g)) {
                                        throw new RuntimeException("Group size " + groupIndex + " should be integer but found:" + g);
                                    } else groupSize = (int)Math.round(g);
                                }
                            }
                        }

                        // insert zero-length coalescent intervals
                        int diff = intervals.getCoalescentEvents(j)-1;
                        if (diff > 0) throw new RuntimeException("Don't handle multifurcations!");
                    }

                    state += 1;
                    current += 1;
                }

                Variate[] bins = new Variate[binCount];
                for (int k = 0; k < binCount; k++) {
                    bins[k] = new Variate.Double();

                    double time = (double)k * maxHeight / (double)binCount;

                    for (state = 0; state < stateCount; state++) {

                        if (isLinear) {
                            double lastGroupTime = 0.0;

                            int index = 0;
                            while (index < groupTimes[state].length && groupTimes[state][index] < time) {
                                lastGroupTime = groupTimes[state][index];
                                index += 1;
                            }

                            if (index < groupTimes[state].length - 1) {
                                double t = (time - lastGroupTime) / (groupTimes[state][index] - lastGroupTime);
                                double p1 = getPopSize(index, state);
                                double p2 = getPopSize(index + 1, state);
                                double popsize = p1 + ((p2 - p1) * t);
                                bins[k].add(popsize);
                            }
                        } else {
                            int index = 0;
                            while (index < groupTimes[state].length && groupTimes[state][index] < time) {
                                index += 1;
                            }

                            if (index < groupTimes[state].length) {
                                bins[k].add(getPopSize(index, state));
                            } else {
                                // Do we really want to do this?
//                                bins[k].add(getPopSize(popSizeCount - 1,state));
                            }
                        }
                    }
                    current += 1;
                }

                Variate xData = new Variate.Double();
                Variate yDataMean = new Variate.Double();
                Variate yDataMedian = new Variate.Double();
                Variate yDataUpper = new Variate.Double();
                Variate yDataLower = new Variate.Double();

                double t = 0.0;
                for (int i = 0; i < bins.length; i++) {
                    xData.add(t);
                    yDataMean.add(bins[i].getMean());
                    yDataMedian.add(bins[i].getQuantile(0.5));
                    yDataLower.add(bins[i].getQuantile(0.025));
                    yDataUpper.add(bins[i].getQuantile(0.975));

                    t += delta;
                }

                frame.setupDemographic("Bayesian Skyline", xData,
                                        yDataMean, yDataMedian,
                                        yDataUpper, yDataLower,
                                        timeMean, timeMedian,
                                        timeUpper, timeLower);

            } catch (Importer.ImportException ie) {
                JOptionPane.showMessageDialog(frame, "Error parsing file: " + ie.getMessage(),
                        "Error parsing file",
                        JOptionPane.ERROR_MESSAGE);
            } catch (java.io.IOException ioe) {
                JOptionPane.showMessageDialog(frame, "Error reading file: " + ioe.getMessage(),
                                                            "Error reading file",
                                                            JOptionPane.ERROR_MESSAGE);
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(frame, "Fatal exception (email the authors):" + ex.getMessage(),
                                                            "Fatal exception",
                                                            JOptionPane.ERROR_MESSAGE);
                ex.printStackTrace(System.out);
            }

			return null;
		}

        private double getPopSize(int index, int state) {
            return popSizes[index][state];
        }
	};

}