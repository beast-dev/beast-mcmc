package dr.app.tracer.analysis;

import jebl.evolution.io.*;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.coalescent.Intervals;
import jebl.evolution.coalescent.IntervalList;
import dr.util.Variate;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.framework.DocumentFrame;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.util.LongTask;
import dr.app.tracer.traces.TraceDistribution;
import dr.app.tracer.traces.TraceList;

import javax.swing.*;
import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
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
	private JComboBox maxHeightCombo = new JComboBox(new String[] {
			"Lower 95% HPD", "Median", "Mean", "Upper 95% HPD"});
	private JComboBox rootHeightCombo;
	private JCheckBox manualRangeCheckBox;
	private RealNumberField minTimeField;
	private RealNumberField maxTimeField;
	private JComboBox stepwiseLinearCombo = new JComboBox(new String[] {"Stepwise (Constant)", "Linear Change"});
	private String rootHeightTrace = "None selected";

	private RealNumberField ageOfYoungestField = new RealNumberField();

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

		manualRangeCheckBox = new JCheckBox("Use manual range for bins:");

		maxTimeField = new RealNumberField(0.0, Double.MAX_VALUE);
		maxTimeField.setColumns(12);

		minTimeField = new RealNumberField(0.0, Double.MAX_VALUE);
		minTimeField.setColumns(12);

		ageOfYoungestField.setValue(0.0);
		ageOfYoungestField.setColumns(12);

		optionPanel = new OptionsPanel(12, 12);
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

	public int showDialog(TraceList traceList, TemporalAnalysisFrame temporalAnalysisFrame) {

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

		setArguments(temporalAnalysisFrame);

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
			dialog.setVisible(true);

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

	private void setArguments(TemporalAnalysisFrame temporalAnalysisFrame) {
		optionPanel.removeAll();

		JLabel label = new JLabel(
				"<html>Warning! This analysis should only be run on traces where<br>" +
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

		optionPanel.addComponentWithLabel("Maximum time is the root height's:", maxHeightCombo);

		optionPanel.addComponentWithLabel("Select the trace of the root height:", rootHeightCombo);

		if (temporalAnalysisFrame == null) {
			optionPanel.addSeparator();
			optionPanel.addComponentWithLabel("Number of bins:", binCountField);

			optionPanel.addSpanningComponent(manualRangeCheckBox);
			final JLabel label1 = optionPanel.addComponentWithLabel("Minimum time:", minTimeField);
			final JLabel label2 = optionPanel.addComponentWithLabel("Maximum time:", maxTimeField);

			if (manualRangeCheckBox.isSelected()) {
				label1.setEnabled(true);
				minTimeField.setEnabled(true);
				label2.setEnabled(true);
				maxTimeField.setEnabled(true);
			} else {
				label1.setEnabled(false);
				minTimeField.setEnabled(false);
				label2.setEnabled(false);
				maxTimeField.setEnabled(false);
			}

			manualRangeCheckBox.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent changeEvent) {
					if (manualRangeCheckBox.isSelected()) {
						label1.setEnabled(true);
						minTimeField.setEnabled(true);
						label2.setEnabled(true);
						maxTimeField.setEnabled(true);
					} else {
						label1.setEnabled(false);
						minTimeField.setEnabled(false);
						label2.setEnabled(false);
						maxTimeField.setEnabled(false);
					}
				}
			});
		}

		optionPanel.addSeparator();
		optionPanel.addComponentWithLabel("Age of youngest tip:", ageOfYoungestField);
		JLabel label3 = new JLabel(
				"<html>You can set the age of sampling of the most recent tip in<br>" +
						"the tree. If this is set to zero then the plot is shown going<br>" +
						"backwards in time, otherwise forwards in time.</html>");
		label3.setFont(label3.getFont().deriveFont(((float)label3.getFont().getSize() - 2)));
		optionPanel.addSpanningComponent(label3);
	}

	javax.swing.Timer timer = null;

	public void createBayesianSkylineFrame(TraceList traceList, DocumentFrame parent) {

		TemporalAnalysisFrame frame;

		int binCount = binCountField.getValue();
		double minTime = 0.0;
		double maxTime = 0.0;
		boolean manualRange = manualRangeCheckBox.isSelected();
		if (manualRange) {
			minTime = minTimeField.getValue();
			maxTime = maxTimeField.getValue();

			if (minTime >= maxTime) {
				JOptionPane.showMessageDialog(parent,
						"The minimum time value should be less than the maximum.",
						"Error creating Bayesian skyline",
						JOptionPane.ERROR_MESSAGE);
				return;
			}

			frame = new TemporalAnalysisFrame(parent, "", binCount, minTime, maxTime);
		} else {
			frame = new TemporalAnalysisFrame(parent, "", binCount);
		}
		frame.initialize();
		addToTemporalAnalysis(traceList, frame);
	}

	public void addToTemporalAnalysis(TraceList traceList, TemporalAnalysisFrame frame) {

		int firstPopSize = traceList.getTraceIndex(argumentTraces[0]);
		int popSizeCount = getTraceRange(traceList, firstPopSize);
		int firstGroupSize = traceList.getTraceIndex(argumentTraces[1]);
		int groupSizeCount = getTraceRange(traceList, firstGroupSize);

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
				frame);

		final ProgressMonitor progressMonitor = new ProgressMonitor(frame,
				"Analysing Bayesian Skyline",
				"", 0, analyseTask.getLengthOfTask());
		progressMonitor.setMillisToPopup(0);
		progressMonitor.setMillisToDecideToPopup(0);

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
		TemporalAnalysisFrame frame;
		File treeFile;
		int firstPopSize;
		int firstGroupSize;
		int popSizeCount;
		int groupSizeCount;
		int binCount;
		boolean rangeSet;
		double minTime;
		double maxTime;
		double ageOfYoungest;

		int stateCount;

		double[][] popSizes;
		double[][] groupSizes;

		private int lengthOfTask = 0;
		private int current = 0;
		private String message;
		private boolean isLinear;

		public AnalyseBayesianSkylineTask(TraceList traceList, File treeFile, int firstPopSize, int popSizeCount,
		                                  int firstGroupSize, int groupSizeCount,
		                                  TemporalAnalysisFrame frame) {
			this.traceList = traceList;
			this.frame = frame;
			this.treeFile = treeFile;
			this.firstPopSize = firstPopSize;
			this.firstGroupSize = firstGroupSize;
			this.popSizeCount = popSizeCount;
			this.groupSizeCount = groupSizeCount;

			this.binCount = frame.getBinCount();
			this.rangeSet = frame.isRangeSet();

			isLinear = stepwiseLinearCombo.getSelectedIndex() > 0;
			ageOfYoungest = ageOfYoungestField.getValue();

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
				case 0: maxHeight = timeLower; break;
				case 1: maxHeight = timeMedian; break;
				case 2: maxHeight = timeMean; break;
				case 3: maxHeight = timeUpper; break;
			}

			if (rangeSet) {
				minTime = frame.getMinTime();
				maxTime = frame.getMaxTime();
			} else {
				if (ageOfYoungest > 0.0) {
					minTime = ageOfYoungest - maxHeight;
					maxTime = ageOfYoungest;
				} else {
					minTime = 0.0;
					maxTime = maxHeight - ageOfYoungest;
				}

				frame.setRange(minTime, maxTime);
			}

			if (ageOfYoungest > 0.0) {
				// reverse them if ageOfYoungest is set positive
				timeMean = ageOfYoungest - timeMean;
				timeMedian = ageOfYoungest - timeMedian;
				timeUpper = ageOfYoungest - timeUpper;
				timeLower = ageOfYoungest - timeLower;

				// setting a timeXXXX to -1 means that it won't be displayed...
				if (minTime >= timeLower) timeLower = -1;
				if (minTime >= timeMean) timeMean = -1;
				if (minTime >= timeMedian) timeMedian = -1;
				if (minTime >= timeUpper) timeUpper = -1;
			} else {
				// otherwise use use ageOfYoungest as an offset
				timeMean = timeMean - ageOfYoungest;
				timeMedian = timeMedian - ageOfYoungest;
				timeUpper = timeUpper - ageOfYoungest;
				timeLower = timeLower - ageOfYoungest;

				// setting a timeXXXX to -1 means that it won't be displayed...
				if (maxTime <= timeLower) timeLower = -1;
				if (maxTime <= timeMean) timeMean = -1;
				if (maxTime <= timeMedian) timeMedian = -1;
				if (maxTime <= timeUpper) timeUpper = -1;
			}

			double delta = (maxTime - minTime) / (binCount - 1);

			try {
				BufferedReader reader = new BufferedReader(new FileReader(treeFile));

				String line = reader.readLine();

				TreeImporter importer = null;
				if (line.toUpperCase().startsWith("#NEXUS")) {
					importer = new NexusImporter(reader);
				} else {
					importer = new NewickImporter(reader, false);
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

				try {
					while (importer.hasTree()) {
						RootedTree tree = (RootedTree)importer.importNextTree();
						IntervalList intervals = new Intervals(tree);
						int intervalCount = intervals.getIntervalCount();
						tips = tree.getExternalNodes().size();

						// get the coalescent intervales only
						groupTimes[state] = new double[groupSizeCount];
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

							if (intervals.getIntervalType(j) == IntervalList.IntervalType.COALESCENT) {
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

				} catch (ImportException ie) {
					JOptionPane.showMessageDialog(frame, "Error parsing file: " + ie.getMessage(),
							"Error parsing file",
							JOptionPane.ERROR_MESSAGE);
				} catch (Exception ex) {
					JOptionPane.showMessageDialog(frame, "Fatal exception (email the authors):" + ex.getMessage(),
							"Fatal exception",
							JOptionPane.ERROR_MESSAGE);
					ex.printStackTrace(System.out);
				}

				Variate[] bins = new Variate[binCount];
				double height = 0.0;
				if (ageOfYoungest > 0.0) {
					height = ageOfYoungest - maxTime;
				} else {
					height = ageOfYoungest;
				}


				for (int k = 0; k < binCount; k++) {
					bins[k] = new Variate.Double();

					if (height >= 0.0 && height <= maxHeight) {
						for (state = 0; state < stateCount; state++) {

							if (isLinear) {
								double lastGroupTime = 0.0;

								int index = 0;
								while (index < groupTimes[state].length && groupTimes[state][index] < height) {
									lastGroupTime = groupTimes[state][index];
									index += 1;
								}

								if (index < groupTimes[state].length - 1) {
									double t = (height - lastGroupTime) / (groupTimes[state][index] - lastGroupTime);
									double p1 = getPopSize(index, state);
									double p2 = getPopSize(index + 1, state);
									double popsize = p1 + ((p2 - p1) * t);
									bins[k].add(popsize);
								}
							} else {
								int index = 0;
								while (index < groupTimes[state].length && groupTimes[state][index] < height) {
									index += 1;
								}

								if (index < groupTimes[state].length) {
									double popSize = getPopSize(index, state);
									if (popSize == 0.0) {
										throw new RuntimeException("Zero pop size");
									}

									bins[k].add(popSize);
								} else {
									// Do we really want to do this?
//                                bins[k].add(getPopSize(popSizeCount - 1,state));
								}
							}
						}
					}
					height += delta;
					current += 1;
				}

				Variate xData = new Variate.Double();
				Variate yDataMean = new Variate.Double();
				Variate yDataMedian = new Variate.Double();
				Variate yDataUpper = new Variate.Double();
				Variate yDataLower = new Variate.Double();

				double t;
				if (ageOfYoungest > 0.0) {
					t = maxTime;
				} else {
					t = minTime;
				}
				for (int i = 0; i < bins.length; i++) {
					xData.add(t);
					if (bins[i].getCount() > 0) {
						yDataMean.add(bins[i].getMean());
						yDataMedian.add(bins[i].getQuantile(0.5));
						yDataLower.add(bins[i].getQuantile(0.025));
						yDataUpper.add(bins[i].getQuantile(0.975));
					} else {
						yDataMean.add(Double.NaN);
						yDataMedian.add(Double.NaN);
						yDataLower.add(Double.NaN);
						yDataUpper.add(Double.NaN);
					}
					if (ageOfYoungest > 0.0) {
						t -= delta;
					} else {
						t += delta;
					}
				}

				frame.addDemographic("Bayesian Skyline: " +  traceList.getName(), xData,
						yDataMean, yDataMedian,
						yDataUpper, yDataLower,
						timeMean, timeMedian,
						timeUpper, timeLower);

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