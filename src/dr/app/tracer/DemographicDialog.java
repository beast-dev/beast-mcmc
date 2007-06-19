/*
 * DemographicDialog.java
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

import dr.evolution.coalescent.*;
import dr.evolution.util.Units;
import dr.util.Variate;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.util.LongTask;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/** 
 * DemographicDialog.java
 *
 * Title:			Tracer
 * Description:		An application for analysing MCMC trace files.
 * @author			Andrew Rambaut	
 * @author			Alexei Drummond	
 * @version			$Id: DemographicDialog.java,v 1.12 2005/07/11 14:07:25 rambaut Exp $
 */
public class DemographicDialog {
	
	private JFrame frame;
	
	private JComboBox demographicCombo; 
	private WholeNumberField binCountField;

    public static String[] demographicModels = { "Constant Population",
                                                 "Exponential Growth (Growth Rate)",
                                                 "Exponential Growth (Doubling Time)",
                                                 "Logistic Growth (Growth Rate)",
                                                 "Logistic Growth (Doubling Time)",
                                                 "Expansion (Growth Rate)",
                                                 "Expansion (Doubling Time)",
                                                 "Boom-Bust"};
		
	private String[][] argumentGuesses = { 
			{ "populationsize", "population", "popsize", "n0", "size", "pop" },
			{ "ancestralsize", "ancestralproportion", "proportion", "ancestral", "n1" },
			{ "growthrate", "growth", "rate", "r" },
            { "doublingtime", "doubling", "time", "t" },
			{ "logisticshape", "halflife", "t50", "time50", "shape" },
			{ "spikefactor", "spike", "factor", "f" },
			{ "cataclysmtime", "cataclysm", "time", "t" } };

	private String[] argumentNames = new String[] {
		"Population Size", "Ancestral Proportion", "Growth Rate", "Doubling Time", "Logistic Shape", "Spike Factor", "Spike Time"
	};
	
	private int[][] argumentIndices = { {0}, {0, 2}, {0, 3}, {0, 2, 4}, {0, 3, 4}, {0, 1, 2}, {0, 1, 3}, {0, 2, 5, 6} };
	
	private String[] argumentTraces = new String[argumentNames.length];
	private JComboBox[] argumentCombos = new JComboBox[argumentNames.length]; 
    private JComboBox maxHeightCombo = new JComboBox(new String[] {"Lower 95% HPD", "Median", "Mean", "Upper 95% HPD"});
	private JComboBox rootHeightCombo; 
	private String rootHeightTrace = "None selected";
	private OptionsPanel optionPanel;
	
	public DemographicDialog(JFrame frame) {
		this.frame = frame;

		demographicCombo = new JComboBox(demographicModels);
		
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
	
	public int getDemographicModel() { return demographicCombo.getSelectedIndex(); }

	private int findArgument(JComboBox comboBox, String argument) {
		for (int i = 0; i < comboBox.getItemCount(); i++) {
			String item = ((String)comboBox.getItemAt(i)).toLowerCase();
			if (item.indexOf(argument) != -1) return i;
		}
		return -1;
	}

 	public int showDialog(TraceList traceList) {
 		
		for (int i = 0; i < argumentCombos.length; i++) {
			argumentCombos[i].removeAllItems();
			for (int j = 0; j < traceList.getTraceCount(); j++) {
				String statistic = traceList.getTraceName(j);
				argumentCombos[i].addItem(statistic);
			}
		
			int index = findArgument(argumentCombos[i], argumentTraces[i]);
			
			for (int j = 0; j < argumentGuesses[i].length; j++) {
				if (index != -1) break;
				
				index = findArgument(argumentCombos[i], argumentGuesses[i][j]);
			}
			if (index == -1) index = 0;
			
			argumentCombos[i].setSelectedIndex(index);
		}

		setDemographicArguments();

		for (int j = 0; j < traceList.getTraceCount(); j++) {
			String statistic = traceList.getTraceName(j);
			rootHeightCombo.addItem(statistic);
		}
		int index = findArgument(rootHeightCombo, rootHeightTrace);
		if (index == -1) index = findArgument(rootHeightCombo, "root");
		if (index == -1) index = findArgument(rootHeightCombo, "height");
		if (index == -1) index = 0;
		rootHeightCombo.setSelectedIndex(index);

		JOptionPane optionPane = new JOptionPane(optionPanel,
                                    JOptionPane.QUESTION_MESSAGE,
                                    JOptionPane.OK_CANCEL_OPTION,
                                    null,
                                    null,
                                    null);
		optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

		final JDialog dialog = optionPane.createDialog(frame, "Demographic Analysis");	
   		dialog.pack();

		demographicCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				setDemographicArguments();
   				dialog.pack();
			}});

		dialog.show();
		
		int result = JOptionPane.CANCEL_OPTION;
		Integer value = (Integer)optionPane.getValue();
		if (value != null && value.intValue() != -1) {
			result = value.intValue();
		}
		
		if (result == JOptionPane.OK_OPTION) {
			for (int i = 0; i < argumentCombos.length; i++) {
				argumentTraces[i] = (String)argumentCombos[i].getSelectedItem();
			}
			rootHeightTrace = (String)rootHeightCombo.getSelectedItem();
		}
		
		return result;
	}
	
	private void setDemographicArguments() {
 		optionPanel.removeAll();
 		
		optionPanel.addComponents(new JLabel("Demographic Model:"), demographicCombo);
		
		JLabel label = new JLabel("<html>Warning! Do not select a model other than that which was<br>" +
 									"specified in BEAST to generate the trace being analysed.<br>" +
 									"<em>Any other model will produce meaningless results.</em></html>");
		label.setFont(label.getFont().deriveFont(((float)label.getFont().getSize() - 2)));
 		optionPanel.addSpanningComponent(label);
		optionPanel.addSeparator();
			
 		optionPanel.addLabel("Select the traces to use for the arguments:");

 		int demo = demographicCombo.getSelectedIndex();
 		
 		for (int i = 0; i < argumentIndices[demo].length; i++) {
 			int k = argumentIndices[demo][i];
			optionPanel.addComponentWithLabel(argumentNames[k] + ":", 
													argumentCombos[k]);
		}
		
		optionPanel.addSeparator();

        optionPanel.addComponentWithLabel("Maximum time obtained from root height's:", maxHeightCombo);

		optionPanel.addComponentWithLabel("Select the trace of the root height:", rootHeightCombo);

		optionPanel.addSeparator();
		optionPanel.addComponentWithLabel("Number of bins:", binCountField);
	}
		
	javax.swing.Timer timer = null;
	
    public void createDemographicFrame(TraceList traceList, TracerFrame parent) {
    
		DemographicFrame frame = new DemographicFrame(parent);
		frame.initialize();

		int binCount = binCountField.getValue().intValue();
		final AnalyseDemographicTask analyseTask = new AnalyseDemographicTask(traceList, frame, binCount);
	
		final ProgressMonitor progressMonitor = new ProgressMonitor(parent,
												"Analysing Demographic Model",
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
	
	class AnalyseDemographicTask extends LongTask {
			
		TraceList traceList;
		DemographicFrame frame;
		int binCount;
		
		private int lengthOfTask = 0;
		private int current = 0;
		private String message;
		
		public AnalyseDemographicTask(TraceList traceList, DemographicFrame frame, int binCount) {
			this.traceList = traceList;
			this.frame = frame;
			this.binCount = binCount;
			
			lengthOfTask = traceList.getStateCount() * 2;

		}
		
		public int getCurrent() { return current; }
		public int getLengthOfTask() { return lengthOfTask; }

		public String getDescription() { return "Calculating demographic reconstruction..."; }
		public String getMessage() { return message; }
		
		public Object doWork() {
				

			int n = traceList.getStateCount();

			current = 0;
			
			int[] argIndices = argumentIndices[demographicCombo.getSelectedIndex()];
			double[][] values = new double[argIndices.length][n];
			
			Variate[] bins = new Variate[binCount];
			for (int k = 0; k < binCount; k++) {
				bins[k] = new Variate.Double();
			}	
			
			int index = traceList.getTraceIndex(rootHeightTrace);
			double[] heights = new double[n];
			traceList.getValues(index, heights);
			
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
			
			for (int j = 0; j < argIndices.length; j++) {
				index = traceList.getTraceIndex(argumentTraces[argIndices[j]]);
				traceList.getValues(index, values[j]);
			}
			
			String title = "";
			
			if (demographicCombo.getSelectedIndex() == 0) { // Constant Size
				title = "Constant Population Size";
				ConstantPopulation demo = new ConstantPopulation(Units.YEARS);
				for (int i = 0; i < n; i++) {
					demo.setN0(values[0][i]);
				
					addDemographic(bins, binCount, delta, demo);
					current++;
				}
				
			} else if (demographicCombo.getSelectedIndex() == 1) { // Exponential Growth (Growth Rate)
				title = "Exponential Growth";
				ExponentialGrowth demo = new ExponentialGrowth(Units.YEARS);
				for (int i = 0; i < n; i++) {
					demo.setN0(values[0][i]);
					demo.setGrowthRate(values[1][i]);
				
					addDemographic(bins, binCount, delta, demo);
					current++;
				}

            } else if (demographicCombo.getSelectedIndex() == 2) { // Exponential Growth (Doubling Time)
                title = "Exponential Growth";
                ExponentialGrowth demo = new ExponentialGrowth(Units.YEARS);
                for (int i = 0; i < n; i++) {
                    demo.setN0(values[0][i]);
                    demo.setDoublingTime(values[1][i]);
				
					addDemographic(bins, binCount, delta, demo);
					current++;
				}
				
			} else if (demographicCombo.getSelectedIndex() == 3) { // Logistic Growth (Growth Rate)
				title = "Logistic Growth";
				LogisticGrowth demo = new LogisticGrowth(Units.YEARS);
				for (int i = 0; i < n; i++) {
					demo.setN0(values[0][i]);
					demo.setGrowthRate(values[1][i]);
					demo.setTime50(values[2][i]);
				
					addDemographic(bins, binCount, delta, demo);
					current++;
				}
				
            } else if (demographicCombo.getSelectedIndex() == 4) { // Logistic Growth (Doubling Time)
                title = "Logistic Growth";
                LogisticGrowth demo = new LogisticGrowth(Units.YEARS);
                for (int i = 0; i < n; i++) {
                    demo.setN0(values[0][i]);
                    demo.setDoublingTime(values[1][i]);
                    demo.setTime50(values[2][i]);

                    addDemographic(bins, binCount, delta, demo);
                    current++;
                }

			} else if (demographicCombo.getSelectedIndex() == 5) { // Expansion (Growth Rate)
				title = "Expansion";
				Expansion demo = new Expansion(Units.YEARS);
				for (int i = 0; i < n; i++) {
					demo.setN0(values[0][i]);
					demo.setProportion(values[1][i]);
					demo.setGrowthRate(values[2][i]);
				
					addDemographic(bins, binCount, delta, demo);
					current++;
				}
				
            } else if (demographicCombo.getSelectedIndex() == 6) { // Expansion (Doubling Time)
                title = "Expansion";
                Expansion demo = new Expansion(Units.YEARS);
				for (int i = 0; i < n; i++) {
					demo.setN0(values[0][i]);
                    demo.setProportion(values[1][i]);
                    demo.setDoublingTime(values[2][i]);
				
					addDemographic(bins, binCount, delta, demo);
					current++;
				}
				
			} else if (demographicCombo.getSelectedIndex() == 7) { // Cataclysm
				title = "Boom-Bust";
				CataclysmicDemographic demo = new CataclysmicDemographic(Units.YEARS);
				for (int i = 0; i < n; i++) {
					demo.setN0(values[0][i]);
					demo.setGrowthRate(values[1][i]);
					demo.setCataclysmTime(values[3][i]);
					demo.setSpikeFactor(values[2][i]);
				
					addDemographic(bins, binCount, delta, demo);
					current++;
				}
				
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
	
				current += n / binCount;
				
				t += delta;
			}
				
            frame.setupDemographic(title, xData,
                                    yDataMean, yDataMedian,
                                    yDataUpper, yDataLower,
                                    timeMean, timeMedian,
                                    timeUpper, timeLower);
			
			return null;
		}
		
		private void addDemographic(Variate[] bins, int binCount, double delta, DemographicFunction demo) {
			double t = 0;
			for (int k = 0; k < binCount; k++) {
				bins[k].add(demo.getDemographic(t));
				t += delta;
//					if (t > heights[i]) break;
			}	
			current++;
		}
	};
		
}