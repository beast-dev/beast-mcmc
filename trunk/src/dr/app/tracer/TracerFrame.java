/*
 * TracerFrame.java
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

import org.virion.jam.framework.DocumentFrame;
import org.virion.jam.panels.ActionPanel;
import org.virion.jam.util.LongTask;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class TracerFrame extends DocumentFrame implements AnalysisMenuHandler {

	private TracePanel tracePanel = null;
	
	private JTable traceTable = null;
	private TraceTableModel traceTableModel = null;
	private JSplitPane splitPane1 = null;
	private JCheckBox combinedCheckBox = null;
	private JPanel topPanel = null;
	
	private JTable statisticTable = null;
	private StatisticTableModel statisticTableModel = null;
	
	private JScrollPane scrollPane1 = null;
	private JScrollPane scrollPane2 = null;
	private JPanel leftPanel = null;
	
	private JLabel progressLabel;
	private JProgressBar progressBar;
	
	private ArrayList traceLists = new ArrayList();
	private CombinedTraces combinedTraces = null;
	private TraceList currentTraceList = null;
	
	private DemographicDialog demographicDialog = null;
    private BayesianSkylineDialog bayesianSkylineDialog = null;
	
	private int dividerLocation = -1;

	private int storedSelectedTrace = 0;
	
	public TracerFrame(String title) {
    	super();
    	
		setTitle(title);
		
        getOpenAction().setEnabled(false);
		getSaveAction().setEnabled(false);
		getSaveAsAction().setEnabled(false);

		getCutAction().setEnabled(false);
		getCopyAction().setEnabled(true);
		getPasteAction().setEnabled(false);
		getDeleteAction().setEnabled(false);
		getSelectAllAction().setEnabled(false);
		getFindAction().setEnabled(false);

		getZoomWindowAction().setEnabled(false);

		getDemographicAction().setEnabled(false);
        getBayesianSkylineAction().setEnabled(false);

        setImportAction(importAction);
        setExportAction(exportAction);
	}

	public void initializeComponents() {

		setSize(new java.awt.Dimension(800, 600));

		tracePanel = new TracePanel(this);
		tracePanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 6, 12, 12)));

        traceTableModel = new TraceTableModel();
        traceTable = new JTable(traceTableModel);
		TableRenderer renderer = new TableRenderer(SwingConstants.LEFT, 4); 
 		traceTable.getColumnModel().getColumn(0).setCellRenderer(renderer); 
		traceTable.getColumnModel().getColumn(1).setPreferredWidth(50); 
		traceTable.getColumnModel().getColumn(1).setCellRenderer(renderer); 
		traceTable.getColumnModel().getColumn(2).setPreferredWidth(50);
		traceTable.getColumnModel().getColumn(2).setCellRenderer(renderer); 
		traceTable.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        
		traceTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent evt) { traceTableSelectionChanged(); }
		});

		scrollPane1 = new JScrollPane(traceTable, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, 
												JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

        ActionPanel actionPanel = new ActionPanel(false);
        actionPanel.setAddAction(getImportAction());
        actionPanel.setRemoveAction(getRemoveTraceAction());
        getRemoveTraceAction().setEnabled(false);

		combinedCheckBox = new JCheckBox("Combined");
		combinedCheckBox.setEnabled(false);
		combinedCheckBox.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent evt) { setCombinedChanged(); }
		});
        JPanel controlPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		controlPanel.add(actionPanel);
        controlPanel.add(combinedCheckBox);
		
		topPanel = new JPanel(new BorderLayout(0,0));
		topPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(0, 0, 6, 0)));
		topPanel.add(new JLabel("Trace Files:"), BorderLayout.NORTH);
		topPanel.add(scrollPane1, BorderLayout.CENTER);
		topPanel.add(controlPanel, BorderLayout.SOUTH);
		
		statisticTableModel = new StatisticTableModel();
		statisticTable = new JTable(statisticTableModel);
 		statisticTable.getColumnModel().getColumn(0).setCellRenderer(renderer); 
		statisticTable.getColumnModel().getColumn(1).setPreferredWidth(50); 
		statisticTable.getColumnModel().getColumn(1).setCellRenderer(renderer); 
		statisticTable.getColumnModel().getColumn(2).setPreferredWidth(50);
		statisticTable.getColumnModel().getColumn(2).setCellRenderer(renderer); 
		statisticTable.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		
		statisticTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent evt) { statisticTableSelectionChanged(); }
		});
		
		scrollPane2 = new JScrollPane(statisticTable, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, 
												JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

		JPanel bottomPanel = new JPanel(new BorderLayout(0,0));
		bottomPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(6, 0, 0, 0)));
		bottomPanel.add(new JLabel("Traces:"), BorderLayout.NORTH);
		bottomPanel.add(scrollPane2, BorderLayout.CENTER);

		leftPanel = new JPanel(new BorderLayout(0, 0));
        leftPanel.setPreferredSize(new Dimension(400,300));
		splitPane1 = new JSplitPane(JSplitPane.VERTICAL_SPLIT, true, topPanel, bottomPanel);
		splitPane1.setBorder(null);
		
		JPanel progressPanel = new JPanel(new BorderLayout(0,0));
		progressLabel = new JLabel("");
		progressBar = new JProgressBar();
		progressPanel.add(progressLabel, BorderLayout.NORTH);
		progressPanel.add(progressBar, BorderLayout.CENTER);
		progressPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(6, 0, 0, 0)));
		
		leftPanel.add(splitPane1, BorderLayout.CENTER);
		leftPanel.add(progressPanel, BorderLayout.SOUTH);
		leftPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 6)));

		JSplitPane splitPane2 = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, true, leftPanel, tracePanel);
		splitPane2.setBorder(null);
		splitPane2.setDividerLocation(300);

		getContentPane().setLayout(new java.awt.BorderLayout(0, 0));
		getContentPane().add(splitPane2, BorderLayout.CENTER);
		
		splitPane1.setDividerLocation(0.5);

	}
	
	private void setupDividerLocation() {
		
		if (dividerLocation == -1 || dividerLocation == splitPane1.getDividerLocation()) {
		   	int h0 = topPanel.getHeight();
		   	int h1 = scrollPane1.getViewport().getHeight();
		   	int h2 = traceTable.getPreferredSize().height;
   			dividerLocation = h2 + h0 - h1;

//		   	int h0 = topPanel.getHeight() - scrollPane1.getViewport().getHeight();
// 			dividerLocation = traceTable.getPreferredSize().height + h0;
   		
   			if (dividerLocation > 400) dividerLocation = 400;
			splitPane1.setDividerLocation(dividerLocation);
		}
	}
	
	public void addTraceList(TraceList traceList) {
		traceLists.add(traceList);
		analyseTraceList(traceList);
		
		if (traceLists.size() > 1) {
			combinedCheckBox.setEnabled(true);
            updateCombinedTraces();
			setCombinedChanged();
			}
			
		getCopyAction().setEnabled(true);
		getDemographicAction().setEnabled(true);
        getBayesianSkylineAction().setEnabled(true);

		traceTableModel.fireTableDataChanged();
			if (combinedCheckBox.isSelected()) {
            traceTable.getSelectionModel().setSelectionInterval(0, traceTable.getRowCount());
        } else {
            int row = traceLists.size() - 1;
            traceTable.getSelectionModel().setSelectionInterval(row, row);
        }
		setupDividerLocation();
	}

    private void removeTraceList() {
        int row = traceTable.getSelectedRow();
				if (row < 0) row = 0;
        TraceList traceList = (TraceList)traceLists.get(row);
        traceLists.remove(traceList);

        if (traceLists.size() < 2) {
            combinedCheckBox.setEnabled(false);
            combinedTraces = null;
        } else {
            updateCombinedTraces();
            setCombinedChanged();
		}
		
        if (traceLists.size() < 1) {
            getRemoveTraceAction().setEnabled(false);
            getCopyAction().setEnabled(false);
            getDemographicAction().setEnabled(false);
            getBayesianSkylineAction().setEnabled(false);
            currentTraceList = null;
            statisticTableModel.fireTableDataChanged();
        }

		traceTableModel.fireTableDataChanged();
        if (traceLists.size() > 0) {
            if (row == traceLists.size()) {
                row--;
            }
		traceTable.getSelectionModel().setSelectionInterval(row, row);
		setupDividerLocation();
        }
	}
	
	public void setBurnIn(int index, int burnIn) {
		Traces traces = (Traces)traceLists.get(index);
		traces.setBurnIn(burnIn);
		statisticTableModel.fireTableDataChanged();
		analyseTraceList((TraceList)traceLists.get(index));
        updateCombinedTraces();
	}
	
    public void updateCombinedTraces() {
        TraceList[] traces = new TraceList[traceLists.size()];
        traceLists.toArray(traces);
        try {
            combinedTraces = new CombinedTraces("Combined", traces);

            analyseTraceList(combinedTraces);
        } catch (TraceException te) {
            // do nothing
        }
    }

	public void setCombinedChanged() {
		int row = statisticTable.getSelectedRow();
		if (row < 0) row = 0;
		

		if (combinedCheckBox.isSelected()) {
			storedSelectedTrace = traceTable.getSelectedRow();

			traceTable.setEnabled(false);

			currentTraceList = combinedTraces;

			traceTable.getSelectionModel().setSelectionInterval(0, traceTable.getRowCount());
		} else {
			traceTable.setEnabled(true);

			currentTraceList = (TraceList)traceLists.get(storedSelectedTrace);
			
			traceTable.getSelectionModel().setSelectionInterval(storedSelectedTrace, storedSelectedTrace);
		}
		
		statisticTableModel.fireTableDataChanged();
		statisticTable.getSelectionModel().setSelectionInterval(row, row);
	}
	
	public void traceTableSelectionChanged() {
		if (!combinedCheckBox.isSelected()) {
			int row = traceTable.getSelectedRow();
			if (row < 0) {
                getRemoveTraceAction().setEnabled(false);
                return;
            }

            getRemoveTraceAction().setEnabled(true);
            getBayesianSkylineAction().setEnabled(true);
			
			currentTraceList = (TraceList)traceLists.get(row);
		 
			int[] rows = statisticTable.getSelectedRows();		
			statisticTableModel.fireTableDataChanged();
			
			if (rows != null && rows.length > 0) {
				for (int i = 0; i < rows.length; i++) {
					statisticTable.getSelectionModel().addSelectionInterval(rows[i], rows[i]);
				}
			} else {
				statisticTable.getSelectionModel().setSelectionInterval(0, 0);
			}

		} else {
            getRemoveTraceAction().setEnabled(false);
            getBayesianSkylineAction().setEnabled(false);
		}
	}

	public void statisticTableSelectionChanged() {
		
		if (combinedCheckBox.isSelected()) {
			int[] selRows = statisticTable.getSelectedRows();
			tracePanel.setCombinedTraces(combinedTraces, selRows);
		} else {
			int[] selRows = statisticTable.getSelectedRows();
			tracePanel.setTraces(currentTraceList, selRows);
		}
	}

	public void analyseTraceList(TraceList traceList) {
		
		if (analyseTask == null) {
			analyseTask = new AnalyseTraceTask();
	
			javax.swing.Timer timer = new javax.swing.Timer(1000, new ActionListener() {
					public void actionPerformed(ActionEvent evt) {
						progressBar.setMaximum(analyseTask.getLengthOfTask());
						progressBar.setValue(analyseTask.getCurrent());
					}
				});
	        
			analyseTask.go();
			timer.start();
		}

		analyseTask.addTraceList(traceList);
	}
	
	AnalyseTraceTask analyseTask = null;
	
	class AnalyseTraceTask extends LongTask {
	
		class AnalysisStack {
		    private ArrayList traceLists = new ArrayList();

		    public synchronized void add(TraceList traceList) {
				traceLists.add(traceList);
		    }

		    public synchronized int getCount() {
				 return traceLists.size();
		    }

		    public synchronized TraceList get(int index) {
		        return (TraceList)traceLists.get(index);
		    }
		    
		    public synchronized void remove(int index) {
		        traceLists.remove(index);
		    }
		};

		private AnalysisStack analysisStack = new AnalysisStack();
		
		public AnalyseTraceTask() {
		}
		
		public void addTraceList(TraceList traceList) {
			analysisStack.add(traceList);
			current = 0;
		}
		
		public int getCurrent() { return current; }
		public int getLengthOfTask() {
			int count = 0;
			for (int i = 0; i < analysisStack.getCount(); i++) {
				count += analysisStack.get(i).getTraceCount();
			}
			return count;
		}

		public void stop() { }
		public boolean done() { return false; }

		public String getDescription() { return "Calculating summary statistics..."; }
		public String getMessage() { return message; }
		
		public Object doWork() {
				
			current = 0;
			boolean textCleared = true;

			do {				
				if (analysisStack.getCount() > 0) {
					TraceList tl = analysisStack.get(0);
					
					for (int i = 0; i < tl.getTraceCount(); i++) {
						progressLabel.setText("Analysing " + tl.getName() + ":");
						textCleared = false;
						tl.analyseTrace(i);
						repaint();
						current += 1;
					}
					
					analysisStack.remove(0);
				} else {
					if (!textCleared) {
						progressLabel.setText("");
						textCleared = true;
					}
					try {
						Thread.sleep(500);
					} catch (InterruptedException ie) { 
						// do nothing 
					}
				}
			} while (true);
		}
		
		private int lengthOfTask = 0;
		private int current = 0;
		private String message;
	};

    public final void doExport() {

        FileDialog dialog = new FileDialog(this,
                                       "Export Data...",
                                       FileDialog.SAVE);

        dialog.setVisible(true);
        if (dialog.getFile() != null) {
            File file = new File(dialog.getDirectory(), dialog.getFile());

            try {
                exportToFile(file);


            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Unable to write file: " + ioe,
                                                            "Unable to write file",
                                                            JOptionPane.ERROR_MESSAGE);
            }
        }

    }

    protected void exportToFile(File file) throws IOException {
        FileWriter writer = new FileWriter(file);
        writer.write(tracePanel.getExportText());
        writer.close();
    }
		
	public final void doImport() {
		
		FileDialog dialog = new FileDialog(this,
	                                   "Import Trace File...",
	                                   FileDialog.LOAD);

		dialog.setVisible(true);
		if (dialog.getFile() != null) {
		    File file = new File(dialog.getDirectory(), dialog.getFile());

			try {
				importFromFile(file);


			} catch (FileNotFoundException fnfe) {
				JOptionPane.showMessageDialog(this, "Unable to open file: File not found", 
															"Unable to open file",
															JOptionPane.ERROR_MESSAGE);
			} catch (IOException ioe) {
				JOptionPane.showMessageDialog(this, "Unable to read file: " + ioe, 
															"Unable to read file",
															JOptionPane.ERROR_MESSAGE);
			}
		}

	}
    
	protected boolean importFromFile(File file) throws FileNotFoundException, IOException {
	
		try {
			final String fileName = file.getName();
			final ProgressMonitorInputStream in = new ProgressMonitorInputStream(
                                  this,
                                  "Reading " + fileName,
                                  new FileInputStream(file)); 	
			
			final Reader reader = new InputStreamReader(in);
			final JFrame frame = this;
			
			// the monitored activity must be in a new thread. 
			Thread readThread = new Thread () {
				public void run() {
					try {
						final TraceList traceList = Traces.loadTraces(reader, fileName);

						EventQueue.invokeLater ( 
							new Runnable () {
								public void run () {
									addTraceList(traceList);
								}
							}); 
			
					} catch (final TraceException tex) {
						EventQueue.invokeLater ( 
							new Runnable () {
								public void run () {
									JOptionPane.showMessageDialog(frame, "Error reading trace file: " + tex, 
																	"Error reading trace file", 
																	JOptionPane.ERROR_MESSAGE);
								}
							}); 
					} catch (final InterruptedIOException iioex) {
						// The cancel dialog button was pressed - do nothing
					} catch (final IOException ioex) {
						EventQueue.invokeLater ( 
							new Runnable () {
								public void run () {
									JOptionPane.showMessageDialog(frame, "File I/O Error: " + ioex, 
																	"File I/O Error", 
																	JOptionPane.ERROR_MESSAGE);
								}
							}); 
					} catch (final Exception ex) {
						EventQueue.invokeLater ( 
							new Runnable () {
								public void run () {
									JOptionPane.showMessageDialog(frame, "Fatal exception: " + ex, 
																	"Error reading file",
																	JOptionPane.ERROR_MESSAGE);
								}
							}); 
					}

				}
			};
			readThread.start(); 			
			
		} catch (IOException ioex) {
			JOptionPane.showMessageDialog(this, "File I/O Error: " + ioex, 
												"File I/O Error", 
												JOptionPane.ERROR_MESSAGE);
			return false;
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this, "Fatal exception: " + ex, 
												"Error reading file",
												JOptionPane.ERROR_MESSAGE);
			return false;
		}
		return true;
	}
	
	protected boolean readFromFile(File file) throws FileNotFoundException, IOException {
		throw new RuntimeException("Cannot read file - use import instead");
	}

	protected boolean writeToFile(File file) {
		throw new RuntimeException("Cannot write file - this is a read-only application");
	}

	public void doCopy() {
		tracePanel.doCopy();
	}
    
 	public void doDemographic() {
		if (demographicDialog == null) {
			demographicDialog = new DemographicDialog(this);
		} 
		
		TraceList traceList = currentTraceList;
		if (combinedCheckBox.isSelected()) {
			traceList = combinedTraces;
		}

		if (demographicDialog.showDialog(traceList) == JOptionPane.CANCEL_OPTION) {
			return;
		}
		
		demographicDialog.createDemographicFrame(traceList, this);
	}

    public void doBayesianSkyline() {
       if (bayesianSkylineDialog == null) {
           bayesianSkylineDialog = new BayesianSkylineDialog(this);
       }

       TraceList traceList = currentTraceList;
       if (combinedCheckBox.isSelected()) {
           traceList = combinedTraces;
       }

       if (bayesianSkylineDialog.showDialog(traceList) == JOptionPane.CANCEL_OPTION) {
           return;
       }

       bayesianSkylineDialog.createBayesianSkylineFrame(traceList, this);
	}
    
   public JComponent getExportableComponent() {
    		
		return tracePanel.getExportableComponent();
	} 	

	class TraceTableModel extends AbstractTableModel {
        final String[] columnNames = {"Log File", "States", "Burn-In"};

		public int getColumnCount() { return columnNames.length; }
		public int getRowCount() {
			return traceLists.size(); 
		}
		public String getColumnName(int col) { return columnNames[col]; }

		public Object getValueAt(int row, int col) {
			Traces traces = (Traces)traceLists.get(row);
			switch (col) {
				case 0: return traces.getName();
				case 1: return new Integer(traces.getMaxState());
				case 2: return new Integer(traces.getBurnIn());
			}
			return null;
		}
		
		public void setValueAt(Object value, int row, int col) {
			if (col == 2) {
				setBurnIn(row, ((Integer)value).intValue());
			}
		}

        public Class getColumnClass(int c) { 
        	return getValueAt(0, c).getClass();
        }
        
        public boolean isCellEditable(int row, int col) {
            if (col < 2) { 
                return false;
            } else {
                return true;
            }
        }
    };

    class StatisticTableModel extends AbstractTableModel {
        final String[] columnNames = {"Statistic", "Mean", "ESS"};

		private DecimalFormat formatter = new DecimalFormat("0.###E0");
		private DecimalFormat formatter2 = new DecimalFormat("####0.###");

		public int getColumnCount() { return columnNames.length; }
		public int getRowCount() { 
			if (currentTraceList == null) return 0;
			return currentTraceList.getTraceCount();
		}
		public String getColumnName(int col) { return columnNames[col]; }

		public Object getValueAt(int row, int col) {
			if (col == 0) return currentTraceList.getTraceName(row);
			
			TraceDistribution td = currentTraceList.getDistributionStatistics(row);
			if (td == null) return "-";
			
			double value = 0.0;
			boolean warning = false;
			switch (col) {
				case 1: value = td.getMean(); break;
				case 2:
					if (!td.isValid()) return "-";
					value = td.getESS(); 
					if (value < 100.0) warning = true;
				break;
			}
			
			String string;
			if (value > -100 && (Math.abs(value) < 0.1 || Math.abs(value) >= 10000.0)) {
				string = formatter.format(value);
			} else {
                string = formatter2.format(value);
            }
			
			if (warning) {
				return "<html><font color=\"#EE0000\">" + string + "</font></html> ";
			}
			
			return string;
		}

        public Class getColumnClass(int c) { 
        	return String.class;
        }
    };

    private AbstractAction importAction = new AbstractAction("Import...") {
        public void actionPerformed(ActionEvent ae) {
            doImport();
        }
    };

    private AbstractAction exportAction = new AbstractAction("Export...") {
        public void actionPerformed(ActionEvent ae) {
            doExport();
        }
    };

    public Action getRemoveTraceAction() { return removeTraceAction; }

     private AbstractAction removeTraceAction = new AbstractAction() {
          public void actionPerformed(ActionEvent ae) {
              removeTraceList();
        }
    };

	public Action getDemographicAction() { return demographicAction; }
	
 	private AbstractAction demographicAction = new AbstractAction(AnalysisMenuFactory.DEMOGRAPHIC_RECONSTRUCTION) {
  		public void actionPerformed(ActionEvent ae) {
  			doDemographic();
  		}
  	};
   	
    public Action getBayesianSkylineAction() { return bayesianSkylineAction; }

     private AbstractAction bayesianSkylineAction = new AbstractAction(AnalysisMenuFactory.BAYESIAN_SKYLINE_RECONSTRUCTION) {
          public void actionPerformed(ActionEvent ae) {
              doBayesianSkyline();
          }
      };
}