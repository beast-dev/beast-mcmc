package dr.app.tracer.application;

import org.virion.jam.framework.DocumentFrame;
import org.virion.jam.panels.ActionPanel;
import org.virion.jam.table.TableRenderer;
import org.virion.jam.util.LongTask;
import dr.app.tracer.analysis.*;
import dr.app.tracer.traces.*;

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
	private JPanel topPanel = null;

	private JTable statisticTable = null;
	private StatisticTableModel statisticTableModel = null;

	private JScrollPane scrollPane1 = null;

	private JLabel progressLabel;
	private JProgressBar progressBar;

	private java.util.List<LogFileTraces> traceLists = new ArrayList<LogFileTraces>();
	private java.util.List<TraceList> currentTraceLists = new ArrayList<TraceList>();
	private CombinedTraces combinedTraces = null;

	private int dividerLocation = -1;

	private DemographicDialog demographicDialog = null;
	private BayesianSkylineDialog bayesianSkylineDialog = null;
	private NewTemporalAnalysisDialog createTemporalAnalysisDialog = null;

	public TracerFrame(String title) {
		super();

		setTitle(title);

		getOpenAction().setEnabled(false);
		getSaveAction().setEnabled(false);
		getSaveAsAction().setEnabled(false);

		getCutAction().setEnabled(false);
		getCopyAction().setEnabled(false);
		getPasteAction().setEnabled(false);
		getDeleteAction().setEnabled(false);
		getSelectAllAction().setEnabled(false);
		getFindAction().setEnabled(false);

		getZoomWindowAction().setEnabled(false);

		setImportAction(importAction);
		setExportAction(exportAction);

		setAnalysesEnabled(false);
	}

	public void initializeComponents() {

		setSize(new java.awt.Dimension(1000, 700));

		tracePanel = new TracePanel(this);
		tracePanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 6, 12, 12)));

		traceTableModel = new TraceTableModel();
		traceTable = new JTable(traceTableModel);
		TableRenderer renderer = new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4));
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

		ActionPanel actionPanel1 = new ActionPanel(false);
		actionPanel1.setAddAction(getImportAction());
		actionPanel1.setRemoveAction(getRemoveTraceAction());
		getRemoveTraceAction().setEnabled(false);

		JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
		controlPanel1.add(actionPanel1);

		topPanel = new JPanel(new BorderLayout(0,0));
		topPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(0, 0, 6, 0)));
		topPanel.add(new JLabel("Trace Files:"), BorderLayout.NORTH);
		topPanel.add(scrollPane1, BorderLayout.CENTER);
		topPanel.add(controlPanel1, BorderLayout.SOUTH);

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

		JScrollPane scrollPane2 = new JScrollPane(statisticTable, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

		JPanel bottomPanel = new JPanel(new BorderLayout(0,0));
		bottomPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(6, 0, 0, 0)));
		bottomPanel.add(new JLabel("Traces:"), BorderLayout.NORTH);
		bottomPanel.add(scrollPane2, BorderLayout.CENTER);

		JPanel leftPanel = new JPanel(new BorderLayout(0, 0));
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

		splitPane1.setDividerLocation(2000);
	}

	public void setVisible(boolean b) {
		super.setVisible(b);
		setupDividerLocation();
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

	public void setAnalysesEnabled(boolean enabled) {
		getDemographicAction().setEnabled(enabled);
		getBayesianSkylineAction().setEnabled(enabled);
		getExportAction().setEnabled(enabled);
		getCopyAction().setEnabled(true);
	}

	public void addTraceList(LogFileTraces traceList) {

		int[] selRows = traceTable.getSelectedRows();

		traceLists.add(traceList);

		updateCombinedTraces();

		setAnalysesEnabled(true);

		traceTableModel.fireTableDataChanged();

		int newRow = traceLists.size() - 1;
		traceTable.getSelectionModel().setSelectionInterval(newRow, newRow);
		if (selRows.length > 1) {
			for (int row : selRows) {
				if (row == traceLists.size() - 1) {
					row = traceLists.size();
				}
				traceTable.getSelectionModel().addSelectionInterval(row, row);
			}
		}

		setupDividerLocation();
	}

	private void removeTraceList() {
		int[] selRows = traceTable.getSelectedRows();
		LogFileTraces[] tls = new LogFileTraces[selRows.length];
		int i = 0;
		for (int row : selRows) {
			tls[i] = traceLists.get(row);
			i++;
		}
		for (LogFileTraces tl : tls) {
			traceLists.remove(tl);
		}
		updateCombinedTraces();

		if (traceLists.size() == 0) {
			getRemoveTraceAction().setEnabled(false);

			setAnalysesEnabled(false);

			currentTraceLists.clear();
			statisticTableModel.fireTableDataChanged();
		}

		traceTableModel.fireTableDataChanged();
		if (traceLists.size() > 0) {
			int row = selRows[0];
			if (row >= traceLists.size()) {
				row = traceLists.size() - 1;
			}
			traceTable.getSelectionModel().addSelectionInterval(row, row);
		}
		setupDividerLocation();
	}

	public void setBurnIn(int index, int burnIn) {
		LogFileTraces trace = (LogFileTraces)traceLists.get(index);
		trace.setBurnIn(burnIn);
		analyseTraceList(trace);
		updateCombinedTraces();
		statisticTableModel.fireTableDataChanged();
	}

	public void updateCombinedTraces() {
		if (traceLists.size() > 1) {
			TraceList[] traces = new TraceList[traceLists.size()];
			traceLists.toArray(traces);
			try {
				combinedTraces = new CombinedTraces("Combined", traces);

				analyseTraceList(combinedTraces);
			} catch (TraceException te) {
				// do nothing
			}
		} else {
			combinedTraces = null;
		}
		traceTableModel.fireTableDataChanged();
		statisticTableModel.fireTableDataChanged();
	}

	public void traceTableSelectionChanged() {
		int[] selRows = traceTable.getSelectedRows();

		if (selRows.length == 0) {
			getRemoveTraceAction().setEnabled(false);
			setAnalysesEnabled(false);
			return;
		}

		setAnalysesEnabled(true);

		getRemoveTraceAction().setEnabled(true);

		currentTraceLists.clear();

		for (int row : selRows) {
			if (row == traceLists.size()) {
				// Combined is include in the selection so disable remove
				getRemoveTraceAction().setEnabled(false);
				currentTraceLists.add(combinedTraces);
			}
		}

		for (int row : selRows) {
			if (row < traceLists.size()) {
				currentTraceLists.add(traceLists.get(row));
			}
		}

		int[] rows = statisticTable.getSelectedRows();
		statisticTableModel.fireTableDataChanged();

		if (rows.length > 0) {
			for (int i = 0; i < rows.length; i++) {
				statisticTable.getSelectionModel().addSelectionInterval(rows[i], rows[i]);
			}
		} else {
			statisticTable.getSelectionModel().setSelectionInterval(0, 0);
		}

	}

	public void statisticTableSelectionChanged() {

		int[] selRows = statisticTable.getSelectedRows();

		boolean isIncomplete = false;
		for (TraceList tl : currentTraceLists) {
			if (tl.getTraceCount() == 0 || tl.getStateCount() == 0) isIncomplete = true;
		}

		if (currentTraceLists.size() == 0 || isIncomplete) {
			tracePanel.setTraces(null, selRows);
		} else {
			TraceList[] tl = new TraceList[currentTraceLists.size()];
			currentTraceLists.toArray(tl);
			tracePanel.setTraces(tl, selRows);
		}
	}

	public void analyseTraceList(TraceList job) {

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

		analyseTask.add(job);
	}

	AnalyseTraceTask analyseTask = null;

	class AnalyseTraceTask extends LongTask {

		class AnalysisStack<T> {
			private java.util.List<T> jobs = new ArrayList<T>();

			public synchronized void add(T job) {
				jobs.add(job);
			}

			public synchronized int getCount() {
				return jobs.size();
			}

			public synchronized T get(int index) {
				return jobs.get(index);
			}

			public synchronized void remove(int index) {
				jobs.remove(index);
			}
		};

		private AnalysisStack<TraceList> analysisStack = new AnalysisStack<TraceList>();

		public AnalyseTraceTask() {
		}

		public void add(TraceList job) {
			analysisStack.add(job);
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

		public String getDescription() { return "Analysing Trace File..."; }
		public String getMessage() { return message; }

		public Object doWork() {

			current = 0;
			boolean textCleared = true;

			do {
				if (analysisStack.getCount() > 0) {
					Object job = analysisStack.get(0);
					TraceList tl = (TraceList)job;

					try {
						for (int i = 0; i < tl.getTraceCount(); i++) {
							progressLabel.setText("Analysing " + tl.getName() + ":");
							textCleared = false;
							tl.analyseTrace(i);
							repaint();
							current += 1;
						}
					} catch (final Exception ex) {
                        // do nothing. An exception is sometimes fired when burnin is changed whilst in the
                        // middle of an analysis. This doesn't seem to matter as the analysis is restarted.

                        ex.printStackTrace();
//                        EventQueue.invokeLater (
//								new Runnable () {
//									public void run () {
//										JOptionPane.showMessageDialog(TracerFrame.this, "Fatal exception: " + ex.getMessage(),
//												"Error reading file",
//												JOptionPane.ERROR_MESSAGE);
//									}
//								});
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

			LogFileTraces traces = new LogFileTraces(dialog.getFile(), file);

			processTraces(traces);

		}

	}

	protected void processTraces(final LogFileTraces traces) {

		try {
			final String fileName = traces.getName();
			final ProgressMonitorInputStream in = new ProgressMonitorInputStream(
					this,
					"Reading " + fileName,
					new FileInputStream(traces.getFile()));

			final Reader reader = new InputStreamReader(in);
			final JFrame frame = this;

			// the monitored activity must be in a new thread.
			Thread readThread = new Thread () {
				public void run() {
					try {
						traces.loadTraces(reader);

						EventQueue.invokeLater (
								new Runnable () {
									public void run () {
										analyseTraceList(traces);
										addTraceList(traces);
									}
								});

					} catch (final TraceException te) {
						EventQueue.invokeLater (
								new Runnable () {
									public void run () {
										JOptionPane.showMessageDialog(frame, "Problem with trace file: " + te.getMessage(),
												"Problem with tree file",
												JOptionPane.ERROR_MESSAGE);
									}
								});
					} catch (final InterruptedIOException iioex) {
						// The cancel dialog button was pressed - do nothing
					} catch (final IOException ioex) {
						EventQueue.invokeLater (
								new Runnable () {
									public void run () {
										JOptionPane.showMessageDialog(frame, "File I/O Error: " + ioex.getMessage(),
												"File I/O Error",
												JOptionPane.ERROR_MESSAGE);
									}
								});
//                    } catch (final Exception ex) {
//                        EventQueue.invokeLater (
//                                new Runnable () {
//                                    public void run () {
//                                        JOptionPane.showMessageDialog(frame, "Fatal exception: " + ex.getMessage(),
//                                                "Error reading file",
//                                                JOptionPane.ERROR_MESSAGE);
//                                    }
//                                });
					}

				}
			};
			readThread.start();

		} catch (FileNotFoundException fnfe) {
			JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
					"Unable to open file",
					JOptionPane.ERROR_MESSAGE);
		} catch (IOException ioex) {
			JOptionPane.showMessageDialog(this, "File I/O Error: " + ioex,
					"File I/O Error",
					JOptionPane.ERROR_MESSAGE);
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
					"Error reading file",
					JOptionPane.ERROR_MESSAGE);
		}
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

	private TemporalAnalysisFrame temporalAnalysisFrame = null;

	private void doCreateTemporalAnalysis() {
		if (createTemporalAnalysisDialog == null) {
			createTemporalAnalysisDialog = new NewTemporalAnalysisDialog(this);
		}

		if (createTemporalAnalysisDialog.showDialog() == JOptionPane.CANCEL_OPTION) {
			return;
		}

		temporalAnalysisFrame = createTemporalAnalysisDialog.createTemporalAnalysisFrame(this);
	}

	public void doDemographic(boolean add) {
		if (demographicDialog == null) {
			demographicDialog = new DemographicDialog(this);
		}

		if (currentTraceLists.size() != 1) {
			JOptionPane.showMessageDialog(this, "Please select exactly one trace to do\n" +
					"this analysis on, or select the Combined trace.",
					"Unable to perform analysis",
					JOptionPane.INFORMATION_MESSAGE);
		}


		if (add) {
			if (demographicDialog.showDialog(currentTraceLists.get(0), temporalAnalysisFrame) == JOptionPane.CANCEL_OPTION) {
				return;
			}

			demographicDialog.addToTemporalAnalysis(currentTraceLists.get(0), temporalAnalysisFrame);
		} else {
			if (demographicDialog.showDialog(currentTraceLists.get(0), null) == JOptionPane.CANCEL_OPTION) {
				return;
			}

			demographicDialog.createDemographicFrame(currentTraceLists.get(0), this);
		}
	}

	public void doBayesianSkyline(boolean add) {
		if (bayesianSkylineDialog == null) {
			bayesianSkylineDialog = new BayesianSkylineDialog(this);
		}

		if (currentTraceLists.size() != 1) {
			JOptionPane.showMessageDialog(this, "Please select exactly one trace to do\n" +
					"this analysis on, (but not the Combined trace).",
					"Unable to perform analysis",
					JOptionPane.INFORMATION_MESSAGE);
		}

		if (add) {
			if (bayesianSkylineDialog.showDialog(currentTraceLists.get(0), temporalAnalysisFrame) == JOptionPane.CANCEL_OPTION) {
				return;
			}

			bayesianSkylineDialog.addToTemporalAnalysis(currentTraceLists.get(0), temporalAnalysisFrame);
		} else {
			if (bayesianSkylineDialog.showDialog(currentTraceLists.get(0), null) == JOptionPane.CANCEL_OPTION) {
				return;
			}

			bayesianSkylineDialog.createBayesianSkylineFrame(currentTraceLists.get(0), this);
		}
	}

	private void doAddTimeDensity() {
	}

	public JComponent getExportableComponent() {

		return tracePanel.getExportableComponent();
	}

	class TraceTableModel extends AbstractTableModel {
		final String[] columnNames = {"Tree File", "States", "Burn-In"};

		public int getColumnCount() { return columnNames.length; }
		public int getRowCount() {
			int n = traceLists.size();
			if (n == 0 || combinedTraces != null) n++;
			return n;
		}
		public String getColumnName(int col) { return columnNames[col]; }

		public Object getValueAt(int row, int col) {
			TraceList traceList = null;

			if (traceLists.size() == 0) {
				switch (col) {
					case 0: return "No files loaded";
					case 1: return "";
					case 2: return "";
				}
			} else if (row == traceLists.size()) {
				traceList = combinedTraces;
				switch (col) {
					case 0: return traceList.getName();
					case 1: return new Integer(traceList.getMaxState());
					case 2: return "-";
				}
			} else {
				traceList = traceLists.get(row);
				switch (col) {
					case 0: return traceList.getName();
					case 1: return new Integer(traceList.getMaxState());
					case 2: return new Integer(traceList.getBurnIn());
				}
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
			if (col == 2 && row < traceLists.size()) {
				return true;
			}
			return false;
		}
	};

	class StatisticTableModel extends AbstractTableModel {
		final String[] columnNames = {"Statistic", "Mean", "ESS"};

		private DecimalFormat formatter = new DecimalFormat("0.###E0");
		private DecimalFormat formatter2 = new DecimalFormat("####0.###");

		public int getColumnCount() { return columnNames.length; }
		public int getRowCount() {
			if (currentTraceLists.size() == 0) return 0;
			return currentTraceLists.get(0).getTraceCount();
		}
		public String getColumnName(int col) { return columnNames[col]; }

		public Object getValueAt(int row, int col) {
			if (col == 0) return currentTraceLists.get(0).getTraceName(row);

			TraceDistribution td = currentTraceLists.get(0).getDistributionStatistics(row);
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
			if (Math.abs(value) < 0.1 || Math.abs(value) >= 10000.0) {
				string = formatter.format(value);
			} else string = formatter2.format(value);

			if (warning) {
				return "<html><font color=\"#EE0000\">" + string + "</font></html> ";
			}

			return string;
		}

		public Class getColumnClass(int c) {
			return String.class;
		}
	};


	public Action getRemoveTraceAction() { return removeTraceAction; }

	public Action getDemographicAction() { return demographicAction; }

	public Action getBayesianSkylineAction() { return bayesianSkylineAction; }

	public Action getCreateTemporalAnalysisAction() { return createTemporalAnalysisAction; }

	public Action getAddDemographicAction() { return addDemographicAction; }

	public Action getAddBayesianSkylineAction() { return addBayesianSkylineAction; }

	public Action getAddTimeDensityAction() { return addTimeDensity; }

	private AbstractAction demographicAction = new AbstractAction(AnalysisMenuFactory.DEMOGRAPHIC_RECONSTRUCTION) {
		 public void actionPerformed(ActionEvent ae) {
			 doDemographic(false);
		 }
	 };

    private AbstractAction bayesianSkylineAction = new AbstractAction(AnalysisMenuFactory.BAYESIAN_SKYLINE_RECONSTRUCTION) {
         public void actionPerformed(ActionEvent ae) {
             doBayesianSkyline(false);
         }
     };

	private AbstractAction createTemporalAnalysisAction = new AbstractAction(AnalysisMenuFactory.CREATE_TEMPORAL_ANALYSIS) {
		 public void actionPerformed(ActionEvent ae) {
			 doCreateTemporalAnalysis();
		 }
	 };

	private AbstractAction addDemographicAction = new AbstractAction(AnalysisMenuFactory.ADD_DEMOGRAPHIC_RECONSTRUCTION) {
		 public void actionPerformed(ActionEvent ae) {
			 doDemographic(true);
		 }
	 };

	private AbstractAction addBayesianSkylineAction = new AbstractAction(AnalysisMenuFactory.ADD_BAYESIAN_SKYLINE_RECONSTRUCTION) {
		 public void actionPerformed(ActionEvent ae) {
			 doBayesianSkyline(true);
		 }
	 };

	private AbstractAction addTimeDensity = new AbstractAction(AnalysisMenuFactory.ADD_TIME_DENSITY) {
		 public void actionPerformed(ActionEvent ae) {
			 doAddTimeDensity();
		 }
	 };

	private AbstractAction removeTraceAction = new AbstractAction() {
		public void actionPerformed(ActionEvent ae) {
			removeTraceList();
		}
	};

	private AbstractAction importAction = new AbstractAction("Import Trace File...") {
		public void actionPerformed(ActionEvent ae) {
			doImport();
		}
	};

	private AbstractAction exportAction = new AbstractAction("Export...") {
		public void actionPerformed(ActionEvent ae) {
			doExport();
		}
	};

}