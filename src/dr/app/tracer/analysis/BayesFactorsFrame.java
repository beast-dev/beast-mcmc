package dr.app.tracer.analysis;

import dr.inference.trace.MarginalLikelihoodAnalysis;
import org.virion.jam.framework.AuxilaryFrame;
import org.virion.jam.framework.DocumentFrame;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.text.DecimalFormat;
import java.util.List;
import java.util.ArrayList;

public class BayesFactorsFrame extends AuxilaryFrame {

	private List<MarginalLikelihoodAnalysis> marginalLikelihoods = new ArrayList<MarginalLikelihoodAnalysis>();
	private JPanel contentPanel;

	private BayesFactorsModel bayesFactorsModel;

	public BayesFactorsFrame(DocumentFrame frame, String title) {

		super(frame);

		setTitle(title);

		bayesFactorsModel = new BayesFactorsModel();
		JTable bayesFactorsTable = new JTable(bayesFactorsModel);

		JScrollPane scrollPane1 = new JScrollPane(bayesFactorsTable,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

		contentPanel = new JPanel(new BorderLayout(0, 0));
		contentPanel.setOpaque(false);
		contentPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(
				new java.awt.Insets(0, 0, 6, 0)));
		contentPanel.add(scrollPane1, BorderLayout.CENTER);


		setContentsPanel(contentPanel);

		getSaveAction().setEnabled(false);
		getSaveAsAction().setEnabled(false);

		getCutAction().setEnabled(false);
		getCopyAction().setEnabled(true);
		getPasteAction().setEnabled(false);
		getDeleteAction().setEnabled(false);
		getSelectAllAction().setEnabled(false);
		getFindAction().setEnabled(false);

		getZoomWindowAction().setEnabled(false);
	}

	public void addMarginalLikelihood(MarginalLikelihoodAnalysis marginalLikelihood) {
		this.marginalLikelihoods.add(marginalLikelihood);
		bayesFactorsModel.fireTableStructureChanged();
	}

	public void initializeComponents() {

		setSize(new Dimension(640, 480));
	}

	public boolean useExportAction() { return true; }

	public JComponent getExportableComponent() {
		return contentPanel;
	}

	public void doCopy() {
		java.awt.datatransfer.Clipboard clipboard =
				Toolkit.getDefaultToolkit().getSystemClipboard();

		java.awt.datatransfer.StringSelection selection =
				new java.awt.datatransfer.StringSelection(bayesFactorsModel.toString());

		clipboard.setContents(selection, selection);
	}



	class BayesFactorsModel extends AbstractTableModel {

		String[] columnNames = {"Trace", "log Marginal Likelihood"};

		private DecimalFormat formatter = new DecimalFormat("0.###E0");
		private DecimalFormat formatter2 = new DecimalFormat("####0.###");

		public BayesFactorsModel() {
		}

		public int getColumnCount() {
			return columnNames.length  + marginalLikelihoods.size();
		}

		public int getRowCount() {
			return marginalLikelihoods.size();
		}

		public Object getValueAt(int row, int col) {

			if (col == 0) {
				return marginalLikelihoods.get(row).getTraceName();
			}

			if (col == 1) {
				return formatter2.format(marginalLikelihoods.get(row).getLogMarginalLikelihood());
			} else {
				if (col - 2 > row) {
					double lnML1 = marginalLikelihoods.get(row).getLogMarginalLikelihood();
					double lnML2 = marginalLikelihoods.get(col - 2).getLogMarginalLikelihood();
					return formatter.format(lnML1 - lnML2);
				} else {
					return "-";
				}
			}
		}

		public String getColumnName(int column) {
			if (column < columnNames.length) {
				return columnNames[column];
			}
			return marginalLikelihoods.get(column - 2).getTraceName();
		}

		public Class getColumnClass(int c) {
			return getValueAt(0, c).getClass();
		}

		public String toString() {
			StringBuffer buffer = new StringBuffer();

			buffer.append(getColumnName(0));
			for (int j = 1; j < getColumnCount(); j++) {
				buffer.append("\t");
				buffer.append(getColumnName(j));
			}
			buffer.append("\n");

			for (int i = 0; i < getRowCount(); i++) {
				buffer.append(getValueAt(i, 0));
				for (int j = 1; j < getColumnCount(); j++) {
					buffer.append("\t");
					buffer.append(getValueAt(i, j));
				}
				buffer.append("\n");
			}

			return buffer.toString();
		}
	}
}