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
    private JTable bayesFactorsTable;


    public BayesFactorsFrame(DocumentFrame frame, String title, String info, boolean hasErrors) {

		super(frame);

		setTitle(title);

        bayesFactorsModel = new BayesFactorsModel(hasErrors);
		bayesFactorsTable = new JTable(bayesFactorsModel);
        bayesFactorsTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        
        JScrollPane scrollPane1 = new JScrollPane(bayesFactorsTable,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

		contentPanel = new JPanel(new BorderLayout(0, 0));
		contentPanel.setOpaque(false);
		contentPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(
				new java.awt.Insets(0, 0, 6, 0)));
		contentPanel.add(scrollPane1, BorderLayout.CENTER);

        JLabel label = new JLabel(info);
        label.setFont(UIManager.getFont("SmallSystemFont"));
        contentPanel.add(label, BorderLayout.NORTH);

        JLabel label2 = new JLabel("<html>If you use these results we would encourage you to cite the following paper:<br>" +
                "Suchard MA, Weiss RE and Sinsheimer JS (2001) Bayesian Selection of Continuous-Time Markov Chain Evolutionary Models." +
                "<i>Mol Biol and Evol</i> <b>18</b>: 1001-1013</html>");
        label2.setFont(UIManager.getFont("SmallSystemFont"));
        contentPanel.add(label2, BorderLayout.SOUTH);

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
        bayesFactorsTable.repaint();
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

		String[] columnNames = {"Trace", "log Marginal Likelihood", "S.E."};

		private DecimalFormat formatter = new DecimalFormat("0.###E0");
		private DecimalFormat formatter2 = new DecimalFormat("####0.###");

        private int columnCount;

        public BayesFactorsModel(boolean hasErrors) {
            this.columnCount = (hasErrors ? 3 : 2);
		}

		public int getColumnCount() {
            return columnCount + marginalLikelihoods.size();
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
			} else if (columnCount > 2 && col == 2 ) {
                return " +/- " + formatter2.format(marginalLikelihoods.get(row).getBootstrappedSE());
			} else {
				if (col - columnCount != row) {
					double lnML1 = marginalLikelihoods.get(row).getLogMarginalLikelihood();
					double lnML2 = marginalLikelihoods.get(col - columnCount).getLogMarginalLikelihood();
					return formatter2.format(lnML1 - lnML2);
				} else {
					return "-";
				}
			}
		}

		public String getColumnName(int column) {
			if (column < columnCount) {
				return columnNames[column];
			}
			return marginalLikelihoods.get(column - columnCount).getTraceName();
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