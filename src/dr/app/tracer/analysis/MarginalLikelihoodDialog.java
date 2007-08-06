package dr.app.tracer.analysis;

import dr.inference.trace.TraceList;
import org.virion.jam.framework.DocumentFrame;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.util.LongTask;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class MarginalLikelihoodDialog {

    private JFrame frame;

    private JComboBox likelihoodCombo;

    private String[] likelihoodGuesses = {
            "likelihood", "treelikelihood", "lnl", "lik"
    };

    private String likelihoodTrace = "None selected";

    private OptionsPanel optionPanel;

    public MarginalLikelihoodDialog(JFrame frame) {
        this.frame = frame;

        likelihoodCombo = new JComboBox();

        optionPanel = new OptionsPanel(12, 12);
    }

    private int findArgument(JComboBox comboBox, String argument) {
        for (int i = 0; i < comboBox.getItemCount(); i++) {
            String item = ((String) comboBox.getItemAt(i)).toLowerCase();
            if (item.indexOf(argument) != -1) return i;
        }
        return -1;
    }

    public String getLikelihoodTrace() {
        return likelihoodTrace;
    }

    public int showDialog(TraceList traceList) {

        setArguments();

        for (int j = 0; j < traceList.getTraceCount(); j++) {
            String statistic = traceList.getTraceName(j);
            likelihoodCombo.addItem(statistic);
        }
        int index = -1;
        for (String guess : likelihoodGuesses) {
            if (index != -1) break;

            index = findArgument(likelihoodCombo, guess);
        }
        if (index == -1) index = 0;
        likelihoodCombo.setSelectedIndex(index);

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Marginal Likelihood Analysis");
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        if (result == JOptionPane.OK_OPTION) {
            likelihoodTrace = (String) likelihoodCombo.getSelectedItem();
        }

        return result;
    }

    private void setArguments() {
        optionPanel.removeAll();

        optionPanel.addLabel("Select the trace to analyse:");

        optionPanel.addComponents(new JLabel("Likelihood trace:"), likelihoodCombo);

        optionPanel.addSeparator();

    }

    Timer timer = null;

    public void createMarginalLikelihoodFrame(TraceList traceList, DocumentFrame parent) {

        final MarginalLikelihoodTask analyseTask = new MarginalLikelihoodTask(traceList);

        final ProgressMonitor progressMonitor = new ProgressMonitor(frame,
                "Estimating Marginal Likelihood",
                "", 0, analyseTask.getLengthOfTask());
        progressMonitor.setMillisToPopup(0);
        progressMonitor.setMillisToDecideToPopup(0);

        timer = new Timer(1000, new ActionListener() {
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

    class MarginalLikelihoodTask extends LongTask {

        TraceList traceList;

        private int lengthOfTask = 0;
        private int current = 0;
        private String message;

        public MarginalLikelihoodTask(TraceList traceList) {
            this.traceList = traceList;

        }

        public int getCurrent() {
            return current;
        }

        public int getLengthOfTask() {
            return lengthOfTask;
        }

        public String getDescription() {
            return "Estimating marginal likelihood...";
        }

        public String getMessage() {
            return message;
        }

        public Object doWork() {

            current = 0;

            int index = traceList.getTraceIndex(likelihoodTrace);
            double[] likelihoods = new double[traceList.getStateCount()];
            traceList.getValues(index, likelihoods);

            throw new UnsupportedOperationException("Not implemented yet...");

            // return null;
        }

    }

}