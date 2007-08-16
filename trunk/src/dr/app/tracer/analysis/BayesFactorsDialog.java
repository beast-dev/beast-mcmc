package dr.app.tracer.analysis;

import dr.inference.trace.TraceList;
import dr.inference.trace.MarginalLikelihoodAnalysis;
import dr.util.TaskListener;
import org.virion.jam.framework.DocumentFrame;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.util.LongTask;
import org.virion.jam.components.WholeNumberField;

import javax.swing.*;
import javax.swing.Timer;
import javax.swing.border.EmptyBorder;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.*;
import java.util.*;
import java.util.List;

public class BayesFactorsDialog {

    private JFrame frame;

    private JComboBox likelihoodCombo;
    private JCheckBox harmonicOnlyCheck;
    private WholeNumberField bootstrapCountField;

    private String[] likelihoodGuesses = {
            "likelihood", "treelikelihood", "lnl", "lik"
    };

    private String likelihoodTrace = "None selected";

    private OptionsPanel optionPanel;

    public BayesFactorsDialog(JFrame frame) {
        this.frame = frame;

        likelihoodCombo = new JComboBox();
        harmonicOnlyCheck = new JCheckBox("Calculate harmonic mean only (no smoothing)");
        bootstrapCountField = new WholeNumberField(0, Integer.MAX_VALUE);
        bootstrapCountField.setValue(1000);
        
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

    public int showDialog(List<TraceList> traceLists) {

        setArguments();

        List<String> statistics = new ArrayList<String>();
        TraceList tl = traceLists.get(0);
        for (int j = 0; j < tl.getTraceCount(); j++) {
            statistics.add(tl.getTraceName(j));
        }

        for (int i = 1; i < traceLists.size(); i++) {
            tl = traceLists.get(i);
            Set<String> statistics2 = new HashSet<String>();
            for (int j = 0; j < tl.getTraceCount(); j++) {
                statistics2.add(tl.getTraceName(j));
            }

            statistics.retainAll(statistics2);
        }

        if (statistics.size() == 0) {
            JOptionPane.showMessageDialog(frame,
                    "These trace files don't seem to have any traces in common.",
                    "Unable to perform analysis",
                    JOptionPane.ERROR_MESSAGE);
            return JOptionPane.CANCEL_OPTION;
        }

        for (String statistic : statistics) {
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

        optionPanel.addComponent(harmonicOnlyCheck);
        bootstrapCountField.setColumns(12);
        optionPanel.addComponents(new JLabel("Bootstrap replicates:"), bootstrapCountField);

        optionPanel.addSeparator();
    }

    Timer timer = null;

    public void createBayesFactorsFrame(List<TraceList> traceLists, DocumentFrame parent) {

        boolean harmonicOnly = harmonicOnlyCheck.isSelected();
        int bootstrapLength = bootstrapCountField.getValue();

        String info = "log P(" + likelihoodTrace + "|Data): ";
        if (harmonicOnly) {
            info += "harmonic mean";
        } else {
            info += "smoothed estimate";
        }
        if (bootstrapLength > 1) {
            info += " (S.E. estimated using " + bootstrapLength + " bootstrap replicates)";
        }
        
        BayesFactorsFrame frame = new BayesFactorsFrame(parent, "Bayes Factors", info, bootstrapLength > 1);
        frame.initialize();
        frame.setVisible(true);

        final MarginalLikelihoodTask analyseTask = new MarginalLikelihoodTask(frame, traceLists, harmonicOnly, bootstrapLength);

        final ProgressMonitor progressMonitor = new ProgressMonitor(frame,
                "Estimating Marginal Likelihoods",
                "", 0, analyseTask.getLengthOfTask());
        progressMonitor.setMillisToPopup(0);
        progressMonitor.setMillisToDecideToPopup(0);

        timer = new Timer(100, new ActionListener() {
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

        List<TraceList> traceLists;
        BayesFactorsFrame frame;
        boolean harmonicOnly;
        int bootstrapLength;

        private int lengthOfTask = 0;
        private int current = 0;
        private String message;


        public MarginalLikelihoodTask(BayesFactorsFrame frame, List<TraceList> traceLists, boolean harmonicOnly, int bootstrapLength) {
            this.traceLists = traceLists;
            this.harmonicOnly = harmonicOnly;
            this.bootstrapLength = bootstrapLength;
            this.frame = frame;
            lengthOfTask = traceLists.size() * 100;
        }

        public int getCurrent() {
            return current;
        }

        public int getLengthOfTask() {
            return lengthOfTask;
        }

        public String getDescription() {
            return "Estimating marginal likelihoods...";
        }

        public String getMessage() {
            return message;
        }

        public Object doWork() {

            current = 0;
            int i = 0;

            for (TraceList traceList : traceLists) {
                final int offset = i * 100;

                int index = traceList.getTraceIndex(likelihoodTrace);
                double[] likelihoods = new double[traceList.getStateCount()];
                traceList.getValues(index, likelihoods);

                final MarginalLikelihoodAnalysis analysis = new MarginalLikelihoodAnalysis(
                        likelihoods,
                        traceList.getName(), traceList.getBurnIn(),
                        harmonicOnly, bootstrapLength);

                analysis.setTaskListener(new TaskListener() {
                    public void progress(double progress) {
                        current = offset + (int)(progress * 100);
                    }
                });

                // call to force calculation...
                analysis.calculate();

                EventQueue.invokeLater(
                        new Runnable() {
                            public void run() {
                                frame.addMarginalLikelihood(analysis);
                            }
                        });

                i++;
            }

            return null;
        }

    }

}