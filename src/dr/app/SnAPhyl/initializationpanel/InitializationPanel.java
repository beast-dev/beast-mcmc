package dr.app.SnAPhyl.initializationpanel;

import dr.app.SnAPhyl.SnAPhylFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.enumTypes.FixRateType;
import dr.app.beauti.enumTypes.StartingTreeType;
import dr.app.beauti.options.BeautiOptions;
import dr.evolution.tree.Tree;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/**
 */
public class InitializationPanel extends BeautiPanel {

	private JComboBox startingTreeCombo = new JComboBox(StartingTreeType.values());
	private JComboBox userTreeCombo = new JComboBox();

	private RealNumberField initRootHeightField = new RealNumberField(Double.MIN_VALUE, Double.MAX_VALUE);

	private SnAPhylFrame frame = null;
	private BeautiOptions options = null;
    private OptionsPanel optionsPanel;

    public InitializationPanel(SnAPhylFrame parent) {
		this.frame = parent;
        
        optionsPanel = new OptionsPanel(12, 18);

        startingTreeCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
//                    	partitionTreeModel.setStartingTreeType( (StartingTreeType) startingTreeCombo.getSelectedItem());

                setupPanel();
            }
        });

        userTreeCombo.addItemListener( new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
                    	fireUserTreeChanged();
                    }
                }
        );

        setupPanel();

        setOpaque(false);
        setLayout(new BorderLayout(0, 0));
        setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));

        add(optionsPanel, BorderLayout.CENTER);
	}

    private void setupPanel() {

		optionsPanel.removeAll();

//		if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN
//    			|| options.clockModelOptions.getRateOptionClockModel() == FixRateType.RELATIVE_TO) {
//			initRootHeightField.setValue(partitionTreeModel.getInitialRootHeight());
//			initRootHeightField.setColumns(10);
//			initRootHeightField.setEnabled(false);
//			addComponentWithLabel("The Estimated Initial Root Height:", initRootHeightField);
//		}


        optionsPanel.addComponentWithLabel("Starting Tree:", startingTreeCombo);

        if (startingTreeCombo.getSelectedItem() == StartingTreeType.USER) {
        	optionsPanel.addComponentWithLabel("Select Tree:", userTreeCombo);
        }

        if (options != null) {
            userTreeCombo.removeAllItems();
            if (options.userTrees.size() == 0) {
                userTreeCombo.addItem("no trees loaded");
                userTreeCombo.setEnabled(false);
            } else {
                for (Tree tree : options.userTrees) {
                    userTreeCombo.addItem(tree.getId());
                }
                userTreeCombo.setEnabled(true);
            }
        }

//		generateTreeAction.setEnabled(options != null && options.dataPartitions.size() > 0);

		validate();
		repaint();
	}


	 private void fireUserTreeChanged() {
//		 partitionTreeModel.setUserStartingTree(getSelectedUserTree(options));
	}

    public void setOptions(BeautiOptions options) {
        this.options = options;

//        if (partitionTreeModel == null) {
//            return;
//        }
//
//        startingTreeCombo.setSelectedItem(partitionTreeModel.getStartingTreeType());
//
//        if (partitionTreeModel.getUserStartingTree() != null) {
//        	userTreeCombo.setSelectedItem(partitionTreeModel.getUserStartingTree().getId());
//        }

        setupPanel();
    }

    public void getOptions(BeautiOptions options) {

    }

    private Tree getSelectedUserTree(BeautiOptions options) {
        String treeId = (String) userTreeCombo.getSelectedItem();
        for (Tree tree : options.userTrees) {
            if (tree.getId().equals(treeId)) {
                return tree;
            }
        }
        return null;
    }

    public JComponent getExportableComponent() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}