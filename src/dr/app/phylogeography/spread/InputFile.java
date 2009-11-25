package dr.app.phylogeography.spread;

import dr.evolution.tree.Tree;
import dr.util.NumberFormatter;

import javax.swing.*;
import java.io.File;

import jam.util.IconUtils;

public class InputFile implements MultiLineTableCellContent {
    private final Icon treeIcon = IconUtils.getIcon(this.getClass(), "images/tree.png");
    private final Icon treesIcon = IconUtils.getIcon(this.getClass(), "images/small_trees.png");

    public enum Type {
        POSTERIOR_TREES,
        MODAL_TREE,
        LOG_FILE
    }

    InputFile(File file, Type type) {
        if (file == null) {
            throw new IllegalArgumentException("File argument to InputFile cannot be null");
        }
        this.file = file;
        this.type = type;
    }

    public File getFile() {
        return file;
    }

    public Type getType() {
        return type;
    }

    public Tree getTree() {
        return tree;
    }

    public void setTree(final Tree tree) {
        this.tree = tree;
    }

    public double getMostRecentSampleDate() {
        return mostRecentSampleDate;
    }

    public void setMostRecentSampleDate(double mostRecentSampleDate) {
        this.mostRecentSampleDate = mostRecentSampleDate;
    }

    public Icon getTableCellIcon() {
        return treesIcon;
    }

    private static final NumberFormatter nf = new NumberFormatter(6);

    public String getTableCellContent() {
        StringBuilder sb = new StringBuilder();
        sb.append("<html>");
        sb.append("<b>").append(file.getName()).append(": </b>").append(type.name()).append("<br>");
        sb.append("<small>Tip count: ").append(tree.getExternalNodeCount());
        if (mostRecentSampleDate != 0.0) {
            sb.append(" | Root TMRCA: ").append(nf.format(mostRecentSampleDate - tree.getNodeHeight(tree.getRoot())));
            sb.append(" | Most recent tip: ").append(mostRecentSampleDate);
        } else {
            sb.append(" | Root height: ").append(nf.format(tree.getNodeHeight(tree.getRoot())));
        }
        sb.append("</small><br>");
        sb.append("</html>");
        return sb.toString();
    }

    public String getToolTipContent() {
        return getTableCellContent();
    }

    @Override
    public String toString() {

        return file.getName();
    }


    private final Type type;
    private final File file;

    private Tree tree = null;
    private double mostRecentSampleDate = 0.0;
}