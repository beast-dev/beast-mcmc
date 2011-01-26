package dr.app.treespace;

import dr.util.NumberFormatter;
import dr.app.gui.table.MultiLineTableCellContent;

import javax.swing.*;
import java.io.File;

import jam.util.IconUtils;
import jebl.evolution.trees.RootedTree;

public class InputFile implements MultiLineTableCellContent {
    private final Icon treeIcon = IconUtils.getIcon(this.getClass(), "images/tree.png");
    private final Icon treesIcon = IconUtils.getIcon(this.getClass(), "images/small_trees.png");

    public enum Type {
        POSTERIOR_TREES("Posterior distribution of trees"),
        MODAL_TREE("Single tree"),
        LOG_FILE("Posterior distribution of parameters");

        Type(final String name) {
            this.name = name;
        }

        @Override
        public String toString() {
            return name;
        }

        private final String name;
    }

    InputFile(File file) {
        if (file == null) {
            throw new IllegalArgumentException("File argument to InputFile cannot be null");
        }
        this.file = file;
        this.tree = null;
        this.treeCount = 0;
        this.type = Type.LOG_FILE;
    }

    InputFile(File file, RootedTree tree) {
        if (file == null) {
            throw new IllegalArgumentException("File argument to InputFile cannot be null");
        }
        this.file = file;
        this.tree = tree;
        this.treeCount = 0;
        this.type = Type.MODAL_TREE;
    }

    InputFile(File file, RootedTree tree, int treeCount) {
        if (file == null) {
            throw new IllegalArgumentException("File argument to InputFile cannot be null");
        }
        this.file = file;
        this.tree = tree;
        this.treeCount = treeCount;
        this.type = Type.POSTERIOR_TREES;
    }

    public File getFile() {
        return file;
    }

    public Type getType() {
        return type;
    }

    public RootedTree getTree() {
        return tree;
    }

    public int getTreeCount() {
        return treeCount;
    }

    public void setTreeCount(final int treeCount) {
        this.treeCount = treeCount;
    }

    public int getBurnin() {
        return burnin;
    }

    public void setBurnin(final int burnin) {
        this.burnin = burnin;
    }

    public double getMostRecentSampleDate() {
        return mostRecentSampleDate;
    }

    public void setMostRecentSampleDate(double mostRecentSampleDate) {
        this.mostRecentSampleDate = mostRecentSampleDate;
    }

    public Icon getTableCellIcon() {
        switch (type) {
            case LOG_FILE:
                return null;
            case MODAL_TREE:
                return treeIcon;
            case POSTERIOR_TREES:
                return treesIcon;
        }
        return null;
    }

    private static final NumberFormatter nf = new NumberFormatter(6);

    public String getTableCellContent() {
        StringBuilder sb = new StringBuilder();
        sb.append("<html>");
        sb.append("<b>").append(file.getName()).append(": </b>").append(type.toString()).append("<br>");
        if (type != Type.LOG_FILE) {
            sb.append("<small>Tip count: ").append(tree.getExternalNodes().size());
            if (mostRecentSampleDate != 0.0) {
                if (type == Type.POSTERIOR_TREES) {
                    sb.append(" | Tree count: ").append(treeCount < 0 ? " counting..." : treeCount);
                } else {
                    sb.append(" | Root TMRCA: ").append(nf.format(mostRecentSampleDate - tree.getHeight(tree.getRootNode())));
                }
                sb.append(" | Most recent tip: ").append(mostRecentSampleDate);
            } else {
                if (type == Type.POSTERIOR_TREES) {
                    sb.append(" | Tree count: ").append(treeCount < 0 ? " counting..." : treeCount);
                } else {
                    sb.append(" | Root height: ").append(nf.format(tree.getHeight(tree.getRootNode())));
                }
            }
            if (type == Type.POSTERIOR_TREES) {
                sb.append(" | Burn-in: ").append(burnin);
            }
            sb.append("</small><br>");
        }
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
    private int treeCount;
    private int burnin = 0;

    private RootedTree tree = null;
    private double mostRecentSampleDate = 0.0;
}