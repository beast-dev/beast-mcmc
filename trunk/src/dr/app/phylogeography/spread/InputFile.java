package dr.app.phylogeography.spread;

import dr.evolution.tree.Tree;

import javax.swing.*;
import java.io.File;

import jam.util.IconUtils;

public class InputFile implements MultiLineTableCellContent {
    private final Icon treeIcon = IconUtils.getIcon(this.getClass(), "images/tree.png");

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

    public Icon getTableCellIcon() {
        return treeIcon;
    }

    public String getTableCellContent() {
        StringBuilder sb = new StringBuilder();
        sb.append("<html>");
        sb.append("<b>").append(file.getName()).append(": </b>").append(type.name()).append("<br>");
        sb.append("<small>Tip count: 200 | Most recent tip: ").append("2008.9").append("</small><br>");
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

    private Tree tree;
}