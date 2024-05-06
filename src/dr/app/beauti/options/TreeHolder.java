/*
 * Copyright (c) 2024. Lorem ipsum dolor sit amet, consectetur adipiscing elit.
 * Morbi non lorem porttitor neque feugiat blandit. Ut vitae ipsum eget quam lacinia accumsan.
 * Etiam sed turpis ac ipsum condimentum fringilla. Maecenas magna.
 * Proin dapibus sapien vel ante. Aliquam erat volutpat. Pellentesque sagittis ligula eget metus.
 * Vestibulum commodo. Ut rhoncus gravida arcu.
 */

package dr.app.beauti.options;

import dr.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

public class TreeHolder {
    public TreeHolder(Tree tree, String fileNameStem) {
        this.trees = Collections.singletonList(tree);
        this.fileNameStem = fileNameStem;
    }

    public TreeHolder(Collection<Tree> trees, String fileNameStem) {
        this.trees = new ArrayList<>(trees);
        this.fileNameStem = fileNameStem;
    }

    public List<Tree> getTrees() {
        return trees;
    }

    public String getFileNameStem() {
        return fileNameStem;
    }

    @Override
    public String toString() {
        return fileNameStem;
    }

    private final List<Tree> trees;
    private final String fileNameStem;

}
