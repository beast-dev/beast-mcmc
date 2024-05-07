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
    public TreeHolder(Tree tree, String name, String fileName) {
        this.trees = Collections.singletonList(tree);
        this.name = name;
        this.fileName = fileName;
    }

    public TreeHolder(Collection<Tree> trees, String name, String fileName) {
        this.trees = new ArrayList<>(trees);
        this.name = name;
        this.fileName = fileName;
    }

    public List<Tree> getTrees() {
        return trees;
    }

    public String getName() {
        return name;
    }

    public String getFileName() {
        return fileName;
    }

    @Override
    public String toString() {
        return name;
    }

    private final List<Tree> trees;
    private final String name;
    private final String fileName;

}
