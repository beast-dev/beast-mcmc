/*
 * BetaTreeDiversityStatistic.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.treestat.statistics;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;

import java.io.IOException;
import java.util.Set;

/**
 * @author Alexei Drummond
 * @version $Id: BetaTreeDiversityStatistic.java,v 1.3 2006/05/09 10:24:27 rambaut Exp $
 */
public class BetaTreeDiversityStatistic extends AbstractTreeSummaryStatistic {

    private BetaTreeDiversityStatistic() {
        this.taxonList = null;
    }

    public void setTaxonList(TaxonList taxonList) {
        this.taxonList = taxonList;
    }

    public double[] getSummaryStatistic(Tree tree) {
        if (taxonList == null) {
            return new double[]{1.0};
        }

        double TL = TreeUtils.getTreeLength(tree, tree.getRoot());

        double betaDiversity = (TL - getUniqueBranches(tree, tree.getRoot())) / TL;

        return new double[]{betaDiversity};

    }

    private double getUniqueBranches(Tree tree, NodeRef node) {

        if (tree.isExternal(node)) {
            return tree.getBranchLength(node);
        } else {
            double length = 0;
            if (isUnique(taxonList, tree, node)) {
                length = tree.getBranchLength(node);
                //System.out.println("length = " + length);
            }
            for (int i = 0; i < tree.getChildCount(node); i++) {
                length += getUniqueBranches(tree, tree.getChild(node, i));
            }
            //System.out.println("length of node " + node + " = " + length);
            return length;
        }
    }

    private boolean isUnique(TaxonList taxonList, Tree tree, NodeRef node) {
        Set<String> taxa = TreeUtils.getDescendantLeaves(tree, node);
        int count = 0;
        for (String taxon : taxa) {
            count += (taxonList.getTaxonIndex(taxon) >= 0 ? 1 : 0);
            //System.out.print(taxon + "\t");
        }
        boolean unique = (count == 0) || (count == taxa.size());

        //System.out.println(count + "\t" + unique);

        return unique;
    }

    public String getSummaryStatisticName() {
        if (taxonList != null) {
            return "betaDiversity(" + taxonList.getId() + ")";
        } else {
            return "betaDiversity(null)";
        }
    }

    public String getSummaryStatisticDescription() {
        return "The beta diversity of the given taxonList assuming the remaining taxa are from the second location.";
    }

    public String getSummaryStatisticReference() {
        return FACTORY.getSummaryStatisticReference();
    }

    public boolean allowsPolytomies() {
        return FACTORY.allowsPolytomies();
    }

    public boolean allowsNonultrametricTrees() {
        return FACTORY.allowsNonultrametricTrees();
    }

    public boolean allowsUnrootedTrees() {
        return FACTORY.allowsUnrootedTrees();
    }

    public Category getCategory() {
        return FACTORY.getCategory();
    }

    public static final Factory FACTORY = new Factory() {

        public TreeSummaryStatistic createStatistic() {
            return new BetaTreeDiversityStatistic();
        }

        public String getSummaryStatisticName() {
            return "Beta Diversity";
        }

        public String getSummaryStatisticDescription() {
            return "The ratio of shared diversity to total diversity (tree length) between given taxa and the remainder.";
        }

        public String getSummaryStatisticReference() {
            return "-";
        }

        public boolean allowsPolytomies() {
            return true;
        }

        public boolean allowsNonultrametricTrees() {
            return true;
        }

        public boolean allowsUnrootedTrees() {
            return false;
        }

        public Category getCategory() {
            return Category.GENERAL;
        }

        public boolean allowsWholeTree() {
            return true;
        }

        public boolean allowsCharacter() {
            return false;
        }

        public boolean allowsCharacterState() {
            return false;
        }

        public boolean allowsTaxonList() {
            return true;
        }
    };

    private TaxonList taxonList = null;

    public static void main(String[] arg) throws IOException, Importer.ImportException {
        Tree tree = (new NewickImporter("((A:1,B:1):1,(C:1, D:1):1);")).importNextTree();

        BetaTreeDiversityStatistic statistic = (BetaTreeDiversityStatistic) BetaTreeDiversityStatistic.FACTORY.createStatistic();

        Taxa taxa = new Taxa();
        taxa.addTaxon(new Taxon("A"));
        taxa.addTaxon(new Taxon("C"));

        statistic.setTaxonList(taxa);

        System.out.println(statistic.getSummaryStatistic(tree)[0]);
    }
}

