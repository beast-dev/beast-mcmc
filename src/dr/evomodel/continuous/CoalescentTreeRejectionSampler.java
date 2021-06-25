package dr.evomodel.continuous;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.coalescent.demographicmodel.DemographicModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * An independent coalescent sampler that rejects samples if they do not comply with any specified monophyletic constraints, based on the coalescent simulator code.
 * This sampler draws from the coalescent prior and accepts each draw with probability 1, and hence should not be used for posterior inference.
 *
 * @author Kanika Nahata (knahata15@gmail.com)
 *
 */

public class CoalescentTreeRejectionSampler extends SimpleMCMCOperator implements GibbsOperator {
    public static final String OPERATOR_NAME = "IterativeGibbsIndependentCoalescentOperator";
    public static final String HEIGHT = "height";
    private DefaultTreeModel treeModel;
    private DemographicModel demoModel;
    private TaxonList allTaxa;
    private List<Set> subtree_nodes;
    private final CoalescentSimulator simulator = new CoalescentSimulator();

    public CoalescentTreeRejectionSampler(TaxonList allTaxa, List<Set> subtree_nodes, DefaultTreeModel treeModel, DemographicModel demoModel, double weight) {

        this.allTaxa = allTaxa;
        this.subtree_nodes = subtree_nodes;
        this.treeModel = treeModel;
        this.demoModel = demoModel;
        setWeight(weight);

    }

    @Override
    public void setPathParameter(double beta) {
        //do nothing
    }
    public String getPerformanceSuggestion() {
        return "";
    }
    @Override
    public String getOperatorName() {
        return "CoalescentTreeRejectionSampler";
    }
    public int getStepCount() {
        return 1;
    }

    @Override
    public double doOperation() {
        Tree newTree = simulateRecursiveTree();
        treeModel.beginTreeEdit();
        treeModel.adoptTreeStructure(newTree);
        treeModel.endTreeEdit();
        return 0;
    }

    private Tree simulateRecursiveTree(){
        // simulate tree with all taxa
        Tree newTree = simulator.simulateTree(allTaxa, demoModel);
        // check if all taxa in subtree are also together in tree
        // recursive
        for (int i = 0; i < subtree_nodes.size(); i++) {
            if (!TreeUtils.isMonophyletic(newTree, subtree_nodes.get(i))){
                return simulateRecursiveTree();
            }
        }
        return newTree;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return OPERATOR_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
            DefaultTreeModel treeModel = (DefaultTreeModel) xo.getChild(TreeModel.class);
            DemographicModel demoModel = (DemographicModel) xo.getChild(DemographicModel.class);

            List<TaxonList> allTaxa = new ArrayList<TaxonList>();
            List<Tree> subtrees = new ArrayList<Tree>();

            // parse xo to get subtrees and the TaxonList
            for (int i = 0; i < xo.getChildCount(); i++) {
                final Object child = xo.getChild(i);
                    if (child instanceof TreeModel){
                        // added to make sure the TreeModel is not added in subtree
                        }
                    else if (child instanceof Tree) {
                        // CoalescentSimulatorParser returns a tree
                        subtrees.add((Tree) child);
                    }
                    else if (child instanceof TaxonList){
                        allTaxa.add((TaxonList) child);
                        // TODO: don't use this else-if loop and get the taxa (TaxonList) from the treeModel outside the for loop
                        // TODO: convert the tree with all taxa List<Taxon> from treemodel asList() to TaxonList
                    }
            }

            // get leaf nodes of the subtrees to test monophyletic constraints in the simulated tree
            List<Set> subtree_nodes = new ArrayList<Set>();
            for (int i = 0; i < subtrees.size(); i++) {
                subtree_nodes.add(TreeUtils.getLeafSet(subtrees.get(i)));
                System.out.println(subtree_nodes.get(i));
            }

            return new CoalescentTreeRejectionSampler(allTaxa.get(0), subtree_nodes, treeModel, demoModel, weight);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(TreeModel.class),
                new ElementRule(DemographicModel.class),
                new ElementRule(Tree.class, 0, Integer.MAX_VALUE),
                new ElementRule(TaxonList.class, 0, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "This element returns an iterative coalescent sampler used in cases of monophyletic constraints on a tree, disguised as a Gibbs operator, from a demographic model.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

    };


}
