package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodelxml.branchratemodel.MixtureModelBranchRatesParser;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.*;

/**
 * @author Wai Lok Sibon Li
 * @version $Id: DiscretizedBranchRates.java,v 1.11 2009/12/01 17:44:30 rambaut Exp $
 */
public class MixtureModelBranchRates extends AbstractBranchRateModel {

    private final ParametricDistributionModel[] distributionModels;

    // The rate categories of each branch
    final TreeParameterModel rateCategoryQuantiles;

    private Parameter distributionIndexParameter;  

    private final double[] rates;
    private boolean useQuantilesForRates = true;
    private boolean normalize = false;
    private double normalizeBranchRateTo = Double.NaN;
    private double scaleFactor = 1.0;
    private TreeModel treeModel;
    private Tree tree;

    public MixtureModelBranchRates(
            TreeModel tree,
            Parameter rateCategoryQuantilesParameter,
            ParametricDistributionModel[] models,
            Parameter distributionIndexParameter) {
        this(tree, rateCategoryQuantilesParameter, models, distributionIndexParameter, true, false, Double.NaN);
    }

    public MixtureModelBranchRates(
            TreeModel tree,
            Parameter rateCategoryQuantilesParameter,
            ParametricDistributionModel[] models,
            Parameter distributionIndexParameter,
            boolean normalize,
            double normalizeBranchRateTo) {
        this(tree, rateCategoryQuantilesParameter, models, distributionIndexParameter, true, normalize, normalizeBranchRateTo);
    }

    public MixtureModelBranchRates(
            TreeModel tree,
            Parameter rateCategoryQuantilesParameter,
            ParametricDistributionModel[] models,
            Parameter distributionIndexParameter,
            boolean useQuantilesForRates) {
        this(tree, rateCategoryQuantilesParameter, models, distributionIndexParameter, useQuantilesForRates, false, Double.NaN);
    }

    public MixtureModelBranchRates(
            TreeModel tree,
            Parameter rateCategoryQuantilesParameter,
            ParametricDistributionModel[] models,
            Parameter distributionIndexParameter,
            boolean useQuantilesForRates,
            boolean normalize,
            double normalizeBranchRateTo) {
        super(MixtureModelBranchRatesParser.MIXTURE_MODEL_BRANCH_RATES);

        this.useQuantilesForRates = useQuantilesForRates;

        this.rateCategoryQuantiles = new TreeParameterModel(tree, rateCategoryQuantilesParameter, false);

        rates = new double[tree.getNodeCount()];

        this.normalize = normalize;

        this.treeModel = tree;
        this.distributionModels = models;
        this.normalizeBranchRateTo = normalizeBranchRateTo;

        this.tree = new SimpleTree(tree);

        this.distributionIndexParameter = distributionIndexParameter;
        addVariable(this.distributionIndexParameter);

        //Force the boundaries of rateCategoryParameter to match the category count
        //d Parameter.DefaultBounds bound = new Parameter.DefaultBounds(categoryCount - 1, 0, rateCategoryParameter.getDimension());
        //d rateCategoryParameter.addBounds(bound);
        //rateCategoryQuantilesParameter.;



        Parameter.DefaultBounds bound = new Parameter.DefaultBounds(1.0, 0.0, rateCategoryQuantilesParameter.getDimension());
        rateCategoryQuantilesParameter.addBounds(bound);

        
        Parameter.DefaultBounds bound2 = new Parameter.DefaultBounds(models.length - 1, 0.0, 1);
        distributionIndexParameter.addBounds(bound2);
        distributionIndexParameter.setParameterValue(0, 0);

        //Parameter distributionIndexParameter;


        for (ParametricDistributionModel distributionModel : distributionModels) {
            addModel(distributionModel);
        }
        // AR - commented out: changes to the tree are handled by model changed events fired by rateCategories
//        addModel(tree);
        //d addModel(rateCategories);

        addModel(rateCategoryQuantiles);

        //addModel(treeModel); // Maybe
        // AR - commented out: changes to rateCategoryParameter are handled by model changed events fired by rateCategories
//        addVariable(rateCategoryParameter);

        if (normalize) {
            tree.addModelListener(new ModelListener() {

                public void modelChangedEvent(Model model, Object object, int index) {
                    computeFactor();
                }

                public void modelRestored(Model model) {
                    computeFactor();
                }
            });
        }

        setupRates();
    }

    // compute scale factor
    private void computeFactor() {

        //scale mean rate to 1.0 or separate parameter

        double treeRate = 0.0;
        double treeTime = 0.0;

        //normalizeBranchRateTo = 1.0;
        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);
            if (!treeModel.isRoot(node)) {
//d                int rateCategory = (int) Math.round(rateCategories.getNodeValue(treeModel, node));
//d                 treeRate += rates[rateCategory] * treeModel.getBranchLength(node);
                treeTime += treeModel.getBranchLength(node);

// d              System.out.println("rates and time\t" + rates[rateCategory] + "\t" + treeModel.getBranchLength(node));
            }
        }
        //treeRate /= treeTime;

        scaleFactor = normalizeBranchRateTo / (treeRate / treeTime);
        System.out.println("scaleFactor\t\t\t\t\t" + scaleFactor);
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        //System.out.println("if you dont know: " + model.getClass().getName());
        for (ParametricDistributionModel distributionModel : distributionModels) {
            if (model == distributionModel) {
                setupRates();
                fireModelChanged();
            } //else if (model == rateCategories) {
            // AR - commented out: if just the rate categories have changed the rates will be the same
            //            setupRates();
            //  fireModelChanged(null, index);
            //}
        }
        if (model == rateCategoryQuantiles) {
            setupRates();   // Maybe
            //rateCategories.fireModelChanged();
            fireModelChanged(null, index);
        }
        /*else if(model == distributionIndexParameter) { // Not a model
            setupRates();
        }*/
        /*else if (model == treeModel) {
            setupRates(); // Maybe
        }*/
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if(variable==distributionIndexParameter) {
            //System.out.print(Math.round(distributionIndexParameter.getValue(0)));
            setupRates();
            fireModelChanged();
        }
        // AR - commented out: changes to rateCategoryParameter are handled by model changed events
        //setupRates();   // Maybe
    }

    protected void storeState() {
        //setupRates();   // Maybe
    }

    protected void restoreState() {
        setupRates();
    }

    protected void acceptState() {
        //setupRates();   // Maybe
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {
        assert !tree.isRoot(node) : "root node doesn't have a rate!";
        //int rateCategory = (int) Math.round(rateCategories.getNodeValue(tree, node));
        return rates[node.getNumber()] * scaleFactor;
    }

    /**
     * Calculates the actual rates corresponding to the category indices.
     */
    protected void setupRates() {


        //System.out.println("BRRRTTZZZ " + distributionIndexParameter.getValue(0));
        for (int i = 0; i < tree.getNodeCount(); i++) {
            //rates[i] = distributionModel.quantile(rateCategoryQuantiles.getNodeValue(
            // rateCategoryQuantiles.getTreeModel(), rateCategoryQuantiles.getTreeModel().getNode(i) ));
            if (!tree.isRoot(tree.getNode(i))) {

                if(useQuantilesForRates) {  /* Using quantiles to represent rates */
                    rates[tree.getNode(i).getNumber()] = distributionModels[(int) Math.round(distributionIndexParameter.getValue(0))]
                            .quantile(rateCategoryQuantiles.getNodeValue(tree, tree.getNode(i)));
                }
                else { /* Not using quantiles to represent rates. This is practically useless for anything else other than simulation */
                    rates[tree.getNode(i).getNumber()] = rateCategoryQuantiles.getNodeValue(tree, tree.getNode(i));
                }
            }
        }
        /*System.out.print(distributionModels[(int) Math.round(distributionIndexParameter.getValue(0))].getClass().getName() + "\t" + (int) Math.round(distributionIndexParameter.getValue(0)) + "\t" + rates[1] + "\t" + rateCategoryQuantiles.getNodeValue(tree, tree.getNode(1)));// + "\t" + distributionModels[(int) Math.round(distributionIndexParameter.getValue(0))].);
        if(distributionModels[(int) Math.round(distributionIndexParameter.getValue(0))].getClass().getName().equals("dr.inference.distribution.LogNormalDistributionModel")) {
            LogNormalDistributionModel lndm = (LogNormalDistributionModel) distributionModels[(int) Math.round(distributionIndexParameter.getValue(0))];
            System.out.println("\t" + lndm.getS());
        }
        else if (distributionModels[(int) Math.round(distributionIndexParameter.getValue(0))].getClass().getName().equals("dr.inference.distribution.InverseGaussianDistributionModel")) {
            InverseGaussianDistributionModel lndm = (InverseGaussianDistributionModel) distributionModels[(int) Math.round(distributionIndexParameter.getValue(0))];
            System.out.println("\t" + lndm.getS());
        }*/
        if (normalize) computeFactor();
    }

}
