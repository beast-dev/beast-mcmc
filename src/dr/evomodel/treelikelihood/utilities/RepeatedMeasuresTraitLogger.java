package dr.evomodel.treelikelihood.utilities;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitSimulator;
import dr.evomodel.treedatalikelihood.continuous.WishartStatisticsWrapper;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.CompoundParameter;
import dr.xml.*;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;


public class RepeatedMeasuresTraitLogger implements Loggable, Reportable {

    private final RepeatedMeasuresTraitSimulator traitSimulator;
    private final RepeatedMeasuresTraitDataModel traitModel;
    private final CompoundParameter traitParameter;
    private final TreeDataLikelihood traitLikelihood;
    private final TreeTrait tipTrait;
    private double[] missingTaxa;

    RepeatedMeasuresTraitLogger(RepeatedMeasuresTraitDataModel traitModel, TreeDataLikelihood traitLikelihood) {

        this.traitSimulator = new RepeatedMeasuresTraitSimulator(traitModel, traitLikelihood);
        this.traitModel = traitModel;
        this.traitParameter = traitModel.getParameter();
        this.traitLikelihood = traitLikelihood;
        this.tipTrait = traitLikelihood.getTreeTrait(REALIZED_TIP_TRAIT + "." + traitModel.getTraitName());

    }

    @Override
    public LogColumn[] getColumns() {


        int m = traitModel.getMissingIndices().size();
        int n = traitParameter.getParameterCount();
        int dim = traitModel.getTraitDimension();

        LogColumn[] columns = new LogColumn[m];

        for (int i = 0; i < m; i++) {

            int index = traitModel.getMissingIndices().get(i);

            final int taxon = index / dim;
            final int trait = index - taxon * dim;
            final int currentInd = i;

            columns[i] = new LogColumn.Abstract(traitParameter.getParameterName() + "." +
                    traitParameter.getParameter(taxon).getParameterName() + "." +
                    Integer.toString(trait)) {
                @Override
                protected String getFormattedValue() {

                    if (currentInd == 0) {
                        resample();
                    }

                    return Double.toString(traitParameter.getParameterValue(trait, taxon));

                }
            };
        }

        return columns;
    }

    private void resample() {

        ContinuousDataLikelihoodDelegate delegate = (ContinuousDataLikelihoodDelegate)
                traitLikelihood.getDataLikelihoodDelegate();

        delegate.fireModelChanged(); //Force new sample

        double[] tips = (double[]) tipTrait.getTrait(traitLikelihood.getTree(), null);

        traitSimulator.simulateMissingData(tips);
    }

    @Override
    public String getReport() {
        return null;
    }

    public static final String RM_TRAIT_LOGGER = "repeatedMeasuresTraitLogger";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeDataLikelihood dataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
            RepeatedMeasuresTraitDataModel traitModel =
                    (RepeatedMeasuresTraitDataModel) xo.getChild(RepeatedMeasuresTraitDataModel.class);

            return new RepeatedMeasuresTraitLogger(traitModel, dataLikelihood);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(RepeatedMeasuresTraitDataModel.class),
                new ElementRule(TreeDataLikelihood.class)
        };

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return RepeatedMeasuresTraitLogger.class;
        }

        @Override
        public String getParserName() {
            return RM_TRAIT_LOGGER;
        }
    };
}
