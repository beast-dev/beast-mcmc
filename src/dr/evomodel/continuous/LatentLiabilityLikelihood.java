package dr.evomodel.continuous;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.TwoStates;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.xml.*;

import java.util.List;
import java.util.logging.Logger;


/**
 * A class to model multivariate {0,1} data as realizations from a latent (liability) multivariate Brownian diffusion
 *
 * @author Marc A. Suchard
 * @author Joe Felsenstein
 * @version $Id$
 */

public class LatentLiabilityLikelihood extends AbstractModelLikelihood {

    public final static String LATENT_LIABILITY_LIKELIHOOD = "latentLiabilityLikelihood";

    public LatentLiabilityLikelihood(TreeModel treeModel, PatternList patternList, CompoundParameter tipTraitParameter) {
        super(LATENT_LIABILITY_LIKELIHOOD);
        this.treeModel = treeModel; // TODO Do we really need the tree?
        this.patternList = patternList;
        this.tipTraitParameter = tipTraitParameter;

        addVariable(tipTraitParameter);

        setTipDataValuesForAllNodes();
    }

    private void setTipDataValuesForAllNodes() {
        if (tipData == null) {
            tipData = new boolean[treeModel.getExternalNodeCount()][patternList.getPatternCount()];
        }
        
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            NodeRef node = treeModel.getExternalNode(i);
            String id = treeModel.getTaxonId(i);
            int index = patternList.getTaxonIndex(id);
            setTipDataValuesForNode(node, index);
        }
    }

    private void setTipDataValuesForNode(NodeRef node, int indexFromPatternList) {
        // Set tip data values
        int index = node.getNumber();
        if (index != indexFromPatternList) {
            throw new RuntimeException("Need to figure out the indexing");
        }

//        Parameter oneTipTraitParameter = tipTraitParameter.getParameter(index);

        for (int datum = 0; datum < patternList.getPatternCount(); ++datum) {
            tipData[index][datum] = patternList.getPattern(datum)[index]  == 1;
//            System.err.println("Data = " + tipData[index][datum] + " : " + oneTipTraitParameter.getParameterValue(datum));
        }
//        double[] traitValue = traitParameter.getParameter(index).getParameterValues();
//        if (traitValue.length < dim) {
//            throw new RuntimeException("The trait parameter for the tip with index, " + index + ", is too short");
//        }
//        System.arraycopy(traitValue, 0, meanCache, dim * index, dim);
    }

//    private void ensureTipTraitsAreDataCompatible() {
//        patternList.getPattern()
//    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = true;
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    public void makeDirty() {
        likelihoodKnown = false;
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = computeLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public String toString() {
        return getClass().getName() + "(" + getLogLikelihood() + ")";
    }    

    // This function can be overwritten to implement other sampling densities, i.e. discrete ranks
    protected double computeLogLikelihood() {
        boolean valid = true;
        for (int tip = 0; tip < tipData.length && valid; ++tip) {
            Parameter oneTipTraitParameter = tipTraitParameter.getParameter(tip);
            boolean[] data = tipData[tip];
            for (int index = 0; index < data.length && valid; ++index) {
                boolean datum = data[index];
                double trait = oneTipTraitParameter.getParameterValue(index);
                boolean positive = trait > 0.0;
                valid =  datum == positive;
//                if (!valid) {
//                    System.err.println("Invalid: " + data[index] + " and " + trait + " at tip " + tip + " entry " + datum);
//                }
            }
        }
        if (valid) {
            return 0.0;
        } else {
            return Double.NEGATIVE_INFINITY;
        }
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
//        public final static String FILE_NAME = "fileName";
        public final static String TIP_TRAIT = "tipTrait";
//        public final static String LOCATIONS = "locations";
//        public static final String MDS_DIMENSION = "mdsDimension";
//        public static final String MDS_PRECISION = "mdsPrecision";

        public String getParserName() {
            return LATENT_LIABILITY_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            AbstractMultivariateTraitLikelihood traitLikelihood = (AbstractMultivariateTraitLikelihood)
                    xo.getChild(AbstractMultivariateTraitLikelihood.class);
            PatternList patternList = (PatternList) xo.getChild(PatternList.class);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            CompoundParameter tipTraitParameter = (CompoundParameter) xo.getElementFirstChild(TIP_TRAIT);

            int numTaxa = treeModel.getTaxonCount();
            int numData = traitLikelihood.getNumData();
            int dimTrait = traitLikelihood.getDimTrait();

            if (tipTraitParameter.getDimension() != numTaxa * numData * dimTrait) {
                throw new XMLParseException("Tip trait parameter is wrong dimension in latent liability model");
            }

            if (!(patternList.getDataType() instanceof TwoStates)) {
                throw new XMLParseException("Latent liability model currently only works for binary data");
            }

            if (patternList.getPatternCount() != numData * dimTrait) {
                throw new XMLParseException("Binary data is wrong dimension in latent liability model");
            }

            return new LatentLiabilityLikelihood(treeModel, patternList, tipTraitParameter);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Provides the likelihood of a latent liability model on multivariate-binary trait data";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(AbstractMultivariateTraitLikelihood.class, "The model for the latent random variables"),                                
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree"),
                new ElementRule(PatternList.class, "The binary tip data"),
                new ElementRule(TreeModel.class, "The tree model"),
        };

        public Class getReturnType() {
            return LatentLiabilityLikelihood.class;
        }
    };


    private TreeModel treeModel;
    private PatternList patternList;
    private CompoundParameter tipTraitParameter;

    private boolean[][] tipData;

    private int numTaxa;
    
    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

}