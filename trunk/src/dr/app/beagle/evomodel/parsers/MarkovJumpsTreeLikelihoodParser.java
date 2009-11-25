package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.sitemodel.BranchSiteModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.app.beagle.evomodel.treelikelihood.MarkovJumpsBeagleTreeLikelihood;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;
import dr.inference.model.Parameter;
import dr.inference.markovjumps.MarkovJumpsCore;

/**
 * @author Marc Suchard
 */

public class MarkovJumpsTreeLikelihoodParser extends AncestralStateTreeLikelihoodParser {

    public static final String MARKOV_JUMP_TREE_LIKELIHOOD = "markovJumpsTreeLikelihood";
    public static final String JUMP_TAG = "jumps";
    public static final String JUMP_TAG_NAME = "jumpTagName";


    public String getParserName() {
        return MARKOV_JUMP_TREE_LIKELIHOOD;
    }

    protected BeagleTreeLikelihood createTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                                        BranchSiteModel branchSiteModel, GammaSiteRateModel siteRateModel,
                                                        BranchRateModel branchRateModel,
                                                        boolean useAmbiguities, PartialsRescalingScheme scalingScheme,
                                                        XMLObject xo) throws XMLParseException {

        SubstitutionModel substModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);

        DataType dataType = substModel.getDataType();

        String stateTag = xo.getAttribute(RECONSTRUCTION_TAG_NAME,RECONSTRUCTION_TAG);
        String jumpTag = xo.getAttribute(JUMP_TAG_NAME, JUMP_TAG);

        Parameter registerMatrixParameter = (Parameter) xo.getChild(Parameter.class);

        if (registerMatrixParameter != null) {
            if (registerMatrixParameter.getDimension() != dataType.getStateCount() * dataType.getStateCount() ) {
                throw new XMLParseException("Matrix "+registerMatrixParameter.getId()+" is of the wrong dimension");
            }
        } else { // Some default values for testing
            int from = 1;
            int to = 2;
            double[] registration = new double[dataType.getStateCount()*dataType.getStateCount()];
            MarkovJumpsCore.fillRegistrationMatrix(registration,from,to,dataType.getStateCount());
            registerMatrixParameter = new Parameter.Default(registration);
        }

        return new MarkovJumpsBeagleTreeLikelihood(
                patternList,
                treeModel,
                branchSiteModel,
                siteRateModel,
                branchRateModel,
                useAmbiguities,
                scalingScheme,
                dataType,
                stateTag,
                substModel,
                registerMatrixParameter,
                jumpTag
        );
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
            AttributeRule.newBooleanRule(TreeLikelihoodParser.USE_AMBIGUITIES, true),
            AttributeRule.newStringRule(RECONSTRUCTION_TAG_NAME, true),
            AttributeRule.newStringRule(JUMP_TAG_NAME, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(GammaSiteRateModel.class),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(SubstitutionModel.class),
            AttributeRule.newStringRule(TreeLikelihoodParser.SCALING_SCHEME, true),
            new ElementRule(Parameter.class,true),
        };
    }
}