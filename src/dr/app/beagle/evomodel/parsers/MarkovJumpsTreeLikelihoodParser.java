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
import dr.inference.markovjumps.MarkovJumpsType;

/**
 * @author Marc Suchard
 */

public class MarkovJumpsTreeLikelihoodParser extends AncestralStateTreeLikelihoodParser {

    public static final String MARKOV_JUMP_TREE_LIKELIHOOD = "markovJumpsTreeLikelihood";
    public static final String JUMP_TAG = "jumps";
    public static final String JUMP_TAG_NAME = "jumpTagName";
    public static final String COUNTS = MarkovJumpsType.COUNTS.getText();
    public static final String REWARDS = MarkovJumpsType.REWARDS.getText();
    public static final String SCALE_REWARDS = "scaleRewardsByTime";


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

        boolean scaleRewards = xo.getAttribute(SCALE_REWARDS,true);

        MarkovJumpsBeagleTreeLikelihood treeLikelihood = new MarkovJumpsBeagleTreeLikelihood(
                patternList,
                treeModel,
                branchSiteModel,
                siteRateModel,
                branchRateModel,
                useAmbiguities,
                scalingScheme,
                dataType,
                stateTag,
                substModel
        );

        int registersFound = parseAllChildren(xo, treeLikelihood, dataType.getStateCount(), jumpTag,
                        MarkovJumpsType.COUNTS, false); // For backwards compatibility

        XMLObject cxo = xo.getChild(COUNTS);
        if (cxo != null) {
            registersFound += parseAllChildren(cxo, treeLikelihood, dataType.getStateCount(), jumpTag,
                        MarkovJumpsType.COUNTS, false);
        }

        cxo = xo.getChild(REWARDS);
        if (cxo != null) {
            registersFound += parseAllChildren(cxo, treeLikelihood, dataType.getStateCount(), jumpTag,
                        MarkovJumpsType.REWARDS, scaleRewards);
        }
        
        if (registersFound == 0) { // Some default values for testing
            double[] registration = new double[dataType.getStateCount()*dataType.getStateCount()];
            MarkovJumpsCore.fillRegistrationMatrix(registration,dataType.getStateCount()); // Count all transitions
            Parameter registerParameter = new Parameter.Default(registration);
            registerParameter.setId(jumpTag);
            treeLikelihood.addRegister(registerParameter,
                                       MarkovJumpsType.COUNTS,
                                       false);
        }

        return treeLikelihood;
    }

    private int parseAllChildren(XMLObject xo,
                                 MarkovJumpsBeagleTreeLikelihood treeLikelihood,
                                 int stateCount,
                                 String jumpTag,
                                 MarkovJumpsType type,
                                 boolean scaleRewards) throws XMLParseException {
        int registersFound = 0;
        for(int i = 0; i < xo.getChildCount(); i++) {
            Object obj = xo.getChild(i);
            if (obj instanceof Parameter) {
                Parameter registerParameter = (Parameter) obj;
                if ((type == MarkovJumpsType.COUNTS &&
                     registerParameter.getDimension() != stateCount * stateCount) ||
                    (type == MarkovJumpsType.REWARDS &&
                     registerParameter.getDimension() != stateCount)
                   ) {
                    throw new XMLParseException("Register parameter " + registerParameter.getId() + " is of the wrong dimension");
                }
                if (registerParameter.getId() == null) {
                    registerParameter.setId(jumpTag+(registersFound+1));
                }
                treeLikelihood.addRegister(registerParameter, type, scaleRewards);
                registersFound++;
            }
        }
        return registersFound;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
            AttributeRule.newBooleanRule(TreeLikelihoodParser.USE_AMBIGUITIES, true),
            AttributeRule.newStringRule(RECONSTRUCTION_TAG_NAME, true),
            AttributeRule.newStringRule(JUMP_TAG_NAME, true),
            AttributeRule.newBooleanRule(SCALE_REWARDS,true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(GammaSiteRateModel.class),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(SubstitutionModel.class),
            AttributeRule.newStringRule(TreeLikelihoodParser.SCALING_SCHEME, true),
            new ElementRule(Parameter.class,0,Integer.MAX_VALUE), // For backwards compatibility
            new ElementRule(COUNTS,
                    new XMLSyntaxRule[] {
                            new ElementRule(Parameter.class,0,Integer.MAX_VALUE)
                    },true),
            new ElementRule(REWARDS,
                    new XMLSyntaxRule[] {
                            new ElementRule(Parameter.class,0,Integer.MAX_VALUE)
                    },true),
        };
    }
}