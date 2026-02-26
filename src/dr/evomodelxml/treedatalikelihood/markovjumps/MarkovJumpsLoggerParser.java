package dr.evomodelxml.treedatalikelihood.markovjumps;

import dr.evolution.tree.Tree;
import dr.evomodel.coalescent.basta.BastaLikelihood;
import dr.evomodel.treedatalikelihood.markovjumps.CompleteHistoryAddOn;
import dr.evomodel.treedatalikelihood.markovjumps.MarkovJumpRewardAddOn;
import dr.evomodel.treedatalikelihood.preorder.AbstractRealizedDiscreteTraitDelegate;
import dr.evomodel.treelikelihood.MarkovJumpsTraitProvider;
import dr.inference.model.Parameter;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;

import java.util.List;

public class MarkovJumpsParserUtils {

    private static final String COUNTS = "counts";
    private static final String REWARDS = "rewards";

    private static boolean isTotalCounts(Parameter p, int stateCount) {
        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                double expectedValue = i == j ? 0.0 : 1.0;
                if (p.getParameterValue(i * stateCount + j) != expectedValue) {
                    return false;
                }
            }
        }
        return true;
    }

    public static void parseXMLObject(XMLObject xo,
                                      Tree tree,
                                      BastaLikelihood likelihood,
                                      String tag) throws XMLParseException {

        int stateCount = likelihood.getPatternList().getDataType().getStateCount();
        AbstractRealizedDiscreteTraitDelegate traitDelegate = likelihood.getRealizedTraitDelegate();
        CompleteHistoryAddOn completeHistoryAddOn = null;

        XMLObject countXo = xo.getChild(COUNTS);
        if (countXo != null) {

            List<Parameter> countParameters = countXo.getAllChildren(Parameter.class);

            int indexOfTotalCounts = -1;
            for (int i = 0; i < countParameters.size(); ++i) {
                Parameter p = countParameters.get(i);
                if (p.getDimension() != stateCount * stateCount) {
                    throw new XMLParseException("Markov jump parameter has wrong dimension");
                }
                if (isTotalCounts(p, stateCount)) {
                    if (indexOfTotalCounts != -1) {
                        throw new XMLParseException("Can only provide one total-count parameter");
                    }
                    indexOfTotalCounts = i;
                }
            }
            if (indexOfTotalCounts == -1) {
                throw new XMLParseException("Must provide a total-count parameter");
            }

            countParameters.remove(indexOfTotalCounts);

            completeHistoryAddOn = new CompleteHistoryAddOn(tag, tree,
                    likelihood.getEvolutionaryProcessDelegate().getBranchSubstitutionModel().getSubstitutionModels(),
                    likelihood.getSiteRateModel(), traitDelegate,
                    MarkovJumpsTraitProvider.ValueScaling.RAW);

            traitDelegate.registerAddOn(completeHistoryAddOn);

            for (Parameter register : countParameters) {

                if (register.getDimension() != stateCount * stateCount) {
                    throw new XMLParseException("Invalid Markov jump register dimension");
                }

                MarkovJumpRewardAddOn markovJumpRewardAddOn = new MarkovJumpRewardAddOn(tag, tree,
                        traitDelegate, completeHistoryAddOn, register, MarkovJumpsTraitProvider.ValueScaling.RAW);

                traitDelegate.registerAddOn(markovJumpRewardAddOn);
            }
        }

        XMLObject rewardXo = xo.getChild(REWARDS);
        if (rewardXo != null) {
            if (completeHistoryAddOn == null) {
                throw new XMLParseException("Must provide a total-count parameter");
            }

            List<Parameter> rewardParameters = rewardXo.getAllChildren(Parameter.class);
            for (Parameter register : rewardParameters) {

                if (register.getDimension() != stateCount) {
                    throw new XMLParseException("Invalid Markov reward register dimension");
                }

                MarkovJumpRewardAddOn markovJumpRewardAddOn = new MarkovJumpRewardAddOn(tag, tree,
                        traitDelegate, completeHistoryAddOn, register, MarkovJumpsTraitProvider.ValueScaling.RAW);

                traitDelegate.registerAddOn(markovJumpRewardAddOn);
            }
        }
    }
}
