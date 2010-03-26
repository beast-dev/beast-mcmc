package dr.app.beagle.tools.parsers;

import dr.app.beagle.evomodel.parsers.MarkovJumpsTreeLikelihoodParser;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.CodonLabeling;
import dr.app.beagle.tools.CompleteHistorySimulator;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.inference.markovjumps.MarkovJumpsType;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */

public class CompleteHistorySimulatorParser extends AbstractXMLObjectParser {

    /* standard xml parser stuff follows */
    public static final String HISTORY_SIMULATOR = "completeHistorySimulator";
    public static final String TREE = "tree";
    public static final String REPLICATIONS = "replications";
    public static final String COUNTS = MarkovJumpsTreeLikelihoodParser.COUNTS;
    public static final String REWARDS = MarkovJumpsTreeLikelihoodParser.REWARDS;
    public static final String JUMP_TAG_NAME = MarkovJumpsTreeLikelihoodParser.JUMP_TAG_NAME;
    public static final String JUMP_TAG = MarkovJumpsTreeLikelihoodParser.JUMP_TAG;

    public static final String SYN_JUMPS = "reportSynonymousMutations";
    public static final String NON_SYN_JUMPS = "reportNonSynonymousMutations";
    public static final String SUM_SITES = "sumAcrossSites";

    public String getParserName() {
        return HISTORY_SIMULATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int nReplications = xo.getIntegerAttribute(REPLICATIONS);

        Tree tree = (Tree) xo.getChild(Tree.class);
        GammaSiteRateModel siteModel = (GammaSiteRateModel) xo.getChild(GammaSiteRateModel.class);
        BranchRateModel rateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        if (rateModel == null)
            rateModel = new DefaultBranchRateModel();

        DataType dataType = siteModel.getSubstitutionModel().getDataType();

        String jumpTag = xo.getAttribute(JUMP_TAG_NAME, JUMP_TAG);

        boolean sumAcrossSites = xo.getAttribute(SUM_SITES, false);

        CompleteHistorySimulator history = new CompleteHistorySimulator(tree, siteModel, rateModel, nReplications,
                sumAcrossSites);

        XMLObject cxo = xo.getChild(COUNTS);
        if (cxo != null) {
            MarkovJumpsTreeLikelihoodParser.parseAllChildren(cxo, history, dataType.getStateCount(),
                    jumpTag, MarkovJumpsType.COUNTS, false);
        }

        cxo = xo.getChild(REWARDS);
        if (cxo != null) {
            MarkovJumpsTreeLikelihoodParser.parseAllChildren(cxo, history, dataType.getStateCount(),
                    jumpTag, MarkovJumpsType.REWARDS, false);
        }

        if (dataType instanceof Codons) {
            Codons codons = (Codons) dataType;
            if (xo.getAttribute(SYN_JUMPS, false)) {
                double[] synRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.SYN, codons, false); // use base 61
                Parameter registerParameter = new Parameter.Default(synRegMatrix);
                registerParameter.setId("S");
                history.addRegister(registerParameter, MarkovJumpsType.COUNTS, false);
            }
            if (xo.getAttribute(NON_SYN_JUMPS, false)) {
                double[] nonSynRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.NON_SYN, codons, false); // use base 61
                Parameter registerParameter = new Parameter.Default(nonSynRegMatrix);
                registerParameter.setId("N");
                history.addRegister(registerParameter, MarkovJumpsType.COUNTS, false);
            }
        }

        history.simulate();
        return history;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A SequenceSimulator that generates random sequences for a given tree, sitemodel and branch rate model";
    }

    public Class getReturnType() {
        return Alignment.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Tree.class),
            new ElementRule(GammaSiteRateModel.class),
            new ElementRule(BranchRateModel.class, true),
            AttributeRule.newIntegerRule(REPLICATIONS),
            AttributeRule.newBooleanRule(SYN_JUMPS, true),
            AttributeRule.newBooleanRule(NON_SYN_JUMPS, true),
            AttributeRule.newBooleanRule(SUM_SITES, true),
    };
}