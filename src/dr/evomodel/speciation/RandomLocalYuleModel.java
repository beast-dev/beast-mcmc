package dr.evomodel.speciation;

import dr.evolution.tree.NodeAttributeProvider;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.randomlocalmodel.RandomLocalTreeVariable;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.text.NumberFormat;
import java.util.Locale;

/**
 * This class contains methods that describe a Yule speciation model whose rate of birth changes
 * at different points in the tree.
 *
 * @author Alexei Drummond
 */
public class RandomLocalYuleModel extends SpeciationModel implements NodeAttributeProvider, RandomLocalTreeVariable {

    public static final String YULE_MODEL = "randomLocalYuleModel";
    public static String MEAN_RATE = "meanRate";
    public static String BIRTH_RATE = "birthRates";
    public static String BIRTH_RATE_INDICATORS = "indicators";
    public static String RATES_AS_MULTIPLIERS = "ratesAsMultipliers";

    public RandomLocalYuleModel(Parameter birthRates, Parameter indicators, Parameter meanRate,
                                boolean ratesAsMultipliers, Type units, int dp) {

        super(RandomLocalYuleModel.YULE_MODEL, units);

        this.birthRates = birthRates;

        addParameter(birthRates);
        birthRates.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, birthRates.getDimension()));

        this.indicators = indicators;

        for (int i = 0; i < indicators.getDimension(); i++) {
            indicators.setParameterValueQuietly(i, 0.0);
        }

        addParameter(indicators);

        this.meanRate = meanRate;
        addParameter(meanRate);

        birthRatesAreMultipliers = ratesAsMultipliers;

        format.setMaximumFractionDigits(dp);
    }

    public final String getVariableName() {
        return "birthRate";
    }

    public final String getSelectionIndicatorName() {
        return "birthRateIndicator";
    }

    public double getVariable(TreeModel tree, NodeRef node) {
        return tree.getNodeTrait(node, getVariableName());
    }

    public boolean isVariableSelected(TreeModel tree, NodeRef node) {
        return (int) Math.round(tree.getNodeTrait(node, getSelectionIndicatorName())) == 1;
    }

    /**
     * @param tree the tree
     * @param node the node to retrieve the birth rate of
     * @return the birth rate of the given node;
     */
    double getBirthRate(TreeModel tree, NodeRef node) {

        double birthRate;
        if (!tree.isRoot(node)) {

            double parentRate = getBirthRate(tree, tree.getParent(node));
            if (isVariableSelected(tree, node)) {
                birthRate = getVariable(tree, node);
                if (birthRatesAreMultipliers) {
                    birthRate *= parentRate;
                } else {
                    throw new RuntimeException("Rates must be multipliers in current implementation! " +
                            "Otherwise root rate might be ignored");
                }
            } else {
                birthRate = parentRate;
            }
        } else {
            birthRate = meanRate.getParameterValue(0);
        }
        return birthRate;
    }

    public String[] getNodeAttributeLabel() {
        return new String[]{"I", "b"};
    }

    public String[] getAttributeForNode(Tree tree, NodeRef node) {

        String rateString = format.format(getBirthRate((TreeModel) tree, node));

        if (tree.isRoot(node)) {
            return new String[]{"0", rateString};
        }

        return new String[]{(isVariableSelected((TreeModel) tree, node) ? "1" : "0"), rateString};
    }

    //
    // functions that define a speciation model
    //
    public double logTreeProbability(int taxonCount) {
        return 0.0;
    }

    //
    // functions that define a speciation model
    //
    public double logNodeProbability(Tree tree, NodeRef node) {


        if (tree.isRoot(node)) {
            return 0.0;
        } else {

            double lambda = getBirthRate((TreeModel) tree, node);

            double branchLength = tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node);

            double logP = -lambda * branchLength;

            if (tree.isExternal(node)) logP += Math.log(lambda);

            return logP;
        }
    }

    public boolean includeExternalNodesInLikelihoodCalculation() {
        return true;
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public org.w3c.dom.Element createElement(org.w3c.dom.Document d) {
        throw new RuntimeException("createElement not implemented");
    }

    /**
     * Parses an element from an DOM document into a SpeciationModel. Recognises
     * RandomLocalYuleModel.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return RandomLocalYuleModel.YULE_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Type units = XMLParser.Utils.getUnitsAttr(xo);

            XMLObject cxo = (XMLObject) xo.getChild(RandomLocalYuleModel.BIRTH_RATE);
            Parameter brParameter = (Parameter) cxo.getChild(Parameter.class);

            cxo = (XMLObject) xo.getChild(RandomLocalYuleModel.BIRTH_RATE_INDICATORS);
            Parameter indicatorsParameter = (Parameter) cxo.getChild(Parameter.class);

            Parameter meanRate = (Parameter) xo.getElementFirstChild(RandomLocalYuleModel.MEAN_RATE);

            boolean ratesAsMultipliers = xo.getBooleanAttribute(RATES_AS_MULTIPLIERS);

            int dp = xo.getAttribute("dp", 4);

            return new RandomLocalYuleModel(brParameter, indicatorsParameter, meanRate, ratesAsMultipliers, units, dp);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A speciation model of a Yule process whose rate can evolve down the tree.";
        }

        public Class getReturnType() {
            return RandomLocalYuleModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(RandomLocalYuleModel.BIRTH_RATE,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(RandomLocalYuleModel.BIRTH_RATE_INDICATORS,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(RandomLocalYuleModel.MEAN_RATE,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                AttributeRule.newBooleanRule(RATES_AS_MULTIPLIERS),
                AttributeRule.newIntegerRule("dp", true),
                XMLUnits.SYNTAX_RULES[0]
        };
    };

    //Protected stuff
    Parameter birthRates;
    Parameter indicators;

    Parameter meanRate;

    boolean birthRatesAreMultipliers = false;

    NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);
}
