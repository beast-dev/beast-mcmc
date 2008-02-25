package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Birth Death model based on  S. Nee equation (21) in his 1994 paper
 * "The Reconstructed Evolutionary Process"
 *
 * The inference is directly on b-d (positive) and d/b (constrained between 0 and 1)
 *
 * Not fully tested yet but seem to work on the two examples in examples/treePriors.
 *
 * @author joseph
 *         Date: 24/02/2008
 */
public class BirthDeathNee94Model extends SpeciationModel {

    public static final String BIRTH_DEATH_MODEL = "Nee94birthDeathModel";
    public static String BIRTHDIFF_RATE = "birthMinusDeathRate";
    public static String RELATIVE_DEATH_RATE = "relativeDeathRate";

    private Parameter relativeDeathRateParameter;
    private Parameter birthDiffRateParameter;


    public BirthDeathNee94Model(Parameter birthDiffRateParameter, Parameter relativeDeathRateParameter, Type units) {

        super(BIRTH_DEATH_MODEL, units);

        this.birthDiffRateParameter = birthDiffRateParameter;
        addParameter(birthDiffRateParameter);
        birthDiffRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.relativeDeathRateParameter = relativeDeathRateParameter;
        addParameter(relativeDeathRateParameter);
        relativeDeathRateParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));
    }

    public double getR() {
        return birthDiffRateParameter.getParameterValue(0);
    }

    public double getA() {
        return relativeDeathRateParameter.getParameterValue(0);
    }

    public double logTreeProbability(int taxonCount) {
        /* leave log(gamma(n)) assuming the likelihood is used inside the MCMC on the "same"
        * tree, i.e. one which has fixed taxa  */
        return /*loggamma(taxonXount)  + */
        (taxonCount-2) * Math.log(getR()) + taxonCount * Math.log(1 - getA());
    }

    public double logNodeProbability(Tree tree, NodeRef node) {
        final double height = tree.getNodeHeight(node);
        final double x = getR() * height;  
        final double z = Math.exp(x) - getA();
        double l = -2*Math.log(z);

        if( tree.getRoot() != node ) {
            l += x;
        }
        return l;
    }

    public boolean includeExternalNodesInLikelihoodCalculation() {
        return false;
    }

    /**
     * Parses an element from an DOM document into a SpeciationModel. Recognises
     * birthDeathModel.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return BIRTH_DEATH_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Type units = XMLParser.Utils.getUnitsAttr(xo);

            Parameter birthParameter = (Parameter) xo.getSocketChild(BIRTHDIFF_RATE);
            Parameter deathParameter = (Parameter) xo.getSocketChild(RELATIVE_DEATH_RATE);

            return new BirthDeathNee94Model(birthParameter, deathParameter, units);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Nee (1994) model of speciation (equation 21).";
        }

        public Class getReturnType() {
            return BirthDeathNee94Model.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(BIRTHDIFF_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(RELATIVE_DEATH_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                XMLUnits.SYNTAX_RULES[0]
        };
    };
}
