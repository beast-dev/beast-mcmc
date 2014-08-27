package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;
import org.apache.commons.math.special.Gamma;

import java.util.logging.Logger;

/**
 * Birth Death model based on Gernhard 2008  "The conditioned reconstructed process"
 * doi:10.1016/j.jtbi.2008.04.005 (http://dx.doi.org/10.1016/j.jtbi.2008.04.005)
 * <p/>
 * This derivation conditions directly on fixed N taxa.
 * <p/>
 * The inference is directly on b-d (strictly positive) and d/b (constrained in [0,1))
 * <p/>
 * Vefified using simulated trees generated by Klass tree sample. (http://www.klaashartmann.com/treesample/)
 *
 * @author joseph
 *         Date: 24/02/2008
 */
public class BirthDeathGernhard08Model extends SpeciationModel {

    public static final String BIRTH_DEATH_MODEL = "birthDeathModel";
    public static String BIRTHDIFF_RATE = "birthMinusDeathRate";
    public static String RELATIVE_DEATH_RATE = "relativeDeathRate";

    public static String BIRTHDIFF_RATE_PARAM_NAME = "birthDeath.BminusDRate";
    public static String RELATIVE_DEATH_RATE_PARAM_NAME = "birthDeath.DoverB";

    private Parameter relativeDeathRateParameter;
    private Parameter birthDiffRateParameter;


    public BirthDeathGernhard08Model(Parameter birthDiffRateParameter, Parameter relativeDeathRateParameter, int units) {

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
        return Gamma.logGamma(taxonCount + 1) +
                (taxonCount - 1) * Math.log(getR()) + taxonCount * Math.log(1 - getA());
    }

    public double logNodeProbability(Tree tree, NodeRef node) {
        final double height = tree.getNodeHeight(node);
        final double mrh = -getR() * height;
        final double z = Math.log(1 - getA() * Math.exp(mrh));
        double l = -2 * z + mrh;

        if (tree.getRoot() == node) {
            l += mrh - z;
        }
        return l;
    }

    public boolean includeExternalNodesInLikelihoodCalculation() {
        return false;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return BIRTH_DEATH_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            int units = XMLParser.Utils.getUnitsAttr(xo);

            Parameter birthParameter = (Parameter) xo.getSocketChild(BIRTHDIFF_RATE);
            Parameter deathParameter = (Parameter) xo.getSocketChild(RELATIVE_DEATH_RATE);

            Logger.getLogger("dr.evomodel").info("Using Gernhard 2008 birth-death model: Gernhard T (2008) J Theor Biol, In press");
            return new BirthDeathGernhard08Model(birthParameter, deathParameter, units);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Gernhard (2008) model of speciation (equation at bottom of page 19 of draft).";
        }

        public Class getReturnType() {
            return BirthDeathGernhard08Model.class;
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