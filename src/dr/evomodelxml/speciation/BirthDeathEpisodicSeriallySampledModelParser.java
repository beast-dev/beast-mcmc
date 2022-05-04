package dr.evomodelxml.speciation;

import dr.evolution.util.Units;
import dr.evomodel.speciation.BirthDeathEpisodicSeriallySampledModel;
import dr.evomodel.speciation.BirthDeathSerialSamplingModel;
import dr.evomodel.speciation.BirthDeathSerialSkylineModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

public class BirthDeathEpisodicSeriallySampledModelParser {

    public static final String ESSBDM = "episodicSeriallySampledBirthDeathModel";
    public static final String NUM_GRID_POINTS = "numGridPoints";
    public static final String CUT_OFF = "cutOff";
    public static final String LAMBDA = "birthRate";
    public static final String MU = "deathRate";
    public static final String PSI = "samplingRate";
    public static final String RHO = "samplingProbability";
    public static final String R = "treatmentProbability";
    public static final String ORIGIN = "origin";
    public static final String TREE_TYPE = "type";
    public static final String CONDITION = "conditionOnSurvival";

    public String getParserName() {
        return ESSBDM;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String modelName = xo.getId();
        final Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        final Parameter cutoff = (Parameter) xo.getElementFirstChild(CUT_OFF);
        final Parameter numGridPoints = (Parameter) xo.getElementFirstChild(NUM_GRID_POINTS);

        final Parameter lambda = (Parameter) xo.getElementFirstChild(LAMBDA);
        final Parameter mu = (Parameter) xo.getElementFirstChild(MU);
        final Parameter psi = (Parameter) xo.getElementFirstChild(PSI);

        final Parameter rho    = xo.hasChildNamed(RHO) ? (Parameter) xo.getElementFirstChild(RHO) : new Parameter.Default(0.0);
        final Parameter r      = xo.hasChildNamed(R) ? (Parameter) xo.getElementFirstChild(R) : new Parameter.Default(1.0);

        final Parameter origin = (Parameter) xo.getElementFirstChild(ORIGIN);;

        return new BirthDeathEpisodicSeriallySampledModel(modelName, lambda, mu, psi, r, rho, origin, (int)(numGridPoints.getParameterValue(0)), cutoff.getParameterValue(0), units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public static String getCitation() {
        // Eventually we can cite Yucai's work too
        return "Gavryushkina et al, Bayesian inference of sampled ancestor trees for epidemiology and fossil calibration, " +
                "PLoS Comp. Biol., doi: 10.1371/journal.pcbi.1003919, 2014";
    }

    public String getParserDescription() {
        return "Gavryushkina et al, Bayesian inference of sampled ancestor trees for epidemiology and fossil calibration, " +
                "PLoS Comp. Biol., doi: 10.1371/journal.pcbi.1003919, 2014";
    }

    public Class getReturnType() {
        return BirthDeathSerialSamplingModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TREE_TYPE, true),
            new ElementRule(ORIGIN, Parameter.class, "The origin of the infection, x0 > tree.rootHeight", false),
            new ElementRule(LAMBDA, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(MU, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(PSI, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(R, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(RHO, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(CONDITION, new XMLSyntaxRule[]{new ElementRule(boolean.class)}, true),
            XMLUnits.SYNTAX_RULES[0]
    };
}
