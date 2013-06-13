package dr.evomodelxml.coalescent;

import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.EmergingEpidemicModel;
import dr.evomodel.coalescent.ExponentialGrowthModel;
import dr.evomodel.tree.TreeModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parser for the EmergingEpidemic exponential model.
 */
public class EmergingEpidemicModelParser extends AbstractXMLObjectParser {

    public static String EMERGING_EPIDEMIC_MODEL = "emergingEpidemic";

    // the parameters
    public static String GROWTH_RATE = "growthRate";
    public static String GENERATION_TIME = "generationTime";
    public static String GENERATION_DISTRIBUTION_SHAPE = "generationShape";
    public static String OFFSPRING_DISPERSION = "offspringDispersion";

    public static String TREE = "epidemicTree";


    public String getParserName() {
        return EMERGING_EPIDEMIC_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        Parameter growthRateParameter =  (Parameter)xo.getElementFirstChild(GROWTH_RATE);
        Parameter generationTimeParameter = (Parameter)xo.getElementFirstChild(GENERATION_TIME);
        Parameter generationShapeParameter = (Parameter)xo.getElementFirstChild(GENERATION_DISTRIBUTION_SHAPE);
        Parameter offspringDispersionParameter = (Parameter)xo.getElementFirstChild(OFFSPRING_DISPERSION);
        TreeModel tree = (TreeModel)xo.getElementFirstChild(TREE);

        return new EmergingEpidemicModel(growthRateParameter, generationTimeParameter, generationShapeParameter, offspringDispersionParameter, tree, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A self-consistent model of emerging epidemics.";
    }

    public Class getReturnType() {
        return EmergingEpidemicModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(GROWTH_RATE, Parameter.class),
            new ElementRule(GENERATION_TIME, Parameter.class),
            new ElementRule(GENERATION_DISTRIBUTION_SHAPE, Parameter.class),
            new ElementRule(OFFSPRING_DISPERSION, Parameter.class),
            new ElementRule(TREE, Tree.class),

            XMLUnits.SYNTAX_RULES[0]
    };


}
