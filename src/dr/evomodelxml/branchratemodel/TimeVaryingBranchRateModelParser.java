package dr.evomodelxml.branchratemodel;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.TimeVaryingBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.inference.model.Parameter;
import dr.xml.*;

import static dr.evomodelxml.coalescent.smooth.SmoothSkygridLikelihoodParser.getGridPoints;

public class TimeVaryingBranchRateModelParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "timeVaryingRates";

    private static final String RATES = StrictClockBranchRatesParser.RATE;
    private static final String GRID_POINTS = GMRFSkyrideLikelihoodParser.GRID_POINTS;
    private static final String NUM_GRID_POINTS = GMRFSkyrideLikelihoodParser.NUM_GRID_POINTS;
    private static final String CUT_OFF = GMRFSkyrideLikelihoodParser.CUT_OFF;

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree) xo.getChild(Tree.class);
        Parameter rates = (Parameter) xo.getElementFirstChild(RATES);
        Parameter gridPoints = getGridPoints(xo);

        if (rates.getDimension() != gridPoints.getDimension() + 1) {
            throw new XMLParseException("Rates dimension != gridPoints dimension + 1");
        }

        return new TimeVaryingBranchRateModel(tree, rates, gridPoints);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return TimeVaryingBranchRateModelParser.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Tree.class),
            new ElementRule(RATES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new XORRule(
                    new ElementRule(GRID_POINTS, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
                    new AndRule(
                            new ElementRule(CUT_OFF, new XMLSyntaxRule[]{
                                    new ElementRule(Parameter.class)
                            }),
                            new ElementRule(NUM_GRID_POINTS, new XMLSyntaxRule[]{
                                    new ElementRule(Parameter.class)
                            })
                    )
            ),
    };
}
