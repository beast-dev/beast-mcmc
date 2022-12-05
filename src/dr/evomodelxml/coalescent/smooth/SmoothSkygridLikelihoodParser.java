package dr.evomodelxml.coalescent.smooth;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.TreeIntervals;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.smooth.OldSmoothSkygridLikelihood;
import dr.evomodel.coalescent.smooth.SmoothSkygridLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class SmoothSkygridLikelihoodParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "smoothSkygridLikelihood";
    private static final String POPULATION_PARAMETER = GMRFSkyrideLikelihoodParser.POPULATION_PARAMETER;
    private static final String PRECISION_PARAMETER = GMRFSkyrideLikelihoodParser.PRECISION_PARAMETER;
    private static final String INTERVALS = GMRFSkyrideLikelihoodParser.INTERVALS;
    private static final String POPULATION_TREE = GMRFSkyrideLikelihoodParser.POPULATION_TREE;
    private static final String GRID_POINTS = GMRFSkyrideLikelihoodParser.GRID_POINTS;

    private static final String NUM_GRID_POINTS = GMRFSkyrideLikelihoodParser.NUM_GRID_POINTS;

    private static final String CUT_OFF = GMRFSkyrideLikelihoodParser.CUT_OFF;
    private static final String OLD = "old";
    private static final String SMOOTH_RATE = "smoothRate";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter logPopSizes = (Parameter) xo.getElementFirstChild(POPULATION_PARAMETER);
        Parameter gridPoints = getGridPoints(xo);

        boolean isOld = xo.getAttribute(OLD, false);

        if (isOld) {
            List<IntervalList> intervalList = new ArrayList<>();
            List<TreeIntervals> debugIntervalList = new ArrayList<>();
            if (xo.hasChildNamed(INTERVALS)) {
                XMLObject cxo = xo.getChild(INTERVALS);
                for (int i = 0; i < cxo.getChildCount(); ++i) {
                    intervalList.add((IntervalList) cxo.getChild(i));
                }
            } else {
                XMLObject cxo = xo.getChild(POPULATION_TREE);
                for (int i = 0; i < cxo.getChildCount(); ++i) {
                    TreeModel tree = (TreeModel) cxo.getChild(i);
                    intervalList.add(new BigFastTreeIntervals(tree));
                    debugIntervalList.add(new TreeIntervals(tree));
                }
            }

            if (!OldSmoothSkygridLikelihood.checkValidParameters(logPopSizes, gridPoints)) {
                throw new XMLParseException("Invalid initial parameters");
            }

            OldSmoothSkygridLikelihood likelihood = new OldSmoothSkygridLikelihood(xo.getId(), intervalList, logPopSizes, gridPoints);
            likelihood.setDebugIntervalList(debugIntervalList);
            return likelihood;
        } else {
            List<TreeModel> trees = new ArrayList<>();
            XMLObject cxo = xo.getChild(POPULATION_TREE);
            for (int i = 0; i < cxo.getChildCount(); i++) {
                trees.add((TreeModel) cxo.getChild(i));
            }
            Parameter smoothRate = (Parameter) xo.getElementFirstChild(SMOOTH_RATE);
            SmoothSkygridLikelihood likelihood = new SmoothSkygridLikelihood(xo.getId(), trees, logPopSizes, gridPoints, smoothRate);
            return likelihood;
        }
    }

    public static Parameter getGridPoints(XMLObject xo) throws XMLParseException {
        Parameter gridPoints;
        if (xo.getChild(GRID_POINTS) != null) {
            gridPoints = (Parameter) xo.getElementFirstChild(GRID_POINTS);
        } else {
            int numGridPoints = (int) ((Parameter) xo.getElementFirstChild(NUM_GRID_POINTS)).getParameterValue(0);
            double cutOff = ((Parameter) xo.getElementFirstChild(CUT_OFF)).getParameterValue(0);
            double[] gridLocations = new double[numGridPoints];

            for (int pt = 0; pt < gridLocations.length; pt++) {
                gridLocations[pt] = (pt + 1) * (cutOff / numGridPoints);
            }
            gridPoints = new Parameter.Default(gridLocations);
        }
        return gridPoints;
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
        return OldSmoothSkygridLikelihood.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(POPULATION_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(PRECISION_PARAMETER, new XMLSyntaxRule[]{
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

            new XORRule(
                    new ElementRule(INTERVALS, new XMLSyntaxRule[]{
                            new ElementRule(IntervalList.class, 1, Integer.MAX_VALUE)
                    }),

                    new AndRule(
                            new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
                                    new ElementRule(TreeModel.class, 1, Integer.MAX_VALUE)
                            }),
                            new ElementRule(SMOOTH_RATE, new XMLSyntaxRule[]{
                                    new ElementRule(Parameter.class, "Smooth rate for sigmoid functions")
                            },true)
                    )
            )
    };
}
