package dr.evomodelxml.coalescent;

import dr.evolution.tree.Tree;
import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 *
 */
public class GMRFSkyrideLikelihoodParser extends AbstractXMLObjectParser {

	public static final String SKYLINE_LIKELIHOOD = "gmrfSkyrideLikelihood";
    public static final String SKYRIDE_LIKELIHOOD = "skyrideLikelihood";

	public static final String POPULATION_PARAMETER = "populationSizes";
	public static final String GROUP_SIZES = "groupSizes";
	public static final String PRECISION_PARAMETER = "precisionParameter";
	public static final String POPULATION_TREE = "populationTree";
	public static final String LAMBDA_PARAMETER = "lambdaParameter";
	public static final String BETA_PARAMETER = "betaParameter";
	public static final String COVARIATE_MATRIX = "covariateMatrix";
	public static final String RANDOMIZE_TREE = "randomizeTree";
	public static final String TIME_AWARE_SMOOTHING = "timeAwareSmoothing";

    public static final String RESCALE_BY_ROOT_ISSUE = "rescaleByRootHeight";
    public static final String GRID_POINTS = "gridPoints";
    public static final String OLD_SKYRIDE = "oldSkyride";
    public static final String NUM_GRID_POINTS = "numGridPoints";
    public static final String CUT_OFF = "cutOff";
    public static final String PHI_PARAMETER = "phiParameter";
    

    public String getParserName() {
        return SKYLINE_LIKELIHOOD;
    }

    public String[] getParserNames(){
        return new String[]{getParserName(), SKYRIDE_LIKELIHOOD}; // cannot duplicate 
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(POPULATION_PARAMETER);
        Parameter popParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(PRECISION_PARAMETER);
        Parameter precParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(POPULATION_TREE);

        List<Tree> treeList = new ArrayList<Tree>();
        for (int i = 0; i < cxo.getChildCount(); i++) {
            Object testObject = cxo.getChild(i);
            if (testObject instanceof Tree) {
                treeList.add((TreeModel) testObject);
            }
        }

//        TreeModel treeModel = (TreeModel) cxo.getChild(TreeModel.class);

        cxo = xo.getChild(GROUP_SIZES);
        Parameter groupParameter = null;
        if (cxo != null) {
            groupParameter = (Parameter) cxo.getChild(Parameter.class);
        
            if (popParameter.getDimension() != groupParameter.getDimension())
                throw new XMLParseException("Population and group size parameters must have the same length");
        }

        Parameter lambda;
        if (xo.getChild(LAMBDA_PARAMETER) != null) {
            cxo = xo.getChild(LAMBDA_PARAMETER);
            lambda = (Parameter) cxo.getChild(Parameter.class);
        } else {
            lambda = new Parameter.Default(1.0);
        }
        /*
        Parameter gridPoints = null;
        if (xo.getChild(GRID_POINTS) != null) {
            cxo = xo.getChild(GRID_POINTS);
            gridPoints = (Parameter) cxo.getChild(Parameter.class);
        }
        */
        Parameter numGridPoints = null;
        if (xo.getChild(NUM_GRID_POINTS) != null) {
            cxo = xo.getChild(NUM_GRID_POINTS);
            numGridPoints = (Parameter) cxo.getChild(Parameter.class);
        }

        Parameter cutOff = null;
        if (xo.getChild(CUT_OFF) != null) {
            cxo = xo.getChild(CUT_OFF);
            cutOff = (Parameter) cxo.getChild(Parameter.class);
        }

        Parameter phi = null;
        if (xo.getChild(PHI_PARAMETER) != null) {
            cxo = xo.getChild(PHI_PARAMETER);
            phi = (Parameter) cxo.getChild(Parameter.class);
        }


        Parameter beta = null;
        if (xo.getChild(BETA_PARAMETER) != null) {
            cxo = xo.getChild(BETA_PARAMETER);
            beta = (Parameter) cxo.getChild(Parameter.class);
        }

        MatrixParameter dMatrix = null;
        if (xo.getChild(COVARIATE_MATRIX) != null) {
            cxo = xo.getChild(COVARIATE_MATRIX);
            dMatrix = (MatrixParameter) cxo.getChild(MatrixParameter.class);
        }

        boolean timeAwareSmoothing = GMRFSkyrideLikelihood.TIME_AWARE_IS_ON_BY_DEFAULT;
        if (xo.hasAttribute(TIME_AWARE_SMOOTHING)) {
            timeAwareSmoothing = xo.getBooleanAttribute(TIME_AWARE_SMOOTHING);
        }

        if ((dMatrix != null && beta == null) || (dMatrix == null && beta != null))
            throw new XMLParseException("Must specify both a set of regression coefficients and a design matrix.");

        if (dMatrix != null) {
            if (dMatrix.getRowDimension() != popParameter.getDimension())
                throw new XMLParseException("Design matrix row dimension must equal the population parameter length.");
            if (dMatrix.getColumnDimension() != beta.getDimension())
                throw new XMLParseException("Design matrix column dimension must equal the regression coefficient length.");
        }

        if (xo.getAttribute(RANDOMIZE_TREE, false)) {
            for (Tree tree : treeList) {
                if (tree instanceof TreeModel) {
                    GMRFSkyrideLikelihood.checkTree((TreeModel) tree);
                } else {
                    throw new XMLParseException("Can not randomize a fixed tree");
                }
            }
        }

        boolean rescaleByRootHeight = xo.getAttribute(RESCALE_BY_ROOT_ISSUE, true);

        Logger.getLogger("dr.evomodel").info("The " + SKYLINE_LIKELIHOOD + " has " +
                (timeAwareSmoothing ? "time aware smoothing" : "uniform smoothing"));

        if (xo.getAttribute(OLD_SKYRIDE, true)) {

             return new GMRFSkyrideLikelihood(treeList, popParameter, groupParameter, precParameter,
                lambda, beta, dMatrix, timeAwareSmoothing, rescaleByRootHeight);

        } else {

             return new GMRFMultilocusSkyrideLikelihood(treeList, popParameter, groupParameter, precParameter,
                lambda, beta, dMatrix, timeAwareSmoothing,cutOff.getParameterValue(0),(int) numGridPoints.getParameterValue(0), phi);

        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of the tree given the population size vector.";
    }

    public Class getReturnType() {
        return GMRFSkyrideLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(POPULATION_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(PRECISION_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(PHI_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true), // Optional
            new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
                    new ElementRule(TreeModel.class, 1, Integer.MAX_VALUE)
            }),
            new ElementRule(GROUP_SIZES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            AttributeRule.newBooleanRule(RESCALE_BY_ROOT_ISSUE, true),
            AttributeRule.newBooleanRule(RANDOMIZE_TREE, true),
            AttributeRule.newBooleanRule(TIME_AWARE_SMOOTHING, true),
            AttributeRule.newBooleanRule(OLD_SKYRIDE, true)
    };

}
