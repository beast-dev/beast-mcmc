package dr.evomodel.coalescent;

import dr.evolution.coalescent.MultiLociTreeSet;
import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Joseph Heled
 * @version $Id$
 */
public class VariableDemographicModel extends DemographicModel implements MultiLociTreeSet {
    static final String MODEL_NAME = "variableDemographic";
    public static final String POPULATION_SIZES = "populationSizes";
    public static final String INDICATOR_PARAMETER = "indicators";
    public static final String POPULATION_TREES = "trees";
    private static final String PLOIDY = "ploidy";
    public static String POP_TREE = "ptree";

    public static final String LOG_SPACE = "logUnits";
    public static final String USE_MIDPOINTS = "useMidpoints";

    public static final String TYPE = "type";
    public static final String STEPWISE = "stepwise";
    public static final String LINEAR = "linear";
    public static final String EXPONENTIAL = "exponential";

    public static final String demoElementName = "demographic";

    private final Parameter popSizeParameter;
    private final Parameter indicatorParameter;

    public Type getType() {
        return type;
    }

    private final Type type;
    private final boolean logSpace;
    private final boolean mid;
    private final TreeModel[] trees;
    private VDdemographicFunction demoFunction = null;
    private VDdemographicFunction savedDemoFunction = null;
    private final double[] populationFactors;

    public Parameter getPopulationValues() {
        return popSizeParameter;
    }

    public enum Type {
        STEPWISE,
        LINEAR,
        EXPONENTIAL
    }

    public VariableDemographicModel(TreeModel[] trees, double[] popFactors,
                                    Parameter popSizeParameter, Parameter indicatorParameter,
                                    Type type, boolean logSpace, boolean mid) {
        super(MODEL_NAME);

        this.popSizeParameter = popSizeParameter;
        this.indicatorParameter = indicatorParameter;

        this.populationFactors = popFactors;

        int events = 0;
        for (Tree t : trees) {
            // number of coalescent envents
            events += t.getExternalNodeCount() - 1;
            // we will have to handle this I guess
            assert t.getUnits() == trees[0].getUnits();
        }
        // all trees share time 0, need fixing for serial data

        events += type == Type.STEPWISE ? 0 : 1;

        final int popSizes = popSizeParameter.getDimension();
        final int nIndicators = indicatorParameter.getDimension();
        this.type = type;
        this.logSpace = logSpace;
        this.mid = mid;

        if (popSizes != events) {
            throw new IllegalArgumentException("Dimension of population parameter (" + popSizes +
                    ") must be the same as the number of internal nodes in the tree. (" + events + ")");
        }

        if (nIndicators != popSizes - 1) {
            throw new IllegalArgumentException("Dimension of indicator parameter must one less than the number of internal nodes in the tree. ("
                    + nIndicators + " != " + (events - 1) + ")");
        }

        this.trees = trees;

        for (TreeModel t : trees) {
            addModel(t);
        }

        addParameter(indicatorParameter);
        addParameter(popSizeParameter);
    }

    public int nLoci() {
        return trees.length;
    }

    public Tree getTree(int k) {
        return trees[k];
    }

    public TreeIntervals getTreeIntervals(int nt) {
        return getDemographicFunction().getTreeIntervals(nt);
    }

    public double getPopulationFactor(int nt) {
        return populationFactors[nt];
    }

    public void storeTheState() {
        // as a demographic model store/restore is already taken care of 
    }

    public void restoreTheState() {
        // as a demographic model store/restore is already taken care of
    }

    public VDdemographicFunction getDemographicFunction() {
        if (demoFunction == null) {
            demoFunction = new VDdemographicFunction(trees, type,
                    indicatorParameter.getParameterValues(), popSizeParameter.getParameterValues(), logSpace, mid);
        } else {
            demoFunction.setup(trees, indicatorParameter.getParameterValues(), popSizeParameter.getParameterValues(),
                    logSpace, mid);
        }
        return demoFunction;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // tree has changed
        //System.out.println("model changed: " + model);
        if (demoFunction != null) {
            if (demoFunction == savedDemoFunction) {
                demoFunction = new VDdemographicFunction(demoFunction);
            }
            for (int k = 0; k < trees.length; ++k) {
                if (model == trees[k]) {
                    demoFunction.treeChanged(k);
                    //System.out.println("tree changed: " + k + " " + Arrays.toString(demoFunction.dirtyTrees)
                    //       + " " + demoFunction.dirtyTrees);
                    break;
                }
                assert k + 1 < trees.length;
            }
        }
        super.handleModelChangedEvent(model, object, index);
        fireModelChanged(this);
    }

    protected final void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
        //System.out.println("parm changed: " + parameter);
        super.handleParameterChangedEvent(parameter, index, type);
        if (demoFunction != null) {
            if (demoFunction == savedDemoFunction) {
                demoFunction = new VDdemographicFunction(demoFunction);
            }
            demoFunction.setDirty();
        }
        fireModelChanged(this);
    }

    protected void storeState() {
        savedDemoFunction = demoFunction;
    }

    protected void restoreState() {
        //System.out.println("restore");
        demoFunction = savedDemoFunction;
        savedDemoFunction = null;
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return VariableDemographicModel.MODEL_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = (XMLObject) xo.getChild(VariableSkylineLikelihood.POPULATION_SIZES);
            Parameter popParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = (XMLObject) xo.getChild(VariableSkylineLikelihood.INDICATOR_PARAMETER);
            Parameter indicatorParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = (XMLObject) xo.getChild(POPULATION_TREES);

            final int nc = cxo.getChildCount();
            TreeModel[] treeModels = new TreeModel[nc];
            double[] populationFactor = new double[nc];

            for (int k = 0; k < treeModels.length; ++k) {
                final XMLObject child = (XMLObject) cxo.getChild(k);
                populationFactor[k] = child.hasAttribute(PLOIDY) ? child.getDoubleAttribute(PLOIDY) : 1.0;

                treeModels[k] = (TreeModel) child.getChild(TreeModel.class);
            }

            Type type = Type.STEPWISE;

            if (xo.hasAttribute(TYPE)) {
                final String s = xo.getStringAttribute(TYPE);
                if (s.equalsIgnoreCase(STEPWISE)) {
                    type = Type.STEPWISE;
                } else if (s.equalsIgnoreCase(LINEAR)) {
                    type = Type.LINEAR;
                } else if (s.equalsIgnoreCase(EXPONENTIAL)) {
                    type = Type.EXPONENTIAL;
                } else {
                    throw new XMLParseException("Unknown Bayesian Skyline type: " + s);
                }
            }

            final boolean logSpace = xo.getAttribute(LOG_SPACE, false);
            final boolean useMid = xo.getAttribute(USE_MIDPOINTS, false);

            Logger.getLogger("dr.evomodel").info("Variable demographic: " + type.toString() + " control points");

            return new VariableDemographicModel(treeModels, populationFactor, popParam, indicatorParam, type,
                    logSpace, useMid);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of the tree given the population size vector.";
        }

        public Class getReturnType() {
            return DemographicModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(VariableSkylineLikelihood.TYPE, true),
                AttributeRule.newBooleanRule(LOG_SPACE, true),
                AttributeRule.newBooleanRule(USE_MIDPOINTS, true),

                new ElementRule(VariableSkylineLikelihood.POPULATION_SIZES, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(VariableSkylineLikelihood.INDICATOR_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(POPULATION_TREES, new XMLSyntaxRule[]{
                        new ElementRule(POP_TREE, new XMLSyntaxRule[]{
                                AttributeRule.newDoubleRule(PLOIDY, true),
                                new ElementRule(TreeModel.class),
                        }, 1, Integer.MAX_VALUE)
                })
        };
    };
}
