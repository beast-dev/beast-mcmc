package dr.evomodel.continuous;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeAttributeProvider;
import dr.inference.model.*;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * @author Marc Suchard
 */


public class MultivariateDiffusionModel extends AbstractModel implements TreeAttributeProvider {

    public static final String DIFFUSION_PROCESS = "multivariateDiffusionModel";
    public static final String DIFFUSION_CONSTANT = "precisionMatrix";
    public static final String PRECISION_TREE_ATTRIBUTE = "precision";

    public static final double LOG2PI = Math.log(2*Math.PI);

    /**
     * Construct a diffusion model.
     */

    public MultivariateDiffusionModel(MatrixParameter diffusionPrecisionMatrixParameter) {

        super(DIFFUSION_PROCESS);

        this.diffusionPrecisionMatrixParameter = diffusionPrecisionMatrixParameter;
        calculatePrecisionInfo();
        addVariable(diffusionPrecisionMatrixParameter);

    }

    public MultivariateDiffusionModel() {
        super(DIFFUSION_PROCESS);
    }


//    public void randomize(Parameter trait) {
//    }

    public void check(Parameter trait) throws XMLParseException {
        assert trait != null;
    }

    public Parameter getPrecisionParameter() {checkVariableChanged();

        return diffusionPrecisionMatrixParameter;
    }

    public double[][] getPrecisionmatrix() {
        if (diffusionPrecisionMatrixParameter != null) {
            checkVariableChanged();
            return diffusionPrecisionMatrixParameter.getParameterAsMatrix();
        }
        return null;
    }

    public double getDeterminantPrecisionMatrix() {  checkVariableChanged();
        return determinatePrecisionMatrix; }

    /**
     * @return the log likelihood of going from start to stop in the given time
     */
    public double getLogLikelihood(double[] start, double[] stop, double time) {

        if (time == 0) {
            boolean equal = true;
            for(int i=0; i<start.length; i++) {
                if( start[i] != stop[i] ) {
                    equal = false;
                    break;
                }
            }
            if (equal)
                return 0.0;
            return Double.NEGATIVE_INFINITY;
        }

        return calculateLogDensity(start, stop, time);
    }

    protected void checkVariableChanged(){
        if(variableChanged){
            calculatePrecisionInfo();
            variableChanged=false;
        }
    }

    protected double calculateLogDensity(double[] start, double[] stop, double time) {
        checkVariableChanged();
        final double logDet = Math.log(determinatePrecisionMatrix);
        return MultivariateNormalDistribution.logPdf(stop, start, diffusionPrecisionMatrix, logDet, time);
    }

    // todo should be a test, no?
    public static void main(String[] args) {
        double[] start = {1, 2};
        double[] stop = {0, 0};
        double[][] precision = {{2, 0.5}, {0.5, 1}};
        double scale = 0.2;
        MatrixParameter precMatrix = new MatrixParameter("Hello");
        precMatrix.addParameter(new Parameter.Default(precision[0]));
        precMatrix.addParameter(new Parameter.Default(precision[1]));
        MultivariateDiffusionModel model = new MultivariateDiffusionModel(precMatrix);
        System.err.println("logPDF = " + model.calculateLogDensity(start, stop, scale));
        System.err.println("Should be -19.948");
    }

    protected void calculatePrecisionInfo() {
        diffusionPrecisionMatrix = diffusionPrecisionMatrixParameter.getParameterAsMatrix();
        determinatePrecisionMatrix =
                MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(
                        diffusionPrecisionMatrix);
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************
    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        variableChanged=true;
//        calculatePrecisionInfo();
    }

    protected void storeState() {
        savedDeterminatePrecisionMatrix = determinatePrecisionMatrix;
        savedDiffusionPrecisionMatrix = diffusionPrecisionMatrix;
        storedVariableChanged=variableChanged;
    }

    protected void restoreState() {
        determinatePrecisionMatrix = savedDeterminatePrecisionMatrix;
        diffusionPrecisionMatrix = savedDiffusionPrecisionMatrix;
        variableChanged=storedVariableChanged;
    }

    protected void acceptState() {
    } // no additional state needs accepting

    public String[] getTreeAttributeLabel() {
        return new String[] {PRECISION_TREE_ATTRIBUTE};
    }

    public String[] getAttributeForTree(Tree tree) {
        if (diffusionPrecisionMatrixParameter != null) {
            return new String[] {diffusionPrecisionMatrixParameter.toSymmetricString()};
        }
        return new String[] { "null" };
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document document) {
        throw new RuntimeException("Not implemented!");
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DIFFUSION_PROCESS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = xo.getChild(DIFFUSION_CONSTANT);
            MatrixParameter diffusionParam = (MatrixParameter) cxo.getChild(MatrixParameter.class);

            return new MultivariateDiffusionModel(diffusionParam);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Describes a multivariate normal diffusion process.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(DIFFUSION_CONSTANT,
                        new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
        };

        public Class getReturnType() {
            return MultivariateDiffusionModel.class;
        }
    };

    // **************************************************************
    // Private instance variables
    // **************************************************************

    protected MatrixParameter diffusionPrecisionMatrixParameter;
    private double determinatePrecisionMatrix;
    private double savedDeterminatePrecisionMatrix;
    private double[][] diffusionPrecisionMatrix;
    private double[][] savedDiffusionPrecisionMatrix;

    private boolean variableChanged=true;
    private boolean storedVariableChanged;

}

