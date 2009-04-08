package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */

public class SubstitutionEpochModel extends AbstractSubstitutionModel {

    public static final String SUBSTITUTION_EPOCH_MODEL = "substitutionEpochModel";
    public static final String MODELS = "models";
    public static final String TRANSITION_TIMES = "transitionTimes";

    public static final boolean DEBUG = false;

    SubstitutionEpochModel(String name,
                           List<SubstitutionModel> modelList,
                           Parameter transitionTimes,
                           DataType dataType, FrequencyModel freqModel) {

        super(name, dataType, freqModel);

        this.modelList = modelList;
        this.transitionTimesParameter = transitionTimes;
        this.transitionTimes = transitionTimesParameter.getParameterValues();

        addParameter(transitionTimes);

        for (SubstitutionModel model : modelList)
            addModel(model);

        numberModels = modelList.size();
        weight = new double[numberModels];
        stateCount = dataType.getStateCount();
        stepMatrix = new double[stateCount * stateCount];
        productMatrix = new double[stateCount * stateCount];
        resultMatrix = new double[stateCount * stateCount];
    }


    protected void frequenciesChanged() {
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == freqModel)
            frequenciesChanged();
        else // This is an epoch model and I need to pass the info on
            fireModelChanged(object, index);
    }

    protected void ratesChanged() {
    }

    protected void setupRelativeRates() {
    }

    public void getTransitionProbabilities(double startTime, double endTime, double distance, double[] matrix) {
        int matrixCount = 0;
        boolean oneMatrix = (getEpochWeights(startTime, endTime, weight) == 1);
        for (int m = 0; m < numberModels; m++) {
            if (weight[m] > 0) {
                SubstitutionModel model = modelList.get(m);
                if (matrixCount == 0) {
                    if (oneMatrix) {
                        model.getTransitionProbabilities(distance, matrix);
                        break;
                    } else
                        model.getTransitionProbabilities(distance * weight[m], resultMatrix);
                    matrixCount++;
                } else {
                    model.getTransitionProbabilities(distance * weight[m], stepMatrix);
                    // Sum over unobserved state
                    int index = 0;
                    for (int i = 0; i < stateCount; i++) {
                        for (int j = 0; j < stateCount; j++) {
                            productMatrix[index] = 0;
                            for (int k = 0; k < stateCount; k++) {
                                productMatrix[index] += resultMatrix[i * stateCount + k] * stepMatrix[k * stateCount + j];
                            }
                            index++;
                        }
                    }
                    // Swap pointers
                    double[] tmpMatrix = resultMatrix;
                    resultMatrix = productMatrix;
                    productMatrix = tmpMatrix;
                }
            }
        }
        if (!oneMatrix)
            System.arraycopy(resultMatrix, 0, matrix, 0, stateCount * stateCount);
    }

    private int getEpochWeights(double startTime, double endTime, double[] weights) {

        int matrixCount = 0;
        final double lengthTime = endTime - startTime;
        final int lastTime = numberModels - 2;

        // model 0, 1, 2, ..., K-2, K-1
        // times   0, 1,  ...,   K-2,
        // where K = numberModels

        // First epoch: 0 -> transitionTimes[0];
        if (startTime <= transitionTimes[0]) {
            if (endTime <= transitionTimes[0])
                weights[0] = 1;
            else
                weights[0] = (transitionTimes[0] - startTime) / lengthTime;
            matrixCount++;
        } else
            weights[0] = 0;

        // Middle epoches:
        for (int i = 1; i <= lastTime; i++) {
            if (startTime <= transitionTimes[i]) {
                double start = Math.max(startTime, transitionTimes[i - 1]);
                double end = Math.min(endTime, transitionTimes[i]);
                weights[i] = (end - start) / lengthTime;
                matrixCount++;
            } else
                weights[i] = 0;
        }

        // Last epoch: transitionTimes[K-2] -> Infinity
        if (lastTime >= 0) {
            if (endTime > transitionTimes[lastTime]) {
                double start = Math.max(startTime, transitionTimes[lastTime]);
                weights[lastTime + 1] = (endTime - start) / lengthTime;
                matrixCount++;
            } else
                weights[lastTime + 1] = 0;
        }

        if (DEBUG) {
            double totalWeight = 0;
            for (int i = 0; i < numberModels; i++)
                totalWeight += weights[i];
            System.err.println("Start: " + startTime + " End: " + endTime + " Count: " + matrixCount + " Weight: " + totalWeight + " - " + new Vector(weights));
            if (totalWeight > 1.001) System.exit(-1);
            if (totalWeight < 0.999) System.exit(-1);
        }

        return matrixCount;
    }

    public void getTransitionProbabilities(double distance, double[] matrix) {
        throw new RuntimeException("Should not get here in a substitution epoch model.");
    }

    protected final void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
        super.handleParameterChangedEvent(parameter, index, type);

        if (parameter == transitionTimesParameter) {
            transitionTimes = transitionTimesParameter.getParameterValues();
            fireModelChanged(parameter, index);
        }
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            DataType dataType = null;
            FrequencyModel freqModel = null;
            List<SubstitutionModel> modelList = new ArrayList<SubstitutionModel>();
            XMLObject cxo = (XMLObject) xo.getChild(MODELS);
            for (int i = 0; i < cxo.getChildCount(); i++) {
                SubstitutionModel model = (SubstitutionModel) cxo.getChild(i);

                if (dataType == null) {
                    dataType = model.getDataType();
                } else if (dataType != model.getDataType())
                    throw new XMLParseException("Substitution models across epoches must use the same data type.");

                if (freqModel == null) {
                    freqModel = model.getFrequencyModel();
                } else if (freqModel != model.getFrequencyModel())
                    throw new XMLParseException("Substitution models across epoches must currently use the same frequency model.\nHarass Marc to fix this.");

                modelList.add(model);
            }

            Parameter transitionTimes = (Parameter) xo.getChild(Parameter.class);

            if (transitionTimes.getDimension() != modelList.size() - 1) {
                throw new XMLParseException("# of transition times must equal # of substitution models - 1\n" + transitionTimes.getDimension() + "\n" + modelList.size());
            }

            return new SubstitutionEpochModel(SUBSTITUTION_EPOCH_MODEL,
                    modelList, transitionTimes,
                    dataType, freqModel);
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(MODELS,
                            new XMLSyntaxRule[]{
                                    new ElementRule(AbstractSubstitutionModel.class, 1, Integer.MAX_VALUE),
                            }),
                    new ElementRule(Parameter.class),
            };
        }

        public String getParserDescription() {
            return null;
        }

        public Class getReturnType() {
            return SubstitutionEpochModel.class;
        }

        public String getParserName() {
            return SUBSTITUTION_EPOCH_MODEL;
        }
    };

    private List<SubstitutionModel> modelList;
    private Parameter transitionTimesParameter;
    private double[] transitionTimes;
    private double[] weight;
    private double[] stepMatrix;
    private double[] productMatrix;
    private double[] resultMatrix;
    private int numberModels;
    private int stateCount;

}
