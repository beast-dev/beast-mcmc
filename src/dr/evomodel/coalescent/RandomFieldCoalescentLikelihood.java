package dr.evomodel.coalescent;

import dr.inference.distribution.RandomField;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.Reportable;

/**
 * @author Filippo Monti
 * @author Xiang Ji
 * @author Marc A. Suchard
 */

public class RandomFieldCoalescentLikelihood extends AbstractModel implements Reportable {

    private final MultilocusNonparametricCoalescentLikelihood coalescentLikelihood;
    private final RandomField fieldLikelihood;

    private boolean coalescentLikelihoodKnown;
    private boolean fieldLikelihoodKnown;
    private boolean logLikelihoodKnown;

    private double logFieldLikelihood;
    private double logCoalescentLikelihood;
    private double logLikelihood;

    public RandomFieldCoalescentLikelihood(MultilocusNonparametricCoalescentLikelihood coalescentLikelihood,
                                           RandomField fieldLikelihood) {
        super("RandomFieldCoalescentLikelihood");

        this.coalescentLikelihood = coalescentLikelihood;
        this.fieldLikelihood = fieldLikelihood;

        this.coalescentLikelihoodKnown = false;
        this.fieldLikelihoodKnown = false;
        this.logLikelihoodKnown = false;

        addModel(coalescentLikelihood);
        addModel(fieldLikelihood);
    }

    public double getLogLikelihood() {
        if (!coalescentLikelihoodKnown) {
            logCoalescentLikelihood = coalescentLikelihood.getLogLikelihood();
            coalescentLikelihoodKnown = true;
            logLikelihoodKnown = false;
        }
        if (!fieldLikelihoodKnown) {
            logFieldLikelihood = fieldLikelihood.getLogLikelihood();
            fieldLikelihoodKnown = true;
            logLikelihoodKnown = false;
        }
        if (!logLikelihoodKnown) {
            logLikelihood = logCoalescentLikelihood + logFieldLikelihood;
            logLikelihoodKnown = true;
        }
        return logLikelihood;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == coalescentLikelihood) {
            coalescentLikelihoodKnown = false;
            fireModelChanged();
        } else if (model == fieldLikelihood) {
            fieldLikelihoodKnown = false;
            fireModelChanged();
        } else {
            throw new RuntimeException("Unknown model changed event");
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }


    // TODO adapt this commeted out code? It should be already partially taken care of downstream
//    protected void storeState() {
//        super.storeState();
//        System.arraycopy(coalescentIntervals, 0, storedCoalescentIntervals, 0, coalescentIntervals.length);
//        System.arraycopy(sufficientStatistics, 0, storedSufficientStatistics, 0, sufficientStatistics.length);
//        storedWeightMatrix = weightMatrix.copy();
//        storedLogFieldLikelihood = logFieldLikelihood;
//    }
//
//    protected void restoreState() {
//        super.restoreState();
//        // TODO Just swap pointers
//        System.arraycopy(storedCoalescentIntervals, 0, coalescentIntervals, 0, storedCoalescentIntervals.length);
//        System.arraycopy(storedSufficientStatistics, 0, sufficientStatistics, 0, storedSufficientStatistics.length);
//        weightMatrix = storedWeightMatrix;
//        logFieldLikelihood = storedLogFieldLikelihood;
//    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    public String getReport() {
        return null;
    }
}



