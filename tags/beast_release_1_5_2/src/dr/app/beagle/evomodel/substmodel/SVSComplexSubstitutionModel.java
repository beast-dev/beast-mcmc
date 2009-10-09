package dr.app.beagle.evomodel.substmodel;

import dr.inference.model.*;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.evolution.datatype.DataType;

/**
 * @author Marc Suchard
 */
public class SVSComplexSubstitutionModel extends ComplexSubstitutionModel implements Likelihood,
        BayesianStochasticSearchVariableSelection {
    
    public SVSComplexSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel,
                                          Parameter ratesParameter, Parameter indicatorsParameter) {
           super(name, dataType, freqModel, ratesParameter);

           this.indicatorsParameter  = indicatorsParameter;
           addVariable(indicatorsParameter);
       }

       @Override
       protected void setupRelativeRates(double[] rates) {
           for (int i = 0; i < rates.length; i++) {
               rates[i] = ratesParameter.getParameterValue(i) * indicatorsParameter.getParameterValue(i);
           }
       }

       public Parameter getIndicators() {
           return indicatorsParameter;
       }

       public boolean validState() {
           return !updateMatrix || connectedAndWellConditioned();
       }

       protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
           if (variable == ratesParameter && indicatorsParameter.getParameterValue(index) == 0)
               return; // Does not affect likelihood
           super.handleVariableChangedEvent(variable,index,type);
       }

       public Model getModel() {
           return this;
       }

       public double getLogLikelihood() {
           if (updateMatrix) {
               if (!connectedAndWellConditioned()) {
                   return Double.NEGATIVE_INFINITY;
               }
           }
           return 0;
       }

//       public double getLogLikelihood() {
//        double logL = super.getLogLikelihood();
//        if (logL == 0 && considerConnectedness && !isStronglyConnected()) { // Also check that graph is connected
//            logL = Double.NEGATIVE_INFINITY;
//        }
//        return logL;
//    }

       private boolean connectedAndWellConditioned() {
           if (probability == null)
               probability = new double[stateCount*stateCount];

           getTransitionProbabilities(1.0,probability);
           return Utils.connectedAndWellConditioned(probability);
       }

       /**
        * Forces a complete recalculation of the likelihood next time getLikelihood is called
        */
       public void makeDirty() {
          updateMatrix = true;
       }

       /**
        * @return A detailed name of likelihood for debugging.
        */
       public String prettyName() {
           return "SVSComplexSubstitutionModel-connectedness";
       }

         // **************************************************************
       // Loggable IMPLEMENTATION
       // **************************************************************

       public LogColumn[] getColumns() {
           return new LogColumn[]{
                   new LikelihoodColumn(getId())
           };
       }

       protected class LikelihoodColumn extends NumberColumn {
           public LikelihoodColumn(String label) {
               super(label);
           }

           public double getDoubleValue() {
               return getLogLikelihood();
           }
       }

       private double[] probability = null;

       private final Parameter indicatorsParameter;
}
