package dr.evomodel.epidemiology;

import dr.inference.model.Parameter;
import java.util.List;

public class SIRCompartmentalModel extends CompartmentalModel {

    protected Parameter originTimeNumS;

    public SIRCompartmentalModel(
            List<Parameter> rateParams,
            List<Parameter> compartmentCounts,
            Parameter origin,
            Parameter originTimeNumS,
            int numReactionChannels,
            int numGridPoints,
            double cutOff) {
        super("SIRCompartmentalModel");

        this.rateParameters = rateParams;
        for(int i = 0; i < rateParameters.size(); i++) {
            addVariable(rateParameters.get(i));
        }

        this.compartmentCounts = compartmentCounts;
        for(int i = 0; i < compartmentCounts.size(); i++) {
            addVariable(compartmentCounts.get(i));
        }
        this.originOne = origin;
        addVariable(originOne);
        this.numGridPoints = numGridPoints;
        this.cutOff = cutOff;
        this.numReactionChannels = numReactionChannels;
        this.numSpecies = compartmentCounts.size();
        this.vMatrix = getVMatrix();
        this.originTimeNumS = originTimeNumS;
        addVariable(originTimeNumS);
    }

    protected void setOriginTimeCompartmentCounts(int index){
        // initial S value
        compartmentCounts.get(0).setParameterValue(index, originTimeNumS.getParameterValue(0));
        // initial I value
        compartmentCounts.get(1).setParameterValue(index, 1);
        // initial R value
        compartmentCounts.get(2).setParameterValue(index, 0);
    }

    protected void setDefaultCompartmentCounts(int index){
        // default S value
        compartmentCounts.get(0).setParameterValue(index, originTimeNumS.getParameterValue(0)+1);
        // default I value
        compartmentCounts.get(1).setParameterValue(index, 0);
        // default R value
        compartmentCounts.get(2).setParameterValue(index, 0);
    }

    protected int[] getHighestOrdersOfReactions(){
        int[] gVec = new int[numSpecies];
        // g_i's determined by using the highest order of reaction of species i
        // g_S = 2;
        gVec[0] = 2;
        // g_I = 2;
        gVec[1] = 2;
        // g_R = 1;
        gVec[2] = 1;
        return gVec;
    }

    // v matrix describes how count vector changes with reaction
    // row corresponds to species/particle type, column corresponds to reaction channel
    protected int[][] getVMatrix(){
        int[][] v = new int[][]{
                {-1, 0, 1},  // S
                {1, -1, 0},  // I
                {0, 1, -1}   // R
        };
        return v;
    }

    protected double[] getReactionIntensities(double[] currentCounts){
        double[] rVec = new double[numReactionChannels];
        // currentCounts[0] has number of S for current time step
        // currentCounts[1] has number of I for current time step
        // currentCounts[2] has number of R for current time step
        double transmissionRate = rateParameters.get(0).getParameterValue(0);
        double recoveryRate = rateParameters.get(1).getParameterValue(0);
        double resusRate = rateParameters.get(3).getParameterValue(0);
        // infection
        rVec[0] = transmissionRate*currentCounts[0]*currentCounts[1];
        // recovery
        rVec[1] = recoveryRate*currentCounts[1];
        // re-susceptibility
        rVec[2] = resusRate*currentCounts[2];
        return rVec;
    }

    // Moved to TauLeapingSimulator.java
    /*
    protected double[] getMaxFiringTimes(double[] currentCounts, double[] r){
        double[] returnVal = new double[numReactionChannels];
        for (int c = 0; c < numReactionChannels; c++) {
            returnVal[c] = 0;
            for(int i = 0; i < numSpecies; i++) {
                if(vMatrix[i][c] < 0) {
                    double candidate = currentCounts[i] / Math.abs(vMatrix[i][c]);
                    if (r[i] > 0) {
                        if (returnVal[c] == 0 || returnVal[c] > candidate) {
                            returnVal[c] = candidate;
                        }
                    }
                }
            }
        }
        return returnVal;
    }
    */

    // compute derivative of r_j(x) with respect to time for all reactions j
    protected double[] getTimeDerivatives(double[] currentCounts){
        double[] returnVec = new double[numReactionChannels];
        double S = currentCounts[0];
        double I = currentCounts[1];
        double R = currentCounts[2];
        double transRate = rateParameters.get(0).getParameterValue(0);
        double recovRate = rateParameters.get(1).getParameterValue(0);
        double resusRate = rateParameters.get(3).getParameterValue(0);
        // for infection
        returnVec[0] = transRate*I*(-transRate*S*I + resusRate*R) + transRate*S*(transRate*S*I - recovRate*I);
        // for recovery
        returnVec[1] = recovRate*(transRate*S*I - recovRate*I);
        // for resusceptibility
        returnVec[2] = resusRate*(recovRate*I - resusRate*R);
        return returnVec;
    }

    // moved to TauLeapingSimulator.java
    /*
    protected double[] getTauLeapingPoissonIntensities(double[] currentCounts, double[] reactionInt, double tau){
        double[] returnVal = new double[numReactionChannels];
        // for standard tau leaping
        for(int r = 0; r < numReactionChannels; r++) {
            returnVal[r] = reactionInt[r]*tau;
        }
        return returnVal;
    }
    */

    protected double[] getSALPoissonIntensities(double[] currentCounts, double[] reactionInt, double tau){
        double[] returnVal = new double[numReactionChannels];
        // for standard tau leaping
        for(int r = 0; r < numReactionChannels; r++) {
            returnVal[r] = reactionInt[r]*tau;
        }
        // for SAL algorithm, need to add extra terms to account for linear change in intensity
        double[] timeDerivatives = getTimeDerivatives(currentCounts);
        for(int r = 0; r < numReactionChannels; r++) {
            returnVal[r] = returnVal[r] + timeDerivatives[r]*tau*tau*0.5;
        }
        return returnVal;
    }

    protected double[] getUpdatedCompartmentCounts(double[] currentCounts, double[] countsNew){
        double[] updatedCounts = new double[numSpecies];
        // S(t+tau) = S(t) + newS - newI
        updatedCounts[0] = currentCounts[0] + countsNew[0] - countsNew[1];
        // I(t+tau) = I(t) + newI - newR
        updatedCounts[1] = currentCounts[1] + countsNew[1] - countsNew[2];
        // R(t+tau) = R(t) + newR - newS
        updatedCounts[2] = currentCounts[2] + countsNew[2] - countsNew[0];
        return updatedCounts;
    }

    /*
    protected boolean hasMinimalCounts(int[] counts) {
        // check for negative counts
        for (int i = 0; i < counts.length; i++) {
            if (counts[i] < 0) {
                return false;
            }
        }
        // check if less than 1 susceptible
        if(counts[1] < 1){
            return false;
        }
        return true;
    }
    */
}