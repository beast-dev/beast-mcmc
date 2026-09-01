package dr.evomodel.epidemiology;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.VariableListener;

import java.util.ArrayList;
import java.util.List;

public class TwoPathogenModelBirthRates extends Parameter.Abstract implements VariableListener {

    private int numSkippedIntervals;
    private int numGridIntervals;
    private final int pathogenNumber;
    private final boolean seasonalModel;
    // interval width of every interval for pathogen with most recent sampling date
    // and width of every interval except for possibly the most recent interval for other pathogen
    private final double generalIntervalWidth;
    private final Parameter transmissionRate;
    private final List<Parameter> susCompartmentCounts;
    private final Parameter infectionRateModulationI;
    private final Parameter infectionRateModulationC;
    private final Parameter infectionRateModulationR;
    private final double[] seasonalMultipliers;
    private Bounds bounds = null;

    public TwoPathogenModelBirthRates(int pathogenNumber,
                                      boolean seasonalModel,
                                      Parameter transmissionRate,
                                      Parameter numSS,
                                      Parameter numSI,
                                      Parameter numSC,
                                      Parameter numSR,
                                      Parameter numIS,
                                      Parameter numCS,
                                      Parameter numRS,
                                      Parameter infecRateModI,
                                      Parameter infecRateModC,
                                      Parameter infecRateModR,
                                      Parameter originTimeNumSS,
                                      Parameter seasonalAmplitude,
                                      Parameter seasonalPeakDay,
                                      Parameter seasonalPeriod,
                                      Parameter cutOff,
                                      Parameter numGridPoints,
                                      Parameter mostRecentSamplingDateOne,
                                      Parameter mostRecentSamplingDateTwo
                                       ) {

        this.pathogenNumber = pathogenNumber;
        this.seasonalModel = seasonalModel;
        this.transmissionRate = transmissionRate;
        this.infectionRateModulationI = infecRateModI;
        this.infectionRateModulationC = infecRateModC;
        this.infectionRateModulationR = infecRateModR;

        int numGridPointsVal = (int) numGridPoints.getParameterValue(0);
        double cutOffVal = cutOff.getParameterValue(0);
        this.generalIntervalWidth = cutOffVal/numGridPointsVal;

        setGridDimensions(pathogenNumber, numGridPointsVal, mostRecentSamplingDateOne.getParameterValue(0),
                mostRecentSamplingDateTwo.getParameterValue(0));

        double moreRecentDate = Math.max(mostRecentSamplingDateOne.getParameterValue(0),
                mostRecentSamplingDateTwo.getParameterValue(0));

        this.seasonalMultipliers = getSeasonalMultipliers(moreRecentDate, originTimeNumSS,seasonalAmplitude,
                seasonalPeakDay, seasonalPeriod);

        this.susCompartmentCounts = getSusCompartmentCounts(numSS, numSI, numSC, numSR, numIS, numCS, numRS);

        transmissionRate.addVariableListener(this);
        for (Parameter p : susCompartmentCounts) {
            p.addVariableListener(this);
        }
        infectionRateModulationI.addVariableListener(this);
        infectionRateModulationC.addVariableListener(this);
        infectionRateModulationR.addVariableListener(this);
    }

    private double getTransmissionRate(int dim){
        if(seasonalModel){
            return transmissionRate.getParameterValue(0)*seasonalMultipliers[dim];
        }else{
            return transmissionRate.getParameterValue(0);
        }
    }

    private void setGridDimensions(int pathogenNumber, int numGridPoints, double samplingDateOne,
                                   double samplingDateTwo) {
        boolean isMoreRecent = false;
        if((pathogenNumber == 1 && samplingDateOne >= samplingDateTwo)
                || (pathogenNumber == 2 && samplingDateTwo >= samplingDateOne)) {
            isMoreRecent = true;
        }
        if(isMoreRecent){
            this.numSkippedIntervals = 0;
            this.numGridIntervals = numGridPoints;
        }else{
            double samplingDateDiff = Math.abs(samplingDateOne-samplingDateTwo);
            this.numSkippedIntervals = (int) Math.floor(samplingDateDiff/generalIntervalWidth);
            this.numGridIntervals = numGridPoints - numSkippedIntervals;
        }
    }

    private double[] getSeasonalMultipliers(double moreRecentDate, Parameter originTimeNumSS, Parameter seasonalAmplitude,
                                          Parameter seasonalPeakDay, Parameter seasonalPeriod) {
        if(seasonalModel){
            double[] returnVal = new double[numGridIntervals];
            double origNumSS = originTimeNumSS.getParameterValue(0);
            double seasonalAmp = seasonalAmplitude.getParameterValue(0);
            double seasonalPeak = seasonalPeakDay.getParameterValue(0);
            double period = seasonalPeriod.getParameterValue(0);

            for(int k = 0; k < numGridIntervals; k++){
                double calendarTime = moreRecentDate - (k+1+numSkippedIntervals)*generalIntervalWidth;
                returnVal[k] = origNumSS*(1.0 + seasonalAmp*Math.cos(2.0*Math.PI*(calendarTime-seasonalPeak)/period));
            }
            return returnVal;
        }else{
            return null;
        }
    }

    private List<Parameter> getSusCompartmentCounts(Parameter numSS, Parameter numSI, Parameter numSC, Parameter numSR,
                                                    Parameter numIS, Parameter numCS, Parameter numRS) {
        List<Parameter> susComptCounts = new ArrayList<>();
        susComptCounts.add(numSS);
        if(pathogenNumber == 1){
            susComptCounts.add(numSI);
            susComptCounts.add(numSC);
            susComptCounts.add(numSR);
        }else{
            susComptCounts.add(numIS);
            susComptCounts.add(numCS);
            susComptCounts.add(numRS);
        }
        return susComptCounts;
    }

    public int getDimension() {
        return numGridIntervals;
    }

    protected void storeValues() {
        transmissionRate.storeParameterValues();
        for (Parameter p : susCompartmentCounts) {
            p.storeParameterValues();
        }
        infectionRateModulationI.storeParameterValues();
        infectionRateModulationC.storeParameterValues();
        infectionRateModulationR.storeParameterValues();
    }

    protected void restoreValues() {
        transmissionRate.restoreParameterValues();
        for (Parameter p : susCompartmentCounts) {
            p.restoreParameterValues();
        }
        infectionRateModulationI.restoreParameterValues();
        infectionRateModulationC.restoreParameterValues();
        infectionRateModulationR.restoreParameterValues();
    }

    protected void acceptValues() {
        transmissionRate.acceptParameterValues();
        for (Parameter p : susCompartmentCounts) {
            p.acceptParameterValues();
        }
        infectionRateModulationI.acceptParameterValues();
        infectionRateModulationC.acceptParameterValues();
        infectionRateModulationR.acceptParameterValues();
    }

    protected void adoptValues(Parameter source) {
        throw new RuntimeException("Not implemented");
    }

    public double getParameterValue(int dim) {
        double tRate = getTransmissionRate(dim);
        double iRateModI = infectionRateModulationI.getParameterValue(0);
        double iRateModC = infectionRateModulationC.getParameterValue(0);
        double iRateModR = infectionRateModulationR.getParameterValue(0);

        int adjIndex = dim + numSkippedIntervals;
        double val = susCompartmentCounts.get(0).getParameterValue(adjIndex)
                + iRateModI*susCompartmentCounts.get(1).getParameterValue(adjIndex)
                + iRateModC*susCompartmentCounts.get(2).getParameterValue(adjIndex)
                + iRateModR*susCompartmentCounts.get(3).getParameterValue(adjIndex);
        //System.out.println("weighted compartment counts: " + val);
        //double result = val*tRate;
        //if (Double.isNaN(result) || result <= 0) {
        //    System.out.println("dim=" + dim + " adjIndex=" + adjIndex +
        //            " tRate=" + tRate + " val=" + val + " result=" + result);
        //}
        return val*tRate;
    }

    public void setParameterValue(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    public void setParameterValueNotifyChangedAll(int dim, double value){
        throw new RuntimeException("Not implemented");
    }

    public String getParameterName() {
        return getId();
    }

    public void addBounds(Bounds bounds) {
        this.bounds = bounds;
    }

    public Bounds<Double> getBounds() {
        if(bounds == null){
            return transmissionRate.getBounds();
        }else{
            return bounds;
        }
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not implemented.");
    }

    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        fireParameterChangedEvent(index,type);
    }
}

