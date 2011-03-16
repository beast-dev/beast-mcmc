package dr.app.beauti.options;

import dr.app.beauti.types.FixRateType;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class ClockModelGroup {

    private String name;
    private boolean fixMean = false;
    private double fixMeanRate = 1.0;
    private FixRateType rateTypeOption = FixRateType.RELATIVE_TO;

//    public ClockModelGroup() { }

    public ClockModelGroup(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public boolean isFixMean() {
        return fixMean;
    }

    public void setFixMean(boolean fixMean) {
        this.fixMean = fixMean;
    }

    public FixRateType getRateTypeOption() {
        return rateTypeOption;
    }

    public void setRateTypeOption(FixRateType rateTypeOption) {
        this.rateTypeOption = rateTypeOption;
    }

    public double getFixMeanRate() {

        return fixMeanRate;
    }

    public void setFixMeanRate(double fixMeanRate) {
        this.fixMeanRate = fixMeanRate;
    }


}
