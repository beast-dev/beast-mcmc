package dr.app.beauti.options;

import dr.app.beauti.types.FixRateType;
import dr.evolution.datatype.DataType;

import java.io.Serializable;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class ClockModelGroup implements Serializable {

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
        setFixMean(rateTypeOption == FixRateType.FIX_MEAN);
    }

    public double getFixMeanRate() {
        return fixMeanRate;
    }

    public void setFixMeanRate(double fixMeanRate, BeautiOptions options) {
        this.fixMeanRate = fixMeanRate;
        for (PartitionClockModel model : options.getPartitionClockModels(this)) {
            model.setRate(fixMeanRate, false);
        }
    }

    public boolean contain(DataType dataType, BeautiOptions options) {
        for (AbstractPartitionData pd : options.getDataPartitions(this)) {
           if (pd.getDataType().getType() == dataType.getType()) {
                return true;
            }
        }
        return false;
    }

}
