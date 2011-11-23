package dr.inference.trace;

import dr.app.gui.chart.Axis;
import dr.app.gui.chart.LinearAxis;
import dr.stats.Variate;
import dr.util.FrequencyDistribution;

/**
 * Creates a simple histogram-based density estimate
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class DensityEstimate {

    public DensityEstimate(Double[] samples, int minimumBinCount) {
        this.samples = new Variate.D(samples);

        calculateDensity(this.samples, minimumBinCount);
    }

    protected void calculateDensity(Variate.D data, int minimumBinCount) {
        frequencyDistribution = calculateFrequencies(samples, minimumBinCount);

        xCoordinates = new Variate.D();
        yCoordinates = new Variate.D();

        double x = frequencyDistribution.getLowerBound() - frequencyDistribution.getBinSize();

        xCoordinates.add(x + (frequencyDistribution.getBinSize() / 2.0));
        yCoordinates.add(0.0);
        x += frequencyDistribution.getBinSize();

        for (int i = 0; i < frequencyDistribution.getBinCount(); i++) {
            xCoordinates.add(x + (frequencyDistribution.getBinSize() / 2.0));
            double density = frequencyDistribution.getFrequency(i) / frequencyDistribution.getBinSize() / data.getCount();
            yCoordinates.add(density);
            x += frequencyDistribution.getBinSize();
        }

        xCoordinates.add(x + (frequencyDistribution.getBinSize() / 2.0));
        yCoordinates.add(0.0);
    }

    protected FrequencyDistribution calculateFrequencies(Variate data, int minimumBinCount) {
        double min = (Double) data.getMin();
        double max = (Double) data.getMax();

        if (min == max) {
            if (min == 0) {
                min = -1.0;
            } else {
                min -= Math.abs(min / 10.0);
            }
            if (max == 0) {
                max = 1.0;
            } else {
                max += Math.abs(max / 10.0);
            }
        }

        Axis axis = new LinearAxis(Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK);
        axis.setRange(min, max);

        int majorTickCount = axis.getMajorTickCount();
        axis.setPrefNumTicks(majorTickCount, 4);

        double binSize = axis.getMinorTickSpacing();
        int binCount = (int) ((axis.getMaxAxis() - axis.getMinAxis()) / binSize) + 2;

        if (minimumBinCount > 0) {
            while (binCount < minimumBinCount) {
                majorTickCount++;
                axis.setPrefNumTicks(majorTickCount, 4);

                binSize = axis.getMinorTickSpacing();
                binCount = (int) ((axis.getMaxAxis() - axis.getMinAxis()) / binSize) + 2; // should +2, otherwise the last bar will lose
            }
        }

        FrequencyDistribution frequency = new FrequencyDistribution(axis.getMinAxis(), binCount, binSize);

        for (int i = 0; i < data.getCount(); i++) {
            frequency.addValue((Double) data.get(i));
        }

        return frequency;
    }


    public Variate.D getXCoordinates() {
        return xCoordinates;
    }

    public Variate.D getYCoordinates() {
        return yCoordinates;
    }

    public FrequencyDistribution getFrequencyDistribution() {
        return frequencyDistribution;
    }

    protected final Variate.D samples;
    protected Variate.D xCoordinates;
    protected Variate.D yCoordinates;
    protected FrequencyDistribution frequencyDistribution;
}