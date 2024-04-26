package dr.evomodel.substmodel.codon;

/**
 * @author Marc A. Suchard
 */
public class CodonOptions {

    final boolean isParameterTotalRate;

    public CodonOptions() {
        this.isParameterTotalRate = true;
    }

    public CodonOptions(boolean isParameterTotalRate) {
        this.isParameterTotalRate = isParameterTotalRate;
    }
}