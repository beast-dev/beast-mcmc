package dr.evolution.sequence;

import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;

/**
 * @author Marc A. Suchard
 */
public class DelimitedSequence extends Sequence {

    public DelimitedSequence(Taxon taxon, String sequence, DataType dataType) {
        super(taxon, sequence);
        this.dataType = dataType;
        codes = sequence.split(dataType.getDelimiter());
        original = sequence;
        sequenceString = null;
    }

    @Override
    public char getChar(int index) {
        throw new RuntimeException("Not single character available");
    }

    @Override
    public String getSequenceString() {
        return original;
    }

    @Override
    public int getLength() {
        return codes.length;
    }

    @Override
    public int getState(int index) {
        return dataType.getState(codes[index]);
    }

    @Override
    public void setState(int index, int state) {
        throw new RuntimeException("Not implemented");
    }

    private final String[] codes;

    private final String original;
}
