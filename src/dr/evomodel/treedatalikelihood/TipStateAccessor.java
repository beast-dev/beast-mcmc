package dr.evomodel.treedatalikelihood;

public interface TipStateAccessor {

    void setTipStates(int tipNum, int[] states);

    void getTipStates(int tipNum, int[] states);

    int getPatternCount();

    int getTipCount();

    void makeDirty();
}
