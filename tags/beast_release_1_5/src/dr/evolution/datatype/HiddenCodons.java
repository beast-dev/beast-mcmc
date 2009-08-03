package dr.evolution.datatype;

/**
 * @author Marc A. Suchard
 */
public class HiddenCodons extends Codons implements HiddenDataType {
		
	public static final String DESCRIPTION = "hiddenCodon";

	public static final HiddenCodons UNIVERSAL_HIDDEN_2 = new HiddenCodons(GeneticCode.UNIVERSAL,2);
	public static final HiddenCodons UNIVERSAL_HIDDEN_3 = new HiddenCodons(GeneticCode.UNIVERSAL,3);
	
	/**
	 * Private constructor - DEFAULT_INSTANCE provides the only instance
	 */
	private HiddenCodons(GeneticCode geneticCode, int hiddenClassCount) {
		super(geneticCode);
		this.hiddenClassCount = hiddenClassCount;
	}

	/**
	 * returns an array containing the non-ambiguous states that this state represents.
	 */
	public boolean[] getStateSet(int state) {

	    boolean[] stateSet = new boolean[stateCount*hiddenClassCount];

	    if (!isAmbiguousState(state)) {
		    for(int h=0; h<hiddenClassCount; h++)
	            stateSet[h*stateCount + state] = true;
	    } else {
	        for (int i = 0; i < stateCount; i++) {
	            stateSet[i] = true;
	        }
	    }

	    return stateSet;
	}

	public int getStateCount() {
        return stateCount*hiddenClassCount;
    }

	private int hiddenClassCount;

    public int getHiddenClassCount() {
        return hiddenClassCount;
    }
}
