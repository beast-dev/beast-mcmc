package dr.app.beauti.options;

/**
 * @author Alexei Drummond
 */
public enum TreePrior {

	CONSTANT("Coalescent: Constant Size"),
	EXPONENTIAL("Coalescent: Exponential Growth"),
	LOGISTIC("Coalescent: Logistic Growth"),
	EXPANSION("Coalescent: Expansion Growth"),
	SKYLINE("Coalescent: Bayesian Skyline"),
	EXTENDED_SKYLINE("Coalescent: Extended Bayesian Skyline"),
	GMRF_SKYRIDE("Coalescent: GMRF Bayesian Skyride"),
	YULE("Speciation: Yule Process"),
	BIRTH_DEATH("Speciation: Birth-Death Process");

	TreePrior(String name) {
		this.name = name;
	}

	public String toString() {
		return name;
	}

	private final String name;
}
