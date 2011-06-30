package dr.app.beauti.components.dnds;

import java.util.List;

import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ComponentOptions;
import dr.app.beauti.options.ModelOptions;
import dr.app.beauti.options.Operator;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.types.OperatorType;

/**
 * @author Filip Bielejec
 * @version $Id$
 */

public class DnDsComponentOptions implements ComponentOptions {

	final private BeautiOptions options;
	static public final String CODON_PARTITIONED_ROBUST_COUNTING = "codon.partitioned.robust.counting";

	public DnDsComponentOptions(final BeautiOptions options) {
		this.options = options;
	}

	public void createParameters(ModelOptions modelOptions) {

		modelOptions.createOperator("CODON_PARTITIONED_ROBUST_COUNTING",
				OperatorType.BITFLIP, -1.0, -1.0);

	}

	public void selectOperators(ModelOptions modelOptions, List<Operator> ops) {

		ops.add(modelOptions.getOperator(CODON_PARTITIONED_ROBUST_COUNTING));

	}

	public void selectParameters(ModelOptions modelOptions,
			List<Parameter> params) {

		params
				.add(modelOptions
						.getParameter(CODON_PARTITIONED_ROBUST_COUNTING));

	}

	public void selectStatistics(ModelOptions modelOptions,
			List<Parameter> stats) {
		// TODO Auto-generated method stub

	}

}
