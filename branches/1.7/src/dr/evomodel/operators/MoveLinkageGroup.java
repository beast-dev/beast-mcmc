package dr.evomodel.operators;

import dr.evolution.util.Taxon;
import dr.evomodel.tree.HiddenLinkageModel;
import dr.evomodelxml.operators.MoveLinkageGroupParser;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

/**
 * @author Aaron Darling
 */
public class MoveLinkageGroup extends SimpleMCMCOperator {

	HiddenLinkageModel hlm;	
	int readCount;
	int groupCount;
	
	public MoveLinkageGroup(HiddenLinkageModel hlm, double weight){
		this.hlm = hlm;
		readCount = hlm.getData().getReadsTaxa().getTaxonCount();
		groupCount = hlm.getLinkageGroupCount();
        setWeight(weight);
	}

	public double doOperation() throws OperatorFailedException {
		// pick a read uniformly at random, add it to a linkage group uniformly at random
		int r = MathUtils.nextInt(readCount);		
		Taxon read = hlm.getData().getReadsTaxa().getTaxon(r);
		int g=hlm.getLinkageGroupId(read);

		int newGroup = MathUtils.nextInt(groupCount);
		hlm.moveReadGroup(read, g, newGroup);

		// moves are symmetric -- same forward and backward proposal probabilities
		double logHastings = 0.0;
		return logHastings;
	}

	public String getOperatorName() {
		return MoveLinkageGroupParser.MOVE_LINKAGE_GROUP + "(" + hlm.getId() + ")";
	}

	public String getPerformanceSuggestion() {
		return "Ask Aaron Darling to write a better operator";
	}

}
