package dr.evomodel.operators;

import java.util.ArrayList;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.ARGModel;
import dr.evomodel.tree.ARGModel.Node;
import dr.inference.model.Parameter;
import dr.inference.model.VariableSizeCompoundParameter;
import dr.inference.model.VariableSizeParameter;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.NormalDistribution;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class NewerARGEventOperator extends SimpleMCMCOperator implements CoercableMCMCOperator{

	public static final String NEWER_ARG_EVENT = "newerARGEvent";
	public static final String SWAP_RATES = "swapRates";
	public static final String SWAP_TRAITS = "swapTraits";
	public static final String MAX_VALUE = "maxTips";
	public static final String SINGLE_PARTITION = "singlePartitionProbability";
	public static final String IS_RECOMBINATION = "isRecombination";

	public static final String JUST_INTERNAL = "justInternalNodes";
	public static final String INTERNAL_AND_ROOT = "internalAndRootNodes";
	public static final String NODE_RATES = "nodeRates";
	public static final double LOG_TWO = 0.693147181;
	
	private double size = 0.0;
	private double singlePartitionProbability = 0.0;
	

	private boolean isRecombination = false;
	private int mode = CoercableMCMCOperator.DEFAULT;
	private ARGModel arg;
	private VariableSizeCompoundParameter internalNodeParameters;
	private VariableSizeCompoundParameter internalAndRootNodeParameters;
	private VariableSizeCompoundParameter nodeRates;

	
	public NewerARGEventOperator(ARGModel arg, int weight, double size, boolean gaussian,
            boolean swapRates, boolean swapTraits, int mode,
            VariableSizeCompoundParameter param1,
            VariableSizeCompoundParameter param2,
            VariableSizeCompoundParameter param3,
            double singlePartitionProbability, boolean isRecombination) {

			setWeight(weight);

			this.arg = arg;
			this.size = size;

			this.internalNodeParameters = param1;
			this.internalAndRootNodeParameters = param2;
			this.nodeRates = param3;
			this.singlePartitionProbability = singlePartitionProbability;
			this.isRecombination = isRecombination;
			this.mode = mode;
	}
	
	public double doOperation() throws OperatorFailedException {
		double logHastings = 0;
		try {
			if (MathUtils.nextGaussian() < size){
				logHastings = addOperation() + getAddLogQ();
			}else{
				logHastings = removeOperation() - getAddLogQ();
			}
		}catch (NoReassortmentNodeException ofe) {
					
		}catch (OperatorFailedException ofe){
			System.err.println("Catch: " + ofe.getMessage());
			System.exit(-1);
		}

		return logHastings;
	}
	
	/**
	 * Fills <code>list</code> with the <code>NodeRef</code> of nodes
	 * that can be removed.  As a by-product, this method also returns
	 * the size of <code>list</code>.  If <code>list</code> is <code>null</code>,
	 * the method will still return the number of nodes, it just
	 * will not fill <code>list</code> as it's <code>null</code>
	 * @param list an <code>ArrayList<NodeRef></code> that will be filled.
	 * @return The number of nodes that can be removed.
	 */
	private int findPotentialNodesToRemove(ArrayList<NodeRef> list) {
		int count = 0;
		int totalNumberOfNodes = arg.getNodeCount();

		for (int i = 0; i < totalNumberOfNodes; i++) {
			NodeRef node = arg.getNode(i);
			NodeRef leftParent = arg.getParent(node,ARGModel.LEFT);
			NodeRef rightParent = arg.getParent(node,ARGModel.RIGHT);
			if(  arg.isReassortment(node) && 
			    (arg.isBifurcation(leftParent) || arg.isBifurcation(rightParent))){
				if (list != null){
					list.add(node);
				}
				count++;
			}
		}
		return count;
	}
	
	
	/**
	 * Fills <code>list</code> with possible <code>NodeRef</code> that can 
	 * become possible attachment points for a new arg-event.  
	 * @param time the point on the tree where the arg-event will occur
	 * @param list the <code>ArrayList<NodeRef></code> that get's filled
	 * @return the number of possible attachment points
	 */
	private int findPotentialAttachmentPoints(double time, ArrayList<NodeRef> list) {
		if(time >= arg.getNodeHeight(arg.getRoot())){
			list.add(arg.getRoot());
			return 1;
		}
		int count = 0;
		for (int i = 0, n = arg.getNodeCount(); i < n; i++) {
			NodeRef nr = arg.getNode(i);
			if (!arg.isRoot(nr) && arg.getNodeHeight(nr) < time) {
				if (arg.getNodeHeight(arg.getParent(nr, ARGModel.LEFT)) > time) {
					if (list != null)
						list.add(nr);
					count++;
				}
				if (arg.isReassortment(nr) && arg.getNodeHeight(arg.getParent(nr, ARGModel.RIGHT)) > time) {
					if (list != null)
						list.add(nr);
					count++;
				}
			}
		}
		return count;
	}
	
	private static double nextExponential(double mean){
		return -1.0*Math.log(1 - MathUtils.nextDouble())*mean;
	}
	
	/**
	 * Attempts to add a re-assortment event to <code>ARGModel</code>
	 * @return the hastings ratio
	 * @throws OperatorFailedException if something goes totally wrong
	 */
	private double addOperation() throws OperatorFailedException {
		ArrayList<NodeRef> potentialBifurcationChildren = new ArrayList<NodeRef>();
		ArrayList<NodeRef> potentialReassortmentChildren = new ArrayList<NodeRef>();
		
		double logq = 0;
		double treeHeight = arg.getNodeHeight(arg.getRoot());
		double probabilityBeyondRoot = 0;
		double additionalHeightParameter = treeHeight / 4; 
		
		////Start of operator/////
		
		//This value always needs to fall below the current treeHeight.
		double newReassortmentHeight = treeHeight * MathUtils.nextDouble();
//		logq += Math.log(treeHeight);
				
		double newBifurcationHeight;
				
		if(MathUtils.nextDouble() < probabilityBeyondRoot){
			newBifurcationHeight = nextExponential(additionalHeightParameter) +	treeHeight;
//			logq += additionalHeightParameter*(newBifurcationHeight) 
//				- Math.log(probabilityBeyondRoot*additionalHeightParameter);
		}else{
			newBifurcationHeight = treeHeight * MathUtils.nextDouble();
//			logq += Math.log(treeHeight/(1-probabilityBeyondRoot));
		}

		
		if (newReassortmentHeight > newBifurcationHeight) {
			double temp = newReassortmentHeight;
			newReassortmentHeight = newBifurcationHeight;
			newBifurcationHeight = temp;
		}
				
		int totalPotentialBifurcationChildren = findPotentialAttachmentPoints(
				newBifurcationHeight, potentialBifurcationChildren);
		int totalPotentialReassortmentChildren = findPotentialAttachmentPoints(
				newReassortmentHeight, potentialReassortmentChildren);
	
		
		
		assert totalPotentialBifurcationChildren > 0 && totalPotentialReassortmentChildren > 0; 

		Node reassortChild = (Node) potentialReassortmentChildren
			.get(MathUtils.nextInt(totalPotentialReassortmentChildren));
		NodeRef reassortLeftParent = arg.getParent(reassortChild, 0);
		NodeRef reassortRightParent =  arg.getParent(reassortChild, 1);
		NodeRef reassortParent = reassortLeftParent;
//	    logq += Math.log(totalPotentialReassortmentChildren);

		

		if (reassortLeftParent != reassortRightParent){
			if (arg.getNodeHeight(reassortLeftParent) < newReassortmentHeight){
				reassortParent = reassortRightParent;
			}else if (arg.getNodeHeight(reassortRightParent) > newReassortmentHeight){
				
				if(MathUtils.nextDouble() > 0.5)	
					reassortParent = reassortRightParent;
//				logq += Math.log(2);
			}
		}

		Node bifurcateChild = (Node) potentialBifurcationChildren.get(MathUtils
				.nextInt(totalPotentialBifurcationChildren));
//	    logq += Math.log(totalPotentialBifurcationChildren);

		
		System.out.println(arg.toARGSummary());
		System.out.println(bifurcateChild);
		System.out.println(reassortLeftParent);
		System.out.println(reassortRightParent);
		System.exit(-1);
	   
		Node bifurcateLeftParent = bifurcateChild.leftParent;
		Node bifurcateRightParent = bifurcateChild.rightParent;
		Node bifurcateParent = bifurcateLeftParent;
		if (bifurcateLeftParent != bifurcateRightParent) {
			if (arg.getNodeHeight(bifurcateLeftParent) < newBifurcationHeight)
				bifurcateParent = bifurcateRightParent;
			else if (arg.getNodeHeight(bifurcateRightParent) > newBifurcationHeight
					&& MathUtils.nextDouble() > 0.5){
				bifurcateParent = bifurcateRightParent;
				logq += Math.log(2);
				
			}
		}

		
		double newBifurcationRateCategory = 1.0;
		double newReassortmentRateCategory = 1.0;

		Node newBifurcation = arg.new Node();
		newBifurcation.heightParameter = new Parameter.Default(newBifurcationHeight);
		newBifurcation.rateParameter = new Parameter.Default(newBifurcationRateCategory);
		newBifurcation.number = arg.getNodeCount();
		newBifurcation.setupHeightBounds();

		Node newReassortment = arg.new Node();
		newReassortment.bifurcation = false;
		newReassortment.heightParameter = new Parameter.Default(newReassortmentHeight);
		newReassortment.rateParameter = new Parameter.Default(newReassortmentRateCategory);
		newReassortment.number = arg.getNodeCount() + 1;
		newReassortment.setupHeightBounds();
		
		
		//Begin mutating the tree.
		arg.beginTreeEdit(); 
		if(bifurcateParent != null){
			arg.removeChild(bifurcateParent, bifurcateChild);
			arg.addChild(bifurcateParent,newBifurcation);
		}
		
		if (bifurcateChild != reassortChild){
			arg.removeChild(reassortParent, reassortChild);
			arg.singleAddChild(newBifurcation, bifurcateChild);
		}
		arg.doubleAddChild(newReassortment, reassortChild);
		
		VariableSizeParameter partitioning = new VariableSizeParameter(arg.getNumberOfPartitions());

        logq += drawRandomPartitioning(partitioning);


		if (bifurcateChild != reassortChild) {
			arg.addChildAsRecombinant(newBifurcation, reassortParent,
					newReassortment, partitioning);
		} else {
			arg.addChildAsRecombinant(newBifurcation, newBifurcation,
					newReassortment, partitioning);
		}

		if(bifurcateParent == null){
			arg.setRoot(newBifurcation);
		}
		arg.expandARGWithRecombinant(newBifurcation, newReassortment,
				internalNodeParameters,
				internalAndRootNodeParameters,
				nodeRates);
		
		arg.pushTreeSizeChangedEvent(); 

		try {
			arg.endTreeEdit();
		} catch (MutableTree.InvalidTreeException ite) {
			throw new RuntimeException(ite.toString() + "\n" + arg.toString()
					+ "\n" + Tree.Utils.uniqueNewick(arg, arg.getRoot()));
		}
		//End mutating the tree.
		
		
        logq -= Math.log(findPotentialNodesToRemove(null));
        if (!(reassortChild.leftParent.isBifurcation() && reassortChild.rightParent.isRoot()) &&
                !(reassortChild.rightParent.isBifurcation() && reassortChild.leftParent.isRoot())){
            logq -= Math.log(2.0);
        }
		return logq;
	}
	
	private double drawRandomPartitioning(Parameter partitioning) {
		double logq = 0;
		int len = arg.getNumberOfPartitions();
		if (len == 2) {
//            boolean first = MathUtils.nextBoolean();
			if (partitioning != null) {
				if (MathUtils.nextBoolean())
					partitioning.setParameterValueQuietly(0, 1.0);
				else
					partitioning.setParameterValueQuietly(1, 1.0);
			}
			return Math.log(2);
		}
		if (isRecombination) {
			logq += drawRandomRecombination(partitioning);
		} else {
			logq += drawRandomReassortment(partitioning);
		}
		return logq;
	}
	
	private double drawRandomReassortment(Parameter partitioning) {
		int len = arg.getNumberOfPartitions();
		double logq = 0;
		if (MathUtils.nextDouble() < singlePartitionProbability) {
			if (partitioning != null)
				partitioning.setParameterValueQuietly(MathUtils.nextInt(len), 1.0);
			return Math.log(len);
		}
		int[] permutation = MathUtils.permuted(len);
		int cut = MathUtils.nextInt(len - 1);
		for (int i = 0; i < len; i++) {
			logq += Math.log(i + 1);
			if (i > cut && partitioning != null)
				partitioning.setParameterValueQuietly(permutation[i], 1.0);
		}
		logq += Math.log(len - 1);
		return logq;
	}
	
	private double drawRandomRecombination(Parameter partitioning) {
		int len = arg.getNumberOfPartitions();
		double logq = 0;
		double leftValue = MathUtils.nextInt(2);
		double rightValue = 1.0 - leftValue;
		logq += Math.log(2);
		if (partitioning != null) {
			int cut = MathUtils.nextInt(len - 1);
			for (int i = 0; i <= cut; i++)
				partitioning.setParameterValueQuietly(i, leftValue);
			for (int i = cut + 1; i < len; i++)
				partitioning.setParameterValueQuietly(i, rightValue);
		}
		logq += Math.log(len - 1);
		return logq;
	}

	/**
	 * Attemps to remove a reassortment node from the <code>ARGModel</code>.
	 * @return The hastings ratio.
	 * @throws OperatorFailedException 
	 * 		if something goes totally wrong
	 * @throws NoReassortmentNodeException
	 * 		if there are no reassortment nodes to remove.
	 */
	private double removeOperation() throws OperatorFailedException,
		NoReassortmentNodeException{
		double logq = 0;

		
		ArrayList<NodeRef> potentialNodes = new ArrayList<NodeRef>();
		int totalPotentials = findPotentialNodesToRemove(potentialNodes);
				
		if (totalPotentials == 0){
			throw new NoReassortmentNodeException("No reassortment nodes to remove.");
		}
		
		arg.beginTreeEdit();
		boolean doneSomething = false;
		
		NodeRef removeNode = potentialNodes.get(MathUtils.nextInt(totalPotentials));
		double reverseReassortmentHeight = arg.getNodeHeight(removeNode);
		double reverseBifurcationHeight = 0;

		NodeRef removeLeftParent = arg.getParent(removeNode, ARGModel.LEFT);
		NodeRef removeLeftChild  = arg.getChild(removeNode, ARGModel.LEFT);
		
		if(arg.getParent(removeNode,ARGModel.LEFT) == arg.getParent(removeNode,ARGModel.RIGHT)){
			assert !arg.isRoot(removeLeftParent);
			
			NodeRef recGrandParent = arg.getParent(removeLeftParent,ARGModel.LEFT);
						
			arg.doubleRemoveChild(recGrandParent, removeLeftParent);
			arg.doubleRemoveChild(removeNode, removeLeftChild);
			arg.addChild(recGrandParent,removeLeftChild);
			
			doneSomething = true;
		}else{ 
			NodeRef keptParent = removeLeftParent;
			NodeRef removeParent = arg.getParent(removeNode,ARGModel.RIGHT);
									
			if(arg.isBifurcation(keptParent) 
					&& arg.isBifurcation(removeParent)){
				if(MathUtils.nextBoolean()){
					removeParent = removeLeftParent;
					keptParent = arg.getParent(removeNode,ARGModel.RIGHT);
				}
				logq += LOG_TWO;
			}else if(arg.isBifurcation(removeParent)){
				removeParent = removeLeftParent;
				keptParent = arg.getParent(removeNode,ARGModel.RIGHT);
			}
			
			if(arg.isRoot(removeParent)){
				NodeRef otherChild = arg.getOtherChild(removeParent, removeNode);
				
				arg.singleRemoveChild(removeParent, otherChild);
				arg.removeChild(keptParent, removeNode);
				arg.doubleRemoveChild(removeNode, removeLeftChild);
				
				//TODO Figure out how to link keptParent and removeLeftChild.
//				arg.addChild(keptParent, child)
				
				arg.setRoot(otherChild);
			}else{
				NodeRef removeGrandParent = arg.getParent(removeParent, ARGModel.LEFT);
				NodeRef otherChild = arg.getOtherChild(removeParent, removeNode);
						
				arg.removeChild(removeGrandParent, removeParent);
				arg.singleRemoveChild(removeParent, otherChild);
				arg.removeChild(keptParent,removeNode);
				arg.doubleRemoveChild(removeNode, removeLeftChild);
			
				
				if(otherChild != removeLeftChild){
					arg.addChild(removeGrandParent, otherChild);
					arg.addChild(keptParent, removeLeftChild);
				}else{
					System.exit(-1);
					arg.addChildWithSingleParent(removeGrandParent, otherChild);
					arg.addChildWithSingleParent(keptParent, removeLeftChild);
				}
			}
			doneSomething = true;
			removeLeftParent = removeParent;
		}
		
		
		assert sanityCheck();

		reverseBifurcationHeight = arg.getNodeHeight(removeLeftParent);

		if (doneSomething) {
			try {
				arg.contractARGWithRecombinant((Node)removeLeftParent, (Node)removeNode,
						internalNodeParameters, internalAndRootNodeParameters, nodeRates);
			} catch (Exception e) {
				System.err.println("here");
				System.err.println(e);

			}
		}

		arg.pushTreeSizeChangedEvent();
		try {
			arg.endTreeEdit();
		} catch (MutableTree.InvalidTreeException ite) {
			throw new RuntimeException(ite.toString() + "\n" + arg.toString()
					+ "\n" + Tree.Utils.uniqueNewick(arg, arg.getRoot()));
		}

		
//	    double d1 = findPotentialAttachmentPoints(reverseBifurcationHeight,null);
//	    double d2 = findPotentialAttachmentPoints(reverseReassortmentHeight,null);
//	    System.err.printf("d1 = %5.4f, d2 = %5.4f\n",d1,d2);

	    logq -= Math.log(findPotentialAttachmentPoints(reverseBifurcationHeight,null));
	    logq -= Math.log(findPotentialAttachmentPoints(reverseReassortmentHeight,null));

//        int nodes = arg.getInternalNodeCount() - 1;
//
        logq += 3*Math.log(2);
//
//        logq = 100;

//	    System.err.println("logq remove = "+logq);
//		logq = 0;
//	    logq -= 10;
	    logq += Math.log(totalPotentials);
		return logq; // 1 / total potentials * 1 / 2 (if valid) * length1 * length2 * attachmentSisters
	}
	/**
	 * Checks whether an ARGmodel makes sense after an operator has occured
	 * @return If the model makes sense.
	 */
	public boolean sanityCheck() {
		int len = arg.getNodeCount();
		for (int i = 0; i < len; i++) {
			Node node = (Node) arg.getNode(i);
			if (node.bifurcation) {
				boolean equalChild = (node.leftChild == node.rightChild);
				if ((equalChild && node.leftChild != null)) {
					if (!node.leftChild.bifurcation && ((node.leftChild).leftParent == node))
						;
					else {
						return false;
					}
				}
			} else {
				if ((node.leftChild != node.rightChild)) {
					return false;
				}
			}
			if (!node.isRoot()) {
				double d;
				d = node.getHeight();
			}
		}
		return true;
	}
	
	public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

		public String getParserName() {
			return NEWER_ARG_EVENT;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			boolean swapRates = false;
			boolean swapTraits = false;
			boolean isRecombination = false;

			double singlePartitionProbability = 0.0;

			//Get Attributes
			int mode = CoercableMCMCOperator.DEFAULT;
			if (xo.hasAttribute(AUTO_OPTIMIZE)) {
				if (xo.getBooleanAttribute(AUTO_OPTIMIZE)) {
					mode = CoercableMCMCOperator.COERCION_ON;
				} else {
					mode = CoercableMCMCOperator.COERCION_OFF;
				}
			}

			if (xo.hasAttribute(SINGLE_PARTITION)) {
				singlePartitionProbability = xo.getDoubleAttribute(SINGLE_PARTITION);
			}

			if (xo.hasAttribute(IS_RECOMBINATION)) {
				isRecombination = xo.getBooleanAttribute(IS_RECOMBINATION);

			}

			if (xo.hasAttribute(SWAP_RATES)) {
				swapRates = xo.getBooleanAttribute(SWAP_RATES);
			}
			if (xo.hasAttribute(SWAP_TRAITS)) {
				swapTraits = xo.getBooleanAttribute(SWAP_TRAITS);
			}
						
			int weight = xo.getIntegerAttribute("weight");
			double size = xo.getDoubleAttribute("size");
			boolean gaussian = xo.getBooleanAttribute("gaussian");

			//Get the arg model
			Object obj = xo.getChild(ARGModel.class);
			
			ARGModel treeModel = null;
			
			if (obj instanceof ARGModel) {
				treeModel = (ARGModel) obj;
			} 
			
			VariableSizeCompoundParameter parameter1 = (VariableSizeCompoundParameter) xo.getSocketChild(JUST_INTERNAL);
			VariableSizeCompoundParameter parameter2 = (VariableSizeCompoundParameter) xo.getSocketChild(INTERNAL_AND_ROOT);
			VariableSizeCompoundParameter parameter3 = (VariableSizeCompoundParameter) xo.getSocketChild(NODE_RATES);

			return new NewerARGEventOperator(treeModel, weight, size, gaussian, swapRates, swapTraits,
					mode, parameter1, parameter2, parameter3, singlePartitionProbability, isRecombination);
		}

		public String getParserDescription() {
			return "An operator that adds and removes reassortment nodes.";
		}

		public Class getReturnType() {
			return AddRemoveARGEventOperator.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newIntegerRule("weight"),
				AttributeRule.newDoubleRule("size"),
				AttributeRule.newBooleanRule("gaussian"),
				AttributeRule.newBooleanRule(SWAP_RATES, true),
				AttributeRule.newBooleanRule(SWAP_TRAITS, true),
				AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
				new ElementRule(ARGModel.class)
		};
	};

	public String getOperatorName() {
		return NEWER_ARG_EVENT;
	}

	public double getCoercableParameter() {
		return Math.log(size);
	}
	
	private double getAddLogQ(){
		double c = NormalDistribution.cdf(size, 0, 1);
		return Math.log((1-c)/c);
	}
	
	public double getSize() {
		return size;
	}

	public int getMode() {
		return mode;
	}

	public double getRawParameter() {
		return size;
	}

	public void setCoercableParameter(double value) {
		size = Math.exp(value);
	}

	public String getPerformanceSuggestion() {
		// TODO Get this guy working.
		return null;
	}

	private class NoReassortmentNodeException extends OperatorFailedException{
		private static final long serialVersionUID = 1L;

		public NoReassortmentNodeException(String a){
			super(a);
		}
	}

	
}
