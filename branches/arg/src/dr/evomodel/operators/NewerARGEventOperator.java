package dr.evomodel.operators;

import java.util.ArrayList;
import java.util.logging.Logger;

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

		NodeRef root = arg.getRoot();
		for (int i = 0; i < totalNumberOfNodes; i++) {
			NodeRef node = arg.getNode(i);
			NodeRef leftParent = arg.getParent(node,ARGModel.LEFT);
			NodeRef rightParent = arg.getParent(node,ARGModel.RIGHT);
			if (arg.isReassortment(node) && 
					((leftParent != root && arg.isBifurcation(leftParent)) ||
					  (rightParent != root && arg.isBifurcation(rightParent)))){
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
	
	/**
	 * Attempts to add a reassortment event to <code>ARGModel</code>
	 * @return the hastings ratio
	 * @throws OperatorFailedException if something goes totally wrong
	 */
	private double addOperation() throws OperatorFailedException {
		double logq = 0;
		
		ArrayList<NodeRef> potentialBifurcationChildren = new ArrayList<NodeRef>();
		ArrayList<NodeRef> potentialReassortmentChildren = new ArrayList<NodeRef>();

		double treeHeight = arg.getNodeHeight(arg.getRoot());
				
		double probabilityBeyondRoot = 0.05;
		double addHeightParm = 2 / treeHeight; //This is 1/mean of the additional height
		double newBifurcationHeight;
		double newReassortmentHeight;
		
		//Choose locations for the new event.
		if(MathUtils.nextDouble() < probabilityBeyondRoot){
			newBifurcationHeight = - Math.log(1 - MathUtils.nextDouble())
															/addHeightParm;
			logq += addHeightParm*(newBifurcationHeight) 
				- Math.log(probabilityBeyondRoot*addHeightParm);
			newBifurcationHeight += treeHeight;
		}else{
			newBifurcationHeight = treeHeight * MathUtils.nextDouble();
			logq += Math.log(treeHeight/(1-probabilityBeyondRoot));
		}
		if(MathUtils.nextDouble() < probabilityBeyondRoot){
			newReassortmentHeight = -Math.log(1 - MathUtils.nextDouble())
															/addHeightParm;
			logq += addHeightParm*(newReassortmentHeight) 
				- Math.log(probabilityBeyondRoot*addHeightParm);
			newReassortmentHeight += treeHeight;
		}else{
			newReassortmentHeight = treeHeight * MathUtils.nextDouble();
			logq += Math.log(treeHeight/(1-probabilityBeyondRoot));
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
	
		//TODO Marc, should this be an assert?
		assert totalPotentialBifurcationChildren > 0 && totalPotentialReassortmentChildren > 0; 
			
		Node recNode = (Node) potentialReassortmentChildren.get(MathUtils
				.nextInt(totalPotentialReassortmentChildren));
	    logq += Math.log(totalPotentialReassortmentChildren);

	    
		Node recParentL = recNode.leftParent;
		Node recParentR = recNode.rightParent;
		Node recParent = recParentL;
		
		if (recParentL != recParentR){
			if (arg.getNodeHeight(recParentL) < newReassortmentHeight){
				recParent = recParentR;
			}else if (arg.getNodeHeight(recParentR) > newReassortmentHeight
					&& MathUtils.nextDouble() > 0.5){
				recParent = recParentR;
				logq += Math.log(2);
			}
		}


		Node sisNode = (Node) potentialBifurcationChildren.get(MathUtils
				.nextInt(totalPotentialBifurcationChildren));
	    logq += Math.log(totalPotentialBifurcationChildren);

	   
		Node sisParentL = sisNode.leftParent;
		Node sisParentR = sisNode.rightParent;
		Node sisParent = sisParentL;
		if (sisParentL != sisParentR) {
			if (arg.getNodeHeight(sisParentL) < newBifurcationHeight)
				sisParent = sisParentR;
			else if (arg.getNodeHeight(sisParentR) > newBifurcationHeight
					&& MathUtils.nextDouble() > 0.5){
				sisParent = sisParentR;
				logq += Math.log(2);
			}
		}

		double newBifurcationRateCategory = 1.0;
		double newReassortmentRateCategory = 1.0;

		Node newBifurcation = arg.new Node();
		newBifurcation.heightParameter = new Parameter.Default(newBifurcationHeight);
		newBifurcation.rateParameter = new Parameter.Default(newBifurcationRateCategory);
		newBifurcation.setupHeightBounds();

		Node newReassortment = arg.new Node();
		newReassortment.bifurcation = false;
		newReassortment.heightParameter = new Parameter.Default(newReassortmentHeight);
		newReassortment.rateParameter = new Parameter.Default(newReassortmentRateCategory);
		newReassortment.setupHeightBounds();

		//Begin mutating the tree.
		arg.beginTreeEdit(); 
		arg.removeChild(sisParent, sisNode);
		if (sisNode != recNode){
			arg.removeChild(recParent, recNode);
		}
		arg.addChild(sisParent,newBifurcation);
		
		if (sisNode != recNode)
			arg.singleAddChild(newBifurcation, sisNode);
		arg.doubleAddChild(newReassortment, recNode);

		VariableSizeParameter partitioning = new VariableSizeParameter(arg.getNumberOfPartitions());

        logq += drawRandomPartitioning(partitioning);


		if (sisNode != recNode) {
			arg.addChildAsRecombinant(newBifurcation, recParent,
					newReassortment, partitioning);
		} else {
			arg.addChildAsRecombinant(newBifurcation, newBifurcation,
					newReassortment, partitioning);
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
        if (!(recNode.leftParent.isBifurcation() && recNode.rightParent.isRoot()) &&
                !(recNode.rightParent.isBifurcation() && recNode.leftParent.isRoot())){
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

		// 1. Get a list of all the reassortment nodes.
		ArrayList<NodeRef> potentialNodes = new ArrayList<NodeRef>();
		int totalPotentials = findPotentialNodesToRemove(potentialNodes);
				
		//2. If there aren't any, throw a failure exception.
		if (totalPotentials == 0){
			throw new NoReassortmentNodeException("No reassortment nodes to remove.");
		}
		arg.beginTreeEdit();
		boolean doneSomething = false;
		
		//3. Choose a reassortment node randomly
		NodeRef recNode = potentialNodes.get(MathUtils.nextInt(totalPotentials));

		double reverseReassortmentHeight = arg.getNodeHeight(recNode);
		double reverseBifurcationHeight = 0;

		NodeRef recLeftParent = arg.getParent(recNode, ARGModel.LEFT);
		NodeRef recLeftChild  = arg.getChild(recNode, ARGModel.RIGHT);
		
		if(arg.getParent(recNode,ARGModel.LEFT) == arg.getParent(recNode,ARGModel.RIGHT)){
			NodeRef recGrandParent = arg.getParent(recLeftParent,ARGModel.LEFT);
						
			arg.doubleRemoveChild(recGrandParent, recLeftParent);
			arg.doubleRemoveChild(recNode, recLeftChild);
			arg.addChild(recGrandParent,recLeftChild);
			
			doneSomething = true;
		}else{ 
			NodeRef recParent1 = recLeftParent;
			NodeRef recParent2 = arg.getParent(recNode,ARGModel.RIGHT);
						
			//TODO Figure out what these orientations are.
			if(arg.isReassortment(recParent1) || arg.isRoot(recParent2)){ 
				recParent1 = arg.getParent(recNode,ARGModel.RIGHT);
				recParent2 = arg.getParent(recNode,ARGModel.LEFT);
			}else if(arg.isBifurcation(recParent1) && !arg.isRoot(recParent2)){ 
				if (MathUtils.nextDouble() < 0.5) { 
					recParent1 = arg.getParent(recNode,ARGModel.RIGHT);
					recParent2 = arg.getParent(recNode,ARGModel.LEFT);
				}
				logq += Math.log(2);
			}
			NodeRef recGrandParent = arg.getParent(recParent1, ARGModel.LEFT);
			NodeRef otherChild = arg.getChild(recParent1, ARGModel.LEFT);
					
			if (otherChild == recNode){
				otherChild = arg.getChild(recParent1, ARGModel.RIGHT);
			}
			//This part takes care of removing and readjusting the arg model after the deletion.
			arg.removeChild(recGrandParent, recParent1);
			arg.singleRemoveChild(recParent1, otherChild);
			arg.removeChild(recParent2,recNode);
			arg.doubleRemoveChild(recNode, recLeftChild);
			
			if(otherChild != recLeftChild){
				arg.addChild(recGrandParent, otherChild);
				arg.addChild(recParent2, recLeftChild);
			}else{
				arg.addChildWithSingleParent(recGrandParent, otherChild);
				arg.addChildWithSingleParent(recParent2, recLeftChild);
			}

			doneSomething = true;

			recLeftParent = recParent1;
		}
		assert sanityCheck();

		reverseBifurcationHeight = arg.getNodeHeight(recLeftParent);

		if (doneSomething) {
			try {
				arg.contractARGWithRecombinant((Node)recLeftParent, (Node)recNode,
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
