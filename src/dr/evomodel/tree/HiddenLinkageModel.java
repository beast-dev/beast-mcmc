/*
 * HiddenLinkageModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.tree;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import dr.evolution.MetagenomeData;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.oldevomodel.substmodel.SubstitutionModel;
import dr.oldevomodel.treelikelihood.GeneralLikelihoodCore;
import dr.oldevomodel.treelikelihood.LikelihoodCore;
import dr.oldevomodel.treelikelihood.NativeAminoAcidLikelihoodCore;
import dr.oldevomodel.treelikelihood.NativeNucleotideLikelihoodCore;
import dr.evomodel.tipstatesmodel.TipStatesModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;

/**
 * A model of linkage among sets of metagenome sequence reads
 * Reads are assumed to derive from a fixed number of organisms, where
 * the identity of the organisms is unknown, and the assignment of reads
 * to those organisms is also unknown.
 * @author Aaron Darling (koadman)
 *
 */
public class HiddenLinkageModel extends TipStatesModel implements PatternList
{

	int linkageGroupCount = 0;
	ArrayList< HashSet<Taxon> > groups = null;
	MetagenomeData data = null;
	ArrayList<Taxon> alignmentTaxa;
	double[][] tipPartials;
	double[][] storedTipPartials;
	boolean[] dirtyTipPartials;
	
	LikelihoodCore core;
	double blen = 0.001;
	SubstitutionModel substitutionModel;
	
	public HiddenLinkageModel(int linkageGroupCount, MetagenomeData data) {
		super("HiddenLinkageModel", data.getReferenceTaxa(), data.getReadsTaxa());
		this.linkageGroupCount = linkageGroupCount;
		this.data = data;

		// initial state: randomly assign reads to groups
		groups = new ArrayList< HashSet<Taxon> >(linkageGroupCount);
		for(int i=0; i<linkageGroupCount; i++)
			groups.add(new HashSet<Taxon>());
		TaxonList reads = data.getReadsTaxa();
		for(int i=0; i<reads.getTaxonCount(); i++){
			int g = MathUtils.nextInt(linkageGroupCount);
			groups.get(g).add(reads.getTaxon(i));
		}
		
		// create an alignment taxa list with reference + hidden groups
		alignmentTaxa = new ArrayList<Taxon>(data.getReferenceTaxa().asList());
		for(int i=0; i<linkageGroupCount; i++)
			alignmentTaxa.add(new Taxon("LinkageGroup_" + i));
		
		int plen = data.getAlignment().getSiteCount() * data.getAlignment().getStateCount();
		tipPartials = new double[alignmentTaxa.size()][plen];
		storedTipPartials = new double[alignmentTaxa.size()][plen];
		dirtyTipPartials = new boolean[alignmentTaxa.size()];
		
		initCore();
		setupMatrices();

		// compute initial partials
		for(int i=0; i<tipPartials.length; i++)
			computeTipPartials(i);
	}

	double[] tipMatrix;
	double[] internalMatrix;

	@Override
	public boolean areUnique() {
		return false;
	}

	@Override
	public boolean areUncertain() {
		return false;
	}

	/*
         * Initializes a likelihoodCore to calculate likelihoods for
         * the tips
         */
	private void initCore(){
		if(data.getAlignment().getDataType() instanceof dr.evolution.datatype.Nucleotides)
			core = new NativeNucleotideLikelihoodCore();
		if(data.getAlignment().getDataType() instanceof dr.evolution.datatype.AminoAcids)
			core = new NativeAminoAcidLikelihoodCore();
		if(data.getAlignment().getDataType() instanceof dr.evolution.datatype.Codons)
			core = new GeneralLikelihoodCore(data.getAlignment().getStateCount());
		// initialize the likelihood core
		core.initialize(data.getReadsTaxa().getTaxonCount()*2, data.getAlignment().getSiteCount(), 1, false);
		for(int i=0; i<data.getReadsTaxa().getTaxonCount(); i++){
			int index = data.getAlignment().getTaxonIndex(data.getReadsTaxa().getTaxon(i));
			int[] states = new int[data.getAlignment().getSiteCount()];
			for(int j=0; j<states.length; j++)
				states[j]=data.getAlignment().getState(index,j);
			core.setNodeStates(i, states);
		}
		for(int i=0; i<data.getReadsTaxa().getTaxonCount(); i++)
			core.createNodePartials(i+data.getReadsTaxa().getTaxonCount());
	}

	/*
	 * set up transition matrices for internal and external nodes
	 */
	private void setupMatrices(){
		tipMatrix=new double[data.getAlignment().getStateCount()*data.getAlignment().getStateCount()];
		internalMatrix=new double[data.getAlignment().getStateCount()*data.getAlignment().getStateCount()];
		double diag = 1.0-blen;
		double offdiag = blen / (data.getAlignment().getStateCount() - 1);
		double internalDiag = 0.99999999999999;
		double internalOffDiag = (1.0-internalDiag) / (data.getAlignment().getStateCount() - 1);
		for(int i=0; i<tipMatrix.length; i++)
		{
			tipMatrix[i]=offdiag;
			internalMatrix[i]=internalOffDiag;
		}
		for(int i=0; i<data.getAlignment().getStateCount(); i++)
		{
			tipMatrix[i*data.getAlignment().getStateCount() + i] = diag;
			internalMatrix[i*data.getAlignment().getStateCount() + i] = internalDiag;
		}

		for(int i=0; i<data.getReadsTaxa().getTaxonCount(); i++)
			core.setNodeMatrix(i, 0, tipMatrix);
		for(int i=0; i<data.getReadsTaxa().getTaxonCount(); i++)
			core.setNodeMatrix(i+data.getReadsTaxa().getTaxonCount(), 0, internalMatrix);
	}

	public int getLinkageGroupCount() {
		return linkageGroupCount;
	}

	public MetagenomeData getData() {
		return data;
	}


	/**
	 * Returns the group ID to which a particular metagenomic read belongs
	 * @param t
	 * @return
	 */
	public int getLinkageGroupId(Taxon t)
	{
		int i=0;
		for(HashSet<Taxon> h : groups){
			if(h.contains(t))
				break;
			i++;
		}
		return i;
	}
	
	private class Move {
		public Move(Taxon read, int fromGroup, int toGroup){
			this.read = read;
			this.fromGroup = fromGroup;
			this.toGroup = toGroup;
		}
		Taxon read;
		int fromGroup;
		int toGroup;
	}
	ArrayList<Move>	movesMade = new ArrayList<Move>();
	
	@Override
	protected void acceptState() {
		movesMade.clear();
		for(int i=0; i<dirtyTipPartials.length; i++)
			dirtyTipPartials[i]=false;
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// TODO Auto-generated method stub

	}

	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		// TODO Auto-generated method stub

	}

	protected void restoreState() {
		// make all moves in reverse
		for(int i=movesMade.size(); i>0; i--)
		{
			Move m = movesMade.get(i-1);
			groups.get(m.toGroup).remove(m.read);
			groups.get(m.fromGroup).add(m.read);
		}
		movesMade.clear();
		for(int i=0; i<dirtyTipPartials.length; i++)
		{
			if(dirtyTipPartials[i]){
				swapTipPartials(i);
				dirtyTipPartials[i]=false;
			}
		}
	}

	protected void storeState() {
		movesMade.clear();
		for(int i=0; i<dirtyTipPartials.length; i++)
			dirtyTipPartials[i]=false;
	}

	/**
	 * Get the members of the i'th linkage group
	 * @param i
	 * @return Members of the i'th group
	 */
	public Set<Taxon> getGroup(int i){
		return groups.get(i);
	}
	/**
	 * Moves a read from one linkage group to another linkage group
	 * @param read
	 * @param oldGroup
	 * @param newGroup
	 */
	public void moveReadGroup(Taxon read, int oldGroup, int newGroup)
	{
		boolean found = groups.get(oldGroup).remove(read);
		if(!found)
			throw new RuntimeException("Error, could not find read " + read + " in linkage group " + oldGroup);
		groups.get(newGroup).add(read);
		movesMade.add(new Move(read, oldGroup, newGroup));
		computeTipPartials(data.getReferenceTaxa().getTaxonCount() + oldGroup);
		computeTipPartials(data.getReferenceTaxa().getTaxonCount() + newGroup);

//		this.fireModelChanged();
		this.fireModelChanged(alignmentTaxa.get(alignmentTaxa.size() - groups.size() + oldGroup));
		this.fireModelChanged(alignmentTaxa.get(alignmentTaxa.size() - groups.size() + newGroup));
	}

	private void swapTipPartials(int nodeIndex){
		double[] tmp = storedTipPartials[nodeIndex];
		storedTipPartials[nodeIndex] = tipPartials[nodeIndex];
		tipPartials[nodeIndex] = tmp;
	}

	private void computeTipPartials(int nodeIndex){
		if(!dirtyTipPartials[nodeIndex]){
			swapTipPartials(nodeIndex);
			dirtyTipPartials[nodeIndex]=true;
		}
		double[] tipPartials = this.tipPartials[nodeIndex];

		// if this is one of the reference organisms, then return the resolved partials
		// if this is a linkage group, return partials that correspond to probabilities of each nucleotide.
		Alignment aln = data.getAlignment();
		int sc = aln.getStateCount();
		for(int i=0; i<tipPartials.length; i++){
			tipPartials[i]=0.0;
		}
		if(nodeIndex < data.getReferenceTaxa().getTaxonCount()){			
			int j=0;
			for(int i=0; i<aln.getSiteCount(); i++){
				int s = aln.getState(nodeIndex, i);
				if(s>=sc){
					for(int k=0; k<sc; k++)
						tipPartials[j + k] = 1.0;
				}else
					tipPartials[j + s] = 1.0;
				j += sc;
			}
		}else{
			int gI = nodeIndex - data.getReferenceTaxa().getTaxonCount();
			HashSet<Taxon> group = groups.get(gI);
			int internalNum = data.getReadsTaxa().getTaxonCount();
			Taxon firstTax=null;
			boolean peeled = false;
			for(Taxon tax : group){
				if(firstTax==null)
				{
					firstTax = tax;
					continue;
				}
				int c2 = data.getReadsTaxa().getTaxonIndex(tax);
				if(!peeled){
					int c1 = data.getReadsTaxa().getTaxonIndex(firstTax);
					core.setNodePartialsForUpdate(internalNum);
	                core.calculatePartials(c1, c2, internalNum);
				}else{
					core.setNodePartialsForUpdate(internalNum);
	                core.calculatePartials(internalNum-1, c2, internalNum);
				}
				internalNum++;
				peeled = true;
			}
			if(group.size()==0)
			{
				for(int i=0; i<tipPartials.length; i++)
					tipPartials[i]=1.0;
			}else if(!peeled)
				getPartialsForGroupSizeOne(firstTax, tipPartials);
			else
				core.getPartials(internalNum-1, tipPartials);
		}
	}
	private void getPartialsForGroupSizeOne(Taxon tax, double[] tipPartials)
	{
		Alignment aln = data.getAlignment();
		int sc = aln.getStateCount();
		int index = aln.getTaxonIndex(tax);
		int j=0;
		for(int i=0; i<aln.getSiteCount(); i++){
			int s = aln.getState(index, i);
			if(s>=sc){
				for(int k=0; k<sc; k++)
					tipPartials[j + k] = 1.0;
			}else
				System.arraycopy(internalMatrix, s*sc, tipPartials, j, sc);
			j += sc;
		}
	}
/*	
	private void computeTipPartials(int nodeIndex){
		if(!dirtyTipPartials[nodeIndex]){
			swapTipPartials(nodeIndex);
			dirtyTipPartials[nodeIndex]=true;
		}
		double[] tipPartials = this.tipPartials[nodeIndex];

		// if this is one of the reference organisms, then return the resolved partials
		// if this is a linkage group, return partials that correspond to probabilities of each nucleotide.
		Alignment aln = data.getAlignment();
		int sc = aln.getStateCount();
		for(int i=0; i<tipPartials.length; i++){
			tipPartials[i]=0.0;
		}
		if(nodeIndex < data.getReferenceTaxa().getTaxonCount()){			
			int j=0;
			for(int i=0; i<aln.getSiteCount(); i++){
				int s = aln.getState(nodeIndex, i);
				if(s>=sc){
					for(int k=0; k<sc; k++)
						tipPartials[j + k] = 1.0;
				}else
					tipPartials[j + s] = 1.0;
				j += sc;
			}
		}else{
			// average the information for each sequence in the linkage group
			int gI = nodeIndex - data.getReferenceTaxa().getTaxonCount();
			HashSet<Taxon> group = groups.get(gI);
			for( Taxon tax : group){
				int sI = data.getAlignment().getTaxonIndex(tax);

				int j=0;
				for(int i=0; i<aln.getSiteCount(); i++){
					int s = aln.getState(sI, i);
					if(s<sc)
						tipPartials[j + s] += 1.0;
					j += sc;
				}
			}

			// now normalize back to probability distributions
			int j=0;
			int l=tipPartials.length / sc;
			for(int i=0; i<l; i++){
				double max=0;
				for(int k=0; k<sc; k++){
					max = max > tipPartials[j+k] ? max : tipPartials[j+k];
				}
				if(max>0){
					for(int k=0; k<sc; k++){
						tipPartials[j+k]/=max;
					}
				}else{
					// if no sequence had info, set it to unknown (1.0 for all values)
					for(int k=0; k<sc; k++){
						tipPartials[j+k]=1.0;
					}
				}
				j += sc;
			}
		}
	}
*/
	/**
	 * NOT YET IMPLEMENTED. Creates a new read linkage group.
	 * @return the new group id
	 */
	public int newGroup(){
		throw new RuntimeException("Not implemented!");
	}
	/**
	 * NOT YET IMPLEMENTED. Deletes a read linkage group.  Must be empty.
	 */
	public void deleteGroup(){
		throw new RuntimeException("Not implemented!");
	}

    @Override
    public Type getModelType() {
        return Type.PARTIALS;
    }

    @Override
    public void getTipStates(int nodeIndex, int[] tipStates) {
        throw new IllegalArgumentException("This model emits only tip partials");
    }

    @Override
	public void getTipPartials(int nodeIndex, double[] tipPartials) {
		int n = nodeIdToMyTaxaMap[tree.getNode(nodeIndex).getNumber()];
		System.arraycopy(this.tipPartials[n], 0, tipPartials, 0, tipPartials.length);
	}

	int[] nodeIdToMyTaxaMap;
	protected void taxaChanged() {
		nodeIdToMyTaxaMap = new int[tree.getNodeCount()];
		for(int i=0; i<nodeIdToMyTaxaMap.length; i++){
			for(int j=0; j<alignmentTaxa.size(); j++){
				if(tree.getTaxon(i)==null)
					continue;
				if(tree.getTaxon(i)==null || alignmentTaxa.get(j) == null)
					System.err.print("asdgasdg\n");
				else if(tree.getTaxon(i).getId()==null || alignmentTaxa.get(j).getId() == null)
					System.err.print("asdgasdg\n");
				if(tree.getTaxon(i).getId().equalsIgnoreCase(alignmentTaxa.get(j).getId())){
					nodeIdToMyTaxaMap[tree.getExternalNode(i).getNumber()] = j;
					break;
				}
			}
		}
	}

	
	
	
	//
	// BEGIN PatternList Implementation
	// This merely delegates to Alignment for most methods
	//
	
	public DataType getDataType() {
		return data.getAlignment().getDataType();
	}

	public int[] getPattern(int patternIndex) {
		return data.getAlignment().getPattern(patternIndex);
	}

	@Override
	public double[][] getUncertainPattern(int patternIndex) {
		return new double[0][];
	}

	public int getPatternCount() {
		return data.getAlignment().getPatternCount();
	}

	public int getPatternLength() {
		return data.getAlignment().getPatternLength();
	}

	public int getPatternState(int taxonIndex, int patternIndex) {
		if(taxonIndex<data.getReferenceTaxa().getTaxonCount())
			return data.getAlignment().getPatternState(taxonIndex, patternIndex);

		return 0;
	}

	@Override
	public double[] getUncertainPatternState(int taxonIndex, int patternIndex) {
		return new double[0];
	}

	public double getPatternWeight(int patternIndex) {
		return data.getAlignment().getPatternWeight(patternIndex);
	}

	public double[] getPatternWeights() {
		return data.getAlignment().getPatternWeights();
	}

	public int getStateCount() {
		return data.getAlignment().getStateCount();
	}

	public double[] getStateFrequencies() {
		return data.getAlignment().getStateFrequencies();
	}

	public List<Taxon> asList() {
		return alignmentTaxa;
	}

	public Taxon getTaxon(int taxonIndex) {
		return alignmentTaxa.get(taxonIndex);
	}

	public Object getTaxonAttribute(int taxonIndex, String name) {
		return alignmentTaxa.get(taxonIndex).getAttribute(name);
	}

	public int getTaxonCount() {
		return alignmentTaxa.size();
	}

	public String getTaxonId(int taxonIndex) {
		return alignmentTaxa.get(taxonIndex).getId();
	}

	public int getTaxonIndex(String id) {
		for(int i=0; i<alignmentTaxa.size(); i++){
			if(alignmentTaxa.get(i).getId().equals(id))
				return i;
		}
		return -1;
	}

	public int getTaxonIndex(Taxon taxon) {
		for(int i=0; i<alignmentTaxa.size(); i++){
			if(alignmentTaxa.get(i).compareTo(taxon)==0)
				return i;
		}
		return -1;
	}

	public Iterator<Taxon> iterator() {
		return alignmentTaxa.iterator();
	}
}
