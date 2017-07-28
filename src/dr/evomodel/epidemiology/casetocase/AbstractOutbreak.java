/*
 * AbstractOutbreak.java
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

package dr.evomodel.epidemiology.casetocase;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.inference.model.AbstractModel;
import dr.inference.model.Parameter;

import java.util.*;

/**
 * Abstract class for outbreaks. Implements PatternList for ease of compatibility with AbstractTreeLikelihood,
 * but there is one and only one pattern.
 *
 * User: Matthew Hall
 * Date: 14/04/13
 */

public abstract class AbstractOutbreak extends AbstractModel implements PatternList {

    protected GeneralDataType caseDataType;
    protected TaxonList taxa;
    private boolean hasLatentPeriods;
    protected final boolean hasGeography;
    private final String CASE_NAME = "hostID";
    protected ArrayList<AbstractCase> cases;
    protected int infectedSize = 0;

    public AbstractOutbreak(String name, Taxa taxa){
        this(name, taxa, false, true);
    }

    public AbstractOutbreak(String name, Taxa taxa, boolean hasLatentPeriods, boolean hasGeography){
        super(name);
        this.taxa = taxa;
        ArrayList<String> caseNames = new ArrayList<String>();
        for(int i=0; i<taxa.getTaxonCount(); i++){
            caseNames.add((String)taxa.getTaxonAttribute(i, CASE_NAME));
        }
        caseDataType = new GeneralDataType(caseNames);
        this.hasLatentPeriods = hasLatentPeriods;
        this.hasGeography = hasGeography;
    }

    public ArrayList<AbstractCase> getCases(){
        return new ArrayList<AbstractCase>(cases);
    }

    public boolean hasLatentPeriods() {
        return hasLatentPeriods;
    }

    public boolean hasGeography(){
        return hasGeography;
    }

    // todo this should be in terms of arbitary distance functions, not kernels

    public abstract double getLatentPeriod(AbstractCase aCase);

    public double getKernelValue(AbstractCase a, AbstractCase b, SpatialKernel kernel){
        if(!hasGeography){
            return 1;
        } else {
            return kernel.value(getDistance(a,b));
        }
    }

    // all the kernel values going TO case a (this is symmetric, usually, but potentially might not be)

    public double[] getKernelValues(AbstractCase aCase, SpatialKernel kernel){
        double[] out = new double[cases.size()];
        for(int i=0; i<out.length; i++){
            out[i] = kernel.value(getDistance(aCase, cases.get(i)));
        }
        return out;
    }

    public int getCaseIndex(AbstractCase thisCase){
        return cases.indexOf(thisCase);
    }

    public int size(){
        return cases.size();
    }

    public int infectedSize(){
        return infectedSize;
    }

    public abstract double getDistance(AbstractCase a, AbstractCase b);

    public AbstractCase getCase(int i){
        return cases.get(i);
    }

    public AbstractCase getCase(String name){
        for(AbstractCase thisCase: cases){
            if(thisCase.getName().equals(name)){
                return thisCase;
            }
        }
        return null;
    }

    public TaxonList getTaxa(){
        return taxa;
    }


    //************************************************************************
    // PatternList implementation
    //************************************************************************

    // not considering the possibility that we are simultaneously reconstructing more than one transmission tree!

    public int getPatternCount(){
        return 1;
    }

    public int getStateCount(){
        return size();
    }

    public int getPatternLength(){
        return taxa.getTaxonCount();
    }

    // @todo if these are never going to be used, get them to throw exceptions

    public int[] getPattern(int patternIndex){
        int[] out = new int[cases.size()];
        for(int i=0; i<cases.size(); i++){
            out[i] = i;
        }
        return out;
    }

    @Override
    public double[][] getUncertainPattern(int patternIndex) {
        throw new UnsupportedOperationException("uncertain patterns not implemented yet");
    }

    public int getPatternState(int taxonIndex, int patternIndex){
        return taxonIndex;
    }

    @Override
    public double[] getUncertainPatternState(int taxonIndex, int patternIndex) {
        throw new UnsupportedOperationException("uncertain patterns not implemented yet");
    }

    public double getPatternWeight(int patternIndex){
        return 1;
    }

    public double[] getPatternWeights(){
        return new double[]{1};
    }

    public double[] getStateFrequencies(){
        double[] out = new double[cases.size()];
        Arrays.fill(out, 1/cases.size());
        return out;
    }

    @Override
    public boolean areUnique() {
        return false;
    }

    @Override
    public boolean areUncertain() {
        return false;
    }

    public DataType getDataType(){
        return caseDataType;
    }

    //************************************************************************
    // TaxonList implementation
    //************************************************************************

    public int getTaxonCount(){
        return taxa.getTaxonCount();
    }

    public Taxon getTaxon(int taxonIndex){
        return taxa.getTaxon(taxonIndex);
    }

    public String getTaxonId(int taxonIndex){
        return taxa.getTaxonId(taxonIndex);
    }

    public int getTaxonIndex(String id){
        return taxa.getTaxonIndex(id);
    }

    public int getTaxonIndex(Taxon taxon){
        return taxa.getTaxonIndex(taxon);
    }

    public List<Taxon> asList(){
        return taxa.asList();
    }

    public Object getTaxonAttribute(int taxonIndex, String name){
        return taxa.getTaxonAttribute(taxonIndex, name);
    }

    public Iterator<Taxon> iterator() {
        if (taxa == null) throw new RuntimeException("Patterns has no TaxonList");
        return taxa.iterator();
    }

}
