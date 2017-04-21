/*
 * GaussianProcessSkytrackTreeOperator.java
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

package dr.evomodel.coalescent.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.coalescent.GaussianProcessSkytrackLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.operators.GaussianProcessSkytrackTreeOperatorParser;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorUtils;
import dr.math.MathUtils;
import dr.util.ComparableDouble;
import dr.util.HeapSort;
import no.uib.cipr.matrix.*;

//import javax.sound.midi.SysexMessage;
import java.util.ArrayList;
//import java.util.logging.Logger;

/* A Metropolis-Hastings operator to update the GP proposal when a new tree is proposed
 *
 * @author Julia Palacios
 * @version $Id:  $
 */
public class GaussianProcessSkytrackTreeOperator extends AbstractCoercableOperator {

	private TreeModel tree = null;
	private double scaleFactor;
//    private double lambdaScaleFactor;
//    private int fieldLength;

//    private int maxIterations;
//    private double stopValue;

	private Parameter popSizeParameter;
	private Parameter changePoints;

	private Parameter CoalCounts;
	private Parameter GPtype;
	private Parameter coalfactor;
	private Parameter precisionParameter;
	private Parameter lambdaParameter;
	private GaussianProcessSkytrackBlockUpdateOperator GPOperator;

	GaussianProcessSkytrackLikelihood gpLikelihood;

	private double[] zeros;

	public GaussianProcessSkytrackTreeOperator(TreeModel treeModel,GaussianProcessSkytrackLikelihood gpLikelihood,
											   double weight, double scaleFactor,CoercionMode mode) {
		super(mode);
		GPOperator=new GaussianProcessSkytrackBlockUpdateOperator();
		this.tree=treeModel;
		this.gpLikelihood = gpLikelihood;
		popSizeParameter = gpLikelihood.getPopSizeParameter();
		changePoints=gpLikelihood.getChangePoints();
		GPtype=gpLikelihood.getGPtype();
		CoalCounts=gpLikelihood.getCoalCounts();
		coalfactor=gpLikelihood.getcoalfactor();
		precisionParameter=gpLikelihood.getPrecisionParameter();
//        precisionParameter = gmrfLikelihood.getPrecisionParameter();
//        lambdaParameter = gmrfLikelihood.getLambdaParameter();

		this.scaleFactor = scaleFactor;
//        lambdaScaleFactor = 0.0;
//        fieldLength = popSizeParameter.getDimension();
//
//        this.maxIterations = maxIterations;
//        this.stopValue = stopValue;
		setWeight(weight);



//        zeros = new double[fieldLength];
	}




	private static void collectAllTimes(Tree tree, NodeRef top, NodeRef[] excludeBelow,
										ArrayList<ComparableDouble> times, ArrayList<Integer> childs) {

		times.add(new ComparableDouble(tree.getNodeHeight(top)));
		childs.add(tree.getChildCount(top));

		for (int i = 0; i < tree.getChildCount(top); i++) {
			NodeRef child = tree.getChild(top, i);
			if (excludeBelow == null) {
				collectAllTimes(tree, child, excludeBelow, times, childs);
			} else {
				// check if this subtree is included in the coalescent density
				boolean include = true;
				for (NodeRef anExcludeBelow : excludeBelow) {
					if (anExcludeBelow.getNumber() == child.getNumber()) {
						include = false;
						break;
					}
				}
				if (include)
					collectAllTimes(tree, child, excludeBelow, times, childs);
			}
		}
	}

//    public double [] getupIntervals() {

//        getTreeIntervals(tree, getMRCAOfCoalescent(tree), getExcludedMRCAs(tree), ti);
	private class Trip1GP{
		private double[] array1;
		private int[] array2;
		private int array3;

		public Trip1GP (double[]  array1, int [] array2, int array3){
			this.array1=array1;
			this.array2=array2;
			this.array3=array3;

		}
		public double [] getInterval() {return array1;}
		public int[] getLineages() {return array2;}
		public int getnIntervals() {return array3;}
	}

	private class TripGP{
		private double[] array1;
		private double[] array2;
		private int array3;

		public TripGP (double[]  array1, double [] array2, int array3){
			this.array1=array1;
			this.array2=array2;
			this.array3=array3;

		}
		public double [] getTimes() {return array1;}
		public double[] getFactor() {return array2;}
		public int getnIntervals() {return array3;}
	}


	public Trip1GP getTreeIntervals(Tree tree, NodeRef root, NodeRef[] exclude){
		int maxIntervalCount = tree.getNodeCount();

		double [] intervals = new double[maxIntervalCount];
		int [] lineageCounts = new int[maxIntervalCount];
		int nIntervals;

		double MULTIFURCATION_LIMIT = 1e-9;

//        System.err.println(intervals.length+"intervals");
		ArrayList<ComparableDouble> times = new ArrayList<ComparableDouble>();
		ArrayList<Integer> childs = new ArrayList<Integer>();
		collectAllTimes(tree, root, exclude, times, childs);
//        private static void collectAllTimes(Tree tree, NodeRef top, NodeRef[] excludeBelow,
//         ArrayList<ComparableDouble> times, ArrayList<Integer> childs) {


		int[] indices = new int[times.size()];

		HeapSort.sort(times, indices);



		// start is the time of the first tip
		double start = times.get(indices[0]).doubleValue();

		int numLines = 0;
		int i = 0;
		int intervalCount = 0;
//        System.err.println("times size"+times.size());


		while (i < times.size()) {


			int lineagesRemoved = 0;
			int lineagesAdded = 0;

			final double finish = times.get(indices[i]).doubleValue();
//            System.err.println("i"+i+" finish"+finish);
			double next = finish;

			while (Math.abs(next - finish) < MULTIFURCATION_LIMIT) {
				final int children = childs.get(indices[i]);
//                System.err.println("children"+children);
				if (children == 0) {
//                    System.err.println("adding");
					lineagesAdded += 1;
				} else {
					lineagesRemoved += (children - 1);
//                    System.err.println("deleting");
				}
				i += 1;
//                System.err.println("now i is"+i+"and time.size"+times.size());
				if (i == times.size()) break;

				next = times.get(indices[i]).doubleValue();
//                System.err.println("next is"+next);
			}
			//System.out.println("time = " + finish + " removed = " + lineagesRemoved + " added = " + lineagesAdded);
			if (lineagesAdded > 0) {
//                System.err.println("lineagesAdded>0 and intervalCount"+intervalCount);

				if (intervalCount > 0 || ((finish - start) > MULTIFURCATION_LIMIT)) {
//                    System.err.println("finish-start"+finish+"and"+start);
//                    double val=finish-start;
//                    System.err.println(val);
//                    System.err.println(intervals.length);
					intervals[intervalCount] = finish - start;
//                    System.err.println("here 1");
					lineageCounts[intervalCount] = numLines;
//                    System.err.println("here 2");
					intervalCount += 1;

				}
//                  System.err.println("start"+start);
				start = finish;
//                System.err.println("start after"+start);
			}
			// add sample event
			numLines += lineagesAdded;

//            System.err.println("numLines is"+numLines);
			if (lineagesRemoved > 0) {
//                System.err.println("lineages Removed>0 and int count is"+intervalCount);
//                System.err.println("finish-start"+finish+"start"+start);

				intervals[intervalCount] = finish - start;
//                System.err.println('0');
				lineageCounts[intervalCount] = numLines;
				intervalCount += 1;
//                System.err.println("lineages Removed>0 and int count is"+intervalCount);
				start = finish;
//                System.err.println("start"+start);
			}
			// coalescent event
			numLines -= lineagesRemoved;
		}


		nIntervals = intervalCount;
		return new Trip1GP(intervals,lineageCounts,nIntervals);
	}

    private double[] getoldCoalTimes(int ncoal){
//        System.err.println("number coal"+ncoal);
        double [] coaltimes= new double[ncoal];
        int k=0;
        for (int j=0; j<changePoints.getSize();j++){
            if (GPtype.getParameterValue(j)==1) {
                coaltimes[k]=changePoints.getParameterValue(j);
//                System.err.println("coaltime"+k+" is"+coaltimes[k]);
                k++;
            }

        }
        return coaltimes;
    }

    private int wherePoint(double newval,double [] intervals){
        double prevle=0;
        double length=0;
        int res=0;
        for (int j=0;j<intervals.length;j++){
            length+=intervals[j];
//            System.err.println(j+" j "+length);
            if (newval<length & newval>prevle){
                res=j;
                break;}
        }
        return res;
    }

	private TripGP getNewCoalTimes(double[] intervals,int nIntervals, int[] lineageCounts){
		double [] NewCoalTimes=new double[nIntervals];
		double [] NewCoalFactor=new double[nIntervals];
		int k=0;
		int res=0;
		double length=0.0;
		for (int i=0;i<nIntervals;i++){
		 length+=intervals[i];
		 res=getEventType(i,nIntervals,lineageCounts);
		 if (res > 0) {
			 NewCoalTimes[k]=length;
			 NewCoalFactor[k]=lineageCounts[i]*(lineageCounts[i]-1.0) / 2.0;
			 k++;
		 }
		}
		return new TripGP(NewCoalTimes,NewCoalFactor,k);
	}

	private final int getEventType(int i,int nIntervals,int[] lineageCounts) {

		if (i >= nIntervals) throw new IllegalArgumentException();
		if (i < nIntervals - 1) {
			return lineageCounts[i] - lineageCounts[i + 1];
		} else {
			return lineageCounts[i] - 1;
		}
	}

//        if (numEvents > 0) return CoalescentEventType.COALESCENT;
//        else if (numEvents < 0) return CoalescentEventType.NEW_SAMPLE;
//        else return CoalescentEventType.NOTHING;


	public double doOperation() {

//		System.err.println("Does GPSkytrack Tree Operator"+changePoints.getSize());
        double hRatio = 0;
		Trip1GP intervalInfo= getTreeIntervals(tree,tree.getRoot(),null);
		double [] intervals =intervalInfo.getInterval();
		int [] lineageCount=intervalInfo.getLineages();
		int nIntervals=intervalInfo.getnIntervals();
        boolean rootIssue=false;
        int k=0;
		int count=0;
        double newval;
		int score=0;
        int oldnum=changePoints.getSize();
		TripGP newCoalTimes=getNewCoalTimes(intervals,nIntervals,lineageCount);

		double [] newTimes=newCoalTimes.getTimes();
      	double [] factor=newCoalTimes.getFactor();
        double [] oldTimes=getoldCoalTimes(newCoalTimes.getnIntervals());
//        for (int j=0;j<newCoalTimes.getNumber();j++){
//            System.err.println("newCoaltime"+j+" is"+newTimes[j]);
//        }

//        int [] changeInPoints=new int[newCoalTimes.getnIntervals()];
//        int [] changeInCoal=new int[newCoalTimes.getnIntervals()];


//		System.err.println(changePoints.getSize()+"and k"+newCoalTimes.getnIntervals());
        System.err.println("Makes GP proposal for coalescent times");

        if (changePoints.getParameterValue(changePoints.getSize()-1)!=newTimes[newCoalTimes.getnIntervals()-1]){
            rootIssue=true;
        } else{
            rootIssue=false;
        }
//        System.err.println("rootissue"+rootIssue);
//		if (changePoints.getParameterValue(changePoints.getSize()-1)==newTimes[newCoalTimes.getnIntervals()-1]){
//			System.err.println("It doesn't change the root");

        if (rootIssue==false){
            int mycase=0;
//            System.err.println("does not change the root");
//            System.err.println("before size"+changePoints.getSize()+"last val"+changePoints.getParameterValue(changePoints.getSize()-1));

            for (int j=0;j<changePoints.getSize();j++){
            if (GPtype.getParameterValue(j)==1){
                mycase=0;
//	    		System.err.println("new coal:"+newTimes[k]+" and coal old"+changePoints.getParameterValue(j)+"and k"+k+" and j"+j);
		    	if (changePoints.getParameterValue(j)!=newTimes[k]){
                	if (changePoints.getParameterValue(j)<newTimes[k] & k<= (newCoalTimes.getnIntervals()-1)){
                        mycase=1;
					}
//				    if (changePoints.getParameterValue(j)<newTimes[k] & k==(newCoalTimes.getnIntervals()-1)){
//				        mycase=3;
//				    }
                    if (changePoints.getParameterValue(j)>newTimes[k]){
                        mycase=2;
					}
			    }
//                System.err.println("mycase"+mycase);
			    if (mycase==1){
//                    System.err.println("move up");
                    double [] currentChangePoints = (double []) changePoints.getParameterValues();
                    DenseVector currentGPvalues= new DenseVector(popSizeParameter.getParameterValues());
                    double currentPrecision = precisionParameter.getParameterValue(0);

                    hRatio+=GPOperator.getDensity(currentChangePoints,currentGPvalues,changePoints.getParameterValue(j),currentPrecision,j);
                    changePoints.removeDimension(j+1);
                    GPtype.removeDimension(j+1);
                    popSizeParameter.removeDimension(j+1);
                    coalfactor.removeDimension(j+1);
                    j--;
                    k--;
                }
                if (mycase==2){
//                    System.err.println("move down");
                    double [] currentChangePoints = (double []) changePoints.getParameterValues();
                    DenseVector currentGPvalues= new DenseVector(popSizeParameter.getParameterValues());
                    double currentPrecision = precisionParameter.getParameterValue(0);
                    double[] val=GPOperator.getGPvalue(currentChangePoints,currentGPvalues,newTimes[k],currentPrecision);
//                    System.err.println("j is"+j+" and size is"+oldnum);
//                    THIS SHOULD NOT HAPPEN
//                    if (j==oldnum-1){
//                        hRatio+=GPOperator.getDensity(currentChangePoints,currentGPvalues,changePoints.getParameterValue(j),currentPrecision,j);
//                        if (j!=val[2]) {System.err.println("WARNING"+val[2]+"and j:"+j);}
////                        System.err.println("root was replaced instead of moved down -WRONG need to add MH value");
//                        popSizeParameter.setParameterValue((int) val[2],val[0]);
//                        changePoints.setParameterValue((int) val[2],newTimes[k]);
//                        GPtype.setParameterValue((int) val[2],1);
//                        coalfactor.setParameterValue(j,factor[k]);
//                        j--;
//                    } else{
////                        System.err.println("add in pos"+val[2]);
                      popSizeParameter.addDimension((int) val[2],val[0]);
                      changePoints.addDimension((int) val[2],newTimes[k]);
                      GPtype.addDimension((int) val[2],1);
                      coalfactor.addDimension((int) val[2],factor[k]);
                      j=(int) val[2];
                      hRatio-=val[1];

                }
//                if (mycase==3){
//                    System.err.println("move up");
//                    double [] currentChangePoints = (double []) changePoints.getParameterValues();
//                    DenseVector currentGPvalues= new DenseVector(popSizeParameter.getParameterValues());
//                    double currentPrecision = precisionParameter.getParameterValue(0);
//
//                    hRatio+=GPOperator.getDensity(currentChangePoints,currentGPvalues,changePoints.getParameterValue(j),currentPrecision,j);
//                    changePoints.removeDimension(j+1);
//                    GPtype.removeDimension(j+1);
//                    popSizeParameter.removeDimension(j+1);
//                    coalfactor.removeDimension(j+1);
//                    j--;
//                    k--;
//                }
            k++;
                if (k==newCoalTimes.getnIntervals()){
                    break;}
	    	}
		}
        }

        if (rootIssue==true){
//            System.err.println("old root"+changePoints.getParameterValue(changePoints.getSize()-1));
//            System.err.println("new root"+changePoints.getParameterValue(changePoints.getSize()-1));

//            System.err.println("the root changes - re-arrangement of latent points here "+newCoalTimes.getnIntervals());
            int k2=0;
            int pos=0;
            int last=0;
            int count2=0;
            double oldlow=0.0;
            double oldupp=oldTimes[k2];
            double newlow=0.0;
            double newupp=newTimes[k2];
            int mycase=0;
//            int temp=changePoints.getSize();
              System.err.println("size before"+changePoints.getSize());
            for (int j=0; j<changePoints.getSize();j++){
                mycase=0;
//                System.err.println(j+"point"+changePoints.getParameterValue(j)+" and newupp"+newupp+"and k"+k2+"type"+GPtype.getParameterValue(j)+"oldupp"+oldupp+" oldlow"+oldlow);
                if (changePoints.getParameterValue(j)<newupp & j<(changePoints.getSize()-1)){
                    if (GPtype.getParameterValue(j)==-1){
                        if (changePoints.getParameterValue(j)>Math.max(oldlow,newlow) & changePoints.getParameterValue(j)<Math.min(oldupp,newupp)){
                            count2++;
//                            System.err.println("count"+count2);
                        } else {
//                            System.err.println("deleting the point"+changePoints.getParameterValue(j));
                            mycase=1;
                        }
                    }else {
//                        System.err.println(j+"deleting a coal point"+changePoints.getParameterValue(j));
                        mycase=1;

                    }

                }
                if (changePoints.getParameterValue(j)<newupp & j==(changePoints.getSize()-1)){
                    mycase=2;

                }

                if (changePoints.getParameterValue(j)>=newupp){
//                    System.err.println("End of interval"+k2+" count was"+count2+" and should be"+CoalCounts.getParameterValue(k2));
                    if (count2<CoalCounts.getParameterValue(k2)){
                        mycase=3;
//                        System.err.println("adds points");


                    } else {
                        if (changePoints.getParameterValue(j)>newupp & k2<newCoalTimes.getnIntervals()-1){
                            mycase=4;
                        }
                        if (changePoints.getParameterValue(j)==newupp & k2<newCoalTimes.getnIntervals()-1){
                            mycase=5;
                        }
                        if (changePoints.getParameterValue(j)==newupp & k2==newCoalTimes.getnIntervals()-1){
                            if (j<changePoints.getSize()-1) {
                            mycase=7;
                            }
                            if (j==(changePoints.getSize()-1)){
                            break;
                            }
                        }
                        if (changePoints.getParameterValue(j)>newupp & k2==newCoalTimes.getnIntervals()-1){
                            mycase=6;
                        }

                }
                }
//                System.err.println(j+"my case"+mycase+" and count"+count2);
                if (mycase==1){
                    double [] currentChangePoints = (double []) changePoints.getParameterValues();
                    DenseVector currentGPvalues= new DenseVector(popSizeParameter.getParameterValues());
                    double currentPrecision = precisionParameter.getParameterValue(0);

                    hRatio+=GPOperator.getDensity(currentChangePoints,currentGPvalues,changePoints.getParameterValue(j),currentPrecision,j);
                    hRatio-=Math.log(oldupp-oldlow);
                    changePoints.removeDimension(j+1);
                    GPtype.removeDimension(j+1);
                    popSizeParameter.removeDimension(j+1);
                    coalfactor.removeDimension(j+1);
                    j--;
                }
                if (mycase==2){
                    double [] currentChangePoints = (double []) changePoints.getParameterValues();
                    DenseVector currentGPvalues= new DenseVector(popSizeParameter.getParameterValues());
                    double currentPrecision = precisionParameter.getParameterValue(0);

                    hRatio+=GPOperator.getDensity(currentChangePoints,currentGPvalues,changePoints.getParameterValue(j),currentPrecision,j);

                    double[] val=GPOperator.getGPvalueroot(currentChangePoints,currentGPvalues,newupp,currentPrecision);
//                          System.err.println(val.length+" "+val[0]+" "+val[1]);
//                        System.err.println("root was replaced instead of moved down");
                    popSizeParameter.setParameterValue(j,val[0]);
                    changePoints.setParameterValue(j,newupp);
//                            GPtype.setParameterValue(j,1);
                    coalfactor.setParameterValue(j,factor[k2-1]);
                    hRatio-=val[1];
                    j--;
                }
                if (mycase==3){
                    int should=(int) CoalCounts.getParameterValue(k2)-count2;
                    j=j+(int) CoalCounts.getParameterValue(k2)-count2-1;
                    for (int i=0;i<should;i++){
                        newval=MathUtils.uniform(newlow,newupp);
                        pos=wherePoint(newval,intervals);
//                            System.err.println("pos"+pos);
                        double [] currentChangePoints = (double []) changePoints.getParameterValues();
                        DenseVector currentGPvalues= new DenseVector(popSizeParameter.getParameterValues());
                        double currentPrecision = precisionParameter.getParameterValue(0);
                        double[] val=GPOperator.getGPvalue(currentChangePoints,currentGPvalues,newval,currentPrecision);
//                            System.err.println("pass"+val[2]);
                        popSizeParameter.addDimension((int) val[2],val[0]);
                        changePoints.addDimension((int) val[2],newval);
                        GPtype.addDimension((int) val[2],-1);
                        coalfactor.addDimension((int) val[2],0.5*lineageCount[pos]*(lineageCount[pos]-1));
                        hRatio-=val[1];
                        hRatio+=Math.log(newupp-newlow);
                        count2++;
                    }

                }
                if (mycase==4){
                    double [] currentChangePoints = (double []) changePoints.getParameterValues();
                    DenseVector currentGPvalues= new DenseVector(popSizeParameter.getParameterValues());
                    double currentPrecision = precisionParameter.getParameterValue(0);

                    double[] val=GPOperator.getGPvalue(currentChangePoints,currentGPvalues,newupp,currentPrecision);

                    popSizeParameter.addDimension((int) val[2],val[0]);
                    changePoints.addDimension((int) val[2],newupp);
                    GPtype.addDimension((int) val[2],1);
                    coalfactor.addDimension((int) val[2],factor[k2]);
                    hRatio-=val[1];

                    count2=0;
                    k2++;
                    last=pos;
                    newlow=newupp;
                    oldlow=oldupp;
                    newupp=newTimes[k2];
                    oldupp=oldTimes[k2];
                }
                if (mycase==5){
                    count2=0;
                    k2++;
                    last=pos;
                    newlow=newupp;
                    oldlow=oldupp;
                    newupp=newTimes[k2];
                    oldupp=oldTimes[k2];


                    double [] currentChangePoints = (double []) changePoints.getParameterValues();
                    DenseVector currentGPvalues= new DenseVector(popSizeParameter.getParameterValues());
                    double currentPrecision = precisionParameter.getParameterValue(0);
                    if (j==changePoints.getSize()-1){
                    double[] val=GPOperator.getGPvalueroot(currentChangePoints,currentGPvalues,newupp,currentPrecision);
                        popSizeParameter.addDimension(j+1,val[0]);
                        changePoints.addDimension(j+1,newupp);
                        GPtype.addDimension(j+1,1);
                        coalfactor.addDimension(j+1,factor[k2]);
                        hRatio-=val[1];
                     }else{
//                        double[] val=GPOperator.getGPvalue(currentChangePoints,currentGPvalues,newupp,currentPrecision);
//                        popSizeParameter.addDimension(j+1,val[0]);
//                        changePoints.addDimension(j+1,newupp);
//                        GPtype.addDimension(j+1,1);
//                        coalfactor.addDimension(j+1,factor[k2]);
//                        hRatio-=val[1];
                    }
                }
                if (mycase==6){
                    double [] currentChangePoints = (double []) changePoints.getParameterValues();
                    DenseVector currentGPvalues= new DenseVector(popSizeParameter.getParameterValues());
                    double currentPrecision = precisionParameter.getParameterValue(0);
                    double[] val=GPOperator.getGPvalue(currentChangePoints,currentGPvalues,newupp,currentPrecision);
                    popSizeParameter.addDimension(j,val[0]);
                    changePoints.addDimension(j,newupp);
                    GPtype.addDimension(j,1);
                    coalfactor.addDimension(j,factor[k2]);
                    hRatio-=val[1];
                    j--;
                }
                if (mycase==7){
                    double [] currentChangePoints = (double []) changePoints.getParameterValues();
                    DenseVector currentGPvalues= new DenseVector(popSizeParameter.getParameterValues());
                    double currentPrecision = precisionParameter.getParameterValue(0);

                    hRatio+=GPOperator.getDensity(currentChangePoints,currentGPvalues,changePoints.getParameterValue(j+1),currentPrecision,j+1);
                    hRatio-=Math.log(oldupp-oldlow);
                    if (j<changePoints.getSize()-2){
                    changePoints.removeDimension(j+2);
                    GPtype.removeDimension(j+2);
                    popSizeParameter.removeDimension(j+2);
                    coalfactor.removeDimension(j+2);
                    }else{
                        changePoints.setDimension(j+1);
                        GPtype.setDimension(j+1);
                        popSizeParameter.setDimension(j+1);
                        coalfactor.setDimension(j+1);
                    }

                    j--;
                }
            }
        }

        System.err.println("size after"+changePoints.getSize()+"and ratio"+hRatio);
		return hRatio;
//        return 0.0;
	}

	//MCMCOperator INTERFACE

	public final String getOperatorName() {
		return GaussianProcessSkytrackTreeOperatorParser.BLOCK_UPDATE_OPERATOR;
	}

	public double getCoercableParameter() {
//        return Math.log(scaleFactor);
		return Math.sqrt(scaleFactor - 1);
	}

	public void setCoercableParameter(double value) {
//        scaleFactor = Math.exp(value);
		scaleFactor = 1 + value * value;
	}

	public double getRawParameter() {
		return scaleFactor;
	}

	public double getScaleFactor() {
		return scaleFactor;
	}

	public double getTargetAcceptanceProbability() {
		return 0.234;
	}

	public double getMinimumAcceptanceLevel() {
		return 0.1;
	}

	public double getMaximumAcceptanceLevel() {
		return 0.4;
	}

	public double getMinimumGoodAcceptanceLevel() {
		return 0.20;
	}

	public double getMaximumGoodAcceptanceLevel() {
		return 0.30;
	}

	public final String getPerformanceSuggestion() {

		double prob = Utils.getAcceptanceProbability(this);
		double targetProb = getTargetAcceptanceProbability();
		dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);

		double sf = OperatorUtils.optimizeWindowSize(scaleFactor, prob, targetProb);

		if (prob < getMinimumGoodAcceptanceLevel()) {
			return "Try setting scaleFactor to about " + formatter.format(sf);
		} else if (prob > getMaximumGoodAcceptanceLevel()) {
			return "Try setting scaleFactor to about " + formatter.format(sf);
		} else return "";
	}


//  public DenseVector oldNewtonRaphson(double[] data, DenseVector currentGamma, SymmTridiagMatrix proposedQ) throws OperatorFailedException{
//  return newNewtonRaphson(data, currentGamma, proposedQ, maxIterations, stopValue);
//
//}
//
//public static DenseVector newtonRaphson(double[] data, DenseVector currentGamma, SymmTridiagMatrix proposedQ,
//int maxIterations, double stopValue) {
//
//DenseVector iterateGamma = currentGamma.copy();
//int numberIterations = 0;
//while (gradient(data, iterateGamma, proposedQ).norm(Vector.Norm.Two) > stopValue) {
//inverseJacobian(data, iterateGamma, proposedQ).multAdd(gradient(data, iterateGamma, proposedQ), iterateGamma);
//numberIterations++;
//}
//
//if (numberIterations > maxIterations)
//throw new RuntimeException("Newton Raphson algorithm did not converge within " + maxIterations + " step to a norm less than " + stopValue);
//
//return iterateGamma;
//}
//
//private static DenseMatrix inverseJacobian(double[] data, DenseVector value, SymmTridiagMatrix Q) {
//
//      SPDTridiagMatrix jacobian = new SPDTridiagMatrix(Q, true);
//      for (int i = 0; i < value.size(); i++) {
//          jacobian.set(i, i, jacobian.get(i, i) + Math.exp(-value.get(i)) * data[i]);
//      }
//
//      DenseMatrix inverseJacobian = Matrices.identity(jacobian.numRows());
//      jacobian.solve(Matrices.identity(value.size()), inverseJacobian);
//
//      return inverseJacobian;
//  }
}
