/*
 * NPAntigenicLikelihood.java
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

package dr.evomodel.antigenic;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.GammaFunction;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Gabriela Cybis
 * @author Marc Suchard
 */
public class NPAntigenicLikelihood extends AbstractModelLikelihood {
	
	  public static final String NP_ANTIGENIC_LIKELIHOOD = "NPAntigenicLikelihood";

	
	
	public NPAntigenicLikelihood (TreeModel treeModel,CompoundParameter traitParameter , Parameter assignments, Parameter links,
                                  Parameter chi, Parameter clusterPrec, Parameter priorMean, Parameter priorPrec,
                                  Parameter transformFactor,Parameter means1, Parameter means2, Parameter locationDrift,
                                  Parameter offsets, Boolean hasDrift){
		  super(NP_ANTIGENIC_LIKELIHOOD);


          this.assignments = assignments;
          this.links = links;
		  this.clusterPrec = clusterPrec;
		  this.priorPrec = priorPrec;
		  this.priorMean = priorMean;
	      this.treeModel= treeModel;
		  this.traitParameter= traitParameter;
		  this.transformFactor=transformFactor;
		  this.means1=means1;
          this.means2=means2;
          this.alpha= chi;
          this.locationDrift=locationDrift;
          this.offsets=offsets;
        //  this.hasDrift=hasDrift;
          this.hasDrift=false;



        addVariable(traitParameter);
        addVariable(assignments);
        addVariable(links);
        addModel(treeModel);
        addVariable(chi);
        addVariable(transformFactor);
        addVariable(alpha);
       // addVariable(locationDrift);
        addVariable(offsets);

        numdata = traitParameter.getParameterCount();
        this.allTips= TreeUtils.getExternalNodes(treeModel,treeModel.getRoot());



       setDepMatrix();





	 for (int i=0; i<numdata; i++){
		 assignments.setParameterValue(i, i);
	     links.setParameterValue(i,i);

     }

	
	       
	       this.logLikelihoodsVector = new double[links.getDimension()+1];
           this.logLikelihoodsVectorKnown = new boolean[links.getDimension()+1];
           this.storedLogLikelihoodsVector = new double[links.getDimension()+1];



        this.m = new double[2];
        m[0]= priorMean.getParameterValue(0);
        m[1]= priorMean.getParameterValue(1);


         this.v0 = 2;
       //  double v1 = 3;

        this.k0= priorPrec.getParameterValue(0)/clusterPrec.getParameterValue(0);
       // double k1= k0+1;


        this.T0Inv= new double[2][2];
        T0Inv[0][0]= v0/clusterPrec.getParameterValue(0);
        T0Inv[1][1]= v0/clusterPrec.getParameterValue(0);
        T0Inv[1][0]= 0.0;
        T0Inv[0][1]= 0.0;


          this.logDetT0= -Math.log(T0Inv[0][0]*T0Inv[1][1]);


      /*  for(int i=0;i<logLikelihoodsVector.length-1;i++){

         double[][] T1Inv = new double[2][2];
            T1Inv[0][0]=T0Inv[0][0]+(k0/k1)* data[i][0]*data[i][0];
            T1Inv[0][1]=T0Inv[0][1]+(k0/k1)* data[i][0]*data[i][1];
            T1Inv[1][0]=T0Inv[1][0]+(k0/k1)* data[i][1]*data[i][0];
            T1Inv[1][1]=T0Inv[1][1]+(k0/k1)* data[i][1]*data[i][1];


            double logDetT1=-Math.log(T1Inv[0][0]*T1Inv[1][1]-T1Inv[0][1]*T1Inv[1][0]);


            logLikelihoodsVector[i]= -(1*2/2)*Math.log(Math.PI);
            logLikelihoodsVector[i]+= Math.log(k0) - Math.log(k1);
            logLikelihoodsVector[i]+= (v1/2)*logDetT1 - (v0/2)*logDetT0;
            logLikelihoodsVector[i]+= GammaFunction.lnGamma(v1/2)+ GammaFunction.lnGamma((v1/2)-0.5);
            logLikelihoodsVector[i]+=-GammaFunction.lnGamma(v0/2)- GammaFunction.lnGamma((v0/2)-0.5);


        }
    */




    }


    /*
    public void setInitialAssignmentsToDates(){







        double[] offsetValues = new double[numdata];
        for (int i=0; i<numdata; i++){
            offsetValues[i] =offsets.getParameterValue(findOffsetIndex(i));

        }


        boolean[] assigned = new boolean[numdata];
        for (int i=0; i<numdata;i++ ){
            assigned[i]=false;
        }

            int group = 0;
            for (int i=0; i<numdata; i++){
               if (!assigned[i]){

                   int last = i;
                  for (int j=0;j<numdata;j++){

                      if (offsetValues[j]==offsetValues[i]||offsetValues[j]==offsetValues[i]+1||offsetValues[j]==offsetValues[i]+2||offsetValues[j]==offsetValues[i]+3||offsetValues[j]==offsetValues[i]+4){
                          links.setParameterValue(j,last);
                          assignments.setParameterValue(j, group);
                          assigned[j]=true;
                          last=j;
                      }

                  }
                   links.setParameterValue(i,last);
                   group ++;

               }



        }


        printInformation(links);
        printInformation(assignments);



    }


    private int findOffsetIndex(int traitParameterIndex){
               String NAME = traitParameter.getParameter(traitParameterIndex).getParameterName();
               boolean notFound =true;
               int i=0;


        while (notFound){

            if(offsets.getDimensionName(i).compareTo(NAME)==0){

                    notFound=false;
                }else{
                    i++;
            }
                }

      //  printInformation((double) i);
      //  printInformation(offsets.getDimensionName(i),NAME);

       return i;
    }




   private void setData(){
        dataMatrixKnown=true;
        int dim = traitParameter.getParameter(0).getSize();


        for (int i=0; i<numdata; i++){
            for (int j=1; j<dim; j++){
                data[i][j]= traitParameter.getParameter(i).getParameterValue(j);
            }
            //if (hasDrift){
               // int offsetIndex = findOffsetIndex(i);
              //  Data[i][0] += locationDrift.getParameterValue(0)*offsets.getParameterValue(offsetIndex);
            //}  else{
                data[i][0]= traitParameter.getParameter(i).getParameterValue(0);
            //}
        }



    }


    private void setDatum(int virus){


        int dim = traitParameter.getParameter(0).getSize();


        for (int j=1; j<dim; j++){
                data[virus][j]= traitParameter.getParameter(virus).getParameterValue(j);
            }
            //if (hasDrift){
            // int offsetIndex = findOffsetIndex(i);
            //  Data[i][0] += locationDrift.getParameterValue(0)*offsets.getParameterValue(offsetIndex);
            //}  else{
            data[virus][0]= traitParameter.getParameter(virus).getParameterValue(0);
            //}




    }
   */


    private void setDepMatrix(){

        depMatrixKnown = true;

        depMatrix=new double[numdata][numdata];
        List<NodeRef> childList = new ArrayList<NodeRef>();

        recursion(treeModel.getRoot(),childList);
        logCorrectMatrix(transformFactor.getParameterValue(0));
        logDepMatrix =  new double[numdata][numdata];
        for(int i=0;i<numdata;i++){
            for(int j=0;j<i;j++){
                logDepMatrix[i][j]=Math.log(depMatrix[i][j]);
                logDepMatrix[j][i]=logDepMatrix[j][i];

            }
        }

    }





    public double getLogLikGroup(int groupNumber){
        double L =0.0;


        int ngroup=0;
        for (int i=0;i<assignments.getDimension(); i++){
            if((int) assignments.getParameterValue(i) == groupNumber){
                ngroup++;}}


        if (ngroup != 0){
            double[][] group = new double[ngroup][2];


            double mean[]=new double[2];

            int count = 0;
            for (int i=0;i<assignments.getDimension(); i++){
                if((int) assignments.getParameterValue(i) == groupNumber){
                    group[count][0] = getData(i,0);
                    group[count][1] = getData(i,0);
                    mean[0]+=group[count][0];
                    mean[1]+=group[count][1];
                    count++;}}

            mean[0]/=ngroup;
            mean[1]/=ngroup;



            double kn= k0+ngroup;
            double vn= v0+ngroup;


            double[][] sumdif=new double[2][2];

            for(int i=0;i<ngroup;i++){
                sumdif[0][0]+= (group[i][0]-mean[0])*(group[i][0]-mean[0]);
                sumdif[0][1]+= (group[i][0]-mean[0])*(group[i][1]-mean[1]);
                sumdif[1][0]+= (group[i][0]-mean[0])*(group[i][1]-mean[1]);
                sumdif[1][1]+= (group[i][1]-mean[1])*(group[i][1]-mean[1]);
            }



            double[][] TnInv = new double[2][2];
            TnInv[0][0]=T0Inv[0][0]+ngroup*(k0/kn)*(mean[0]-m[0])*(mean[0]-m[0])+sumdif[0][0];
            TnInv[0][1]=T0Inv[0][1]+ngroup*(k0/kn)*(mean[1]-m[1])*(mean[0]-m[0])+sumdif[0][1];
            TnInv[1][0]=T0Inv[1][0]+ngroup*(k0/kn)* (mean[0]-m[0])*(mean[1]-m[1])+sumdif[1][0];
            TnInv[1][1]=T0Inv[1][1]+ngroup*(k0/kn)* (mean[1]-m[1])*(mean[1]-m[1])+sumdif[1][1];


            double logDetTn=-Math.log(TnInv[0][0]*TnInv[1][1]-TnInv[0][1]*TnInv[1][0]);


            L+= -(ngroup)*Math.log(Math.PI);
            L+= Math.log(k0) - Math.log(kn);
            L+= (vn/2)*logDetTn - (v0/2)*logDetT0;
            L+= GammaFunction.lnGamma(vn/2)+ GammaFunction.lnGamma((vn/2)-0.5);
            L+=-GammaFunction.lnGamma(v0/2)- GammaFunction.lnGamma((v0/2)-0.5);






        }
        logLikelihoodsVectorKnown[groupNumber]=true;
        return L;

    }







    public Model getModel() {
	        return this;
	    }

	  
	  
	  public double[] getLogLikelihoodsVector(){
		  return logLikelihoodsVector;
	  }

    public Parameter getLinks(){
        return links;
    }

    public Parameter getAssignments(){
        return assignments;
    }


    public double getData(int virus, int dim){
		  return traitParameter.getParameter(virus).getParameterValue(dim);
	  }
	  
	  public double[][] getDepMatrix(){
		  return depMatrix;
	  }
	  

	  public double[][] getLogDepMatrix(){
		  return logDepMatrix;
	  }
	  
	  public Parameter getPriorMean(){
		  return priorMean;
	  }
	  public Parameter getPriorPrec(){
		  return priorPrec;
	  }
	  public Parameter getClusterPrec(){
		  return clusterPrec;
	  }
	  
	  public void setLogLikelihoodsVector(int pos, double value){
		  logLikelihoodsVector[pos]=value;
	  }

    public void setAssingments(int pos, double value){
        assignments.setParameterValue(pos,value);
    }

    public void setLinks(int pos, double value){
        links.setParameterValue(pos,value);
    }

    public void setMeans(int pos, double[] value){
        means1.setParameterValue(pos,value[0]);
        means2.setParameterValue(pos,value[1]);
    }


    public double getLogLikelihood() {
        if (!logLikelihoodKnown) {

            logLikelihood = computeLogLikelihood();
        }

        return logLikelihood;
    }



	  
	  public double computeLogLikelihood() {


          if (!depMatrixKnown ){
          setDepMatrix();
          }




          double logL = 0.0;
		  for (int j=0 ; j<logLikelihoodsVector.length;j++){

              if(!logLikelihoodsVectorKnown[j]){
              logLikelihoodsVector[j]=getLogLikGroup(j);
              }
			  	  logL +=logLikelihoodsVector[j];

		  }
		  
		  for (int j=0 ; j<links.getDimension();j++){
		if(links.getParameterValue(j)==j){
			logL += Math.log(alpha.getParameterValue(0));
		}
		else{logL += Math.log(depMatrix[j][(int) links.getParameterValue(j)]);
			
		}
		
		double sumDist=0.0;
		for (int i=0;i<numdata;i++){
			if(i!=j){sumDist += depMatrix[i][j];
			}
			}
		
		  logL-= Math.log(alpha.getParameterValue(0)+sumDist);
		  }
		


          logLikelihoodKnown=true;
         // printInformation(logL);

          return logL;


      }
	 
	  
	  
	  /* Getting matrix from tree*/
	  
	  void recursion( NodeRef node, List childList){
		 
		  
		 List<NodeRef> leftChildTipList = new ArrayList<NodeRef>(); 
		 List<NodeRef> rightChildTipList = new ArrayList<NodeRef>(); 
		  
		 if(!treeModel.isExternal(node)){
			 recursion(treeModel.getChild(node, 0),leftChildTipList);
			 recursion(treeModel.getChild(node, 1),rightChildTipList);
			 
			
			 double lBranch = treeModel.getBranchLength(treeModel.getChild(node, 0));
			 double rBranch = treeModel.getBranchLength(treeModel.getChild(node, 1));
			 

			 Set<NodeRef> notLeftChildList = new HashSet<NodeRef>();
			 notLeftChildList.addAll(allTips);
			 	for (NodeRef i :leftChildTipList){
			 			notLeftChildList.remove(i);			 
}
			 	Set<NodeRef> notRightChildList = new HashSet<NodeRef>();
				 notRightChildList.addAll(allTips);
			 	for (NodeRef i :rightChildTipList){
			 		notRightChildList.remove(i);			 
}
			 
			 for (NodeRef lChild : leftChildTipList){
				  for (NodeRef Child : notLeftChildList){
				depMatrix[Child.getNumber()][lChild.getNumber()] += lBranch;	 
				depMatrix[lChild.getNumber()][Child.getNumber()] += lBranch;	 
				 }
			 }


			 for (NodeRef rChild : rightChildTipList){
					 for (NodeRef Child : notRightChildList){
					depMatrix[Child.getNumber()][rChild.getNumber()] += rBranch;	 
					depMatrix[rChild.getNumber()][Child.getNumber()] += rBranch;	 
					 }
				 }
			 
			 
			 childList.addAll(leftChildTipList);
			 childList.addAll(rightChildTipList);
		 }
		 else{
			 childList.add(node);
		 }
		 
	  }
	  
	  
	  
	  void logCorrectMatrix(double p){
		  for (int i=0; i<numdata; i++){
		    	for (int j=0; j<i; j++){
		    		depMatrix[i][j]=1/Math.pow(depMatrix[i][j],p);	    
		    		depMatrix[j][i]=depMatrix[i][j];	    
				    }}
		    
	  }
	  
	  
	  
	  

	  
	  
	  
	  public double getTreeDist(int i, int j){
		 double dist=0; 
  
		 NodeRef MRCA = findMRCA(i,j);

		 NodeRef Parent = treeModel.getExternalNode(i);
		 while (Parent!=MRCA){
			 dist+=treeModel.getBranchLength(Parent);
			 Parent = treeModel.getParent(Parent);
		 }

		 Parent = treeModel.getExternalNode(j);
		 while (Parent!=MRCA){
			 dist+=treeModel.getBranchLength(Parent);
			 Parent = treeModel.getParent(Parent);
		 }
		  
			 return dist;
	  }
	  
	  
	  private NodeRef findMRCA(int iTip, int jTip) {
	        Set<String> leafNames = new HashSet<String>();
	        leafNames.add(treeModel.getTaxonId(iTip));
	        leafNames.add(treeModel.getTaxonId(jTip));
	        return TreeUtils.getCommonAncestorNode(treeModel, leafNames);
	    }


	  public void printInformtion(double[][] Mat) {
         StringBuffer sb = new StringBuffer("matrix \n");
         for(int i=0;i <numdata; i++){
        	 sb.append(" \n");
        	 for(int j=0; j<numdata; j++){
        		 sb.append(Mat[i][j]+" \t");
        	 }
         }
         
         Logger.getLogger("dr.evomodel").info(sb.toString()); };


      public void printInformation(Parameter Vec) {
        StringBuffer sb = new StringBuffer("Vector \n");

                 for (int i = 0; i<numdata; i++){
                sb.append(Vec.getParameterValue(i)+" \t");
                 }


        Logger.getLogger("dr.evomodel").info(sb.toString()); };




    public void printInformation(int[] Vec) {
        StringBuffer sb = new StringBuffer("Vector \n");
        for(int i=0;i <numdata; i++){

            sb.append(Vec[i]+" \t");

        }

        Logger.getLogger("dr.evomodel").info(sb.toString()); };




    public void printOrder() {
            StringBuffer sb = new StringBuffer("taxa \n");
            for(int i=0;i <numdata; i++){
           	 sb.append(" \n");
           	 
           		 sb.append(treeModel.getTaxonId(i));
           	             }
            
            Logger.getLogger("dr.evomodel").info(sb.toString()); };
           
         
         
         
         
         public void printInformation(double x) {
             StringBuffer sb = new StringBuffer("Info \n");
             		 sb.append(x);
             
             Logger.getLogger("dr.evomodel").info(sb.toString()); };

    public void printInformation(String x) {
        StringBuffer sb = new StringBuffer("Info \n");
        sb.append(x);

        Logger.getLogger("dr.evomodel").info(sb.toString()); };

    public void printInformation(String x, String y) {
        StringBuffer sb = new StringBuffer("Info \n");
        sb.append(x + " and " + y);

        Logger.getLogger("dr.evomodel").info(sb.toString()); };



    @Override
    protected void storeState() {
        System.arraycopy(logLikelihoodsVector, 0, storedLogLikelihoodsVector, 0, logLikelihoodsVector.length);
    }

    @Override
    protected void restoreState() {
        double[] tmp = logLikelihoodsVector;
        logLikelihoodsVector = storedLogLikelihoodsVector;
        storedLogLikelihoodsVector = tmp;

        depMatrixKnown = !proposedChangeDepMatrix;
        proposedChangeDepMatrix =false;

    //    dataMatrixKnown = !proposedChangeDataMatrix;
    //    proposedChangeDataMatrix =false;

        logLikelihoodKnown = false;
    }





    public void makeDirty() {
	    }

	    public void acceptState() {
	        // DO NOTHING
            proposedChangeDepMatrix =false;
            proposedChangeDataMatrix =false;
	    }


	    protected void handleModelChangedEvent(Model model, Object object, int index) {
            if (model == treeModel)
                depMatrixKnown=false;
                logLikelihoodKnown=false;
            return;
        }


	    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
	          logLikelihoodKnown = false;

	           if (variable == transformFactor)   {
                depMatrixKnown=false;
                proposedChangeDepMatrix=true;
	    }

           // if (variable == traitParameter || variable == locationDrift ){
            if (variable == traitParameter){
               // dataMatrixKnown=false;
               // proposedChangeDataMatrix=true;

                int loc= index / 2;
                int changedGroup=(int)assignments.getParameterValue(loc);
                logLikelihoodsVectorKnown[changedGroup]=false;
            }

        }

	Set<NodeRef> allTips;  
	CompoundParameter traitParameter;  
	Parameter alpha;
	Parameter clusterPrec ;
	Parameter priorPrec ;
	Parameter priorMean ;
	Parameter assignments;
    Parameter links;
    Parameter means2;
    Parameter means1;
    Parameter locationDrift;
    Parameter offsets;
    boolean hasDrift;

    private boolean depMatrixKnown= false;
    private boolean[] dataMatrixKnown;
    private boolean logLikelihoodKnown=false;
    private double logLikelihood =0.0;
    private boolean[] logLikelihoodsVectorKnown;


    boolean proposedChangeDepMatrix=false;
    boolean proposedChangeDataMatrix=false;




    TreeModel treeModel;
	String traitName;
	//double[][] data;
	double[][] depMatrix;
	double[][] logDepMatrix;
	double[] logLikelihoodsVector;
    double[] storedLogLikelihoodsVector;
    int numdata;
	Parameter transformFactor;
	double k0;
    double v0;
    double[][] T0Inv;
	double[] m;
    double logDetT0;

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public final static String CLUSTER_PREC = "clusterPrec";
        public final static String PRIOR_PREC = "priorPrec";
        public final static String PRIOR_MEAN = "priorMean";
        public final static String ASSIGNMENTS = "assignments";
        public final static String LINKS = "links";
        public final static String MEANS_1 = "clusterMeans1";
        public final static String MEANS_2 = "clusterMeans2";
        public final static String TRANSFORM_FACTOR = "transformFactor";
        public final static String CHI = "chi";
        public final static String OFFSETS = "offsets";
        public final static String LOCATION_DRIFT = "locationDrift";




        boolean integrate = false;


        public String getParserName() {
            return NP_ANTIGENIC_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        	TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        	//String traitName = (String) xo.getAttribute(TRAIT_NAME);

	        XMLObject cxo = xo.getChild(CLUSTER_PREC);
	        Parameter clusterPrec = (Parameter) cxo.getChild(Parameter.class);

	        cxo = xo.getChild(PRIOR_PREC);
	        Parameter priorPrec = (Parameter) cxo.getChild(Parameter.class);

	        cxo = xo.getChild(PRIOR_MEAN);
	        Parameter priorMean = (Parameter) cxo.getChild(Parameter.class);

	        cxo = xo.getChild(ASSIGNMENTS);
	        Parameter assignments = (Parameter) cxo.getChild(Parameter.class);


            cxo = xo.getChild(LINKS);
            Parameter links = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(MEANS_2);
            Parameter means2 = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(MEANS_1);
            Parameter means1 = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(CHI);
            Parameter chi = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(TRANSFORM_FACTOR);
            Parameter transformFactor = (Parameter) cxo.getChild(Parameter.class);



                cxo=xo.getChild(LOCATION_DRIFT) ;
                Parameter locationDrift= (Parameter) cxo.getChild(Parameter.class);

                cxo=xo.getChild(OFFSETS);
                Parameter offsets =(Parameter) cxo.getChild(Parameter.class);

            boolean hasDrift = false;
            if (offsets.getDimension()>1){
                hasDrift=true;
            }



	        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();
            String traitName = TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;


            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                    utilities.parseTraitsFromTaxonAttributes(xo, traitName, treeModel, integrate);
           // traitName = returnValue.traitName;
            CompoundParameter traitParameter = returnValue.traitParameter;





	        return new NPAntigenicLikelihood(treeModel,traitParameter,  assignments, links, chi,clusterPrec, priorMean,priorPrec,
                    transformFactor, means1,means2,locationDrift,offsets,hasDrift);
	    }

	    //************************************************************************
	    // AbstractXMLObjectParser implementation
	    //************************************************************************

	    public String getParserDescription() {
	        return "conditional likelihood ddCRP";
	    }

	    public Class getReturnType() {
	        return NPAntigenicLikelihood.class;
	    }

	    public XMLSyntaxRule[] getSyntaxRules() {
	        return rules;
	    }

	    private final XMLSyntaxRule[] rules = {
	    		 new StringAttributeRule(TreeTraitParserUtilities.TRAIT_NAME, "The name of the trait for which a likelihood should be calculated"),

	    		 new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
	                        new ElementRule(Parameter.class)
	                }),
	    		 new ElementRule(PRIOR_PREC,
		    	                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
		           new ElementRule(CLUSTER_PREC,
			                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
   			       new ElementRule(PRIOR_MEAN,
		    	                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(ASSIGNMENTS,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(LINKS,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(TRANSFORM_FACTOR,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),

                new ElementRule(MEANS_1,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(MEANS_2,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(CHI, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class),
                }),
                new ElementRule(OFFSETS, Parameter.class),
                new ElementRule(LOCATION_DRIFT, Parameter.class),

                new ElementRule(TreeModel.class),

	    };
    };


























}
