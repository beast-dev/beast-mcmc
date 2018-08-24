package dr.evomodel.antigenic.phyloclustering.misc.obsolete;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import dr.evolution.tree.NodeRef;
//import dr.evomodel.antigenic.driver.OrderDouble;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Charles Cheung
 * @author Trevor Bedford
 */

//Some suggestion to speed up the code from Charles Cheung 
//(please scroll through to places marked by to see the places of changes   //---------- suggestion from cykc----------------

		
public class ClusterViruses extends AbstractModelLikelihood {
	
	//---------- suggestion from cykc----------------
	private double mostRecentTransformedValue = 0;  //keep a copy of the most recent version of transformFactor, to keep track of whether the transformFactor has changed
	//private boolean ShouldUpdateDepMatrix = true;  // this way of flagging for change is not used anymore
	private boolean treeChanged = false; //a flag that becomes true when treeModel changes
	//---------- End of suggestion from cykc----------------		
        
    public static final String CLUSTER_VIRUSES = "ClusterViruses";
    
    
    
    
    //==============================================================================================================================
    //variables
    
    
    double lambda = 10;
	double sigmaSq = 9;
   // double sigmaSq = 100; //when offset is off

    
    //K - number of parameters
    Parameter K;   // for now, there is no move to change K 
    
    //E|K - excision points
    Parameter excisionPoints;
    
    //C
    Parameter clusterLabels;
    
    //mu - means
    MatrixParameter mu;
    
    
    double[] muLabels;
    
    
   MatrixParameter virusLocations = null;

    
    
    
   // Parameter virusOffsetsParameter; //need to read it from AntigenicLIkelihood   (the year ) // use offsets instead
    
    //need to stop the virus from 
    
    
    
    
    
        
 public ClusterViruses (TreeModel treeModel_in, 
		 				Parameter K_in, 
		 				Parameter excisionPoints_in, 
		 				Parameter clusterLabels_in, 
		 				MatrixParameter mu_in, 
		 				Boolean hasDrift, 
		 		//		Parameter locationDrift_in, 
		 				Parameter offsets_in,
		 				MatrixParameter virusLocations_in){
	 
		super(CLUSTER_VIRUSES);
		
		this.treeModel= treeModel_in;
		this.K = K_in;
		this.excisionPoints = excisionPoints_in;
		this.clusterLabels = clusterLabels_in;
		this.mu = mu_in;
		
		  this.hasDrift=hasDrift;
  //      this.locationDrift=locationDrift_in;
        this.offsets=offsets_in;
		//this.hasDrift=false;
        this.virusLocations = virusLocations_in;
        numdata = offsets.getSize();
        
        System.out.println("numdata = " + numdata);
        
        //initialize clusterLabels
        clusterLabels.setDimension(numdata);
         for (int i = 0; i < numdata; i++) {      	
            clusterLabels.setParameterValue(i, 0);
        }
         addVariable(clusterLabels);
         
        
         
         //initialize mu
         mu.setColumnDimension(2);
         //mu.setRowDimension(numdata);  //in reality, only K of them are used
         int K_int =  (int) K.getParameterValue(0);
         mu.setRowDimension(K_int);  //in reality, only K of them are used
         //System.out.println((int) K.getParameterValue(0));
         //System.exit(0);
         
         for(int i=0; i < K_int; i++){
        	 //can I set initial condition to be like this? and then let the algorithm set it properly later?
        	 double zero=0;
        	 mu.getParameter(i).setValue(0, zero);
        	 mu.getParameter(i).setValue(1, zero);
         }
         
         //adding the pre-clustering step.
         preClustering();
         
        /*
        this.abc.setColumnDimension(2);  //set dimension equal to 2
        abc.setRowDimension(strains.size());
        for (int i = 0; i < strains.size(); i++) {
            abc.getParameter(i).setId(strains.get(i));
        }
        */
        
         
		addVariable(virusLocations);
		addModel(treeModel);
		//addVariable(locationDrift);
		addVariable(offsets);
		 addVariable(K);
		 addVariable(excisionPoints);
		 addVariable(mu);
		 System.out.println("Finished loading the constructor for ClusterViruses");

 }    
    
    
 
private void preClustering() {

	int numViruses = offsets.getSize();
	System.out.println("# offsets = " + offsets.getSize());
	//for(int i=0; i < offsets.getSize(); i++){
		//System.out.println(offsets.getParameterValue(i));
	//}
	
	
	List<OrderDouble> list = new ArrayList<OrderDouble>();
	for(int i=0; i < numViruses; i++){ 
		 list.add(new OrderDouble(i,  offsets.getParameterValue(i)));  //offset of 1
	}
	Collections.sort(list, new OrderDouble());	

	int initialEqualBinSize = numViruses/(int)(K.getParameterValue(0) -1);
	System.out.println("initial bin size = " + initialEqualBinSize);
	System.out.println("Initial cluster assignment:");
//	System.out.println("virus index\tOffset\tCluster label");
	for(int i=0; i < numViruses; i++){
//		  System.out.println(list.get(i).getIndex() + "\t" + list.get(i).getValue()  +"\t"+ i/ initialEqualBinSize   );
		  
		  int label = i/initialEqualBinSize;
		  clusterLabels.setParameterValue(list.get(i).getIndex() , label);
	}
	
	
	
	
	/*
	//borrow from getLogLikelihood:
	int K_int = (int) K.getParameterValue(0);

	double[] meanYear = new double[K_int];
	double[] groupCount = new double[K_int];
	for(int i=0; i < numdata; i++){
		int label = (int) clusterLabels.getParameterValue(i);
		double year  = 0;
        if (offsets != null) {
            //	System.out.print("virus Offeset Parameter present"+ ": ");
            //	System.out.print( virusOffsetsParameter.getParameterValue(i) + " ");
            //	System.out.print(" drift= " + drift + " ");
                year = offsets.getParameterValue(i);   //just want year[i]
                		//make sure that it is equivalent to double offset  = year[virusIndex] - firstYear;
            }
            else{
            	System.out.println("virus Offeset Parameter NOT present. We expect one though. Something is wrong.");
            }
		meanYear[ label] += year;
		
		groupCount[ label  ] = groupCount[ label ]  +1; 
	}
	*/
	int maxLabel=0;
	for(int i=0;i< numdata; i++){
		if(maxLabel < (int) clusterLabels.getParameterValue(i)){
			maxLabel = (int) clusterLabels.getParameterValue(i);
		}
	}
	
	
	/*
	for(int i=0; i <= maxLabel; i++){
		meanYear[i] = meanYear[i]/groupCount[i];
		//System.out.println(meanYear[i]);
	}
	*/
	
	//double beta = locationDrift.getParameterValue(0);
	
	
	//now, change the mu..
	for(int i=0; i <= maxLabel; i++){
		//System.out.println(meanYear[i]*beta);
		//mu.getParameter(i).setParameterValue(0, meanYear[i]*beta);//now separate out mu from virusLocation
		mu.getParameter(i).setParameterValue(0, 0);
		mu.getParameter(i).setParameterValue(1, 0);
	}	
	
	//now change the clusterOffsets
	//... is it necessary?
	
	//System.exit(0);
}



public double getLogLikelihood() {
	
	//System.out.println("getLogLikelihood of ClusterViruses");
	
	double logL = 0;

	
	int maxLabel=0;
	for(int i=0;i< numdata; i++){
		if(maxLabel < (int) clusterLabels.getParameterValue(i)){
			maxLabel = (int) clusterLabels.getParameterValue(i);
		}
	}

	//P(K=k)
	int K_int = (int) K.getParameterValue(0);
	//logL += Math.log(K.getParameterValue(0)) - lambda*K.getParameterValue(0) - Math.log( (double) factorial(K_int));
	logL += -lambda + K.getParameterValue(0)*Math.log(lambda) - Math.log( (double) factorial(K_int));

	
	// p(C | K= k)
	logL -= numdata * Math.log(K.getParameterValue(0));
	
	//p(mu_j | C, years) ~ N( theta , sigma^2)
	//logL -= Math.log(2*Math.PI);
	//logL -= 0.5*Math.log( sigmaSq*sigmaSq);
	
	
	logL -= K.getParameterValue(0) * ( Math.log(2)  + Math.log(Math.PI)+ 0.5*Math.log(sigmaSq) +  0.5*Math.log(sigmaSq) );
	
	//double[] meanYear = new double[K_int];
	double[] groupCount = new double[K_int];
	for(int i=0; i < numdata; i++){
		int label = (int) clusterLabels.getParameterValue(i);
		groupCount[ label  ] = groupCount[ label ]  +1; 
	}
	
	//for(int i=0; i <= maxLabel; i++){
		//meanYear[i] = meanYear[i]/groupCount[i];
		//System.out.println(meanYear[i]);
	//}
	
	
	for(int i=0; i <= maxLabel; i++){
		double mu_i0 = mu.getParameter(i).getParameterValue(0);
		double mu_i1 = mu.getParameter(i).getParameterValue(1);
		
		//double beta = locationDrift.getParameterValue(0);
		if( groupCount[i] >0){
			//System.out.println("meanYear = " + meanYear[i]);
			//logL -=	0.5*(  (mu_i0 - beta*meanYear[i] )*(mu_i0 - beta*meanYear[i] ) + ( mu_i1 -0)*( mu_i1 - 0)   )/sigmaSq;
			logL -=	0.5*(  (mu_i0  )*(mu_i0  ) + ( mu_i1 )*( mu_i1 )   )/sigmaSq;
		}
	}
	
	//System.out.println(logL);

//	System.out.println("logL=" + logL);	
	//System.out.println("done getLogLikelihood of ClusterViruses");

	//System.out.println("logL=" + logL);

	
	//double logL = 0; // for testing purpose only
	return(logL);
/*	
// if treeModel changes, compute the depMatrix from scratch
// if only the transformFactor change, go back to the latest copy of the untransformed deptMatrix,
//and transform it, so it doesn't have to go through the treeModel to get the distance.
if(treeChanged==true){
  // setDepMatrix();   //the super slow step          
 if(treeChanged ==true){
	  treeChanged = false;
 }
 
}


 double logL = 0.0;
         for (int j=0 ; j<logLikelihoodsVector.length;j++){
                  logLikelihoodsVector[j]=getLogLikGroup(j);

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
     
        
        return logL;
      */
 }

    
    
//=====================================================================================================================
        
        
public  int factorial(int n) {
    int fact = 1; // this  will be the result
    for (int i = 1; i <= n; i++) {
        fact *= i;
    }
    return fact;
}
	




    public Model getModel() {
                return this;
            }

          




    public void makeDirty() {
            }

            public void acceptState() {
                // DO NOTHING
            }

            public void restoreState() {
                // DO NOTHING
            }

            public void storeState() {
                // DO NOTHING
            }

            

            protected void handleModelChangedEvent(Model model, Object object, int index) {
                //---------- suggestion from cykc----------------
                // I am making an assumption that whenever treeModel changes, the changes get caught here.             	
            	if(model == treeModel){
            		//System.out.println("==========Tree model changes!!!!!!!!=====");
            		treeChanged = true;
            	}
                else{
                }
            	//---------- End of suggestion from cykc----------------      	
            }

            
          //---------- suggestion from cykc----------------
            //I tried to catch the transformedFactor changes through this routine, but it seems that it doesn't catch all,
            //so I abandoned this routine and now use 'mostRecentTransformedValue' to directly test if transformedValue has changed.
            //This ShouldUpdateDepMatrix never gets used and is now an obsolete variable.
            
            //I am noticing that handleVariableChangedEvent doesn't always catch when transformFactor changes.
            //I am observing that sometimes transformFactor can change more than once within a single MCMC sample - I don't know why,
            //if getLogLikelihood() gets called after the transformFactor changes but ShouldUpdateDepMatrix flag doesn't catch it,
            //then the getLogLikelihood() will not be calculated correctly.
            //Hence, this way of catching when transformFactor changes now becomes obsolete
            //
            protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            	/*
            	if (variable == transformFactor) {
            	//System.out.println("TransformFactor gets updated and is caught here!!");
            		ShouldUpdateDepMatrix = true;

                } else {
                	//has to change at another sample instead of setting it to false right after
                	//updating the matrix because the transformFactor value 
                	//can update twice within a sample, so setDepMatrix has to update twice
                	if(ShouldUpdateDepMatrix == true){
                		//System.out.println("ShouldUpdateDepMatrix changes from true to false");
                		ShouldUpdateDepMatrix = false;
                	}
                }
            	 */
            }
            //---------- End of suggestion from cykc----------------

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

    TreeModel treeModel;
        String traitName;
        double[][] data;
        double[][] depMatrix;
        double[][] logDepMatrix;
        double[][] cur_untransformedMatrix; //---------- suggestion from cykc----------------
        
        double[] logLikelihoodsVector;
        int numdata;
        Parameter transformFactor;
        double k0;
    double v0;
    double[][] T0Inv;
        double[] m;
    double logDetT0;
    
    LinkedList<Integer>[] assignmentsLL; 
    int seqLength;
    
    public int getSeqLength() {
		return seqLength;
	}

    char[][] seqData;




    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
    	
    	
    	public final static String TREEMODEL = "treeModel";
    	public final static String K = "k";
    	public final static String EXCISIONPOINTS = "excisionPoints";
    	public final static String CLUSTERLABELS = "clusterLabels";
    	public final static String  MU = "mu";
    	//public final static String HASDRIFT = ??
 //   	public final static String LOCATION_DRIFT = "locationDrift";
    	public final static String OFFSETS = "offsets";
    	public final static String VIRUS_LOCATIONS = "virusLocations";

        boolean integrate = false;
        
        
        public String getParserName() {
            return CLUSTER_VIRUSES;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
     
                TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
                //String traitName = (String) xo.getAttribute(TRAIT_NAME);
                
                XMLObject cxo = xo.getChild(K);
                Parameter k = (Parameter) cxo.getChild(Parameter.class);
                

                cxo = xo.getChild(EXCISIONPOINTS);
                Parameter excisionPoints = (Parameter) cxo.getChild(Parameter.class);

                cxo = xo.getChild(CLUSTERLABELS);
                Parameter clusterLabels = (Parameter) cxo.getChild(Parameter.class);

                cxo = xo.getChild(MU);
                MatrixParameter mu = (MatrixParameter) cxo.getChild(MatrixParameter.class);

                //alternative way to load in the MatrixParameter?
              //  MatrixParameter serumLocationsParameter = null;
               // if (xo.hasChildNamed(SERUM_LOCATIONS)) {
                //    serumLocationsParameter = (MatrixParameter) xo.getElementFirstChild(SERUM_LOCATIONS);
                //}

                
                
      //          cxo=xo.getChild(LOCATION_DRIFT) ;
   //             Parameter locationDrift= (Parameter) cxo.getChild(Parameter.class);

                cxo=xo.getChild(OFFSETS);
                Parameter offsets =(Parameter) cxo.getChild(Parameter.class);

                cxo=xo.getChild(VIRUS_LOCATIONS);
                MatrixParameter virusLocations =(MatrixParameter) cxo.getChild(MatrixParameter.class);
                
            boolean hasDrift = false;
            if (offsets.getDimension()>1){
                hasDrift=true;
            }

            
           // TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();
           // String traitName = TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;
         
            
          //  TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
         //           utilities.parseTraitsFromTaxonAttributes(xo, traitName, treeModel, integrate);
         //  // traitName = returnValue.traitName;
         //   CompoundParameter traitParameter = returnValue.traitParameter;
                
 
			//return new 	ClusterViruses (treeModel,traitParameter ,  K, excisionPoints, clusterLabels, mu,  hasDrift); 
		//	return new ClusterViruses(treeModel, k, excisionPoints, clusterLabels, mu, hasDrift, locationDrift, offsets, virusLocations);
        	return new ClusterViruses(treeModel, k, excisionPoints, clusterLabels, mu, hasDrift, offsets, virusLocations);

					
            }

            //************************************************************************
            // AbstractXMLObjectParser implementation
            //************************************************************************

            public String getParserDescription() {
                return "clustering viruses";
            }

            public Class getReturnType() {
                return ClusterViruses.class;
            }

            public XMLSyntaxRule[] getSyntaxRules() {
                return rules;
            }
            
            
            private final XMLSyntaxRule[] rules = {
                   // new StringAttributeRule(TreeTraitParserUtilities.TRAIT_NAME, "The name of the trait for which a likelihood should be calculated"),
                   // new ElementRule(TREEMODEL, Parameter.class),
                    new ElementRule(K, Parameter.class),
                    new ElementRule(EXCISIONPOINTS, Parameter.class),
                    new ElementRule(CLUSTERLABELS, Parameter.class),
                    new ElementRule(MU, MatrixParameter.class),
      //              new ElementRule(LOCATION_DRIFT, Parameter.class),
                    new ElementRule(OFFSETS, Parameter.class),
                    new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class)
            };
            
    };


    
  

    String Atribute = null;

        
}


