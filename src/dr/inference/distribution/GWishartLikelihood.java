/*
 * GWishartLikelihood.java
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

package dr.inference.distribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.GammaFunction;
import dr.math.MathUtils;
import dr.math.distributions.GWishartStatistics;//--------------Felipe
import dr.math.distributions.MultivariateDistribution;
import dr.math.distributions.WishartStatistics;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;



/**
 * @author Felipe Pinheiro
 */
 


public class GWishartLikelihood extends AbstractModelLikelihood  implements MultivariateDistribution, GWishartStatistics, WishartStatistics {
	public static final String DF = "df";
	public static final String RATE_MATRIX = "rateMatrix";
	public static final String ADJ_INDICATORS = "adjIndicators";
    public final static String GWISHART_LIKELIHOOD = "GWishartLikelihood";
    public static final String GWISHART_PRIOR = "multivariateGWishartPrior";
    public static final String NON_INFORMATIVE = "nonInformative";
    public static final String DATA = "data";

    private double df;
    private int dim;
    private double[][] rateMatrix;
    private double[][] scaleMatrix;
    public final MatrixParameter precisionParameter;
    public final Parameter adjParameter;
    public double[][] adjMatrix; 

    private double[] D;
    private Matrix DMat;     
    private double[] S;
    private Matrix SMat;    
    
    private boolean logNormalizationConstantKnown = false;
    private double storedLogNormalizationConstant;
    private double logNormalizationConstant;

    protected double logLikelihood;
    private double storedLogLikelihood;
    protected boolean logLikelihoodKnown = false;
    
   

   
    /**
     * A G-Wishart distribution class for df degrees of freedom and rate matrix D
     * @param df (or delta)	degrees of freedom or shape \delta 
     * @param rateMatrix	rate matrix
     * @param adjMatrix 	adjMatrix
     */
    
    //======================================================================================================  CONSTRUCTORS
    //Constructor for SparsePrecisionGibbsOperator
    public GWishartLikelihood(MatrixParameter precisionParameter, double df,  double[][] Matrix, Parameter adjParameter, boolean isScale) {
        super(GWISHART_LIKELIHOOD);
    	this.df = df;        
        this.dim = Matrix.length;
        this.precisionParameter = precisionParameter;
        addVariable(precisionParameter);
        
        if (adjParameter == null) {
            this.adjParameter = new Parameter.Default(dim*(dim-1)/2, 1.0);
        } else {
            this.adjParameter = adjParameter;
            addVariable(adjParameter);
        }
        this.adjMatrix = SymmetricMatrix.compoundSymmetricMatrix(new double[dim], adjParameter.getParameterValues(), dim).toComponents();
        
        
        
        
        //System.out.println(">---------------------------------------------------------------------------------->>> GWishartLikelihood Constructor 4");
        //System.out.println("dim="+dim);
        //System.out.println("df="+df);      
        //System.out.println("precisionParameter="+Arrays.deepToString(precisionParameter.getValues()));      

        //System.out.println("adjParameter="+Arrays.deepToString(adjParameter.getValues()));      
        //System.out.println("adjMatrix="+Arrays.deepToString(adjMatrix));       
           
        //Make adjMatrix upper triangular
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j <= i; j++) {
           	 	adjMatrix[i][j]=0.0;
            }
        }
        //System.out.println("adj upper="+Arrays.deepToString(adjMatrix)); 

        if(isScale) {
        	this.scaleMatrix = Matrix;
        	SMat = new Matrix(scaleMatrix);
            double[][] tmp = SMat.toComponents();
            S = new double[dim * dim];
            for (int i = 0; i < dim; i++) {
                System.arraycopy(tmp[i], 0, S, i * dim, dim);
            }
            DMat = SMat.inverse();
            double[][] tmp2 = DMat.toComponents();
            D = new double[dim * dim];
            for (int i = 0; i < dim; i++) {
                System.arraycopy(tmp2[i], 0, D, i * dim, dim);
            }
            this.rateMatrix=SMat.inverse().toComponents();
        }else{
        	this.rateMatrix = Matrix;
        	DMat = new Matrix(rateMatrix);
            double[][] tmp = DMat.toComponents();
            D = new double[dim * dim];
            for (int i = 0; i < dim; i++) {
                System.arraycopy(tmp[i], 0, D, i * dim, dim);
            }
            SMat = DMat.inverse();
            double[][] tmp2 = SMat.toComponents();
            S = new double[dim * dim];
            for (int i = 0; i < dim; i++) {
                System.arraycopy(tmp2[i], 0, S, i * dim, dim);
            }
            this.scaleMatrix=DMat.inverse().toComponents();

        }
       
        //System.out.println("DMat="+Arrays.deepToString(DMat.toComponents()));
        //System.out.println("D="+Arrays.toString(D));
        //System.out.println("SMat="+Arrays.deepToString(SMat.toComponents()));
        //System.out.println("S="+Arrays.toString(S));
        
        
        computeNormalizationConstant();
        //System.out.println("<<<----------------------------------------------------------------------------------< GWishartLikelihood Constructor 4");

    }
     
    public GWishartLikelihood(int dim) { // returns a non-informative (unormalizable) density
        super(GWISHART_LIKELIHOOD);
    	this.df = 0;
        this.scaleMatrix = null;
        this.rateMatrix = null;
        this.adjMatrix = null;
        this.adjParameter = null;
        this.precisionParameter = null;

        this.dim = dim;
        logNormalizationConstant = 0.0;
    }


    
    
    
    
    
    //==============================================================================================  NORMALIZING CONSTANT
    private void computeNormalizationConstant() {
        logNormalizationConstant = computeNormalizationConstant( df, scaleMatrix, adjParameter, dim);
    }

    public static double computeNormalizationConstant( double df, double[][] scaleMatrix, Parameter adjParameter, int dim) {
    	//System.out.println("\n------------------------------------ Calculating the log Normalizing Constant log I_G(delta, D)");
        if (df == 0) {
            return 0.0;
        }
        //System.out.println("adjParameter in computeNormalizationConstant="+ Arrays.toString(adjParameter.getParameterValues()));
        double[][] adjMatrix = SymmetricMatrix.compoundSymmetricMatrix(new double[dim], adjParameter.getParameterValues(), dim).toComponents();
        //System.out.println("adjMatrix="+ Arrays.deepToString(adjMatrix));

        //------------------------------------------- Normalizing Constant for FULL GRAPH (this is the Wishart Normalizing Constant)
        //boolean size = getGraphSize(adjMatrix) == dim*(dim-1)/2;
        //System.out.println("Complete graph? "+ size);
        
        if(getGraphSize(adjMatrix) == dim*(dim-1)/2) {
        	//System.out.println("Full Graph");
       
        	/* CAUTION: Wishart-Family Parametrizations
        	 * 
        	 *	G-Wishart is parametrized with shape parameter or degrees of freedom (\delta) and RATE matrix (D):
        	 *			df = \delta,			where \delta = \nu - dim + 1, as a comparison to the df in Wishart distribution (\nu);
        	 *			rate Matrix = D, 		where D is the inverse of scale matrix S: D=S^-1.
        	 *
        	 *	Wishart is parametrized with shape parameter or degrees of freedom (\nu) and SCALE matrix (S):
        	 *			df = \nu, 				where \nu = \delta + dim - 1;  
        	 *			scale Matrix = S
        	 *
        	 *	The relation between \delta and \nu is: \nu = \delta + dim - 1. 
        	 *
        	 *  So, compared to WishartDistribution Class in BEAST:
        	 *		- the df (\nu) in WishartDistribution Class is replaced here by \nu = df + dim - 1 (\nu in terms of \delta);
        	 *		- the scale matrix (S) in WishartDistribution Class is expressed here in terms of rate matrix D as S=D.inverse().
        	 */
        	
        	//Matrix S = new Matrix(rateMatrix).inverse();
        	Matrix S = new Matrix(scaleMatrix);
        	double nu = df + dim - 1;
        	
            double logNormalizationConstant = 0;
            //System.out.println("log IG(nu,S) = "+logNormalizationConstant);
            
            try {
                logNormalizationConstant = -nu / 2.0 * Math.log(S.determinant());
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }
            
            //System.out.println("log IG(nu,S) = "+logNormalizationConstant);
            logNormalizationConstant -= nu * dim / 2.0 * Math.log(2);
            //System.out.println("log IG(nu,S) = "+logNormalizationConstant);

            logNormalizationConstant -= dim * (dim - 1) / 4.0 * Math.log(Math.PI);
            for (int i = 1; i <= dim; i++) {
                logNormalizationConstant -= GammaFunction.lnGamma((nu + 1 - i) / 2.0);
            }
            
            //System.out.println("log IG( nu = delta + p - 1, S = D^-1 ) = log IG(delta, D) = "+ -logNormalizationConstant);
        	//System.out.println("----------------------------------------------------------------------------------------------");

            return logNormalizationConstant;
        }
       
        
        //--------------------------------------------------------------------------Normalizing Constant for a GENERIC GRAPH
        //System.out.println("------------------------- Calculating C(delta,D)"); 
        double[] vi = rowSums(adjMatrix); 
       	//System.out.println("vi="+Arrays.toString(vi));
        
        double[] ki = colSums(adjMatrix); 
       	//System.out.println("ki="+Arrays.toString(ki));
       	    	
       	//Matrix DMat = new Matrix(rateMatrix);  	
       	//double[][] T = getTMatrix(DMat);
       	double[][] T = getTMatrix(scaleMatrix);
       	//System.out.println("T="+Arrays.deepToString(T));

        //System.out.println("---------------------- Calculating f_D(delta,D)"); 	
       	double[][] H = getHMatrix(T);
       	//System.out.println("H="+Arrays.deepToString(H));

       	int mc = 200;
       	//System.out.println("MC iterations="+mc);
       	
       	//double[] f_D = new double [mc];
       	
       	//double Ef_D = getFD(df, H, adjMatrix, vi, dim, f_D, mc );
       	double Ef_D = getFD(df, H, adjMatrix, vi, dim, mc );
       	//System.out.println("Ef_D="+Ef_D);
       	     	
    	double logC_dfD = 0;
       	for (int i = 0; i < dim; i++) {
           	 	logC_dfD += +( (df + 2*vi[i])/2.0 ) * Math.log(2) + 
           	 			vi[i]/2 * Math.log(Math.PI) + 
           	 			GammaFunction.lnGamma( (df + vi[i]) / 2.0 ) +
           	 			( df + vi[i] + ki[i] )*Math.log( T[i][i] );
        }
       	//System.out.println("log C_dfD="+logC_dfD);
       	
       	//System.out.println("------------------------------------------------"); 
       	if (Double.isInfinite(logC_dfD) || Double.isNaN(logC_dfD)) {//####################################################### -Inf or NaN???
            return Double.NaN;
        }
       	
       	double logIG_dfD = logC_dfD + Math.log(Ef_D);
       	//System.out.println("log I_G(delta,D) = "+logIG_dfD);
        //System.out.println("------------------------------------------------"); 
        
     
    	//System.out.println("----------------------------------------------------------------------------------------------");
        return -logIG_dfD;
    }

    

    public static double computeNormalizationConstant( double df, double[][] scaleMatrix, double[][] adjMatrix, int dim) {
    	//System.out.println("\n------------------------------------ Calculating the log Normalizing Constant log I_G(delta, D)");
        if (df == 0) {
            return 0.0;
        }
        //System.out.println("adjMatrix in computeNormalizationConstant="+ Arrays.toString(adjMatrix));
        //double[][] adjMatrix = SymmetricMatrix.compoundSymmetricMatrix(new double[dim], adjParameter.getParameterValues(), dim).toComponents();
       // System.out.println("adjMatrixin computeNormalizationConstant="+ Arrays.deepToString(adjMatrix));

        //------------------------------------------- Normalizing Constant for FULL GRAPH (this is the Wishart Normalizing Constant)
        //boolean size = getGraphSize(adjMatrix) == dim*(dim-1)/2;
        //System.out.println("Complete graph? "+ size);
        
        if(getGraphSize(adjMatrix) == dim*(dim-1)/2) {
        	//System.out.println("Full Graph");
       
        	/* CAUTION: Wishart-Family Parametrizations
        	 * 
        	 *	G-Wishart is parametrized with shape parameter or degrees of freedom (\delta) and RATE matrix (D):
        	 *			df = \delta,			where \delta = \nu - dim + 1, as a comparison to the df in Wishart distribution (\nu);
        	 *			rate Matrix = D, 		where D is the inverse of scale matrix S: D=S^-1.
        	 *
        	 *	Wishart is parametrized with shape parameter or degrees of freedom (\nu) and SCALE matrix (S):
        	 *			df = \nu, 				where \nu = \delta + dim - 1;  
        	 *			scale Matrix = S
        	 *
        	 *	The relation between \delta and \nu is: \nu = \delta + dim - 1. 
        	 *
        	 *  So, compared to WishartDistribution Class in BEAST:
        	 *		- the df (\nu) in WishartDistribution Class is replaced here by \nu = df + dim - 1 (\nu in terms of \delta);
        	 *		- the scale matrix (S) in WishartDistribution Class is expressed here in terms of rate matrix D as S=D.inverse().
        	 */
        	
        	//Matrix S = new Matrix(rateMatrix).inverse();
        	Matrix S = new Matrix(scaleMatrix);
        	double nu = df + dim - 1;
        	
            double logNormalizationConstant = 0;
            //System.out.println("log IG(nu,S) = "+logNormalizationConstant);
            
            try {
                logNormalizationConstant = -nu / 2.0 * Math.log(S.determinant());
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }
            
            //System.out.println("log IG(nu,S) = "+logNormalizationConstant);
            logNormalizationConstant -= nu * dim / 2.0 * Math.log(2);
            //System.out.println("log IG(nu,S) = "+logNormalizationConstant);

            logNormalizationConstant -= dim * (dim - 1) / 4.0 * Math.log(Math.PI);
            for (int i = 1; i <= dim; i++) {
                logNormalizationConstant -= GammaFunction.lnGamma((nu + 1 - i) / 2.0);
            }
            
            //System.out.println("log IG( nu = delta + p - 1, S = D^-1 ) = log IG(delta, D) = "+ -logNormalizationConstant);
        	//System.out.println("----------------------------------------------------------------------------------------------");

            return logNormalizationConstant;
        }
       
        
        //--------------------------------------------------------------------------Normalizing Constant for a GENERIC GRAPH
        //System.out.println("------------------------- Calculating C(delta,D)"); 
        double[] vi = rowSums(adjMatrix); 
       	//System.out.println("vi="+Arrays.toString(vi));
        
        double[] ki = colSums(adjMatrix); 
       	//System.out.println("ki="+Arrays.toString(ki));
       	    	
       	//Matrix DMat = new Matrix(rateMatrix);  	
       	//double[][] T = getTMatrix(DMat);
       	double[][] T = getTMatrix(scaleMatrix);
       	//System.out.println("T="+Arrays.deepToString(T));

        //System.out.println("---------------------- Calculating f_D(delta,D)"); 	
       	double[][] H = getHMatrix(T);
       	//System.out.println("H="+Arrays.deepToString(H));

       	int mc = 200;
       	//System.out.println("MC iterations="+mc);
       	
       	//double[] f_D = new double [mc];
       	
       	//double Ef_D = getFD(df, H, adjMatrix, vi, dim, f_D, mc );
       	double Ef_D = getFD(df, H, adjMatrix, vi, dim, mc );
       	//System.out.println("Ef_D="+Ef_D);
       	     	
    	double logC_dfD = 0;
       	for (int i = 0; i < dim; i++) {
           	 	logC_dfD += +( (df + 2*vi[i])/2.0 ) * Math.log(2) + 
           	 			vi[i]/2 * Math.log(Math.PI) + 
           	 			GammaFunction.lnGamma( (df + vi[i]) / 2.0 ) +
           	 			( df + vi[i] + ki[i] )*Math.log( T[i][i] );
        }
       	//System.out.println("log C_dfD="+logC_dfD);
       	
       	//System.out.println("------------------------------------------------"); 
       	if (Double.isInfinite(logC_dfD) || Double.isNaN(logC_dfD)) {//It should be -Inf or NaN?
            return Double.NaN;
        }
       	
       	double logIG_dfD = logC_dfD + Math.log(Ef_D);
       	//System.out.println("log I_G(delta,D) = "+logIG_dfD);
        //System.out.println("------------------------------------------------"); 
        
     
    	//System.out.println("----------------------------------------------------------------------------------------------");
        return -logIG_dfD;
    }

    
    
    
     
    //==================================================================================  ABSTRACT MODEL LIKELIHOOD METHODS
    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
    	//System.out.println("-------------------------------handleVariableChangedEvent in GWishartLikelihood");

        if(variable instanceof Parameter) {
        	//System.out.println("variable instanceof Parameter, logNormalizingationConstant=false");
        	logNormalizationConstantKnown = false;
        }
        
        if(variable instanceof MatrixParameter) {
        	//System.out.println("variable instanceof MatrixParameter: loglikeknown=false");
        	logLikelihoodKnown = false;
        }
        

    }

    @Override
    public void storeState() {
    	//System.out.println("----------------------------------storeState() in GWishartLikelihood");

        storedLogNormalizationConstant = logNormalizationConstant;
        storedLogLikelihood = logLikelihood;
    }

    @Override
    public void restoreState() {
    	//System.out.println("---------------------------------restoreState() in GWishartLikelihood");
        logNormalizationConstant = storedLogNormalizationConstant;
        logNormalizationConstantKnown = true;
        
        logLikelihood = storedLogLikelihood;
        logLikelihoodKnown = true;
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    public void makeDirty() {
    	//System.out.println("---------------------------------makeDirty() in GWishartLikelihood");
    	logLikelihoodKnown = false;

    }

    public Model getModel() {
        return this;
    }

    public final double getLogLikelihood() {
    	//System.out.println("----------------------------getLogLikelihood in GWishartLikelihood");
    	
    	//System.out.println("logNormalizationConstantKnown?"+logNormalizationConstantKnown);
    	//System.out.println("storedIG BEFORE="+logNormalizationConstant);

    	 if (!logNormalizationConstantKnown) {
             //logNormalizationConstant = computeNormalizationConstant( df, rateMatrix, adjParameter, dim);
             logNormalizationConstant = computeNormalizationConstant( df, scaleMatrix, adjParameter, dim);
             logNormalizationConstantKnown = true;
         	//System.out.println("storedIG AFTER="+logNormalizationConstant);

             logLikelihood = calculateLogLikelihood();
             logLikelihoodKnown = true; 
             return logLikelihood;
         }
      	//System.out.println("storedIG AFTER="+logNormalizationConstant);

    	//System.out.println("logLikelihoodKnown? "+logLikelihoodKnown);
        if (!logLikelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            logLikelihoodKnown = true;
        }
        return logLikelihood;
    }
    
    protected double calculateLogLikelihood() {
    	//System.out.println("--------------computeLogLikelihood() in GWishartLikelihood");
    	//System.out.println("precisionParameter="+Arrays.toString(precisionParameter.getParameterValues()));
    	
    	Matrix W = new Matrix(precisionParameter.getParameterValues(), dim, dim);
    	
    	return logPdf(W, DMat, df, dim, logNormalizationConstant);
    }
    
    public GWishartLikelihood getDistribution() {
        return this;
    }
    
    
    
    
    
    //==================================================================================================  G-WISHART METHODS
    public String getType() {
        return GWISHART_PRIOR;
    }

    public double[][] getRateMatrix() {
        return rateMatrix;
    }
    
    public double[][] getScaleMatrix() {
        return scaleMatrix;
    }
    
    public double[] getMean() {
        return null;
    }
        
    public double getDF() {
        return df;
    }
    
    public Parameter getAdjParameter() {
        return adjParameter;
    }
    
    public double[][] getAdjMatrix() {
        return SymmetricMatrix.compoundSymmetricMatrix(new double[dim], adjParameter.getParameterValues(), dim).toComponents();
    }
   
    
    
    public static double getGraphSize(double[][] a) {
    	int dim = a.length;
        double sum = 0;
      	 for (int i = 0; i < dim; i++) {
             for (int j = i; j < dim; j++) {
            	 sum += a[i][j]; 
             }
         }
        return sum;
    }
    
    public static double[] rowSums(double[][] a) {
    	int dim = a.length;
        double[] rowsums = new double[dim]; 

      	 for (int i = 0; i < dim; i++) {
      		 double sum = 0;
             for (int j = i; j < dim; j++) {
            	 sum += a[i][j]; 
             }
             rowsums[i]= sum;
         }
        return rowsums;
    }
    
    public static double[] colSums(double[][] a) {
    	int dim = a.length;
    	double[] colsums = new double[dim]; 
     	 for (int j = 0; j < dim; j++) {
      		  int sum = 0;
             for (int i = 0; i < j; i++) {
            	 sum += a[i][j]; 
             }
             colsums[j]= sum;
         }
        return colsums;
    }

    
    public static double[][] getTMatrix(double[][] a) {
    	int dim = a.length;
	
        double[][] cholesky = a;
        //System.out.println("D="+Arrays.deepToString(cholesky));
        
        try {
            cholesky = (new CholeskyDecomposition(cholesky)).getL();
            // caution: this returns the lower triangular form
        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Numerical exception in GWishartLikelihood");
        }
        //System.out.println("T'="+Arrays.deepToString(cholesky));
        
        double[][] T = new double[dim][dim]; 
        for (int i = 0; i < dim; i++) {
 	         for (int j = 0; j <dim; j++) {
 	        	 T[i][j]=cholesky[j][i]; //------------------------T is UPPER TRIANGULAR
 	         }
        }
        //System.out.println("T="+Arrays.deepToString(T));
        return T;
    }  
    
    public static double[][] getHMatrix(double[][] a) {
    	int dim = a.length;

        double[][] H = new double[dim][dim]; 
        for (int i = 0; i < dim; i++) {
 	         for (int j = 0; j < dim; j++) {
 	        	 H[i][j]=a[i][j]/a[j][j]; 
 	         }
        }
        //System.out.println("H="+Arrays.deepToString(H));
        //a and H are upper triangular
        return H;
    }  
    
    private static boolean isDiagonal(double[][] x) {
        for (int i = 0; i < x.length; ++i) {
            for (int j = i + 1; j < x.length; ++j) {
                if (x[i][j] != 0.0) {
                    return false;
                }
            }
        }
        return true;
    }
    
    public static double getFD(double df, double[][] H, double[][] G, double nu[], int dim, int mc ){
    	double sumPsi, sumPsiH, sumPsiHi, sumPsiHj, Ef_D;
    	double max_numeric_limits_ld = Double.MAX_VALUE / 1000;
    	double min_numeric_limits_ld = Double.MIN_VALUE * 1000;
    	double[] f_D = new double [mc];
    	double[][] psi = new double[dim][dim];      


    	if( isDiagonal(H) ){  //---------------- D is diagonal, then Ti=chol(solve(D)) is diagonal and Psi_ij, where Gij==0, will be 0 (because T_{<kj]}=0)
    		//System.out.println("H is diagonal");
    		for( int iter = 0; iter < mc; iter++ ) {
    			//---------------------------------------------------------------------------------------PSI = DIAGONAIS + FREE ELEMENTS IN G
    	        for (int i = 0; i < dim; i++) {
    	            psi[i][i] = Math.sqrt(MathUtils.nextGamma((df + nu[i]) * 0.5, 0.5));   //*****************************************************************************************     
    	        }
    			for( int i = 0; i < dim - 1; i++ ) {
    				for( int j = i + 1; j < dim; j++ ){
    					psi[i][j] = ( G[i][j] == 1 ) ? MathUtils.nextGaussian() : 0.0;				
    				}
    			}
    			//-----------------------------------------------------------------------NON-FREE ELEMENTS BY ROW (Off-diagonals where G_ij==0)
    			for(int i = 0; i < dim - 1; i++ ) { 
    				for(int j = i + 1; j < dim; j++ ){
    					if( G[i][j] == 0 ){
    						
    						psi[i][j] = 0.0; // it's not necessary
    						
    						if( i > 0 ){
    							sumPsi = 0.0;
    							for(int r = 0; r < i; r++ ){
    								if( sumPsi == Double.POSITIVE_INFINITY ) sumPsi = max_numeric_limits_ld;	
    								if( sumPsi == Double.NEGATIVE_INFINITY ) sumPsi = min_numeric_limits_ld;	
    								sumPsi += ( psi[r][i] * psi[r][j] );
    							}
    							if( sumPsi == Double.POSITIVE_INFINITY ) sumPsi = max_numeric_limits_ld;	
								if( sumPsi == Double.NEGATIVE_INFINITY ) sumPsi = min_numeric_limits_ld;
    							psi[i][j] = - sumPsi / psi[i][i];
    						}
    						
    						if( psi[i][j] == Double.POSITIVE_INFINITY ) psi[i][j] = max_numeric_limits_ld;	
    						if( psi[i][j] == Double.NEGATIVE_INFINITY ) psi[i][j] = min_numeric_limits_ld;	
    												
    						f_D[iter] += ( psi[i][j] * psi[i][j] ); 
    					}
    				}
    			}
    		
    			// checking Inf values
    			if( f_D[iter] == Double.POSITIVE_INFINITY ) f_D[iter] = max_numeric_limits_ld;			
    		} 
    	}else{//------------------------------------------------------------------------------D is not diagonal, the we have the Ti's
    		//System.out.println("H is NOT diagonal");

    		for (int iter = 0; iter < mc; iter++ ) {
    			
    			//---------------------------------------------------------------------------------------------- Psi  FREE PARAMETERS 
    			for (int i = 0; i < dim; i++){
    	            psi[i][i] = Math.sqrt(MathUtils.nextGamma((df +nu[i]) * 0.5, 0.5)); //-----------Check gamma parametrizations ********
    			}
    
    			for(int i = 0; i < dim - 1; i++ ) { 
    				for(int j = i + 1; j < dim; j++ ) {
    					psi[i][j] = ( G[i][j] == 1 ) ? MathUtils.nextGaussian() : 0.0;
    				}
    			}
    			
    			//------------------------------------------------------------------------ Psi_ij NON-FREE (off-diagonals where G_ij=0)
    			for(int i = 0; i < dim - 1; i++ ) {
    				for(int j = i + 1; j < dim; j++ ) {
    					
    					if( G[i][j] == 0 )
    					{
    						sumPsiH = 0.0;
    						//------------------------------------------------------------------------------------------------- SUM #1
    						for(int k = i; k < j; k++ )//-------Calculate the first expression in Attay and Massam (2005) for all (ij)
    						{
    							//Double.isInfinite(sumPsi)//-----------------------------------------------------------------to check*************
    							if( sumPsiH == Double.POSITIVE_INFINITY ) sumPsiH = max_numeric_limits_ld;	
    							if( sumPsiH == Double.NEGATIVE_INFINITY ) sumPsiH = min_numeric_limits_ld;	
    							sumPsiH += ( psi[i][k] * H[k][j] ); 
    						}
    						psi[i][j] = - sumPsiH;

    						//-------------------------------------------------------------------------------------------------SUM #2
    						if(i > 0)                     //------------------------Second part is calculated only for next lines
    							for(int r = 0; r < i; r++){
    								
    								sumPsiHi = 0.0;
    								for(int h = r; h < i + 1; h++){ //if this sum goes until i (not i-1) it will include psi_ri*t_<ii]=psi_ri and psi_rj*t_<jj]=psi_rj
    									if( sumPsiHi == Double.POSITIVE_INFINITY ) sumPsiHi = max_numeric_limits_ld;	
    									if( sumPsiHi == Double.NEGATIVE_INFINITY ) sumPsiHi = min_numeric_limits_ld;	
    									sumPsiHi += ( psi[r][h] * H[h][i] );	
    								}
    									
    								sumPsiHj = 0.0;
    								for(int h = r; h < j + 1; h++) {
    									if( sumPsiHj == Double.POSITIVE_INFINITY ) sumPsiHj = max_numeric_limits_ld;	
    									if( sumPsiHj == Double.NEGATIVE_INFINITY ) sumPsiHj = min_numeric_limits_ld;
    									sumPsiHj += ( psi[r][h] * H[h][j] );
    								}
    								
    								psi[i][j] -= ( sumPsiHi * sumPsiHj ) / psi[i][i];
    							}

    						if( psi[i][j] == Double.POSITIVE_INFINITY ) psi[i][j] = max_numeric_limits_ld;	
    						if( psi[i][j] == Double.NEGATIVE_INFINITY ) psi[i][j] = min_numeric_limits_ld;	

    						f_D[iter] += ( psi[i][j] * psi[i][j] ); 
    					}
    				}
    		
    			}
    			// checking Inf values
    			if( f_D[iter] == Double.POSITIVE_INFINITY ) f_D[iter] = max_numeric_limits_ld;			
    		}
    	}
    	
    	//System.out.println("f_D="+Arrays.toString(f_D));
    	Ef_D = 0.0;
    	for (int iter = 0; iter < mc; iter++ ) {
    		Ef_D += Math.exp(-f_D[iter]/2);
    	}
    	return Ef_D/mc;
    } 

    
    
    //=================================================================================================  SAMPLING G-WISHART
    public double[][] nextGWishart() {
    	//System.out.println("nextGWishart");
    	//System.out.println("df="+df);
    	//System.out.println("rateMatrix"+Arrays.deepToString(rateMatrix));
    	//System.out.println("scaleMatrix"+Arrays.deepToString(scaleMatrix));
    	//System.out.println("adjMatrix="+Arrays.deepToString(adjMatrix));
        //return nextGWishart(df, rateMatrix, adjMatrix);
        return nextGWishart(df, scaleMatrix, adjMatrix);
    }

    /**
     * Generates a random draw from a G-Wishart distribution
     * Follows Atay and Massam (2005) 
     * <p/>
     *
     * @param df          degrees of freedom
     * @param rateMatrix rateMatrix
     * @param adjMatrix adjMatrix
     * @return a random draw
     */
    

    public static double[][] nextGWishart(double df, double[][] scaleMatrix,  double[][] adjMatrix) {
    	System.out.println("nextGWishart");

    	int dim = scaleMatrix.length;
        double sumPsi, sumPsiH, sumPsiHi, sumPsiHj;
    	double max_numeric_limits_ld = Double.MAX_VALUE / 1000;
    	double min_numeric_limits_ld = Double.MIN_VALUE * 1000;
    	
    	double[] vi = rowSums(adjMatrix); 
       	//System.out.println("vi="+Arrays.toString(vi));
       	    	   	
       	double[][] T = getTMatrix(scaleMatrix);
       	//System.out.println("T="+Arrays.deepToString(T));

       	double[][] H = getHMatrix(T);
       	//System.out.println("H="+Arrays.deepToString(H));
    	
    	double[][] psi = new double[dim][dim];
    	double[][] phi= new double[dim][dim];      
        double[][] draw = new double[dim][dim];   

        //Calculating psi Matrix
    	if( isDiagonal(H) ){  
    		//System.out.println("H is diagonal");
    
    		for (int i = 0; i < dim; i++) {
    			psi[i][i] = Math.sqrt(MathUtils.nextGamma((df + vi[i]) * 0.5, 0.5));   //*********************************** TO CHECK    
    	    }
    		for( int i = 0; i < dim - 1; i++ ) {
    			for( int j = i + 1; j < dim; j++ ){
    				psi[i][j] = ( adjMatrix[i][j] == 1 ) ? MathUtils.nextGaussian() : 0.0;				
    			}
    		}
    		//-----------------------------------------------------------------------NON-FREE ELEMENTS BY ROW (Off-diagonals where G_ij==0)
    		for(int i = 0; i < dim - 1; i++ ) { 
    			for(int j = i + 1; j < dim; j++ ){
    				if( adjMatrix[i][j] == 0 ){
    						
    					psi[i][j] = 0.0; // it's not necessary
    						
    					if( i > 0 ){
    						sumPsi = 0.0;
    						for(int r = 0; r < i; r++ ){
    							if( sumPsi == Double.POSITIVE_INFINITY ) sumPsi = max_numeric_limits_ld;	
    							if( sumPsi == Double.NEGATIVE_INFINITY ) sumPsi = min_numeric_limits_ld;	
    							sumPsi += ( psi[r][i] * psi[r][j] );
    						}
    						if( sumPsi == Double.POSITIVE_INFINITY ) sumPsi = max_numeric_limits_ld;	
							if( sumPsi == Double.NEGATIVE_INFINITY ) sumPsi = min_numeric_limits_ld;
    						psi[i][j] = - sumPsi / psi[i][i];
   						}
   						
    					if( psi[i][j] == Double.POSITIVE_INFINITY ) psi[i][j] = max_numeric_limits_ld;	
    					if( psi[i][j] == Double.NEGATIVE_INFINITY ) psi[i][j] = min_numeric_limits_ld;	
    				}
    			}
    		}		
    	}else{//------------------------------------------------------------------------------D IS NOT DIAGONAL, we now have these T_i's
    		//System.out.println("H is NOT diagonal");
    		
    		//---------------------------------------------------------------------------------------------- Psi_ii and Psi_ij for FREE PARAMETERS 
    		for (int i = 0; i < dim; i++){
    			psi[i][i] = Math.sqrt(MathUtils.nextGamma((df +vi[i]) * 0.5, 0.5)); //************************************* CHECK THIS GAMMA PARAMETRIZATION
    		}
    
    		for(int i = 0; i < dim - 1; i++ ) { 
    			for(int j = i + 1; j < dim; j++ ) {
    				psi[i][j] = ( adjMatrix[i][j] == 1 ) ? MathUtils.nextGaussian() : 0.0;
    			}
    		}
    			
    		//------------------------------------------------------------------------ Psi_ij NON-FREE (off-diagonals where G_ij=0)
   			for(int i = 0; i < dim - 1; i++ ) {
    			for(int j = i + 1; j < dim; j++ ) {
    			
    				if( adjMatrix[i][j] == 0 ){
    					sumPsiH = 0.0;
    					//------------------------------------------------------------------------------------------------- SUM #1
    					for(int k = i; k < j; k++ ){//-----------------------------------------------------------includes all (ij)
    						//Double.isInfinite(sumPsi)//***************************Teste alternativo CONFERIR SE O ATUAL FUNCIONA
    						if( sumPsiH == Double.POSITIVE_INFINITY ) sumPsiH = max_numeric_limits_ld;	
    						if( sumPsiH == Double.NEGATIVE_INFINITY ) sumPsiH = min_numeric_limits_ld;	
    						sumPsiH += ( psi[i][k] * H[k][j] ); 
    					}
    					psi[i][j] = - sumPsiH;

    					//-------------------------------------------------------------------------------------------------SOMA #2
    					if(i > 0)
    						for(int r = 0; r < i; r++){
    							
    							sumPsiHi = 0.0;
    							for(int h = r; h < i + 1; h++){ //This sum includes psi_ri * t_<ii] = psi_ri and psi_rj * t_<jj] = psi_rj, since t_<ii]=t_<jj]=1
    								if( sumPsiHi == Double.POSITIVE_INFINITY ) sumPsiHi = max_numeric_limits_ld;	
    								if( sumPsiHi == Double.NEGATIVE_INFINITY ) sumPsiHi = min_numeric_limits_ld;	
    								sumPsiHi += ( psi[r][h] * H[h][i] );	
   								}
   
    							sumPsiHj = 0.0;
    							for(int h = r; h < j + 1; h++) {
    								if( sumPsiHj == Double.POSITIVE_INFINITY ) sumPsiHj = max_numeric_limits_ld;	
    								if( sumPsiHj == Double.NEGATIVE_INFINITY ) sumPsiHj = min_numeric_limits_ld;
    								sumPsiHj += ( psi[r][h] * H[h][j] );
    							}
    							
   								psi[i][j] -= ( sumPsiHi * sumPsiHj ) / psi[i][i];
    							}

    					if( psi[i][j] == Double.POSITIVE_INFINITY ) psi[i][j] = max_numeric_limits_ld;	
    					if( psi[i][j] == Double.NEGATIVE_INFINITY ) psi[i][j] = min_numeric_limits_ld;	
    				}
    			}
    		}
    	}
        //System.out.println("G = "+Arrays.deepToString(adjMatrix));     
        //System.out.println("psi = "+Arrays.deepToString(psi));     
        //System.out.println("T = "+Arrays.deepToString(T));     
    	
        //phi = psi * T
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {     
                for (int k = 0; k < dim; k++)    
                    phi[i][j] += psi[i][k] * T[k][j];
            }
        }
        //System.out.println("phi = "+Arrays.deepToString(phi));     

        //K = phi' * phi
        for (int i = 0; i < dim; i++) {           
            for (int j = 0; j < dim; j++) {
                for (int k = 0; k < dim; k++)
                    draw[i][j] += phi[k][i] * phi[k][j];   // transpose of 1st element
            }
        }
        //System.out.println("K = "+Arrays.deepToString(draw));     

         
       /* 
        for (int i = 0; i < dim; i++) {           //---------------------------- K
            for (int j = 0; j < dim; j++) {
                for (int k = 0; k < dim; k++)
                   System.out.println("phi_"+k+""+i+" * phi_"+k+""+j);    // transpose of 2nd element
            }
        }
          
        
    	System.out.println("\nTesting the loops against matrix operation methods");
        Matrix psi2 = new Matrix(psi);
        Matrix T2 = new Matrix(T);
    	System.out.println("psi2="+Arrays.deepToString(psi2.toComponents()));
    	System.out.println("T2="+Arrays.deepToString(T2.toComponents()));

        try {
        	Matrix phi2 = psi2.product(T2);
        	System.out.println("phi2="+Arrays.deepToString(phi2.toComponents()));

        	Matrix draw2 = phi2.transposedProduct(phi2);
        	System.out.println("K="+Arrays.deepToString(draw2.toComponents()));

        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Numerical exception in GWishartLikelihood");
        }
  
         */
        
        return draw;
    }

    
    
    
    //==================================================================================================  LOG PDF G-WISHART 
    
    public double logPdf(double[] x) {
    	//System.out.println("----------------------------------Compute G-Wishart logPdf(x)");

        Matrix W = new Matrix(x, dim, dim);
        
        return logPdf(W, DMat, df, dim, logNormalizationConstant);
    }

    public static double logPdf(Matrix W, Matrix D, double df, int dim, double logNormalizationConstant) {
    	//System.out.println("----------------------------------Compute G-Wishart logPdf(W,D,df,dim,IG)");
    	
    	double logDensity = 0;
        //System.out.println("logDensity="+logDensity);

        try {
//            if (!W.isPD()) { // TODO isPD() does not appear to work
//                return Double.NEGATIVE_INFINITY;
//            }

            logDensity = W.logDeterminant(); // Returns NaN is W is not positive-definite.
            //System.out.println("logDensity=logDet(W)="+logDensity);

            if (Double.isInfinite(logDensity) || Double.isNaN(logDensity)) {
                return Double.NEGATIVE_INFINITY;
            }
            
            logDensity *= 0.5;
            //System.out.println("logDensity*0.5="+logDensity);

            logDensity *= df - 2;
            //System.out.println("logDensity*(df-2)="+logDensity);


            // need only diagonal, no? seems a waste to compute
            // the whole matrix
            if (D != null) {
                Matrix product = D.product(W);
                //System.out.println("D%*%W="+Arrays.deepToString(product.toComponents()));
                for (int i = 0; i < dim; i++)
                    logDensity -= 0.5 * product.component(i, i);
            }
            //System.out.println("0.5*D%*%W="+logDensity);

        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
        logDensity += logNormalizationConstant;
        //System.out.println("logPdf="+logDensity+"\n");
        return logDensity;
    }

    
    
    
    
    //=============================================================================================================== TESTS
    
    public static void testBivariateMethod() {

    	
    	//df
        double delta = 3;
        
   
        //S
        double[][] S = new double[][]{{1, 0, 0, 0, 0}, {0, 1, 0, 0, 0},{0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}};
        //double[][] S = new double[][]{{0.47, -0.21, -0.10},{-0.21, 1.05, 0.15},{-0.10, 0.15, 1.84}};
        //double[][] S = new double[][]{{1, 0, 0}, {0, 2, 0}, {0, 0, 3}};
        //double[][] S = new double[][]{{3,-2,-1},{ -2 , 2, 1},{-1 ,  1 ,   0}};
        //System.out.println("S="+Arrays.deepToString(S));

        
        //D
        Matrix SMat= new Matrix(S);
        Matrix DMat= SMat.inverse();
        double[][] D= DMat.toComponents();
        //System.out.println("D="+Arrays.deepToString(D));
                 
        
        //G
        //double[][] G = new double[][]{{0, 1, 1},{1,  0,  1},{1,  1,  0}};        
        //double[][] G = new double[][]{{0, 0, 1},{0,  0,  1},{1,  1,  0}};
        //System.out.println("G="+Arrays.deepToString(G));

        
        //G-Wishart
        //GWishartLikelihood w = new GWishartLikelihood(delta, D, G );
        
        
        
        //g
        //double[] g = new double[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
        //double[] g = new double[]{0, 1, 1};
        //System.out.println("g="+Arrays.toString(g));
        
        int dim = S.length;
        
        Parameter g = new Parameter.Default(dim*(dim-1)/2, 1.0);
        //System.out.println("g="+Arrays.toString(g.getParameterValues()));

        
        
        

   	 	Parameter[] init = new Parameter[dim];
        for (int i = 0; i < dim; i++) {
            init[i] = new Parameter.Default(dim, 0.0);
            init[i].setParameterValue(i, 1.0);
        }
        MatrixParameter precisionParameter = new MatrixParameter(null, init);
       


        for (int i = 0; i < dim; i++) {
            init[i] = new Parameter.Default(dim, 0.002);
            init[i].setParameterValue(i, 0.05);
        }
        MatrixParameter precisionParameter2 = new MatrixParameter(null, init);
       
        //GWishartLikelihood w = new GWishartLikelihood(delta, D, g );
        GWishartLikelihood w2 = new GWishartLikelihood(precisionParameter, delta, S, g, true);


      
      //W        
        //double[] W = new double[]{ 1.00, -0.05, -0.15, -0.05,  0.81,  0.13, -0.15,  0.13,  0.87};
        double[] W = new double[]{0.05, 0.002, 0.002, 0.002, 0.002, 0.002, 0.05, 0.002, 0.002, 0.002, 0.002, 0.002, 0.05, 0.002, 0.002, 0.002, 0.002, 0.002, 0.05, 0.002, 0.002, 0.002, 0.002, 0.002, 0.05};

        
        //logPdf and logNormalizationConstant
        //System.out.println("\nFast logPdf = " + w.logPdf(W));
        //System.out.println("log 1/IG(delta,D) = " + w.logNormalizationConstant);
        System.out.println("\nFast logPdf = " + w2.logPdf(W));
        System.out.println("\n logPdf = " + w2.calculateLogLikelihood());
        System.out.println("log 1/IG(delta,D) = " + w2.logNormalizationConstant);
        
        GWishartLikelihood w3 = new GWishartLikelihood(precisionParameter2, delta, S, g, true);
        w3.getDistribution();
        
        
        
        
        System.out.println("\nFast logPdf = " + w3.logPdf(W));
        System.out.println("\nlogPdf = " + w3.calculateLogLikelihood());
        System.out.println("log 1/IG(delta,D) = " + w3.logNormalizationConstant);
        
        
        System.out.println(w3.getDistribution().calculateLogLikelihood());

       
        //NextGWishart
        double[] g2 = new double[]{0, 0, 1, 1, 1, 1, 1, 1, 1, 1 };
        
        double[][] draw = GWishartLikelihood.nextGWishart(delta, D, SymmetricMatrix.compoundSymmetricMatrix(new double[dim], g2, dim).toComponents());
        //double[][] draw = GWishartLikelihood.nextGWishart(delta, D, G);
        System.out.println("K="+Arrays.deepToString(draw));

       
        //System.out.println("\nadjIndicators="+Arrays.toString(w.getAdjIndicators()));
        //System.out.println("adjMatrix="+Arrays.deepToString(w.getAdjMatrix()));

    }

    public static void main(String[] argv) {
     /*
    	GWishartLikelihood wd = new GWishartLikelihood(2, new double[][]{{500.0}},new double[][]{{0.0}});
        // The above is just an approximation
        GammaDistribution gd = new GammaDistribution(1.0 / 1000.0, 1000.0);
        double[] x = new double[]{1.0};
        System.out.println("Wishart, df=2, scale = 500, PDF(1.0): " + wd.logPdf(x));
        System.out.println("Gamma, shape = 1/1000, scale = 1000, PDF(1.0): " + gd.logPdf(x[0]));

        wd = new GWishartLikelihood(4, new double[][]{{5.0}}, new double[][]{{0.0}});
        gd = new GammaDistribution(2.0, 10.0);
        x = new double[]{1.0};
        System.out.println("Wishart, df=4, scale = 5, PDF(1.0): " + wd.logPdf(x));
        System.out.println("Gamma, shape = 1/1000, scale = 10, PDF(1.0): " + gd.logPdf(x[0]));
        // These tests show the correspondence between a 1D Wishart and a Gamma

      
        
        wd = new GWishartLikelihood(1);
        x = new double[]{0.1};
        System.out.println("Wishart, uninformative, PDF(0.1): " + wd.logPdf(x));
        x = new double[]{1.0};
        System.out.println("Wishart, uninformative, PDF(1.0): " + wd.logPdf(x));
        x = new double[]{10.0};
        System.out.println("Wishart, uninformative, PDF(10.0): " + wd.logPdf(x));
        // These tests show the correspondence between a 1D Wishart and a Gamma
       */  
        
        
        testBivariateMethod();
        
    }
    
    
    
    
    //============================================================================================================== PARSER

    public static XMLObjectParser GWISHART_PRIOR_PARSER = new AbstractXMLObjectParser() { 

    public String getParserName() {
        return GWISHART_PRIOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

    	GWishartLikelihood likelihood;

        if (xo.hasAttribute(NON_INFORMATIVE) && xo.getBooleanAttribute(NON_INFORMATIVE)) {
            // Make non-informative settings
            XMLObject cxo = xo.getChild(DATA);
            int dim = ((MatrixParameter) cxo.getChild(0)).getColumnDimension();
            likelihood = new GWishartLikelihood(dim);
        } else {
            if (!xo.hasAttribute(DF) || !xo.hasChildNamed(RATE_MATRIX)|| !xo.hasChildNamed(ADJ_INDICATORS)) {
                throw new XMLParseException("Must specify df, rateMatrix and adjIndicators");
            }
            
            double df = xo.getDoubleAttribute(DF);

            XMLObject cxo = xo.getChild(RATE_MATRIX);
            MatrixParameter rateMatrix = (MatrixParameter) cxo.getChild(MatrixParameter.class);
            
            cxo = xo.getChild(ADJ_INDICATORS);
            Parameter adjParameter = (Parameter) cxo.getChild(Parameter.class);
            
            
            cxo = xo.getChild(DATA);
            MatrixParameter precisionParameter = (MatrixParameter) cxo.getChild(Parameter.class);

            
            likelihood = new GWishartLikelihood(precisionParameter, df, rateMatrix.getParameterAsMatrix(), adjParameter, false);

        }

        return likelihood;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules;{
        rules = new XMLSyntaxRule[]{
                AttributeRule.newBooleanRule(NON_INFORMATIVE, true),
                AttributeRule.newDoubleRule(DF, true),
                new ElementRule(RATE_MATRIX,
                        new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}, true),
                new ElementRule(ADJ_INDICATORS,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
                new ElementRule(DATA,
                        new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class, 1, Integer.MAX_VALUE)}
                )
        };
    }

    public String getParserDescription() {
        return "Calculates the likelihood of some data (Precision Matrix) under a G-Wishart distribution.";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }
};
}