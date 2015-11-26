/*
 * ODEDemographicFunction.java
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

/*
 * Base class for demographic models based on numerically-integrated ODEs
 * 
*/
package dr.evomodel.epidemiology;

import dr.evolution.coalescent.*;


/**
 * This interface provides methods that describe a demographic function.
 * @author Daniel Wilson
 */

/*public interface ODEDemographicFunction extends DemographicFunction {
	public abstract class Abstract implements ODEDemographicFunction {
		public Abstract(Type units) {
			setUnits(units);
		}
        public double getLogDemographic(double t) {
            return Math.log(getDemographic(t));
        }
		public double getIntegral(double start, double finish)
		{
			return getIntensity(finish) - getIntensity(start);
		}

	}
}*/

public abstract class ODEDemographicFunction extends DemographicFunction.Abstract {
	public ODEDemographicFunction(Type units) {
		super(units);
	}
	
	// Implement abstract types from base class
	
	/**
     * Default implementation
     * @param t
     * @return log(demographic at t)
     */
    /*public double getLogDemographic(double t) {
        return Math.log(getDemographic(t));
    }*/

    /**
	 * Calculates the integral 1/N(x) dx between start and finish.
	 */
	/*public double getIntegral(double start, double finish)
	{
		return getIntensity(finish) - getIntensity(start);
	}*/

    /**
     * Returns the integral of 1/N(x) between start and finish, calling either the getAnalyticalIntegral or
     * getNumericalIntegral function as appropriate.
     */
	public double getNumericalIntegral(double start, double finish) {
		throw new RuntimeException("not implemented");
	}
	
	public double getDemographic(double t) {
		Evaluate(t);
		if(RKfail) return 0.0;
		return getDemographicFromPrevalence(Ynow,t);
	}
	
	public double getIntensity(double t) {
		Evaluate(t);
		if(RKfail) return Math.log(0.0);
		return Ynow[0];
	}
	
	// Implement the following abstract functions:
	
	/**
	 * Calculate the derivatives and store in dydt
	 */
	abstract void derivs(double t, double[] y, double[] dydt);
	
	/**
	 * Set initial values of y
	 * for(i=0;i<nvar;i++) Y[i][0] = ...;
	 */
	abstract void setInit();
	
	/**
	 * Calculate the effective population size from the prevalence, contained in y
	 * @return
	 */
	abstract double getDemographicFromPrevalence(double[] y, double t);
	
	// Implemented base functions
	
	/**
	 * Evaluate the demographic functions at time Tnow.
	 * Store results in Ynow for immediate use by getDemographic or getIntensity
	 */
	void Evaluate(double t) {
		if(RKfail) return;
		// if (t==Tnow) ???
		if(t<0.0) throw new RuntimeException("t cannot be negative");
		int i;
		if(klast==-1 || t>T[klast]) {
			if(klast+1>=kmax) {
				// 3rd December 2011: Previous behaviour was stupid, now increase storage capacity
				// Assume linear, but flag warning
				//RKwarning = true;
				//for(i=0;i<nvar;i++) Ynow[i] = Y[i][klast];
				//Tnow = t;
				//return;
				RKresize();
			}
			// Continue integration
			try {
				RungeKutta(t);
			} catch(RuntimeException e) {
				System.err.println(e.getMessage());
				RKfail = true;
				return;
			}
			for(i=0;i<nvar;i++) Ynow[i] = Y[i][klast];
			Tnow = t;
			return;
		}
		if(t==T[klast]) {
			for(i=0;i<nvar;i++) Ynow[i] = Y[i][klast];
			Tnow = t;
			return;
		}
		// Linearly interpolate
		int k1 = (int)Math.floor(t/dtsav);
		if(k1>klast-1) k1 = klast-1;
		if(k1<0) k1 = 0;
		int k2 = k1+1;
		if(k2>klast || T[k1]>t) {
			while(k2>klast || T[k1]>t) {
				--k1;
				--k2;
			}
		}
		else if(k1<0 || T[k2]<=t) {
			while(k1<0 || T[k2]<=t) {
				++k1;
				++k2;
			}
		}
		if(T[k1]==t) {
			for(i=0;i<nvar;i++) Ynow[i] = Y[i][k1];
		}
		else {
			// Linearly interpolate
			for(i=0;i<nvar;i++) Ynow[i] = Y[i][k1] + (t-T[k1])*(Y[i][k2]-Y[i][k1])/(T[k2]-T[k1]);
		}
		Tnow = t;		
	}

	/**
	 * Initialize RK integration
	 */
	void RKinit() {
		klast = -1;
		RKwarning = false;
		RKfail = false;
		if(Y==null || Y.length!=nvar || Y[0].length!=kmax) {
			// Should be no memory leak
			Y = new double[nvar][kmax];
			T = new double[kmax];
			Ynow = new double[nvar];
			ak2 = new double[nvar];
			ak3 = new double[nvar];
			ak4 = new double[nvar];
			ak5 = new double[nvar];
			ak6 = new double[nvar];
			Ytemp = new double[nvar];
			Yerr = new double[nvar];
			y = new double[nvar];
			dydt = new double[nvar];
			yscal = new double[nvar];
		}
		nok = nbad = 0;
	}
	
	/**
	 * Increase kmax by kinc on the hoof
	 */
	void RKresize() {
		if(Y==null) throw new RuntimeException("Y not yet allocated");
		if(T==null) throw new RuntimeException("T not yet allocated");
		if(kmax==kabsolutemax) throw new RuntimeException("kabsolutemax exceeded");
		// Store old value of kmax before incrementing it
		int oldkmax = kmax;
		kmax += kinc;
		// Temporary pointers to old arrays, should get garbage collected
		double oldY[][] = Y;
		double oldT[] = T;
		// Allocate new memory for enlarged arrays
		Y = new double[nvar][kmax];
		T = new double[kmax];
		// Copy across old values
		int i,j;
		for(i=0;i<oldkmax;i++) {
			for(j=0;j<nvar;j++) {
				Y[j][i] = oldY[j][i];
			}
			T[i] = oldT[i];
		}
	}

	/**
	 * Flag if kmax is exceeded
	 * @return
	 */
	public boolean RKwarn() {
		return RKwarning;
	}
	
	/**
	 * Driving routine for RungeKutta integration
	 */
	void RungeKutta(double t2) {
//		if(h1<0) throw new RuntimeException("h1 must be positive");
		if(klast==kmax-1) throw new RuntimeException("storage space is exceeded");
		if(klast==-1) {
			setInit();	// virtual function, over-ride in derived class
			T[0] = 0;
			++klast;
		}
		if(t2==0.0) return;

		int i,nstp;
		double t1 = T[klast];		// beginning of time range
		double t = t1; 				// current time
		double tsav = t1;			// time of last saved point
		// Copy initial y values
		for (i=0;i<nvar;i++) y[i]=Y[i][klast];
		// Step size
		double h=hinit; // NB h must be positive
			
		for(nstp=0;nstp<MAXSTP;nstp++) {
			double tmp = nstp;
			if(nstp>MAXSTP/2) {
				tmp = tmp+3;
			}
			// Calculate the derivatives at the present time
			derivs(t,y,dydt);
			// Calculate appropriate scalings for the error tolerance
			for(i=0;i<nvar;i++)
				yscal[i]=Math.abs(y[i])+Math.abs(dydt[i]*h)+TINY;
			// Storage (ensure there's room for final state)
			if(klast < kmax-2 && Math.abs(t-tsav) > Math.abs(dtsav)) {
				++klast;
				for (i=0;i<nvar;i++) Y[i][klast]=y[i];
				T[klast]=t;
				tsav=t;
			}
			// Reduce step size to avoid over-stepping target time t2
			if((t+h-t2)*(t+h-t1) > 0.0) h=t2-t;
			// Perform the adaptive integration and update time
			t = rkqs(y,dydt,t,h,yscal);
			// Was the predicted step size used?
			if(hdid == h) ++nok; else ++nbad;
			// If the target t2 has been reached
			if((t-t2)*(t2-t1) >= 0.0) {
				// Store state at t2
				if (klast < kmax-1) {
					++klast;
					for (i=0;i<nvar;i++) Y[i][klast]=y[i];
					T[klast]=t;
				}
				return;
			}
			if(Math.abs(hnext) <= hmin) {
				throw new RuntimeException("Step size too small in odeint");
			}
			h=hnext;
		}
		throw new RuntimeException("Too many steps in routine odeint");
	}
		
	/**
	 * Take an adaptive step
	 *  
	 *  @return New time
	 */
	// Cannot output primitive scalars by modifying arguments. So hdid and hnext become member variables and t is returned.
	double rkqs(double[] y, double[] dydt, double t, double htry, double[] yscal) {
		int i;
		double errmax,h,htemp,tnew;

		h=htry;
		for(;;) {
			rkck(y,dydt,t,h);
			errmax=0.0;
			for(i=0;i<nvar;i++) errmax=Math.max(errmax,Math.abs(Yerr[i]/yscal[i]));
			errmax /= eps;
			if(errmax <= 1.0) break;
			if(Double.isNaN(errmax)) {
				throw new RuntimeException("errmax NaN");
				// Invalid value of one or more of the variables was chosen. Halve the step size
				//h /= 2.0;
			} else {
				htemp=SAFETY*h*Math.pow(errmax,PSHRNK);
				h=(h >= 0.0 ? Math.max(htemp,0.1*h) : Math.min(htemp,0.1*h));
			}
			tnew=t+h;
			if(tnew == t) {
				throw new RuntimeException("stepsize underflow in rkqs");
			}
		}
		if(errmax > ERRCON) hnext=SAFETY*h*Math.pow(errmax,PGROW);
		else hnext=5.0*h;
		t += (hdid=h);
		for(i=0;i<nvar;i++) y[i]=Ytemp[i];
		return t;
	}

	/**
	 * Take one RK5 step
	 */
	// NB: Arrays such as ak2 and yerr are objects. Pointers to these objects are passed by
	// value into rkck. The whole array is not copied. So this should work...
	void rkck(double[] y, double[] dydt, double t, double h) {
		int i;

		for(i=0;i<nvar;i++)
			Ytemp[i]=y[i]+b21*h*dydt[i];
		derivs(t+a2*h,Ytemp,ak2);
		for(i=0;i<nvar;i++)
			Ytemp[i]=y[i]+h*(b31*dydt[i]+b32*ak2[i]);
		derivs(t+a3*h,Ytemp,ak3);
		for(i=0;i<nvar;i++)
			Ytemp[i]=y[i]+h*(b41*dydt[i]+b42*ak2[i]+b43*ak3[i]);
		derivs(t+a4*h,Ytemp,ak4);
		for(i=0;i<nvar;i++)
			Ytemp[i]=y[i]+h*(b51*dydt[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
		derivs(t+a5*h,Ytemp,ak5);
		for(i=0;i<nvar;i++)
			Ytemp[i]=y[i]+h*(b61*dydt[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
		derivs(t+a6*h,Ytemp,ak6);
		for(i=0;i<nvar;i++)
			Ytemp[i]=y[i]+h*(c1*dydt[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
		for(i=0;i<nvar;i++)
			Yerr[i]=h*(dc1*dydt[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
	}
	
	
	// Member variables
	
	protected int nvar = 0;	//	Default to zero
	
	// Runge-Kutta integration variables
	protected int kmax=200;			//	Maximum storage capacity for integration
	protected int kinc=200;			//	Increment size for kmax when it is exceeded
	protected int kabsolutemax=200000; //	Absolute maximum storage capacity for integration allowable
	protected int klast=-1;			//	Index of the last evaluation point
	protected double hinit=0.1;		//	Initial suggested step size for RungeKutta integration
	protected double[][] Y;		//	Storage for integration results. Y[0] must always contain Lambda, the integrated intensity function
	protected double[] T;			//	Storage for evaluated time points: interpolate between
	protected boolean RKwarning;	//	Warn if kmax is exceeded
	
	// Temporary variables
	protected double[] Ynow;		//	Instead of passing and returning vectors, store immediate value of Y
	protected double Tnow;		//	Immediate value of T
	
	// Static constants and storage used by rkck
	static final double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
	b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42 = -0.9,
	b43=1.2, b51 = -11.0/54.0, b52=2.5, b53 = -70.0/27.0,
	b54=35.0/27.0, b61=1631.0/55296.0, b62=175.0/512.0,
	b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0,
	c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,
	dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,
	dc4=c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6=c6-0.25;
	private double[] ak2, ak3, ak4, ak5, ak6;
	
	// Static constants used by rkqs
	static final double SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;

	// Storage used by rkqs and RungeKutta
	protected double hdid, hnext;
	
	// Storage used by rkqs and rkck
	protected double[] Ytemp, Yerr;
	
	// Storage for RungeKutta
	protected int MAXSTP=10000;
	protected double TINY=1.0e-30;
	protected int nok=0, nbad=0;
	protected double[] y, dydt, yscal;
	protected double dtsav = 0.1;
	protected double hmin = 1.0e-16;
	protected double eps = 1e-4;
	protected boolean RKfail = false;
}

