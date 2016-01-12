/*
 * CoalescentEventsStatistic.java
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

package dr.evomodel.coalescent;

//import dr.evolution.coalescent.IntervalType;
import dr.inference.model.Statistic;

/**
 * @author Guy Baele
 */
public class CoalescentEventsStatistic extends Statistic.Abstract {

    //public static final boolean TEST_NEW_CODE = true;

    public static final boolean DEBUG = false;

    //private static final boolean FULL_FINAL_INTERVAL = true;
    //private static final boolean LOG_COMBINATIONS = true;
    //private static final boolean RETURN_RECIPROCAL = false;

    private final CoalescentIntervalProvider coalescent;
    //treeModel currently only required for debugging purposes
    //private TreeModel treeModel;
    private int coalescentEvents;
    //private double[] coalescentValues;

    //public CoalescentEventsStatistic(CoalescentIntervalProvider coalescent, TreeModel treeModel) {
    public CoalescentEventsStatistic(CoalescentIntervalProvider coalescent) {

        this.coalescent = coalescent;

        //if (TEST_NEW_CODE) {

        this.coalescentEvents = coalescent.getNumberOfCoalescentEvents();

        /*} else {

            this.coalescentEvents = 0;
            if (coalescent instanceof GMRFSkyrideLikelihood) {
                this.coalescentEvents = coalescent.getCoalescentIntervalDimension();
                if (DEBUG) {
                    System.err.println("CoalescentIntervalDimension: " + coalescent.getCoalescentIntervalDimension());
                }
            } else {
                for (int i = 0; i < coalescent.getCoalescentIntervalDimension(); i++) {
                    //Not yet implemented for the skygrid model
                    if (coalescent.getCoalescentIntervalType(i) == IntervalType.COALESCENT) {
                        coalescentEvents++;
                    }
                }
            }
            if (DEBUG) {
                System.err.println("Number of coalescent events: " + this.coalescentEvents);
            }

        }*/

        //this.coalescentValues = new double[coalescentEvents];
        if (DEBUG) {
            System.err.println("CoalescentEventsStatistic dimension: " + this.coalescentEvents);
        }
    }

    public int getDimension() {
        return this.coalescentEvents;
    }

    public double getStatisticValue(int i) {

        //if (TEST_NEW_CODE) {

            return coalescent.getCoalescentEventsStatisticValue(i);

        /*} else {

            //System.err.println(treeModel);
            //i will go from 0 to getDimension()
            //GMRFSkyrideLikelihood
            if (DEBUG) {
                System.err.println("getStatisticValue(int i)");
            }
            if (i == 0) {
                if (DEBUG) {
                    System.err.println("coalescentValues.length = " + coalescentValues.length);
                }
                //reset array of coalescent events
                for (int j = 0; j < coalescentValues.length; j++) {
                    coalescentValues[j] = 0.0;
                }
                //recalculate everything
                int counter = 0;
                if (DEBUG) {
                    System.err.println("coalescent.getCoalescentIntervalDimension() = " + coalescent.getCoalescentIntervalDimension());
                }
                for (int j = 0; j < coalescent.getCoalescentIntervalDimension(); j++) {
                    if (coalescent instanceof GMRFMultilocusSkyrideLikelihood) {
                        if (coalescent.getCoalescentIntervalType(j) == IntervalType.COALESCENT) {
                            if (LOG_COMBINATIONS) {
                                this.coalescentValues[counter] += coalescent.getCoalescentInterval(j)*(coalescent.getCoalescentIntervalLineageCount(j)*(coalescent.getCoalescentIntervalLineageCount(j)-1.0))/2.0;
                            } else {
                                this.coalescentValues[counter] += coalescent.getCoalescentInterval(j);
                            }
                            counter++;
                        } else if (!FULL_FINAL_INTERVAL) {
                            if (coalescent.getCoalescentIntervalType(j) == IntervalType.SAMPLE && counter != 0) {
                                if (LOG_COMBINATIONS) {
                                    this.coalescentValues[counter] += coalescent.getCoalescentInterval(j)*(coalescent.getCoalescentIntervalLineageCount(j)*(coalescent.getCoalescentIntervalLineageCount(j)-1.0))/2.0;
                                } else {
                                    this.coalescentValues[counter] += coalescent.getCoalescentInterval(j);
                                }
                            }
                        } else {
                            if (coalescent.getCoalescentIntervalType(j) == IntervalType.SAMPLE) {
                                if (LOG_COMBINATIONS) {
                                    this.coalescentValues[counter] += coalescent.getCoalescentInterval(j)*(coalescent.getCoalescentIntervalLineageCount(j)*(coalescent.getCoalescentIntervalLineageCount(j)-1.0))/2.0;
                                } else {
                                    this.coalescentValues[counter] += coalescent.getCoalescentInterval(j);
                                }
                            }
                        }
                    } else if (coalescent instanceof GMRFSkyrideLikelihood) {
                        if (DEBUG) {
                            System.err.println("counter = " + counter);
                            System.err.println("((GMRFSkyrideLikelihood)coalescent).getSufficientStatistics()[" + j + "] = " + ((GMRFSkyrideLikelihood)coalescent).getSufficientStatistics()[j]);
                        }
                        this.coalescentValues[counter] = ((GMRFSkyrideLikelihood)coalescent).getSufficientStatistics()[j];
                        counter++;
                    } else {
                        //System.err.println(coalescent.getCoalescentIntervalType(j) + "   " + coalescent.getCoalescentInterval(j));
                        if (coalescent.getCoalescentIntervalType(j) == IntervalType.COALESCENT) {
                            if (LOG_COMBINATIONS) {
                                this.coalescentValues[counter] += coalescent.getCoalescentInterval(j)*(coalescent.getCoalescentIntervalLineageCount(j)*(coalescent.getCoalescentIntervalLineageCount(j)-1.0))/2.0;
                                //System.err.println("interval length: " + coalescent.getCoalescentInterval(j));
                                //System.err.println("lineage count: " + coalescent.getCoalescentIntervalLineageCount(j));
                                //System.err.println("factorial: " + (coalescent.getCoalescentIntervalLineageCount(j)*coalescent.getCoalescentIntervalLineageCount(j)-1.0)/2.0);
                                //System.err.println("counter " + counter + ": " + this.coalescentValues[counter] + "\n");
                                //this.coalescentValues[counter] += coalescent.getCoalescentInterval(j);
                                //this.coalescentValues[counter] = (coalescent.getCoalescentIntervalLineageCount(j)*coalescent.getCoalescentIntervalLineageCount(j)-1.0)/(2.0*this.coalescentValues[counter]);
                            } else {
                                this.coalescentValues[counter] += coalescent.getCoalescentInterval(j);
                            }
                            counter++;
                        } else if (!FULL_FINAL_INTERVAL) {
                            if (coalescent.getCoalescentIntervalType(j) == IntervalType.SAMPLE && counter != 0) {
                                if (LOG_COMBINATIONS) {
                                    this.coalescentValues[counter] += coalescent.getCoalescentInterval(j)*(coalescent.getCoalescentIntervalLineageCount(j)*(coalescent.getCoalescentIntervalLineageCount(j)-1.0))/2.0;
                                    //System.err.println("interval length: " + coalescent.getCoalescentInterval(j));
                                    //System.err.println("lineage count: " + coalescent.getCoalescentIntervalLineageCount(j));
                                    //System.err.println("factorial: " + (coalescent.getCoalescentIntervalLineageCount(j)*coalescent.getCoalescentIntervalLineageCount(j)-1.0)/2.0);
                                    //System.err.println("counter " + counter + ": " + this.coalescentValues[counter] + "\n");
                                    //this.coalescentValues[counter] += coalescent.getCoalescentInterval(j);
                                } else {
                                    this.coalescentValues[counter] += coalescent.getCoalescentInterval(j);
                                }
                            }
                        } else {
                            if (coalescent.getCoalescentIntervalType(j) == IntervalType.SAMPLE) {
                                if (LOG_COMBINATIONS) {
                                    //System.err.println("interval length: " + coalescent.getCoalescentInterval(j));
                                    //System.err.println("lineage count: " + coalescent.getCoalescentIntervalLineageCount(j));
                                    //System.err.println("factorial: " + (coalescent.getCoalescentIntervalLineageCount(j)*coalescent.getCoalescentIntervalLineageCount(j)-1.0)/2.0);
                                    //System.err.println("counter " + counter + ": " + this.coalescentValues[counter] + "\n");
                                    this.coalescentValues[counter] += coalescent.getCoalescentInterval(j)*(coalescent.getCoalescentIntervalLineageCount(j)*(coalescent.getCoalescentIntervalLineageCount(j)-1.0))/2.0;
                                    //this.coalescentValues[counter] += coalescent.getCoalescentInterval(j);
                                } else {
                                    this.coalescentValues[counter] += coalescent.getCoalescentInterval(j);
                                }
                            }
                        }
                    }
                }
            }
            //for (int j = 0; j < this.coalescentEvents; j++) {
            //	System.err.println(this.coalescentValues[j]);
            //}
            //System.exit(0);
            if (RETURN_RECIPROCAL) {
                return 1.0/this.coalescentValues[i];
            } else {
                return this.coalescentValues[i];
            }

        }*/

    }

    public String getStatisticName() {
        return "coalescentEventsStatistic";
    }

}
