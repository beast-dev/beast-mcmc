/*
 * OrderDouble.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.antigenic.phyloclustering.misc.obsolete;

import java.util.Comparator;

public class OrderDouble implements Comparator<OrderDouble>, Comparable<OrderDouble>{
	   private Integer index;
	   private double value;
	   private double value2;
	   
	   OrderDouble(){
	   }

	   OrderDouble(Integer i, double v){
	      index = i;
	      value = v;
	   }
	   
	   OrderDouble(Integer i, double v, double a){
		      index = i;
		      value = v;
		      value2 = a;
		   }

	   public Integer getIndex(){
	      return index;
	   }

	   public double getValue(){
	      return value;
	   }
	   public double getValue2(){
		   return value2;
	   }
	   // Overriding the compareTo method
	   public int compareTo(OrderDouble d){
	      return (this.index).compareTo(d.index);
	   }

	   // Overriding the compare method to sort the value
	   //Returns a negative integer, zero, or a positive integer as the first argument is less than, equal to, or greater than the second.
	   public int compare(OrderDouble d, OrderDouble d1){
		  if(d.value - d1.value > 0){
	      return 1;
		  }
		  else if(d.value - d1.value < 0){
			  return -1;
		  }
		  else{
			  return 0;
		  }
	   }
	}
