/*
 * MyCounter.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.app.simcoal;

/**
 * Created by IntelliJ IDEA.
 * User: adru001
 * Date: Jun 14, 2006
 * Time: 4:58:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class MyCounter{
	private static int y=3;
	private int x;
	public static void main(String[] args){
		MyCounter c1 = new MyCounter();
		MyCounter c2 = new MyCounter();
		c1.y=++y;
		c2.y=y++;

        c1.x=1;
		c2.x=5;
		int x = 0;
		x = x++;
        System.out.println("c1.y = " + c1.y);
        System.out.println("c1.x = " + c1.x);
        System.out.println("c2.y = " + c2.y);
        System.out.println("c2.x = " + c2.x);
        System.out.println("Y: " + y);
        System.out.println("X: " + x);
	}
}
