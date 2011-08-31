/*
 * Assert.java
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

package dr.util;



/**
 * Simple Static class for assertions. This may be replaced in future
 * with Java 1.4's assertion mechanism.
 *
 * @version $Id: Assert.java,v 1.7 2005/07/08 12:15:50 rambaut Exp $
 *
 * @author Andrew Rambaut
 */
//public class Assert {
//
//	/* not sure about exactly what this function should display */
//	private static void fail(String msg)
//	{
//		System.err.println("Assertion failed: " + msg);
//
//		Throwable e = new Throwable();
//		e.printStackTrace();
//	}
//
//	public static void test(boolean b)
//	{
//		if (!b)
//			fail("Statement false");
//	}
//
//	public static void test(long lng)
//	{
//		if (lng == 0L)
//			fail("Zero integer");
//	}
//
//	public static void test(double dbl)
//	{
//		if (dbl == 0.0)
//			fail("Zero double");
//	}
//
//	public static void test(Object ref)
//	{
//		if (ref == null)
//			fail("Null reference");
//	}
//}
