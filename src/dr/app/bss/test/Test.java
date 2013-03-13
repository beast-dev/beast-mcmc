package dr.app.bss.test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import dr.app.bss.PartitionData;
import dr.app.bss.PartitionDataList;
import dr.app.bss.Utils;

public class Test {

	public static void main(String[] args) {

		serialize();
		System.out.println();
		deserialize();

	}// END: main

	private static void serialize() {

		System.out.println("Serializing data:");

		PartitionDataList dataList = new PartitionDataList();

		PartitionData data1 = new PartitionData();
		data1.treeFile = new File("/home/filip/SimTree.figtree");
		dataList.add(data1);
		
		PartitionData data2 = new PartitionData();
		dataList.add(data2);
		
		Utils.printPartitionDataList(dataList);

		try {

			FileOutputStream fileOut = new FileOutputStream("data.ser");
			ObjectOutputStream out = new ObjectOutputStream(fileOut);
			out.writeObject(dataList);
			out.close();
			fileOut.close();

		} catch (IOException e) {
			e.printStackTrace();
		}// END: try-catch block

	}// END: serialize

	private static void deserialize() {

		System.out.println("Deserializing data:");

		PartitionDataList dataList  = null;

		try {

			FileInputStream fileIn = new FileInputStream("data.ser");
			ObjectInputStream in = new ObjectInputStream(fileIn);
			dataList = (PartitionDataList) in.readObject();
			in.close();
			fileIn.close();

			Utils.printPartitionDataList(dataList);

		} catch (IOException e1) {
			e1.printStackTrace();
		} catch (ClassNotFoundException e1) {
			e1.printStackTrace();
		}// END: try-catch block

	}// END: deserialize

}// END: class
