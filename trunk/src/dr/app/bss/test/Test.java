package dr.app.bss.test;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import dr.app.bss.PartitionData;
import dr.app.bss.Utils;

public class Test {

	public static void main(String[] args) {

		serialize();
		deserialize();

	}// END: main

	private static void serialize() {

		System.out.println("Serializing data:");

		PartitionData data = new PartitionData();
		Utils.printPartitionData(data);

		try {

			FileOutputStream fileOut = new FileOutputStream("data.ser");
			ObjectOutputStream out = new ObjectOutputStream(fileOut);
			out.writeObject(data);
			out.close();
			fileOut.close();

		} catch (IOException e) {
			e.printStackTrace();
		}// END: try-catch block

	}// END: serialize

	private static void deserialize() {

		System.out.println("Deserializing data:");

		PartitionData data = null;

		try {

			FileInputStream fileIn = new FileInputStream("data.ser");
			ObjectInputStream in = new ObjectInputStream(fileIn);
			data = (PartitionData) in.readObject();
			in.close();
			fileIn.close();

			Utils.printPartitionData(data);

		} catch (IOException e1) {
			e1.printStackTrace();
		} catch (ClassNotFoundException e1) {
			e1.printStackTrace();
		}// END: try-catch block

	}// END: deserialize

}// END: class
