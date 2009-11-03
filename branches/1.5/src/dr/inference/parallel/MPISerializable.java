package dr.inference.parallel;

/**
 * @author Marc A. Suchard
 */
public interface MPISerializable {

	public void sendState(int toRank);

	public void receiveState(int fromRank);

}
