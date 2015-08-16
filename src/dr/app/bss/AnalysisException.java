package dr.app.bss;

@SuppressWarnings("serial")
public class AnalysisException extends Exception {

	private final String message;
	
	public AnalysisException(String message) {
		
		this.message = message;
		
	}//END: Constructor
	
	@Override
	public String getMessage() {
		return message;
	}
	
}//END: class
