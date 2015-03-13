/**
	Tuple.java: a class which represents a tuple, used in sparse matrix.
**/

public class Tuple<X,Y>
{
	private final X _x;		// row 
	private final Y _y;		// column 

	public Tuple(X x, Y y) {
		this._x = x;
		this._y = y;
	}

	public X getX() {
		return _x;
	}

	public Y getY() {
		return _y;
	}

	public boolean
	equals(Tuple<X,Y> t, Tuple<X,Y> s) {
		return t.getX().equals(s.getX()) && t.getY().equals(s.getY());
	}
}

