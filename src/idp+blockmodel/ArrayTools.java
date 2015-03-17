import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;

public class ArrayTools
{
	public static class ValueComparator<K, V extends Comparable<V>> implements Comparator<K> {
		Map<K, V> map;
 
		public ValueComparator(Map<K, V> map) {
			this.map = map;
		}

		@Override
		public int compare(K keyA, K keyB) {
			V valueA = map.get(keyA);
			Comparable<V> valueB = map.get(keyB);
			if (valueB.equals(valueA)) {
				return 1;					// keys with same value will not be merged
			//	return keyB.compareTo(keyA);			// keys with same value will be merged
			}
			else {
				return valueB.compareTo(valueA);		// sort descending
			}
		}

		public static <K, V extends Comparable<V>> Map<K, V> sortByValue(Map<K, V> unsortedMap) {
			Map<K, V> sortedMap = new TreeMap<K, V>(new ValueComparator<K, V>(unsortedMap));
			sortedMap.putAll(unsortedMap);

			return sortedMap;
		}
	}
} 
