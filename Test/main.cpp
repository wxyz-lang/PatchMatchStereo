#include <vector>
#include <list>
#include <algorithm>

using namespace std; 

int main()
{
	vector<int> vec(5, -1);
	for (int i = 0; i < 5; i++) {
		printf("%d\n", vec[i]);
	}

	return 0;
}