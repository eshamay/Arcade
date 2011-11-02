#include <cstdio>
#include <cstdlib>

int main () {

	FILE * fp = fopen ("wanniers.bin", "rb");

	unsigned int N;
	double vals[3];

	fread(&N, sizeof(unsigned int), 1, fp);
	printf ("N = %d\n", N);

	while (!feof(fp)) {
		for (int i = 0; i < N; i++) {
			fread (vals, sizeof(double), 3, fp);
			printf ("% 15.10lf % 15.10lf % 15.10lf\n", vals[0], vals[1], vals[2]);
		}
	}

	fclose(fp);


	return 0;
}
