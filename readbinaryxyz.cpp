#include <cstdio>
#include <cstdlib>

int main () {

	FILE * fp = fopen ("xyz.bin", "rb");

	unsigned int N;
	double vals[3];
	char name[10];

	fread(&N, sizeof(unsigned int), 1, fp);
	printf ("N = %d\n", N);
	unsigned int len;
	for (int i = 0; i < N; i++) {
		fread (&len, sizeof(unsigned int), 1, fp);
		fread (name, sizeof(char), len+1, fp);
		printf ("%s ", name);
	}
	printf ("\n");

	while (!feof(fp)) {
		for (int i = 0; i < N; i++) {
			fread (vals, sizeof(double), 3, fp);
				printf ("% 15.10lf % 15.10lf % 15.10lf\n", vals[0], vals[1], vals[2]);
		}
	}

	fclose(fp);


	return 0;
}
