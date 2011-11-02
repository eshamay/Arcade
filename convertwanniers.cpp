#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>

int main () {

	FILE * fp = fopen("wanniers", "r");
	FILE * wp = fopen ("wanniers.bin", "wb");

	double d;
	unsigned int N;
	bool NSet = false;
	bool NamesSet = false;
	char dump[1000];
	char name[10];
	double x, y, z;

	fscanf (fp, " %d", &N);
	fwrite (&N, sizeof(unsigned int), 1, wp);
	rewind(fp);

	while (!feof(fp)) {
		fgets(dump, 1000, fp);	 // clear the rest of the line of the header
		fgets(dump, 1000, fp);	// skip the 2nd header line
		for (int i = 0; i < N; i++) {
			fscanf (fp, " %*s %lf %lf %lf %*f ", &x, &y, &z);
			fwrite (&x, sizeof(double), 1, wp);
			fwrite (&y, sizeof(double), 1, wp);
			fwrite (&z, sizeof(double), 1, wp);
		}
	}

	fclose(fp); fclose(wp);
	return 0;
}
