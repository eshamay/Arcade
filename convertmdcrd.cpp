#include <cstdlib>
#include <cstdio>

int main () {

	FILE * fp = fopen ("mdcrd", "r");
	FILE * wp = fopen ("mdcrd.bin", "wb");
	float d;
	while (fscanf (fp, " %f ", &d) != EOF) {
		fwrite(&d, sizeof(float), 1, wp);
	}

	fclose(fp);
	fclose(wp);

	return 0;
}
