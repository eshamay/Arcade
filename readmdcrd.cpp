#include <cstdio>
#include <cstdlib>

void ReadText () {

	float d, total = 0., count = 0.;
	FILE * fp = fopen ("data.dat", "r");

	while (fscanf (fp, " %f ", &d) != EOF) {
		total += d;
		count++;
	}
	fclose(fp);

	printf ("total = %f\ncount = %f\n", total, count);
}

void ReadBinary () {

	float d, total = 0., count = 0.;
	FILE * fp = fopen64 ("mdcrd.bin", "rb");

	while (!feof(fp)) {
		fread(&d, sizeof(float), 1, fp);
		printf ("%f\n", d);
		total += d;
		count++;
	}
	fclose(fp);

	printf ("total = %f\ncount = %f\n", total, count);
}

int main () {

	ReadBinary ();
	//ReadText ();

	return 0;
}
