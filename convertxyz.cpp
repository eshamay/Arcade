#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <iostream>

int main () {

	FILE * fp = fopen("xyz", "r");
	FILE * wp = fopen ("xyz.bin", "wb");
	std::cout << "Covnerting file xyz to xyz.bin" << std::endl;

	double d;
	unsigned int N;
	bool NSet = false;
	bool NamesSet = false;
	char dump[1000];
	char name[10];
	double x, y, z;

	bool set = false;
	int num = 0;
	while (!feof(fp)) {
		fscanf (fp, " %d", &N);

		if (!set) {
			num = N;
			std::cout << "Parsing " << N << " atoms from the frame" << std::endl;
			set = true;
		}

		if (num != N) {
			std::cout << "frame has " << N << " atoms instead of the expected " << num << " atoms" << std::endl;
		}

		// grab the header with the number of atoms to parse
		if (!NSet) {
			fwrite (&N, sizeof(unsigned int), 1, wp);
			fgets(dump, 1000, fp);	 // clear the rest of the line of the header
			fgets(dump, 1000, fp);	// skip the 2nd header line

			// grab all the atom names
			unsigned int len;
			for (int i = 0; i < N; i++) {
				fscanf (fp, " %s ", name);

				len = (unsigned int)strlen(name);
				fwrite (&len, sizeof(unsigned int), 1, wp);

				fwrite (name, sizeof(char), len+1, wp);

				//printf ("len %s = %u\n", name, len);
				fgets(dump, 1000, fp);	// skip the 2nd header line
			}
			rewind(fp);

			NSet = true;
		}

		fgets(dump, 1000, fp);	 // clear the rest of the line of the header
		fgets(dump, 1000, fp);	// skip the 2nd header line
		
		for (int i = 0; i < N; i++) {
			fscanf (fp, " %s %lf %lf %lf", name, &x, &y, &z);
			fwrite (&x, sizeof(double), 1, wp);
			fwrite (&y, sizeof(double), 1, wp);
			fwrite (&z, sizeof(double), 1, wp);
		}
	}

	fclose(fp); fclose(wp);
	return 0;
}
