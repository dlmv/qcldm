#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_scfout.h"

int main(int argc, char *argv[]) 
{

	int ct_AN,h_AN,Gh_AN,i,j,TNO1,TNO2,l,m,n, type;	
	int spin,Rn;
	double x,y,z, a, b, c;

	static FILE *f;

	read_scfout(argv);

	/*check cells*/
	int max_l = 0;
	int min_l = 0;
	int max_m = 0;
	int min_m = 0;
	int max_n = 0;
	int min_n = 0;		
	for (Rn=0; Rn<=TCpyCell; Rn++){
			l = atv_ijk[Rn][1];
			m = atv_ijk[Rn][2];
			n = atv_ijk[Rn][3];
			if (l > max_l) {max_l = l;}
			if (l < min_l) {min_l = l;}
			if (m > max_m) {max_m = m;}
			if (m < min_m) {min_m = m;}
			if (n > max_n) {max_n = n;}
			if (n < min_n) {min_n = n;}
	}
	
	int d_l = max_l - min_l;
	int d_m = max_m - min_m;
	int d_n = max_n - min_n;

	/*init neighbours*/
	
	int *****NB = (int*****)malloc(sizeof(int****)*(atomnum+1));
	
	int AN, AN1;
	for (AN=1; AN<=atomnum; AN++){
		NB[AN] = (int****)malloc(sizeof(int***)*(atomnum+1));
		for (AN1=1; AN1<=atomnum; AN1++){
			NB[AN][AN1] = (int***)malloc(sizeof(int**)*(d_l+1));
			for (l=0; l<=d_l; l++){
				NB[AN][AN1][l] = (int**)malloc(sizeof(int*)*(d_m+1));
				for (m=0; m<=d_m; m++){
					NB[AN][AN1][l][m] = (int*)malloc(sizeof(int)*(d_n+1));
					for (n=0; n<=d_n; n++){
						NB[AN][AN1][l][m][n] = -1;
					}
				}
			}
		}
	}

	int ****allAtoms = (int****)malloc(sizeof(int***)*(atomnum+1));
	for (AN1=1; AN1<=atomnum; AN1++){
		allAtoms[AN1] = (int***)malloc(sizeof(int**)*(d_l+1));
		for (l=0; l<=d_l; l++){
			allAtoms[AN1][l] = (int**)malloc(sizeof(int*)*(d_m+1));
			for (m=0; m<=d_m; m++){
				allAtoms[AN1][l][m] = (int*)malloc(sizeof(int)*(d_n+1));
				for (n=0; n<=d_n; n++){
					allAtoms[AN1][l][m][n] = 0;
				}
			}
		}
	}

	/*load neighbours*/
	for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {
		for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
			Gh_AN = natn[ct_AN][h_AN];
			Rn = ncn[ct_AN][h_AN];
			l = atv_ijk[Rn][1] - min_l;
			m = atv_ijk[Rn][2] - min_l;
			n = atv_ijk[Rn][3] - min_l;
			NB[ct_AN][Gh_AN][l][m][n] = h_AN;
			allAtoms[Gh_AN][l][m][n] = 1;
		}
	}

	/*counting all atoms*/
	int allAtomNum = 0;
	for (AN1=1; AN1<=atomnum; AN1++){
		for (l=0; l<=d_l; l++){
			for (m=0; m<=d_m; m++){
				for (n=0; n<=d_n; n++){
					allAtomNum += allAtoms[AN1][l][m][n];
				}
			}
		}
	}

	/*loading all atoms*/
	int allAtomTypes[allAtomNum], allAtomL[allAtomNum], allAtomM[allAtomNum], allAtomN[allAtomNum];
	i = 0;
	for (AN1=1; AN1<=atomnum; AN1++){
		allAtomTypes[i] = AN1;
		allAtomL[i] = - min_l;
		allAtomM[i] = - min_m;
		allAtomN[i] = - min_n;
		++i;
		for (l=0; l<=d_l; l++){
			for (m=0; m<=d_m; m++){
				for (n=0; n<=d_n; n++){
					if (allAtoms[AN1][l][m][n] && (l != -min_l || m != -min_m || n != -min_n)) {
						allAtomTypes[i] = AN1;
						allAtomL[i] = l;
						allAtomM[i] = m;
						allAtomN[i] = n;
						++i;
					}
				}
			}
		}
	}

	f = fopen("atoms.mat", "w");
	/*print all atoms*/
	fprintf(f, "all atoms: %i\n",allAtomNum);
	for (i=0; i<allAtomNum; ++i){
		type = allAtomTypes[i];
		l = allAtomL[i];
		m = allAtomM[i];
		n = allAtomN[i];
		fprintf(f, "%4i %4i %4i %4i %4i\n", i + 1, type, l + min_l, m + min_m, n + min_n);
	}
	fclose(f);

	/*print dm*/
	int AN2;

	const char * spin_array[] = {"aa_re", "bb_re","ab_re","ab_im"};
	const char * spin_iarray[] = {"aa_im","bb_im"};
	int matrix_switch = (SpinP_switch==3) ? SpinP_switch+2 : SpinP_switch;
	for (spin=0; spin<=matrix_switch; spin++){
		char* dm_index;
		if (spin<=SpinP_switch) {
			asprintf(&dm_index, spin_array[spin]);
		} else {
			asprintf(&dm_index, spin_iarray[spin-SpinP_switch-1]);
		}
		printf("DM_%s\n", dm_index);
		char* fname;
		asprintf(&fname, "dm_%s.mat", dm_index);
		f = fopen(fname, "w");
		free(fname);

		fprintf(f,"\nDensity matrix spin = %s\n", dm_index);
		free(dm_index);
		

		/*print matrix*/
		int GAN1=0;
		for (GAN1=1; GAN1<atomnum+1; GAN1++){
			fprintf(f,"\nDensity matrix part for cell atom %i\n", GAN1);
			fprintf(f,"**********");
			int l1 = 0;
			int m1 = 0;
			int n1 = 0;
			for (AN2=0; AN2<allAtomNum; ++AN2){
				int GAN2 = allAtomTypes[AN2];
				int l2 = allAtomL[AN2] + min_l;
				int m2 = allAtomM[AN2] + min_m;
				int n2 = allAtomN[AN2] + min_n;
				int dl = l2 - l1 - min_l;
				int dm = m2 - m1 - min_m;
				int dn = n2 - n1 - min_n;
				if (dl < d_l && dl >=0 && dm < d_m && dm >=0 && dn < d_n && dn >=0) {
					h_AN = NB[GAN1][GAN2][dl][dm][dn];
				} else {
					h_AN = -1;
				}
				if (h_AN >= 0) {
					TNO2 = Total_NumOrbs[GAN2];
					for (i=0; i<TNO2; i++) {
						fprintf(f, " %12i %2i", AN2 + 1, i + 1);
					}
				}
			}
			TNO1 = Total_NumOrbs[GAN1];
			for (i=0; i<TNO1; i++) {
				fprintf(f, "\n %12i %2i", GAN1, i + 1);
				for (AN2=0; AN2<allAtomNum; AN2++) {
					int GAN2 = allAtomTypes[AN2];
					int l2 = allAtomL[AN2] + min_l;
					int m2 = allAtomM[AN2] + min_m;
					int n2 = allAtomN[AN2] + min_n;
					int dl = l2 - l1 - min_l;
					int dm = m2 - m1 - min_m;
					int dn = n2 - n1 - min_n;
					if (dl < d_l && dl >=0 && dm < d_m && dm >=0 && dn < d_n && dn >=0) {
						h_AN = NB[GAN1][GAN2][dl][dm][dn];
					} else {
						h_AN = -1;
					}
					if (h_AN >= 0) {
						TNO2 = Total_NumOrbs[GAN2];
						for (j=0; j<TNO2; j++){
							double dme;
							if (spin<=SpinP_switch) {
								dme = DM[spin][GAN1][h_AN][i][j];
							} else {
								dme = iDM[spin-SpinP_switch-1][GAN1][h_AN][i][j];
							}
							fprintf(f, " %15.12f", dme);
						}
					}
				}
			}
		}
		fclose(f);
	}

	f = fopen("olp.mat", "w");
	printf("overlap");
	fprintf(f,"\nOverlap matrix\n");

	/*print matrix*/
	int GAN1=0;
	for (GAN1=1; GAN1<atomnum+1; GAN1++){
		fprintf(f,"\nOverlap matrix part for cell atom %i\n", GAN1);
		fprintf(f,"**********");
		int l1 = 0;
		int m1 = 0;
		int n1 = 0;
		for (AN2=0; AN2<allAtomNum; ++AN2){
			int GAN2 = allAtomTypes[AN2];
			int l2 = allAtomL[AN2] + min_l;
			int m2 = allAtomM[AN2] + min_m;
			int n2 = allAtomN[AN2] + min_n;
			int dl = l2 - l1 - min_l;
			int dm = m2 - m1 - min_m;
			int dn = n2 - n1 - min_n;
			if (dl < d_l && dl >=0 && dm < d_m && dm >=0 && dn < d_n && dn >=0) {
				h_AN = NB[GAN1][GAN2][dl][dm][dn];
			} else {
				h_AN = -1;
			}
			if (h_AN >= 0) {
				TNO2 = Total_NumOrbs[GAN2];
				for (i=0; i<TNO2; i++) {
					fprintf(f, " %12i %2i", AN2 + 1, i + 1);
				}
			}
		}
		TNO1 = Total_NumOrbs[GAN1];
		for (i=0; i<TNO1; i++) {
			fprintf(f, "\n %12i %2i", GAN1, i + 1);
			for (AN2=0; AN2<allAtomNum; AN2++) {
				int GAN2 = allAtomTypes[AN2];
				int l2 = allAtomL[AN2] + min_l;
				int m2 = allAtomM[AN2] + min_m;
				int n2 = allAtomN[AN2] + min_n;
				int dl = l2 - l1 - min_l;
				int dm = m2 - m1 - min_m;
				int dn = n2 - n1 - min_n;
				if (dl < d_l && dl >=0 && dm < d_m && dm >=0 && dn < d_n && dn >=0) {
					h_AN = NB[GAN1][GAN2][dl][dm][dn];
				} else {
					h_AN = -1;
				}
				if (h_AN >= 0) {
					TNO2 = Total_NumOrbs[GAN2];
					for (j=0; j<TNO2; j++){
						double olpe = OLP[GAN1][h_AN][i][j];
						fprintf(f, " %15.12f", olpe);
					}
				}
			}
		}
	}
	fclose(f);



}
