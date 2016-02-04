#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct ComplexShort {
	unsigned short re;
	unsigned short im;
} ComplexShort;

ComplexShort ComplexMult(ComplexShort a, ComplexShort b);
ComplexShort ComplexAdd(ComplexShort a, ComplexShort b);
ComplexShort ComplexSub(ComplexShort a, ComplexShort b);

static int CurrentW = 0;
static ComplexShort *W = NULL;
static int *BitRev = NULL;
static double *WindowFunction = NULL;
static int CurrentWindowType = -1;
static void MakeWindowFunction(int DataLength, int wtype);
static void MakeW(int DataLength, int N, int wtype);
static void FFT(unsigned short *data, ComplexShort *spectrum, int DataLength,
	int WindowType);
static void CalcPow(ComplexShort *spectrum, double *pow, int DataLength,
	bool LogScale, int FullScale);


int main() {
	int i;
	const int datalength = 64;
	int chs = 1;

	ComplexShort *spectrum = (ComplexShort *)malloc(
		sizeof(ComplexShort) * datalength * chs);
	double *pow = (double *)malloc(sizeof(double) * datalength * chs);
	memset(spectrum, 0, sizeof(ComplexShort) * datalength * chs);
	memset(pow, 0, sizeof(double) * datalength * chs);
	unsigned short *rdata = (unsigned short *)malloc(sizeof(unsigned short) * datalength * chs);

	int tmp_ch;
	for (i = 0; i < datalength; i++) {
		rdata[i] = 0x7fff * i / datalength;//cos(2 * 3.14159 * i / datalength * 3);
										   //		rdata[i] = 0x7fff * cos(2 * 3.14159 * i / datalength * 3);
	}
	for (i = 0; i<datalength; i++) {
		//map_ch1[i] = (unsigned short)(16384 * sin(2 * M_PI * i / memlen) + 1.0);
		rdata[i] = 0x7fff * sin(2 * 3.141592653589793 * 10 * i / (double)datalength);
	}
	for (i = 0; i<datalength; i++) {
		//map_ch1[i] = (unsigned short)(16384 * sin(2 * M_PI * i / memlen) + 1.0);
		rdata[i] = 0x7fff * sin(0.78539816339744830961566084581988 * i);
	}

	tmp_ch = 0;
	FFT(&rdata[tmp_ch * datalength], &spectrum[tmp_ch * datalength],
		datalength, 0);
	//	CalcPow(&spectrum[tmp_ch * datalength], &pow[tmp_ch * datalength],
	//			datalength / 2, true, 65535);

	for (i = 0; i < datalength; i++) {
		//	for (i = 0; i < datalength / 2; i++) {
		printf("%05d\t", i);
		int re = (short)spectrum[datalength * tmp_ch + i].re;
		int im = (short)spectrum[datalength * tmp_ch + i].im;
		printf("%d\t%d", re, im);
		printf("\n");
	}

	/*
	ComplexShort A,B;
	A.re = 0xc000;
	A.im = 0xc000;
	B.re = 0xc000;
	B.im = 0xc000;

	ComplexShort C = ComplexMult(A,B);
	printf(": %x %x\n",C.re,C.im);
	*/

	free(rdata);
	free(spectrum);
	free(pow);
}

void FFT(unsigned short *data, ComplexShort *spectrum, int DataLength,
	int WindowType) {
	int stage, i, j, k;
	int N;

	switch (DataLength) {
	case 4:
		N = 2;
		break;
	case 8:
		N = 3;
		break;
	case 16:
		N = 4;
		break;
	case 32:
		N = 5;
		break;
	case 64:
		N = 6;
		break;
	case 128:
		N = 7;
		break;
	case 256:
		N = 8;
		break;
	case 512:
		N = 9;
		break;
	case 1024:
		N = 10;
		break;
	case 2048:
		N = 11;
		break;
	case 4096:
		N = 12;
		break;
	case 8192:
		N = 13;
		break;
	case 16384:
		N = 14;
		break;
	case 32768:
		N = 15;
		break;
	case 65536:
		N = 16;
		break;
	case 131072:
		N = 17;
		break;
	case 262144:
		N = 18;
		break;
	case 524288:
		N = 19;
		break;
	case 1048576:
		N = 20;
		break;
	case 2097152:
		N = 21;
		break;
	case 4194304:
		N = 22;
		break;
	case 8388608:
		N = 23;
		break;
	case 16777216:
		N = 24;
		break;
	default:
		return;
	}
	//	MakeWindowFunction(DataLength, WindowType);
	MakeW(DataLength, N, WindowType);


	/*
	for (i = 0; i < DataLength; i++) {
	int v = (int) data[i];
	spectrum[BitRev[i]].re = (double) v / DataLength * WindowFunction[i];
	spectrum[BitRev[i]].im = 0;
	}
	*/
	for (i = 0; i < DataLength; i++) {
		//		spectrum[BitRev[i]].re = data[i];
		//		spectrum[BitRev[i]].im = 0;
		spectrum[i].re = data[i];
		spectrum[i].im = 0;
	}

	/*
	// normal
	int b = 2;
	for (stage = 0; stage < N; stage++) {
	for (i = 0; i < DataLength / b; i++) {
	int offset = b * i;
	for (j = 0; j < b / 2; j++) {
	ComplexShort A = spectrum[offset + j];
	ComplexShort B = ComplexMult(spectrum[offset + j + b / 2],W[j * DataLength / b]);
	spectrum[offset + j] = ComplexAdd(A, B);
	spectrum[offset + j + b / 2] = ComplexSub(A, B);
	printf("stage %d, w %d, bfly %d & %d",stage,j * DataLength / b,offset + j,offset + j + b / 2);
	printf("=> write %d & %d\n",offset + j,offset + j + b / 2);
	}
	}
	b = b * 2;
	}
	*/

	/*
	// normal VFFT
	ComplexShort *work = (ComplexShort *)malloc(sizeof(ComplexShort) * DataLength);
	int mask = ~((2 << (N-2))-1);
	for (stage = 0; stage < N; stage++) {
	int offset = DataLength/2;
	for (j = 0; j < DataLength / 2; j++) {
	ComplexShort A = spectrum[j*2];
	ComplexShort B = ComplexMult(spectrum[j*2+1],W[j & mask]);
	work[j] = ComplexAdd(A, B);
	work[j + DataLength/2] = ComplexSub(A, B);
	printf("stage %d, w %d, bfly %d & %d",stage,j & mask,j*2,j*2+1);
	printf("=> write %d & %d\n",j,j + DataLength/2);
	}
	mask = mask >> 1;
	for (j = 0; j < DataLength; j++) {
	spectrum[j] = work[j];
	}
	}
	free(work);
	*/

	/*
	// 周波数間引き
	int b = 2 << (N-1);
	for (stage = 0; stage < N; stage++) {
	for (i = 0; i < DataLength / b; i++) {
	int offset = b * i;
	for (j = 0; j < b / 2; j++) {
	ComplexShort A = spectrum[offset + j];
	ComplexShort B = spectrum[offset + j + b / 2];
	spectrum[offset + j] = ComplexAdd(A, B);
	spectrum[offset + j + b / 2] = ComplexMult(ComplexSub(A, B),W[j * DataLength / b]);
	printf("stage %d, w %d, bfly %d & %d",stage,j * DataLength / b,offset + j,offset + j + b / 2);
	printf("=> write %d & %d\n",offset + j,offset + j + b / 2);
	}
	}
	b = b / 2;
	}
	*/

	// 周波数間引きVFFT
	ComplexShort *work = (ComplexShort *)malloc(sizeof(ComplexShort) * DataLength);
	int mask = (2 << N) - 1;
	//	for (stage = 0; stage < 4; stage++) {
	for (stage = 0; stage < N; stage++) {
		int offset = DataLength / 2;
		for (j = 0; j < DataLength / 2; j++) {
			ComplexShort A = spectrum[j];
			ComplexShort B = spectrum[j + DataLength / 2];
			work[j * 2] = ComplexAdd(A, B);
			work[j * 2 + 1] = ComplexMult(ComplexSub(A, B), W[j & mask]);
			//			printf("stage %d, w %d, bfly %d & %d",stage,j & mask,j*2,j*2+1);
			//			printf("=> write %d & %d\n",j,j + DataLength/2);
			printf("%d ", j);
			printf("W:%04x+%04xi,", W[j & mask].re, W[j & mask].im);
			printf("A:%04x+%04xi, ", A.re, A.im);
			printf("B:%04x+%04xi => ", B.re, B.im);
			printf("%04x+%04xi , ", work[j * 2].re, work[j * 2].im);
			printf("%04x+%04xi\n", work[j * 2 + 1].re, work[j * 2 + 1].im);
		}
		printf("\n");
		mask = mask << 1;
		for (j = 0; j < DataLength; j++) {
			spectrum[j] = work[j];
		}
	}
	free(work);

	printf("result:\n");
	int s;
	for (s = 11; s <= 11; s++) {
		printf("stage %d\n", s);
		for (i = 0; i<DataLength; i++) {
			printf("%d %04x %04x\n", i, spectrum[i].re, spectrum[i].im);

			//			short re = map_ch1[s*datalen*2+BitRev[i]*2+0];
			//			short im = map_ch1[s*datalen*2+BitRev[i]*2+1];
			//			double pow = sqrt((re / 32768.) * (re / 32768.) +(im / 32768.) * (im / 32768.));

			//			printf("%d\t%f\n",i,pow*20);
			//			printf("%d : %d %d\n",i,map_ch1[s*datalen*2+i*2+0],map_ch1[s*datalen*2+i*2+1]);
		}
	}

}

ComplexShort ComplexMult(ComplexShort a, ComplexShort b) {
	ComplexShort C;
	int ar = (short)a.re;
	int ai = (short)a.im;
	int br = (short)b.re;
	int bi = (short)b.im;

	C.re = ((ar * br) >> 15) - ((ai * bi) >> 15);
	C.im = ((ar * bi) >> 15) + ((ai * br) >> 15);
	return C;
}

ComplexShort ComplexAdd(ComplexShort a, ComplexShort b) {
	ComplexShort C;
	int ar = (short)a.re;
	int ai = (short)a.im;
	int br = (short)b.re;
	int bi = (short)b.im;

	C.re = (ar + br) >> 1;
	C.im = (ai + bi) >> 1;
	return C;
}

ComplexShort ComplexSub(ComplexShort a, ComplexShort b) {
	ComplexShort C;
	int ar = (short)a.re;
	int ai = (short)a.im;
	int br = (short)b.re;
	int bi = (short)b.im;

	C.re = (ar - br) >> 1;
	C.im = (ai - bi) >> 1;
	return C;
}

static void MakeWindowFunction(int DataLength, int wtype) {
	if ((DataLength == CurrentW) && (CurrentWindowType == wtype))
		return;
	CurrentWindowType = wtype;

	// 窓関数を生成
	double pi = 3.1415926535898;

	if (WindowFunction)
		free(WindowFunction);
	WindowFunction = (double *)malloc(sizeof(double) * DataLength);
	int i;
	for (i = 0; i < DataLength; i++) {
		double x = (double)i / DataLength;
		double w;
		switch (wtype) {
		default: // 窓関数なし
			w = 1;
			break;
		case 1: // ハニング窓
			w = 0.5 - 0.5 * cos(2 * pi * x);
			break;
		case 2: // ハミング窓
			w = 0.54 - 0.46 * cos(2 * pi * x);
			break;
		case 3: // ブラックマン窓
			w = 0.42 - 0.5 * cos(2 * pi * x) + 0.08 * cos(4 * pi * x);
			break;
		case 4: // ナットール窓
			w = 0.355768 - 0.487396 * cos(2 * pi * x)
				+ 0.144232 * cos(4 * pi * x) - 0.012604 * cos(6 * pi * x);
			break;
		case 5: // ブラックマン・ハリス窓
			w = 0.35875 - 0.48829 * cos(2 * pi * x) + 0.14128 * cos(4 * pi * x)
				- 0.01168 * cos(6 * pi * x);
			break;
		case 6: // ブラックマン・ナットール窓
			w = 0.3635819 - 0.4891775 * cos(2 * pi * x)
				+ 0.1365995 * cos(4 * pi * x) - 0.0106411 * cos(6 * pi * x);
			break;
		}
		WindowFunction[i] = w;
	}
}

static void MakeW(int DataLength, int N, int wtype) {
	if (DataLength == CurrentW)
		return;
	CurrentW = DataLength;

	// 回転子を生成
	if (W)
		free(W);
	W = (ComplexShort *)malloc(sizeof(ComplexShort) * DataLength);
	int i;
	for (i = 0; i < DataLength; i++) {
		W[i].re = 0x7fff * cos(2 * 3.141592653589793 * (double)i / DataLength);
		W[i].im = 0x7fff * sin(2 * 3.141592653589793 * (double)i / DataLength);
	}

	// ビットリバースを生成
	if (BitRev)
		free(BitRev);
	BitRev = (int *)malloc(sizeof(int) * DataLength);
	for (i = 0; i < DataLength; i++) {
		int br = 0;
		int j;
		for (j = 0; j < N; j++) {
			if (i & (1 << j))
				br |= (1 << N - j - 1);
		}
		BitRev[i] = br;
	}
}

void CalcPow(ComplexShort *spectrum, double *pow, int DataLength,
	bool LogScale, int FullScale) {
	double logfull = log10(FullScale);
	int i;
	for (i = 0; i < DataLength; i++) {
		double p2 = spectrum[i].re * spectrum[i].re
			+ spectrum[i].im * spectrum[i].im;
		if (LogScale) {
			if (p2 == 0)
				pow[i] = -160;
			else
				pow[i] = 10 * (log10(p2) - 2 * logfull);
		}
		else {
			pow[i] = sqrt(p2);
		}
	}
}
