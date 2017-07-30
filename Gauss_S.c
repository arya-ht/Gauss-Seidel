#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// created by Arya HajiTaheri
//#pragma warning(disable:4996)
#define tMax 10
#define dx (M_PI/10)
#define xMAX M_PI
#define x0 0
#define D 1.0
#define F 0.0
#define n 1.0

#define lambda .4
#define dt (lambda*dx*dx/D)
//function declaration
double U0(double);
double exact(double, double);
void Gauss_S_N();

int main(void)
{
    Gauss_S_N();
    return 0;
}

void Gauss_S_N()
{
    // init
    int k;
    printf("Enter k:\n");
    scanf("%d", &k);
    int N = (int)((xMAX - x0) / dx), M = (int)(tMax / dt);
    double u[M + 1][N + 1], error = 0.0;


    int i, j, t;
    FILE *error_f = fopen("Error.csv", "w"), *estimate_f= fopen("Estimate.csv", "w"), *exact_f= fopen("Exact.csv", "w");
    for (t = 0; t <= M; t++)
    {
        for (i = 1; i<N; i++)
        {
            u[t][i] = U0(i*dx);
        }
    }
    for (i = 0; i <= M; i++)
    {
        u[i][0] = 0;
        u[i][N] = 0;
    }

    for (i = 0; i <= M; i++)
    {
        if (u[i][j] == 0) // because error would be infinity (0/0)
        {
            error = 0;
        }

        printf("Estimate: ");
        for (j = 0; j <= N; j++)
        {
            printf("%.8lf ", u[i][j]);
            fprintf(estimate_f, "%lf, ", u[i][j]);
            fprintf(exact_f, "%lf, ", exact(j*dx, i*dt));error = fabs((u[i][j] - exact(j*dx, i*dt)) / u[i][j]) * 100;
            fprintf(error_f, "%lf, ", error);
        }
        printf("\n\nExact:");
        for (j = 0; j <= N; j++)
        {
            printf("%.8lf ", exact(j*dx, i*dt));
        }
        printf("\n");
        printf("\n");
        fprintf(estimate_f, "\n");
        fprintf(exact_f, "\n");
        fprintf(error_f, "\n");

    }
    printf("\n\n\a");

    fclose(exact_f);
    fclose(estimate_f);
    fclose(error_f);
}
double U0(double x)
{
    return  sin(n*x);
}
double exact(double x, double t)
{
    return exp(-n*n*t)*sin(n*x);
}

