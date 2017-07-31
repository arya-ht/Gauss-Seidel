#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
//#pragma warning(disable:4996)
//created by arya HajiTaheri
#define dx (M_PI/10.0)
#define lambda 0.4
#define x_max M_PI
#define t_max 10.0
#define D 1.0
#define N 1.0
#define dt (lambda*dx*dx/D)
#define iteration_K 5000000

#define SIZE_X (int)(ceil(x_max/dx))
#define SIZE_T (int)(ceil(t_max/dt))

int F(double, double);
double exact(double, double);
void Gauss_S();
double error(double, double);
void print(double U[SIZE_T + 1][SIZE_X + 1]);
int main(void)
{
    Gauss_S();
    return 0;
}

int F(double t, double x)
{
    return 0;
}

double exact(double t, double x)
{
    return exp(-N*N*t)*sin(N*x);
}

double error(double exact, double approx)//calculate error %
{
    return (fabs(approx - exact) / fabs(approx)) * 100;
}
void print(double U[SIZE_T + 1][SIZE_X + 1])
{
    int n, j;

    double exact_value = 0, error_value = 0;
    FILE *u, *x, *t, *exact_graph, *error_graph;
    u = fopen("u.csv", "w"), x = fopen("x.csv", "w"), t = fopen("t.csv", "w"),
    exact_graph = fopen("exact.csv", "w"), error_graph = fopen("error.csv", "w");

    for (n = 0; n <= SIZE_T; n++)
    {
        for (j = 0; j <= SIZE_X; j++)
        {
            exact_value = exact(n*dt, j*dx);
            error_value = error(exact_value, U[n][j]);
            if (j == 0 || j == SIZE_X) // cannot divide 0/0
            {
                error_value = 0.0;
            }
            printf("Estimate: %.15lf\tExact: %.15lf\tError: %.15lf%%\n", U[n][j], exact_value, error_value);
            fprintf(u, "%lf\n", U[n][j]);
            fprintf(x, "%d,", j);
            fprintf(t, "%d\n", n);
            fprintf(exact_graph, "%lf\n", exact_value);
            fprintf(error_graph, "%lf\n", error_value);
        }
        printf("\n\n");
    }
    fclose(u), fclose(x), fclose(t),
           fclose(exact_graph), fclose(error_graph);
}
void Gauss_S()
{
    int n, k, i, t, x, tol;

    double u[SIZE_T + 1][SIZE_X + 1];

    //calculating u[][]
    for (k = 0; k <= SIZE_T; k++)
    {
        for (i = 0; i <= SIZE_X; i++)
        {
            u[k][i] = sin(i*dx*N);
        }
    }
    //set right/left bounds
    for (n = 0; n <= SIZE_T; n++)
    {
        u[n][0] = 0;
        u[n][SIZE_X] = 0;
    }

    //init values
    double a = 1.0 + 2.0*lambda, b = -lambda, c = -lambda,
           temp_u[SIZE_T][SIZE_X], tolerance = 0, bool_last_iteration = 0, abort_program = 0;
    for (x = 0; x <= SIZE_X; x++)
    {
        temp_u[0][x] = u[0][x]; // have a u to compare tolerance to
    }

    for (t = 1; t <= SIZE_T; t++) // follow t->k->x loop for accurate process.
    {
        for (k = 1; k<iteration_K; k++)
        {
            bool_last_iteration = 0;
            abort_program = 0;
            for (x = 1; x<SIZE_X; x++)
            {
                u[t][x] = (u[t - 1][x] - b*u[t][x - 1] - c*u[t][x + 1]) / a;
                temp_u[k][x] = u[t][x];
            }
            for (tol = 1; tol<SIZE_X; tol++)
            {
                tolerance = fabs(u[t][tol] - temp_u[k - 1][tol]);
                if (tolerance<0.0000000001) // if tolerance is getting smaller than this, abort.
                {

                    abort_program++;
                    if (abort_program == SIZE_X - 1)
                    {
                        bool_last_iteration = 1;
                        break;
                    }
                }
            }
            if (bool_last_iteration == 1)
            {
                break;
            }
        }
    }
    printf("Done! Now printing...\n\a");
    print(u);
}
