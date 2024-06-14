#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <charm/charm.h>


const double PI = 3.14159265358979323846;
const double PI_2 = 1.5707963267948966;


void fill_vector(double *v, double min, double step, double n){
    for (size_t i = 0; i < n; i++)
    {
        v[i] = min + (double)i * step;
    }
}

void write_vector(const char *name, double *v, int n){

    FILE *file = fopen(name, "w");

    if (file == NULL)
    {
        printf("Failed to open the file.\n");
    }

    for (int i = 0; i < n; i++){
        fprintf(file, "%lf ", v[i]);

        fprintf(file, "\n");
    }

    fclose(file);

}


double deg2rad(double x){
    return x / 180 * PI;
}



void write_array(const char* name, double *a, int rows, int cols){

    FILE *file = fopen(name, "w");

    if (file == NULL)
    {
        printf("Failed to open the file.\n");
    }

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++)
        {
            fprintf(file, "%lf ", a[i * cols + j]);
        }
        fprintf(file, "\n");
    }
    
    fclose(file);

}


void function_array(double *out, double *rows, double *cols, int n_rows, int n_cols){

    for (int i = 0; i < n_rows; i++){
        for (int j = 0; j < n_cols; j++)
        {
            out[i * n_cols + j] = exp(rows[i]) * exp(cols[j]);
        }
    }
    
}



int main(void)

{

    time_t start_t, end_t;
    double diff_t;

    int N = 2160;

    char shcs_in_file[] = "./synthesis_data/Earth2012.RET2012.SHCto2160.dat";


    /* Initialize a "charm_shc" structure */
    charm_shc *coeff = charm_shc_calloc(N, 1.0, 1.0);
    if (coeff == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Initialize a "charm_err" structure */
    charm_err *err = charm_err_init();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_err\" structure.\n");
        exit(CHARM_FAILURE);
    }

    charm_shc_read_tbl(shcs_in_file, N, coeff, err);

    // unsigned long n = 2160;
    // unsigned long m = 2160;
    // printf("C(%3lu,%3lu) = %0.16e\n", n, m, coeff->c[m][n - m]);
    // printf("S(%3lu,%3lu) = %0.16e\n", n, m, coeff->s[m][n - m]);


    // // the grid for performing synthesis
    double sh_lat_max = deg2rad(40.0);
    double sh_lon_min = deg2rad(250);
    int nlat_syn = 1799;
    int nlon_syn = 3599;
    double step = deg2rad(0.003);


    charm_point* grd_pnt = charm_crd_point_calloc(CHARM_CRD_POINT_GRID, nlat_syn, nlon_syn);
    if (grd_pnt == NULL)
    {
        fprintf(stderr,
            "Failed to initialize the \"charm_point\" structure.\n");
        exit(CHARM_FAILURE);
    }
    
    /* define grid points */
    fill_vector(grd_pnt->lat, sh_lat_max, -step, nlat_syn);
    fill_vector(grd_pnt->r, 1, 0, nlat_syn);
    fill_vector(grd_pnt->lon, sh_lon_min, step, nlon_syn);


    /* allocate the memory for the output of the synthesis*/
    double *out = (double *)malloc(nlat_syn * nlon_syn * sizeof(double));
    if (out == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    printf("Starting of the program...\n");
    time(&start_t);

    charm_shs_point(grd_pnt, coeff, N, out, err);

    time(&end_t);
    diff_t = difftime(end_t, start_t);

    printf("Execution time = %f\n", diff_t);

    charm_err_handler(err, 1);

    write_array("./synthesis_data/sph_heights.txt", out, nlat_syn, nlon_syn);


    
}