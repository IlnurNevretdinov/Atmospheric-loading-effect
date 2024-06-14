#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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



// double F(double x, double y)
// {
//     return exp(x) * exp(y);
// }



int main(void)

{

    int N = 300;
    int nlat_ans = N + 1;
    int nlon_ans = 2 * N + 1;

    // where the final coefficients will be stored
    charm_shc *coeff = charm_shc_malloc(N, 1.0, 1.0);
    if (coeff == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_shc\" structure.\n");
        exit(CHARM_FAILURE);
    }

    charm_err *err = charm_err_init();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_err\" structure.\n");
        exit(CHARM_FAILURE);
    }

    /* Initialize a "charm_cell" structure to hold a global grid of cells 1 on 1 degree,
    on which values of function are computed, it will be used at analysis stage */
    charm_cell *grd_cell = charm_crd_cell_malloc(CHARM_CRD_CELL_GRID, nlat_ans, nlon_ans);

    // vector with latutudes of centers of all cells
    double *lat_mean = (double*)malloc(nlat_ans * sizeof(double));
    if (lat_mean == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }
    // vector with longitudes of centers of all cells
    double *lon_mean = (double *)malloc(nlon_ans * sizeof(double));
    if (lon_mean == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }

    /* Define some grid cells */
    for (int i = 0; i < nlat_ans; i++)
    {
        // write the latitudal boundaries for i-th cell into the structure
        grd_cell->latmax[i] = PI_2 -((double)i / (double)nlat_ans) * PI;;
        grd_cell->latmin[i] = PI_2 - ((double)(i + 1) / (double)nlat_ans) * PI;
        // compute the mean latitude of this cell and write it to the corresponding vector
        lat_mean[i] = (grd_cell->latmax[i] + grd_cell->latmin[i]) / 2;
        // the spherical radii have a constant value here
        grd_cell->r[i] = coeff->r;
    }

    for (int j = 0; j < nlon_ans; j++)
    {
        // write the longitudal boundaries for i-th cell into the structure
        grd_cell->lonmin[j] = ((double)j / (double)nlon_ans) * (2.0 * PI);
        grd_cell->lonmax[j] = ((double)(j + 1) / (double)nlon_ans) * (2.0 * PI);
        // compute the mean longitude of this cell and write it to the corresponding vector
        lon_mean[j] = (grd_cell->lonmax[j] + grd_cell->lonmin[j]) / 2;

    }

    /* allocate the memory for the real-valued function */
    double *f = (double *)malloc(nlat_ans * nlon_ans * sizeof(double));
    if (f == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }

    // Fill the allocated memory with values and write them into the file
    function_array(f, lat_mean, lon_mean, nlat_ans, nlon_ans);

    /* perform the analysis using the approximate quadrature */
    charm_sha_cell(grd_cell, f, N, CHARM_SHA_CELL_AQ, coeff, err);

    // the grid for performing synthesis
    double sh_lat_min = deg2rad(35.0);
    double sh_lon_min = deg2rad(200.0);
    int nlat_syn = 21;
    int nlon_syn = 21;
    double step = deg2rad(1.0);

    charm_point* grd_pnt = charm_crd_point_calloc(CHARM_CRD_POINT_GRID, nlat_syn, nlon_syn);
    if (grd_pnt == NULL)
    {
        fprintf(stderr,
            "Failed to initialize the \"charm_point\" structure.\n");
        exit(CHARM_FAILURE);
    }

    // /* define grid points */
    fill_vector(grd_pnt->lat, sh_lat_min, -step, nlat_syn);
    fill_vector(grd_pnt->r, 1, 0, nlat_syn);
    fill_vector(grd_pnt->lon, sh_lon_min, step, nlon_syn);


    write_vector("./synthesis_data/lat_synth.txt", grd_pnt->lat, nlat_syn);
    write_vector("./synthesis_data/r_synth.txt", grd_pnt->r, nlat_syn);
    write_vector("./synthesis_data/lon_synth.txt", grd_pnt->lon, nlon_syn);

    
    /* allocate the memory for the output of the synthesis*/
    double *out = (double *)malloc(nlat_syn * nlon_syn * sizeof(double));
    if (out == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }

    charm_shs_point(grd_pnt, coeff, N, out, err);

    charm_err_handler(err, 1);

    // write into the files all data for perfoming an spherical analysis in Julia
    write_vector("./analysis_data/lat_values.txt", lat_mean, nlat_ans);
    write_vector("./analysis_data/lon_values.txt", lon_mean, nlon_ans);
    write_array("./analysis_data/function_values.txt", f, nlat_ans, nlon_ans);
    charm_shc_write_tbl(coeff, N, "%0.16f", CHARM_SHC_WRITE_M, "./analysis_data/harmonic_coefficients.txt", err);
    // to compare two solutions
    write_array("./synthesis_data/modelled_values.txt", out, nlat_syn, nlon_syn);
    



}
