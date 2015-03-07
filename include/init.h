void init_acou_homog(FILE *filefprintf, double dx, double dz, int nx, int nz, double c0, double *vel, double *vpe);
void init_acou_layer_velocity(FILE *filefprintf, double dx, double dz, int nx, int nz, double c0, double *vel);
void init_acou_layer_density(FILE *filefprintf, double dx, double dz, int nx, int nz, double rho0, double *den);
void init_source_ricker_fwps(FILE *filefprintf, int nt, double dt,int nw, double dw, double *source, double *fren, double *sw);
