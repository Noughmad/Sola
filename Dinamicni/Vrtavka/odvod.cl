#pragma OPENCL EXTENSION cl_intel_printf : enable

__kernel void odvod(__global const char* run,
                    float D,
                    float L,
                    float mg,
                    float a,
                    float R,
                    __global const double* y,
                    __global double* dy)
{
    const int i = get_global_id(0);
    // printf("Processing i=%d, run is %u\n", i, run[i]);

    if (run[i])
    {
        dy[6*i+0] = y[6*i+1] * y[6*i+2] * R + mg * a * y[6*i+4];
        dy[6*i+1] = -y[6*i+2] * y[6*i+0] * R + mg * (L * y[6*i+5] - a * y[6*i+3]);
        dy[6*i+2] = -mg * L * y[6*i+4];

        dy[6*i+3] = y[6*i+1] * y[6*i+5] - y[6*i+2] * y[6*i+4] / D;
        dy[6*i+4] = y[6*i+2] * y[6*i+3] / D - y[6*i+0] * y[6*i+5];
        dy[6*i+5] = y[6*i+0] * y[6*i+4] - y[6*i+1] * y[6*i+3];
    }
    else
    {
        dy[6*i] = 0;
        dy[6*i+1] = 0;
        dy[6*i+2] = 0;
        dy[6*i+3] = 0;
        dy[6*i+4] = 0;
        dy[6*i+5] = 0;
    }
}

__kernel void razdalja(__global const double* y1,
                       __global const double* y2,
                       __global double* d)
{
    const int i = get_global_id(0);

    d[i] = 0;
    for (int k = 0; k < 6; ++k)
    {
        d[i] += (y1[6*i+k] - y2[6*i+k]) * (y1[6*i+k] - y2[6*i+k]);
    }
}