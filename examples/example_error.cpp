//
//  SOFT.cpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 25.05.15.
//
//

#include <pfsoft>

using namespace pfsoft;

void for_back(unsigned int bandwidth, int threads)
{
    // create a grid to fill with values
    grid3D< complex< double > > sample(2 * bandwidth);
    
    // creating fourier coefficients container
    DSOFTFourierCoefficients coef(bandwidth);
    DSOFTFourierCoefficients rec_coef(bandwidth);
    
    uniform_real_distribution< double > ctx;
    ctx.min = -1;
    ctx.max = +1;
    
    rand(coef, ctx);
    
    stopwatch sw = stopwatch::tic();
    FourierTransforms::IDSOFT(coef, sample, threads);
    double time2 = sw.toc();
    
    // perform forward DSOFT transform
    sw = stopwatch::tic();
    FourierTransforms::DSOFT(sample, rec_coef, threads);
    double time = sw.toc();
	
    double max_abs_error = 0.0;
    double max_rel_error = 0.0;
    
    for (int m = 0; m < bandwidth; ++m)
    {
        for (int n = -m; n <= m; ++n)
        {
            for (int k = -m; k <= m; ++k)
            {
                double abs_error = (coef(m,n,k).re - rec_coef(m,n,k).re) * (coef(m,n,k).re - rec_coef(m,n,k).re) + (coef(m,n,k).im - rec_coef(m,n,k).im) * (coef(m,n,k).im - rec_coef(m,n,k).im);
                double rel_error = ((coef(m,n,k).re - rec_coef(m,n,k).re) * (coef(m,n,k).re - rec_coef(m,n,k).re) + (coef(m,n,k).im - rec_coef(m,n,k).im) * (coef(m,n,k).im - rec_coef(m,n,k).im)) / (coef(m,n,k).re * coef(m,n,k).re + coef(m,n,k).im * coef(m,n,k).im);                

                if ( abs_error > max_abs_error ) max_abs_error = abs_error;
                if ( rel_error > max_rel_error ) max_rel_error = rel_error;
            }
        }
    }

    max_abs_error = sqrt (max_abs_error);
    max_rel_error = sqrt (max_rel_error);
    
    printf("Bandbreite:      %d\n", bandwidth);
    printf("DSOFT:           %.6fs\n", time);
    printf("IDSOFT:          %.6fs\n", time2);
    printf("max abs error:   %.2e\n", max_abs_error);
    printf("max rel error:   %.2e\n", max_rel_error);
}

int main(int argc, const char ** argv)
{
    if (argc < 3)
    {
        printf("usage: ./soft_test <Bandwidth> <Threads>\n");
        return 1;
    }
    
    int B = atoi(*(argv + 1));
    int T = atoi(*(argv + 2));
    for_back(B, T);
    
    return 0;
}

