//
//  SOFT.cpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 25.05.15.
//
//

#include <pfsoft>

using namespace pfsoft;

void for_back(unsigned int bandwidth, bool show_coefs)
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
    FourierTransforms::IDSOFT(coef, sample);
    double time2 = sw.toc();
    
    // perform forward DSOFT transform
    sw = stopwatch::tic();
    FourierTransforms::DSOFT(sample, rec_coef);
    double time = sw.toc();
    
    // print Fourier coefficients
    // save outstream flags
    if (show_coefs)
    {
        printf("** Fourier coefficients:\n");
    }
    
    bool equal = true;
    double epsilon = 1e-11;
    int cnt_fc = 0;
	
    double max_abs_error = 0.0;
    double max_rel_error = 0.0;
    
    for (int m = 0; m < bandwidth; ++m)
    {
        for (int n = -m; n <= m; ++n)
        {
            for (int k = -m; k <= m; ++k)
            {
                if (fabs(coef(m,n,k).re - rec_coef(m,n,k).re) > epsilon || fabs(coef(m,n,k).im - rec_coef(m,n,k).im) > epsilon)
                {
                    equal = false;
                }

                double abs_error = (coef(m,n,k).re - rec_coef(m,n,k).re) * (coef(m,n,k).re - rec_coef(m,n,k).re) + (coef(m,n,k).im - rec_coef(m,n,k).im) * (coef(m,n,k).im - rec_coef(m,n,k).im);
                double rel_error = ((coef(m,n,k).re - rec_coef(m,n,k).re) * (coef(m,n,k).re - rec_coef(m,n,k).re) + (coef(m,n,k).im - rec_coef(m,n,k).im) * (coef(m,n,k).im - rec_coef(m,n,k).im)) / (coef(m,n,k).re * coef(m,n,k).re + coef(m,n,k).im * coef(m,n,k).im);                

                if ( abs_error > max_abs_error ) max_abs_error = abs_error;
                if ( rel_error > max_rel_error ) max_rel_error = rel_error;

                cnt_fc++;
                
                if ( show_coefs )
                {
                    // Here the coefficients are printed out on the console
                    printf("l=%4d, M=%4d, M'=%4d: %.4f%s%.4f\n", m,  n, k, coef(m,n,k).re, (coef(m,n,k).im >= 0 ? "+" : ""), coef(m,n,k).im);
                    //printf("%.16f\n", fabs(coef(m,n,k).re - rec_coef(m,n,k).re));
                    
                    if (fabs(coef(m,n,k).re - rec_coef(m,n,k).re) > epsilon || fabs(coef(m,n,k).im - rec_coef(m,n,k).im) > epsilon)
                    {
                        equal = false;
                    }
                }
            }
        }
    }
    if (show_coefs)
    {
        printf("\n");
    }

    max_abs_error = sqrt (max_abs_error);
    max_rel_error = sqrt (max_rel_error);
    
    //printf("#coefficients:   %d\n", cnt_fc);
    printf("Bandbreite:      %d\n", bandwidth);
    printf("DSOFT:           %.6fs\n", time);
    printf("IDSOFT:          %.6fs\n", time2);
    //printf("Correct result:  %s\n", (equal ? "Yes" : "No"));
    printf("max abs error:   %.2e\n", max_abs_error);
    printf("max rel error:   %.2e\n", max_rel_error);
}

int main(int argc, const char ** argv)
{
    if (argc < 2)
    {
        printf("usage: ./soft_test <Bandwidth>\n");
        return 1;
    }
    
    int B = atoi(*(argv + 1));
    for_back(B, false);

//    matrix< double > A(3, 4);
//    
//    int cnt = 1;
//    
//    for (int i = 0; i < A.rows; ++i)
//    {
//        for (int j = 0; j < A.cols; ++j)
//        {
//            A(i, j) = cnt++;
//        }
//    }
//    
//    std::cout << "A = " << A << std::endl;
//    A.transpose();
//    std::cout << "A' = " << A << std::endl;
    
    return 0;
}

