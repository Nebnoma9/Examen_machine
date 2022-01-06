#include <complex>
#include <type_traits>
#include <iostream>
#include <vector>
#include <queue>
#include <cstdint>
#include <algorithm>
#include <chrono>
#include "lodepng/lodepng.h"
#include <omp.h>
#include <mpi.h>

struct sparseSpectralComposition 
{
    std::uint32_t ni, nj;
    std::vector<std::uint32_t> begin_rows;
    std::vector<std::uint32_t> ind_cols;
    std::vector<std::complex<double>> coefficients;
};

std::complex<double>* discretTransformFourier( std::uint32_t width, std::uint32_t height, unsigned char const* pixels )
{
    std::uint32_t firstrow;
    std::uint32_t endrow;
    int rank, size;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    constexpr const double pi = 3.14159265358979324;
    std::uint32_t ni = height;
    std::uint32_t nj = width;
    std::uint32_t niloc = ni / size;
    std::uint32_t niloc_rest = ni%size;
    std::complex<double>* X_loc = new std::complex<double>[niloc*nj];
    std::fill(X_loc, X_loc+niloc*nj, std::complex<double>(0.,0.));
/*
    if (rank==0)
    {
        niloc+=niloc_rest;
    }*/
     
    firstrow = (size-(rank+1)) * niloc;
    endrow = (size-rank) * niloc ; 
    
    
    for( std::uint32_t k1 = firstrow; k1 < endrow; ++k1 )
    {
        for (std::uint32_t k2 = 0; k2 < nj; ++k2)
        {
            for (std::uint32_t n2 = 0; n2 < endrow; ++n2 )
            {
                std::complex<double> exp2(std::cos(-2*pi*n2*k2/height), std::sin(-2*pi*n2*k2/height));
                for (std::uint32_t n1 = 0; n1 < nj; ++n1 )
                {
                    std::complex<double> exp1(std::cos(-2*pi*n1*k1/nj), std::sin(-2*pi*n1*k1/nj));
                    X_loc[k1*nj+k2] += double(pixels[3*(n1+n2*nj)])*exp1*exp2;
                }
            }
        }
    }
    return X_loc;
}

sparseSpectralComposition compressSpectralComposition( std::uint32_t width, std::uint32_t height, const std::complex<double>* plainCoefs, double tauxCompression )
{
    std::uint32_t nb_coefs = std::uint32_t(tauxCompression*width*height);
    struct node 
    {
        node(std::uint32_t t_i, std::uint32_t t_j, double val )
            :   i(t_i),
                j(t_j),
                module(val)
        {}

        bool operator < ( node const& v ) const
        {
            return this->module < v.module;
        }

        bool operator > ( node const& v ) const
        {
            return this->module > v.module;
        }

        bool operator < ( double v ) const
        {
            return this->module < v;
        }

        bool operator > ( double v ) const
        {
            return this->module > v;
        }

        std::uint32_t i,j;
        double module;
    };

    std::priority_queue<node,std::vector<node>,std::greater<node>> queue;
    for ( std::size_t i = 0; i < height; ++i )
        for ( std::size_t j = 0; j < width;  ++j )
        {
            if (queue.size() < nb_coefs)
            {
                queue.emplace(i,j,std::abs(plainCoefs[i*width+j]));
            }
            else
            {
                double val = std::abs(plainCoefs[i*width+j]);
                queue.emplace(i,j,val);
                queue.pop();
            }
        }
    //
    std::vector<std::uint32_t> nbCoefsPerRow(height, 0u);
    auto queue2 = queue;
    while (!queue2.empty())
    {
        std::uint32_t i = queue2.top().i;
        nbCoefsPerRow[i] += 1;
        queue2.pop();       
    }
    //
    sparseSpectralComposition sparse;
    sparse.begin_rows.resize(height+1);
    sparse.begin_rows[0] = 0;
    for (std::uint32_t i = 1; i <= height; ++i )
        sparse.begin_rows[i] = sparse.begin_rows[i-1] + nbCoefsPerRow[i-1];
    std::cout << "Nombre de coefficients de fourier restant : " << sparse.begin_rows[height] << std::endl << std::flush;
    sparse.ind_cols.resize(sparse.begin_rows[height]);
    sparse.coefficients.resize(sparse.begin_rows[height]);
    //
    std::fill(nbCoefsPerRow.begin(), nbCoefsPerRow.end(), 0u);
    //
    while (!queue.empty())
    {
        node n = queue.top();
        queue.pop();
        std::uint32_t ind = sparse.begin_rows[n.i]+nbCoefsPerRow[n.i];
        sparse.ind_cols[ind] = n.j;
        sparse.coefficients[ind] = plainCoefs[n.i*width+n.j];
        nbCoefsPerRow[n.i] += 1;
    }
    sparse.ni = height;
    sparse.nj = width;
    return sparse;
}

unsigned char* inversePartialDiscretTransformFourier( sparseSpectralComposition const& sparse )
{
    constexpr const double pi = 3.14159265358979324;
    unsigned char* x = new unsigned char[3*sparse.ni*sparse.nj];
    std::uint32_t ni = sparse.ni;
    std::uint32_t nj = sparse.nj;
    std::fill(x, x+3*ni*nj, 0);

    for (std::uint32_t n1 = 0; n1 < nj; ++n1 )
    {
        for (std::uint32_t n2 = 0; n2 < ni; ++n2 )
        {
            std::complex<double> val(0.);
            for( std::uint32_t k1 = 0; k1 < ni; ++k1 )
            {
                std::complex<double> exp1(std::cos(+2*pi*n1*k1/ni), std::sin(+2*pi*n1*k1/ni));
                for (std::uint32_t k2 = sparse.begin_rows[k1]; k2 < sparse.begin_rows[k1+1]; ++k2)
                {
                    std::complex<double> exp2(std::cos(+2*pi*n2*sparse.ind_cols[k2]/nj), std::sin(+2*pi*n2*sparse.ind_cols[k2]/nj));
                    val += sparse.coefficients[k2]*exp1*exp2;
                }
            }
            for (std::uint32_t c = 0; c < 3; ++c )
                x[3*(n1+n2*nj)+c] = static_cast<unsigned char>(val.real()/(ni*nj));
        }
    }
    return x;
}

int main(int nargs, char* argv[])
{
    MPI_Init(&nargs, &argv);
    std::chrono::time_point<std::chrono::system_clock> start, start1, start2, end, end1, end2;
    

    std::uint32_t width, height;
    unsigned char* image;
    if (nargs <= 1)

    {
        auto error = lodepng_decode24_file(&image, &width, &height, "data/tiny_lena_gray.png");
        if(error) std::cerr << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    }
    else 
    {
        auto error = lodepng_decode24_file(&image, &width, &height, argv[1]);
        if(error) std::cerr << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    }
    double taux = 0.10; // A changer si on veut pour une image mieux compressée ou de meilleur qualité
    if (nargs > 2)
        taux = std::stod(argv[2]);
    std::cout << "Caracteristique de l'image : " << width << "x" << height << " => " << width*height << " pixels." << std::endl << std::flush;

    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::uint32_t hloc = height / size;
    std::uint32_t niloc_rest = height%size;

    start = std::chrono::system_clock::now();
    std::complex<double>* Xloc = discretTransformFourier( width, height, image );
   
 
   

    if(rank==0)
    {

        //hloc+=niloc_rest;
    
        std::complex<double>* fourier = new std::complex<double>[height*width];
        
        MPI_Gather( Xloc, width*hloc, MPI::DOUBLE_COMPLEX,
                fourier, width*hloc, MPI::DOUBLE_COMPLEX,
                0, // root
                MPI_COMM_WORLD               
               );
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout <<"Temps DFT: " << elapsed_seconds.count() << std::endl;

        start1 = std::chrono::system_clock::now();
        auto sparse = compressSpectralComposition(width, height, fourier, taux);
        end1 = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds1 = end1-start1;
        std::cout <<"Temps p% : " << elapsed_seconds1.count() << std::endl;

        start2 = std::chrono::system_clock::now();
        unsigned char* compressed_img = inversePartialDiscretTransformFourier(sparse);
        end2 = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds2 = end2-start2;
        std::cout <<"Temps reconstitution : " << elapsed_seconds2.count() << std::endl;

        auto error = lodepng_encode24_file("compress.png", compressed_img, width, height);
        if(error) std::cerr << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;

        std::cout << "Fin du programme..." << std::endl << std::flush;
        delete [] compressed_img;
        delete [] Xloc;
        delete [] image;
    }
    else{
        MPI_Gather( Xloc, width*hloc, MPI::DOUBLE_COMPLEX,
                NULL, width*hloc, MPI::DOUBLE_COMPLEX,
                0, // root
                MPI_COMM_WORLD               
               );
               
    }


    MPI_Finalize();
    return EXIT_SUCCESS;
}