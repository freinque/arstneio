/** @file */
#include "arst_neio/arst_neio_gene.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include <tuple>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>

double disso(const Harmo& h1, const Harmo& h2, const Generic& cg)
{
    double w1 = Generic::power(h1.db);
    double w2 = Generic::power(h2.db);
    double w = Generic::power(2*cg.DECI_BASE);
    double one_over_w_times_w1_plus_w2 = 1.0/( w*(w1+w2) );
    double d1 = cg.dist_in_cb(h1,h2);
    double d2 = cg.dist_in_cb(h2,h1);

    double cutoff = 1.0;
    double cutoff_dist = 0.01;
    double a = 1.0e5;

    double disso_out = 0;

    if( d1 < 7 ) 
    {
        //added to kill non-zero gradients near d1=d2=0
        if ( d1<cutoff_dist )
        {
            cutoff = 1 - exp( -a*pow(d1,2) );
        }         
        disso_out = cutoff * 4*w1*w2*one_over_w_times_w1_plus_w2*( w1 * d1 * exp(1-4*d1) + w2 * d2 * exp(1-4*d2) );
    }

    return disso_out;
}


double ddisso_df2(const Harmo& h1, const Harmo& h2, const Generic& cg)
{
    double w1 = Generic::power(h1.db);
    double w2 = Generic::power(h2.db);
    double w = Generic::power(2*cg.DECI_BASE);
    double one_over_w_times_w1_plus_w2 = 1.0/( w*(w1+w2) );
    double d1 = cg.dist_in_cb(h1,h2);
    double d2 = cg.dist_in_cb(h2,h1);
    double exp_one_minus_4d1 = exp(1-4*d1);
    double exp_one_minus_4d2 = exp(1-4*d2);
    double one_over_cb1 = 1.0/cg.cb_fr(h1.fr);
    double one_over_cb2 = 1.0/cg.cb_fr(h2.fr);
    int sign_f2_minus_f1 = 1;
    if ( (h2.fr-h1.fr) < 0)
    {
        sign_f2_minus_f1 = -1;
    }

    double cutoff = 1.0;
    double dcutoff_dd1 = 0.0;
    double cutoff_dist = 0.01;
    double a = 1.0e5;

    double disso_out = 0;
    if( d1<7 )
    {
        //added to kill non-zero gradients near d1=d2=0
        if( d1<cutoff_dist )
        {
            double exp_minus_a_d12 = exp( -a*pow(d1,2) );
            cutoff = 1 - exp_minus_a_d12;
            dcutoff_dd1 = 2*a*d1*exp_minus_a_d12;
        }        
        disso_out =  cutoff * sign_f2_minus_f1 *4*w1*w2*one_over_w_times_w1_plus_w2
            *( 
                    w1 *(1-4*d1)* exp_one_minus_4d1*one_over_cb1 +   
                    w2 *(1-4*d2)* exp_one_minus_4d2*one_over_cb2*( 1 - (h2.fr-h1.fr)*cg.dcb_dfr(h2.fr)*one_over_cb2 ) 
             ) 
            + dcutoff_dd1*sign_f2_minus_f1*one_over_cb1 * 4*w1*w2*one_over_w_times_w1_plus_w2*
            ( w1 * d1 * exp_one_minus_4d1 + w2 * d2 * exp_one_minus_4d2 );  
    } 

    return disso_out;
}


double ddisso_dr(const Harmo& h1, const Harmo& h2, const Generic& cg, double r)
{
    double w1 = Generic::power(h1.db);
    double w2 = Generic::power(h2.db);
    double w = Generic::power(2*cg.DECI_BASE);
    double one_over_w_times_w1_plus_w2 = 1.0/( w*(w1+w2) );
    double d1 = cg.dist_in_cb(h1,h2);
    double d2 = cg.dist_in_cb(h2,h1);
    double one_over_cb1 = 1.0/cg.cb_fr(h1.fr);
    double one_over_cb2 = 1.0/cg.cb_fr(h2.fr);
    double h2fr_minus_h1fr = (h2.fr-h1.fr);
    int sign_f2_minus_f1 = 1;
    if ( (h2.fr-h1.fr) < 0)
    {
        sign_f2_minus_f1 = -1;
    }

    double disso_out=0;

    if( d1<7 )
    {
        disso_out =   sign_f2_minus_f1* 4*w1*w2*one_over_w_times_w1_plus_w2
                *(1.0/r)*h2fr_minus_h1fr*( 
                    w1 *(1-4*d1)* exp(1-4*d1)*one_over_cb1*( 1.0 - cg.dcb_dfr(h1.fr)*h1.fr*one_over_cb1 ) 
                    + w2 *(1-4*d2)* exp(1-4*d2)*one_over_cb2*( 1.0 - cg.dcb_dfr(h2.fr)*h2.fr*one_over_cb2 ) 
                    ); 
    } 

    return disso_out;
}


double disso(const Pitch& p1, const Pitch& p2, const Generic& cg)
{
    double sum_disso=0;

    for (int i=0; i<p1.harmo.size(); i++)
    {
        for (int j=0; j<p2.harmo.size(); j++)
        {
            sum_disso = sum_disso + disso(p1.harmo[i],p2.harmo[j],cg);
        }
    }

    return sum_disso;
}


double ddisso_df2(const Pitch& p1, const Harmo& h2, const Generic& cg)
{
    double sum_disso=0;

    for (int i=0; i<p1.harmo.size(); i++)
    {
        sum_disso = sum_disso + ddisso_df2(p1.harmo[i],h2,cg);
    }

    return sum_disso;
}


double disso(const Chord& c1, const Chord& c2, const Generic& cg)
{
    double sum_disso=0;

    for (int i=0; i<c1.harmo.size(); i++)
    {
        for (int j=0; j<c1.harmo[i].size(); j++)
        {
            for (int k=0; k<c2.harmo.size(); k++)
            {
                for (int l=0; l<c2.harmo[k].size(); l++)
                {
                    double temp = disso(c1.harmo[i][j],c2.harmo[k][l],cg);
                    sum_disso = sum_disso + temp;
                }
            }
        }
    }

    return sum_disso;
}


double ddisso_df2(const Chord& c1, const Harmo& h2, const Generic& cg)
{
    double sum_disso=0;

    for (int i=0; i<c1.harmo.size(); i++)
    {
        for (int j=0; j<c1.harmo[i].size(); j++)
        {
            sum_disso = sum_disso + ddisso_df2(c1.harmo[i][j],h2,cg);
        }
    }

    return sum_disso;
}


double disso(const ChordS& t1, const ChordS& t2, const Generic& cg)
{
    double sum_disso=0;

    for (int i=0; i<t1.harmo.size(); i++)
    {
        for (int j=0; j<t1.harmo[i].size(); j++)
        {
            for (int k=0; k<t1.harmo[i][j].size(); k++) 
            {
                for (int l=0; l<t2.harmo.size(); l++)
                {
                    for (int m=0; m<t2.harmo[l].size(); m++)
                    {
                        for (int n=0; n<t2.harmo[l][m].size(); n++)
                        {
                            sum_disso = sum_disso + disso(t1.harmo[i][j][k],t2.harmo[l][m][n],cg);
                        }
                    }
                }
            }
        }
    }

    return sum_disso;
}


double ddisso_df2(const ChordS& t1, const Harmo& h2, const Generic& cg)
{
    double sum_disso=0;

    for (int i=0; i<t1.harmo.size(); i++)
    {
        for (int j=0; j<t1.harmo[i].size(); j++)
        {
            for (int k=0; k<t1.harmo[i][j].size(); k++)
            {
                sum_disso = sum_disso + ddisso_df2(t1.harmo[i][j][k],h2,cg);
            }
        }
    }

    return sum_disso;
}


double disso(const Scale& s1, const Scale& s2, const Generic& cg)
{
    return disso(s1.chord(cg), s2.chord(cg), cg);
}


double ddisso_df2(const Scale& s1, const Harmo& h2, const Generic& cg)
{
    return ddisso_df2(s1.chord(cg), h2, cg);
}


double i_th(const double x, void* params)
{
    std::tuple<gsl_multimin_function*, int, gsl_vector*>* p = (std::tuple<gsl_multimin_function*, int, gsl_vector*>*) params;

    gsl_vector* x_0_x_i = gsl_vector_alloc( (*(std::get<0>(*p))).n );

    for (int j=0; j<(*(std::get<0>(*p))).n; j++)
    {
        if ( j == std::get<1>(*p) )
        {
            gsl_vector_set(x_0_x_i, j, x );

        } 
        else
        {    
            gsl_vector_set(x_0_x_i, j, gsl_vector_get(std::get<2>(*p),j) );
        }
    }


    return (*(std::get<0>(*p))).f( x_0_x_i, (*(std::get<0>(*p))).params );
}


/*
   gsl_function gsl_i_th (gsl_multimin_function* gsl_f, int i, const gsl_vector* x_0)
   {
   std::cout << " gsl_func called ";
   cin.get();

   gsl_function gsl_f_i;
   gsl_function* gsl_f_i_ptr = &gsl_f_i;

   std::tuple<gsl_multimin_function*, int, const gsl_vector*> params_f_i ( gsl_f, i, x_0 );
   std::tuple<gsl_multimin_function*, int, const gsl_vector*>* params_f_i_ptr = &params_f_i;

   std::cout << " c ";
   cin.get();

   gsl_f_i.function = &i_th;
   gsl_f_i.params = (void *) params_f_i_ptr;

   return gsl_f_i;
   }
*/


void gsl_numerical_gradient(gsl_multimin_function* gsl_f, const gsl_vector* x_0, gsl_vector* h, gsl_vector* grad, gsl_vector* abserr)
{
    for (int i=0; i < (*x_0).size; i++)
    {
        double* grad_i = gsl_vector_ptr(grad,i);
        double* abserr_i = gsl_vector_ptr(abserr,i);

        gsl_function gsl_f_i;
        gsl_function* gsl_f_i_ptr = &gsl_f_i; 

        std::tuple<gsl_multimin_function*, int, const gsl_vector*> params_f_i ( gsl_f, i, x_0 );
        std::tuple<gsl_multimin_function*, int, const gsl_vector*>* params_f_i_ptr = &params_f_i;

        gsl_f_i.function = &i_th;
        gsl_f_i.params = (void *) params_f_i_ptr;

        gsl_deriv_central(gsl_f_i_ptr, gsl_vector_get(x_0,i), gsl_vector_get(h,i),  grad_i, abserr_i); 

        gsl_vector_set(grad, i, *grad_i);
    }
}


double total_disso(const gsl_vector* r, void* params)
{
    std::tuple<double,int,int,double,Generic,gsl_multimin_function*>* p = (std::tuple<double,int,int,double,Generic,gsl_multimin_function*>*) params;

    std::vector<double> std_r;
    std_r.push_back(1);
    for (int i=0; i<(*r).size; i++)
    {
        std_r.push_back( gsl_vector_get(r, i) ); 
    }  
    std_r.push_back(std::get<3>(*p));
    Scale s( std::get<0>(*p) , std::get<1>(*p), std::get<2>(*p), std_r);

    double penalty = 0;
    for (int i=0; i<std_r.size()-1; i++)
    {
        for (int j=i+1; j<std_r.size(); j++)
        {
            if ( std_r[i]-std_r[j] > 0 )
            {
                penalty = penalty + 0.5*1e4*pow(std_r[i]-std_r[j],2);
            }
        }
    }

    return disso(s, s, std::get<4>(*p))  + penalty;
}


void dtotal_disso_dr(const gsl_vector* r, void* params, gsl_vector* grad)
{
    std::tuple<double,int,int,double,Generic,gsl_multimin_function*>* p = (std::tuple<double,int,int,double,Generic,gsl_multimin_function*>*) params;

    std::vector<double> std_r;
    std_r.push_back(1);
    for (int i=0; i<(*r).size; i++)
    {
        std_r.push_back( gsl_vector_get(r, i) ); 
    }  
    std_r.push_back(std::get<3>(*p));
    Scale s( std::get<0>(*p), std::get<1>(*p), std::get<2>(*p), std_r);

    for (int i=1; i<((*r).size+1); i++)
    {
        double dtotal_dr=0;
        Chord s_minus_i = s.chord( std::get<4>(*p));

        for (int k=0; k< s.n_octave*s.n_pitch_octave; k++)
        {
            for (int l=0; l< k; l++)
            {
                dtotal_dr = dtotal_dr + ddisso_dr( s.chord( std::get<4>(*p)).harmo[l][i], s.chord( std::get<4>(*p)).harmo[k][i], std::get<4>(*p), std_r[i] );
            }

            s_minus_i.harmo[k][i].db = -1e10;
        } 


        for (int k=0; k< s.n_octave*s.n_pitch_octave; k++)
        {
            dtotal_dr = dtotal_dr + ddisso_df2( s_minus_i, s.chord( std::get<4>(*p)).harmo[k][i], std::get<4>(*p) ) * s.chord( std::get<4>(*p) ).harmo[k][0].fr;
        }

        double dpenalty_dr = 0;
        for (int j=i+1; j<std_r.size(); j++)
        {
            if ( std_r[i]-std_r[j] > 0 )
            {
                dpenalty_dr = dpenalty_dr + 1e4*(std_r[i]-std_r[j]);
            }
        }
        for (int j=0; j<i; j++)
        {
            if ( std_r[j]-std_r[i] > 0 )
            {
                dpenalty_dr = dpenalty_dr - 1e4*(std_r[j]-std_r[i]);
            }
        }

        dtotal_dr  = 2*dtotal_dr + dpenalty_dr;

        gsl_vector_set(grad, i-1, dtotal_dr);
    }


    //display
    *std::get<4>(*p).logfile << std::fixed;
    *std::get<4>(*p).logfile << std::setprecision(9);

    *std::get<4>(*p).logfile << "df gives ";
    for (int i=0; i < (*grad).size; i++)
    {
        *std::get<4>(*p).logfile << gsl_vector_get(grad,i) << " ";
    }
    *std::get<4>(*p).logfile << " at ";
    for (int i=0; i<(*r).size; i++)
    {
        *std::get<4>(*p).logfile << gsl_vector_get(r,i) << " ";
    }
    *std::get<4>(*p).logfile << std::endl;
}


void dtotal_disso_dr_num(const gsl_vector* r, void* params, gsl_vector* grad)
{
    std::tuple<double,int,int,double,Generic,gsl_multimin_function*>* p = (std::tuple<double,int,int,double,Generic,gsl_multimin_function*>*) params;

    gsl_vector* abserr = gsl_vector_alloc((*r).size);
    gsl_vector* h = gsl_vector_alloc((*r).size);


    for (int i=0; i<(*h).size; i++)
    {
        gsl_vector_set(h,i,1e-7);
    }

    gsl_numerical_gradient(std::get<5>(*p), r, h, grad, abserr);

    //display
    *std::get<4>(*p).logfile << std::fixed;
    *std::get<4>(*p).logfile << std::setprecision(9);

    *std::get<4>(*p).logfile << "df num gives ";
    for (int i=0; i < (*grad).size; i++)
    {
        *std::get<4>(*p).logfile << gsl_vector_get(grad,i) << " "; //<< " pm " << gsl_vector_get(abserr,i) << " ";
    }
    *std::get<4>(*p).logfile << " at ";
    for (int i=0; i<(*r).size; i++)
    {
        *std::get<4>(*p).logfile << gsl_vector_get(r,i) << " ";
    }
    *std::get<4>(*p).logfile << std::endl;
}


void total_disso_dtotal_disso_dr(const gsl_vector* r, void* params, double *f, gsl_vector* grad)
{
    *f = total_disso(r, params);
    dtotal_disso_dr(r, params, grad);
}


void total_disso_dtotal_disso_dr_num(const gsl_vector* r, void* params, double* f, gsl_vector* grad)
{
    *f = total_disso(r, params);
    dtotal_disso_dr_num(r, params, grad);
}


void plomp_levelt(const std::string& o_f_n)
{
    const int N_FR=200000;
    const double DB_BAS=60.0;
    Generic gene(N_FR, DB_BAS); 

    gene.output_file_name = o_f_n;

    std::ofstream logfile;
    gene.logfile = &logfile; 
    gene.logfile_init();
    *gene.logfile << "logfile initiated \n"; //"f1 0 0 1 \"bass_trans_cycl.wav\" 0 0 0\n";

    std::ofstream scorefile;
    gene.scorefile = &scorefile; 
    gene.scorefile_init();
    *gene.scorefile << "f1 0 16384 10 1 0   0 0    0 0    0  \n"; //"f1 0 0 1 \"bass_trans_cycl.wav\" 0 0 0\n";
    *gene.scorefile << "f2 0 16384 10 1 0   0 0    0 0    0  \n"; //"f2 0 0 1 \"bass_trans_cycr.wav\" 0 0 0\n";
    *gene.scorefile << "\n";

    const int N_MINS=5;
    const int ITER_MAX=150;
    const int T_NET=1;
    int BASE_PITCH=200;
    int N_OCTAVE=4;
    int N_PITCH_OCTAVE=8;
    int N_HARMO_PITCH=6;
    double MAX_HARMO=8.0;
    Score score(2*N_MINS-1, 1, N_OCTAVE*N_PITCH_OCTAVE, N_HARMO_PITCH);

        // setting up the gsl_multimin_function to be optimized by gsl
    gsl_multimin_function gsl_total_disso;
    gsl_multimin_function* gsl_total_disso_ptr = &gsl_total_disso;
    std::tuple<double, int, int, double, Generic, gsl_multimin_function*> params (BASE_PITCH, N_OCTAVE, N_PITCH_OCTAVE, MAX_HARMO, gene, gsl_total_disso_ptr);
    std::tuple<double, int, int, double, Generic, gsl_multimin_function*>* params_ptr = &params;
    gsl_total_disso.n = N_HARMO_PITCH-2;
    gsl_total_disso.f = &total_disso;
    gsl_total_disso.params = (void *) params_ptr;

        //initializing Nelder-Mead simplex method
    /*
       const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex2rand;
       gsl_multimin_fminimizer* total_disso_minimizer =  gsl_multimin_fminimizer_alloc(T, N_HARMO_PITCH-2);


       gsl_total_disso.n = N_HARMO_PITCH-2;
       gsl_total_disso.f = &total_disso;
       gsl_total_disso.params = (void *) params_ptr;
    */

        //initializing gradient descent algorithms 
    const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_fr;
    //const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_steepest_descent;
    gsl_multimin_fdfminimizer* total_disso_minimizer =  gsl_multimin_fdfminimizer_alloc(T, N_HARMO_PITCH-2);

    gsl_multimin_function_fdf gsl_total_disso_fdf;
    gsl_multimin_function_fdf* gsl_total_disso_fdf_ptr = &gsl_total_disso_fdf;  
    gsl_total_disso_fdf.n = N_HARMO_PITCH-2;
    gsl_total_disso_fdf.f = &total_disso;
    gsl_total_disso_fdf.df = &dtotal_disso_dr;
    gsl_total_disso_fdf.fdf = &total_disso_dtotal_disso_dr;
    //gsl_total_disso_fdf.df = &dtotal_disso_dr_num; //test against numerical derivatives
    //gsl_total_disso_fdf.fdf = &total_disso_dtotal_disso_dr_num; //test against numerical derivatives
    gsl_total_disso_fdf.params = (void *) params_ptr;
    double step_size = 0.1; //just initial
    double tol = 0.1; //linear optimization tolerance parameter
    double GRAD_NORM_MIN = 1; //norm of gradient staus parameter

    std::vector<std::vector<double>> mins_list;
    for (int n=0; n<N_MINS; n++)
    {
        std::cout << "------- n = " << n << " ---------" << std::endl;

            //either start at a chosen set of ratios
        //std::vector<double> harmo_ratios = {1.0,3.1,3,4.1,5,6.1,8.05,8.0};
        //std::vector<double> harmo_ratios = {1,2.37841423,4,5.656854249,6.727171322,8};
            //or choose ratios randomly
        std::vector<double> harmo_ratios = {1,MAX_HARMO};
        for (int i=0; i<N_HARMO_PITCH-2; i++)
        {
            double r = double(rand())/RAND_MAX;

            harmo_ratios.push_back( 1 + r*(MAX_HARMO-1) );
        } 
        sort(harmo_ratios.begin(),harmo_ratios.end());

        Scale scale(BASE_PITCH, N_OCTAVE, N_PITCH_OCTAVE, harmo_ratios);  
        //score.harmo[0][0] = scale.chord(gene).harmo;

        std::cout << "initial disso : " << disso(scale,scale,gene) << std::endl; 
        *gene.logfile << "initial disso : " << disso(scale,scale,gene) << std::endl; 

        std::cout << "at ";
        *gene.logfile << "at ";
        for (int i=0; i < harmo_ratios.size(); i++)
        {
            std::cout << harmo_ratios[i] << " ";
            *gene.logfile << harmo_ratios[i] << " ";
        }
        std::cout << std::endl;
        *gene.logfile << std::endl;

        gsl_vector* initial_ratios = gsl_vector_alloc(N_HARMO_PITCH-2);
        //gsl_vector* initial_step = gsl_vector_alloc(N_HARMO_PITCH-2); //in case other algos

        for (int i=0; i<N_HARMO_PITCH-2; i++)
        {
            gsl_vector_set(initial_ratios, i, harmo_ratios[i+1]); 
            //gsl_vector_set(initial_step, i, 0.03*harmo_ratios[i+1]);  //in case other algos
        }

        //gsl_multimin_fminimizer_set(total_disso_minimizer, gsl_total_disso_ptr, initial_ratios, initial_step); //in case other algos
        gsl_multimin_fdfminimizer_set(total_disso_minimizer, gsl_total_disso_fdf_ptr, initial_ratios, step_size, tol);


        int status;
        int iter=0;
        //iterations of the minimizing method
        do
        {
            //status = gsl_multimin_fminimizer_iterate(total_disso_minimizer); //in case other algos
            status = gsl_multimin_fdfminimizer_iterate(total_disso_minimizer);

            //std::cout << "after iter " << iter << ", disso = " << gsl_multimin_fdfminimizer_minimum(total_disso_minimizer ) << std::endl;
            *gene.logfile << "after iter " << iter << ", disso = " << gsl_multimin_fdfminimizer_minimum(total_disso_minimizer ) << std::endl;

            std::vector<double> new_harmo_ratios;
            std::vector<double> new_gradient;

            new_harmo_ratios.push_back(1);
            new_gradient.push_back(0);

            for (int i=0; i < total_disso_minimizer->x->size; i++)
            {
                new_harmo_ratios.push_back( gsl_vector_get((*total_disso_minimizer).x, i) );
                new_gradient.push_back( gsl_vector_get((*total_disso_minimizer).gradient, i) );
            }
            new_harmo_ratios.push_back(harmo_ratios[N_HARMO_PITCH-1]);
            new_gradient.push_back(0);

            //std::cout << "new ratios : ";
            *gene.logfile << " ratios = ";
            for (int i=0; i < new_harmo_ratios.size(); i++)
            {
                //std::cout << new_harmo_ratios[i] << " ";
                *gene.logfile << new_harmo_ratios[i] << " ";
            }
            //std::cout << std::endl;
            *gene.logfile << std::endl;

            //std::cout << "-gradient = "; 
            *gene.logfile << "-gradient = "; 
            for (int i=0; i < new_gradient.size(); i++)
            {
                //std::cout << -new_gradient[i] << " ";
                *gene.logfile << -new_gradient[i] << " ";
            }
            //std::cout << std::endl;
            *gene.logfile << std::endl;

            //status = gsl_multimin_test_size(size, 1e-2); //in case other algo
            status = gsl_multimin_test_gradient((*total_disso_minimizer).gradient, GRAD_NORM_MIN);

            if((status == GSL_SUCCESS) || (iter==ITER_MAX) )
            {
                std::cout << "stops after iter " << iter << ", disso = " << gsl_multimin_fdfminimizer_minimum(total_disso_minimizer ) << std::endl;

                if (status==GSL_SUCCESS)
                {
                    *gene.logfile << "minimum found\n";
                    std::cout << "MINIMUM FOUND\n";
                    mins_list.push_back( new_harmo_ratios );
                }
                if (iter==ITER_MAX)
                {
                    *gene.logfile << "max iterations reached\n";
                    std::cout << "max iterations reached\n";
                }

                std::cout << "-gradient is : "; 
                for (int i=0; i < new_gradient.size(); i++)
                {
                    std::cout << -new_gradient[i] << " ";
                }
                std::cout << std::endl;
                std::cout << "at ";
                for (int i=0; i < new_harmo_ratios.size(); i++)
                {
                    std::cout << new_harmo_ratios[i] << " ";
                }
                std::cout << std::endl;


                scale.harmo_ratios = new_harmo_ratios;
                score.harmo[2*n][0] = scale.chord(gene).harmo;

                break;
            }

            iter++;
        }
        while(status == GSL_CONTINUE);
    }

    *gene.logfile << "mins_list : ";
    for (int i=0; i < mins_list.size(); i++)
    { 
        for (int j=0; j < mins_list[i].size(); j++)
        {
            *gene.logfile << mins_list[i][j] << " ";
        }
        *gene.logfile << std::endl;
    }

    gsl_multimin_fdfminimizer_free(total_disso_minimizer);

    for (int i=0; i<score.time.size(); i++)
    {
        score.time[i] = 2+4*i; 
    }

    score.print(gene, T_NET);
    *gene.logfile << "executed" << std::endl;

    (*gene.scorefile).close();
    (*gene.logfile).close();

}
