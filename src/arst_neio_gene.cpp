/** @file */
#include "arst_neio/arst_neio_gene.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

Generic::Generic(int n_fr, double db_bas)
{
    NUMB_FREQ = n_fr;
    DECI_BASE = db_bas;
    INVE_NET=4;
    NET = 0.25;
    
    double MINI_FREQ = 20;
    double STEP = 0.0002442;
    
    fr = fr_init(NUMB_FREQ, MINI_FREQ, STEP);

    cb = cb_init(fr);
    a_we = a_we_init(fr);
    a_we_db = a_we_db_init(a_we);
    
    b_we = b_we_init(fr);
    b_we_db = b_we_db_init(b_we);
    
    ifr_tabl = ifr_tabl_init();

}


double Generic::cb_fr(double fr) const
{
    return 75*pow(1 + 1.4*pow(fr /1000,2),0.69) +25;
}

 
double Generic::dcb_dfr(double fr) const
{
    return (0.0001449*fr)/pow(1.4e-6*pow(fr,2)+1,0.31);
    //return 0.0001449*fr - (6.28866e-11)*pow(fr,3.0) + (5.7667e-17)*pow(fr,5.0);
}


double Generic::power(double db)
{
    return pow(10.0,db/10.0);
}


std::vector<double> Generic::fr_init(int n_f, double m_f, double s)
{	
    std::vector<double> f;
    
	f.push_back( m_f );
    // naivement qu'on fait le pas en fr comme cb
	for( int i=1 ; i<n_f ; i++)
	{
		f.push_back( f.back() + s * cb_fr( f.back() ) );// (75*pow(1 + 1.4*pow(f[i-1]/1000,2),0.69) +25);
	}

	return f;
}


std::vector<double> Generic::cb_init(std::vector<double> f)
{	
    std::vector<double> c;
    
	// on assume naivement que la resolution en freq va comme 
	for( int i=0 ; i<f.size() ; i++)
	{
	    c.push_back( cb_fr( f[i] ) ); //75*pow(1 + 1.4*pow(f[i] /1000,2),0.69) +25;
	}
	return c;
}


std::vector<double> Generic::a_we_init(std::vector<double> f)
{
	std::vector<double> a_w;
	
	// A-weighting
	for( int i=0 ; i<f.size() ; i++)
	{
	a_w.push_back( (pow(12200.0,2.0) * pow(f[i],4.0))/((pow(f[i],2.0)+pow(20.6,2.0))*(pow(f[i],2.0)+pow(12200.0,2.0))*sqrt((pow(f[i],2.0)+pow(107.7,2.0))*(pow(f[i],2.0)+pow(737.9,2.0)))) );
	}
	
	return a_w;
}


std::vector<double> Generic::a_we_db_init(std::vector<double> a_w)
{
	std::vector<double> a_w_d;
	
	for( int i=0 ; i<a_w.size() ; i++)
	{
	a_w_d.push_back( 2 + 20*log10(a_w[i]) );
	}
	
	return a_w_d;
}


std::vector<double> Generic::b_we_init(std::vector<double> f)
{
	std::vector<double> b_w;

	// B-weighting
	for( int i=0 ; i<f.size() ; i++)
	{
	b_w.push_back( (pow(12200,2.0) * pow(f[i],3.0))/((pow(f[i],2.0)+pow(20.6,2.0))*(pow(f[i],2.0)+pow(12200,2.0))*sqrt(pow(f[i],2.0)+pow(158.5,2.0))) );
	}
    return b_w;	
}


std::vector<double> Generic::b_we_db_init(std::vector<double> b_w)
{
	std::vector<double> b_w_d;
	
	for( int i=0 ; i<b_w.size() ; i++)
	{
	b_w_d.push_back( 0.17 + 20*log10(b_w[i]) );
	}
	
	return b_w_d;
}


int Generic::ifr(double freq) const
{
	int ifr_out;

	if ( fr[NUMB_FREQ-1] < freq )
	{
	ifr_out = NUMB_FREQ;
	}
	else
	{
	for( int i= 0; i< NUMB_FREQ; i++)
	{
		if  ( fr[i] >= freq )
			{
			ifr_out = i;
			break;
			}		
	}
	}
	return ifr_out;
}

// to be changed
std::vector<double> Generic::ifr_tabl_init()
{
    std::vector<double> i_t(INVE_NET*20000,99999);
    
	i_t[0] = 0;

	for (int i=1; i<i_t.size(); i++)
	{
	for (int j= i_t[i-1]; j<NUMB_FREQ; j++)
	{
	if  (fr[j] >= NET*i)
	{	
		i_t[i] = j;
		break;
	}
	}
	}
	return i_t;
}


int Generic::ifr_fast(double f) const
{
    int inde_f;
    inde_f = (int)(std::floor(INVE_NET*f+0.5));
    
    return ifr_tabl[inde_f];
}


double Generic::dist_in_cb(Harmo h1, Harmo h2) const
{
    return std::fabs(h2.fr-h1.fr)/cb_fr(h1.fr);
}


int Generic::randd(std::vector<int> a, std::vector<double> b)
{
	if (b.size() != a.size()) {std::cout << "b.size() != a.size() in randd()" << std::endl;}
	
	int randd_out;
	double randd_num;

	// sort le poids total et choisit un point [0,total]
	double sum_b = 0;

	double part_sum_b[b.size()];
	for (int i=0; i<b.size() ; i++)
	{
		sum_b = sum_b + b[i];
		part_sum_b[i] = sum_b;
	}
	if (sum_b == 0) {std::cout << "sum_b=0" << std::endl;}

	randd_num = (double(rand()) / RAND_MAX)*sum_b;
	

	// determine dans quelle section du poids on a tire
	for (int i=0 ;  i<b.size() ; i++)
	{
		if  ( ((sum_b > 0) - (sum_b < 0))*part_sum_b[i] > ((sum_b > 0) - (sum_b < 0))*randd_num )
		{
			randd_out = a[i];
			break;
		}

	}

	return randd_out;
}


void Generic::logfile_init()
{
    std::string log = ".log";
    std::string an = "arst_neio_";
	std::string output_file_name_log = an + output_file_name + log;
    std::cout << "writing log to " << output_file_name_log << std::endl;
    (*logfile).open(output_file_name_log.c_str()); // opens a log file
}


void Generic::scorefile_init()
{
    std::string sco = ".sco";
	std::string an = "arst_neio_";
    std::string output_file_name_sco = an + output_file_name + sco;
    std::cout << "writing score to " << output_file_name_sco << std::endl;
	(*scorefile).open(output_file_name_sco.c_str()); // opens a score file
}


Generic::~Generic()
{    
}

////////////////////////////////////////////////////////////

Harmo::Harmo(int int_freq, double freq, double deci)
{
    f = int_freq;
    fr = freq;
    db = deci;
}


Harmo::~Harmo()
{
}


bool Harmo::equa(const Harmo& a, const Harmo& b)
{
    return (a.f == b.f) && (a.fr == b.fr) && (a.db == b.db);
}
 
//////////////////////////////////////////////////////////////

Pitch::Pitch(int n_harmo_pitch)
{
    harmo = std::vector<Harmo> (n_harmo_pitch);
}


Pitch::~Pitch()
{
}

//////////////////////////////////////////////////////////////

Chord::Chord(int n_pitch_chord, int n_harmo_pitch)
{
    harmo = std::vector<std::vector<Harmo>> (n_pitch_chord, std::vector<Harmo>(n_harmo_pitch) );
}


Chord::~Chord()
{
}


Pitch Chord::pitch(int i)
{
    Pitch p;
    try
    {
    p.harmo = harmo[i];
    }
    catch (const std::bad_alloc& oor) 
    {
    std::cerr << "error Chord::pitch: " << oor.what() << std::endl;
    }
    catch (const std::out_of_range& oor) 
    {
    std::cerr << "error Chord::pitch: " << oor.what() << std::endl;
    }
    catch (...) 
    {
    std::cerr << "unexpected error Chord::pitch: " << std::endl;
    }
    return p;
}

//////////////////////////////////////////////////////////////

ChordS::ChordS(int n_chord, int n_pitch_chord, int n_harmo_pitch)
{
    harmo = std::vector< std::vector< std::vector<Harmo>>> (n_chord, std::vector< std::vector<Harmo>> (n_pitch_chord, std::vector<Harmo> (n_harmo_pitch)));
}


ChordS::~ChordS()
{
}


Chord ChordS::chord(int i)
{
    Chord c;
    try
    {
    c.harmo = harmo[i];
    }
    catch (const std::bad_alloc& oor) 
    {
    std::cerr << "error ChordS::chord: " << oor.what() << std::endl;
    }
    catch (const std::out_of_range& oor) 
    {
    std::cerr << "error ChordS::chord: " << oor.what() << std::endl;
    }
    catch (...) 
    {
    std::cerr << "unexpected error ChordS::chord: " << std::endl;
    }
    return c;
}


Pitch ChordS::pitch(int i,int j)
{
    Pitch p; 
    try
    {
    p.harmo = harmo[i][j];
    }
    catch (const std::bad_alloc& oor) 
    {
    std::cerr << "error ChordS::pitch: " << oor.what() << std::endl;
    }
    catch (const std::out_of_range& oor) 
    {
    std::cerr << "error ChordS::pitch: " << oor.what() << std::endl;
    }
    catch (...) 
    {
    std::cerr << "unexpected error ChordS::pitch: " << std::endl;
    }
    return p;
}

//////////////////////////////////////////////////////////////

Score::Score(int t_fin, int n_chord, int n_pitch_chord, int n_harmo_pitch)
{
    harmo = std::vector<std::vector<std::vector<std::vector<Harmo>>>> ( t_fin, std::vector<std::vector<std::vector<Harmo>>> ( n_chord, std::vector<std::vector<Harmo>> ( n_pitch_chord, std::vector<Harmo> (n_harmo_pitch))));
    time = std::vector<double> (t_fin+1);
}


Score::~Score()
{
}


ChordS Score::t(int i)
{
    ChordS t;
    try
    {
    t.harmo = harmo[i];
    }
    catch (const std::bad_alloc& oor) 
    {
    std::cerr << "error Score::t: " << oor.what() << std::endl;
    }
    catch (const std::out_of_range& oor) 
    {
    std::cerr << "error Score::t: " << oor.what() << std::endl;
    }
    catch (...) 
    {
    std::cerr << "unexpected error Score::t: " << std::endl;
    }
    return t;
}


Chord Score::chord(int i,int j)
{
    Chord c;
    try
    {
    c.harmo = harmo[i][j];
    }
    catch (const std::bad_alloc& oor) 
    {
    std::cerr << "error Score::chord: " << oor.what() << std::endl;
    }
    catch (const std::out_of_range& oor) 
    {
    std::cerr << "error Score::chord: " << oor.what() << std::endl;
    }
    catch (...) 
    {
    std::cerr << "unexpected error Score::chord: " << std::endl;
    }
    return c;
}


Pitch Score::pitch(int i,int j,int k)
{
    Pitch p; 
    try
    {
    p.harmo = harmo[i][j][k];
    }
    catch (const std::bad_alloc& oor) 
    {
    std::cerr << "error Score::pitch: " << oor.what() << std::endl;
    }
    catch (const std::out_of_range& oor) 
    {
    std::cerr << "error Score::pitch: " << oor.what() << std::endl;
    }
    catch (...) 
    {
    std::cerr << "unexpected error Score::pitch: " << std::endl;
    }
    return p;
}


void Score::print(const Generic& cg, int t_net)
{
    int t_fin = harmo.size();    
    int n_chord = harmo[0].size();
    int n_pitch_chord= harmo[0][0].size();
	int n_harmo_pitch = harmo[0][0][0].size();
    //for (int i = 0; i< harmo.size(); i++)
    //{
    //n_chord = std::max(harmo[i].size(), n_chord);
    //for (int j = 0; j< harmo[i].size(); j++)
    //{
    //n_pitch_chord = std::max(harmo[i][j].size(), n_pitch_chord);
    //for (int k = 0; k< t.harmo[i].size(); k++)
    //{
    //n_harmo_pitch = std::max(harmo[i][j][k].size(), n_harmo_pitch);
    //}
    //}
    //}
	
    //switching from the Score to a rescaled resc Score
    int t_fin_over_t_net = int( t_fin/t_net );
    ChordS zero_t(n_chord,n_pitch_chord,n_harmo_pitch);
    Score resc;
    
    resc.harmo.push_back(zero_t.harmo);
    resc.time.push_back(0);
    for (int m = 0; m < t_fin_over_t_net; m++)
	{
	resc.harmo.push_back( harmo[m*t_net] );
	resc.time.push_back( time[m*t_net] );
	}
    resc.harmo.push_back(zero_t.harmo);
    resc.time.push_back(time.back());
    
    std::vector<std::vector<std::vector<std::vector<int>>>> harmo_dura(t_fin_over_t_net, std::vector<std::vector<std::vector<int>>> ( n_chord, std::vector<std::vector<int>> ( n_pitch_chord, std::vector<int> (n_harmo_pitch,0))));
	
	//determines the durations of the Harmo 
    for (int j=0; j<n_chord; j++) {
	for (int k=0; k<n_pitch_chord; k++) {
	for (int l=0; l<n_harmo_pitch; l++) {
	for (int i=1; i<t_fin_over_t_net+1; i++) {
		if (Harmo::equa(resc.harmo[i-1][j][k][l],resc.harmo[i][j][k][l])) {continue;}
		else
		{
        int d = 1;
        while ( Harmo::equa(resc.harmo[i-1+d][j][k][l],resc.harmo[i+d][j][k][l]) )
        {
            d = d+1;
            if (i+d>=t_fin_over_t_net+1) {break;}
        }
		harmo_dura[i-1][j][k][l] = d;
        }
	}
	}
    }
	}

	//writes Harmo to score file separately
    for (int j = 0; j<n_chord; j++) {
	for (int k = 0; k<n_pitch_chord; k++) {
	for (int l = 0; l<n_harmo_pitch; l++) {
	for (int i = 0; i<t_fin_over_t_net; i++) {

	if (harmo_dura[i][j][k][l]==0) {continue;}
	else
	{
		*cg.scorefile << "i"<< j+1 << " " << time[i] +  (k/double(resc.harmo[i+1][j].size()))*(time[i+harmo_dura[i][j][k][l]] - time[i]) << " " << 0.4*(time[i+harmo_dura[i][j][k][l]] - time[i]) << " " << resc.harmo[i+1][j][k][l].fr << " " << resc.harmo[i+1][j][k][l].db << "\n"; //- cg.a_we_db[resc.harmo[i+1][j][k][l].f]
	}
    
    }}}}
}

//////////////////////////////////////////////////////////////

Scale::Scale(double b_p, int n_o, int n_p_o, std::vector<double> r)
{
    base_pitch = b_p;
    n_octave = n_o;
    n_pitch_octave = n_p_o;
    harmo_ratios = r;
}


Scale::~Scale()
{
}


Chord Scale::chord(const Generic& cg) const
{   
    int n_pitch_chord = n_octave*n_pitch_octave;
    Chord c(n_pitch_chord,harmo_ratios.size());
    
    for (int i=0; i<c.harmo.size(); i++) 
    {
        //c.harmo[i][0].f = cg.ifr( base_pitch + base_pitch*(double(i)/n_pitch_octave));
        c.harmo[i][0].fr = base_pitch + base_pitch*(double(i)/n_pitch_octave);
        c.harmo[i][0].db = cg.DECI_BASE;

        double rate = 1.05; //just a constant for db decrease of Harmo of a Pitch
        for (int j=1; j<harmo_ratios.size(); j++)
        {
        //c.harmo[i][j].f = cg.ifr( std::min(harmo_ratios[j]*c.harmo[i][0].fr(cg)  + 0.5 , 19999.0 ) );
        c.harmo[i][j].fr = std::min( harmo_ratios[j]*c.harmo[i][0].fr, 19999.0 );
        c.harmo[i][j].db = cg.DECI_BASE - 10*j*log(rate)/(log(10.0));
        }
	}

    return c;
}


