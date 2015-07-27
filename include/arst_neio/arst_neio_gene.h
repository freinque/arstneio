/** @file */ 
#ifndef HEADER_CSOUND_GENE_H_INCLUDED
#define HEADER_CSOUND_GENE_H_INCLUDED

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

struct Harmo;

/// Contains the basic physiological functions and log/scorefile initializations
struct Generic 
{
    int NUMB_FREQ; ///< Number of freq spectrum subdivisions
    double DECI_BASE; ///< db level of highest Harmo
    //double STEP;
    std::vector<double> fr; ///< Discrete freq spectrum
    std::vector<double> cb; ///< cb values of fr
    std::vector<double> a_we; ///< a_weighting factors of fr
    std::vector<double> a_we_db; ///< a_weighting db factors of fr
    std::vector<double> b_we; ///< b_weighting factors of fr
    std::vector<double> b_we_db; ///< b_weighting db factors of fr
    std::vector<double> ifr_tabl; ///< Discrete inverse of fr

    Generic(int n_fr, double db_bas); ///< Initializes a Generic with NUMB_FREQ=n_fr and DECI_BASE=db_bas 
    ~Generic();
    
    static double power(double db); ///< Returns the power of a db level
    static int randd(std::vector<int> a, std::vector<double> b); ///< Returns a random element of a using weights in b 

    std::vector<double> fr_init(int n_f, double m_f, double s);
    std::vector<double> cb_init(std::vector<double> f);
    std::vector<double> a_we_init(std::vector<double> f);
    std::vector<double> a_we_db_init(std::vector<double> a_w);
    std::vector<double> b_we_init(std::vector<double> sf);
    std::vector<double> b_we_db_init(std::vector<double> b_w);
    std::vector<double> ifr_tabl_init();

    double cb_fr(double fr) const; ///< Returns the critical bandwidth at a frequency
    double dcb_dfr(double fr) const; ///< Returns the derivative of the critical bandwidth at a frequency
    //double dist_in_cb(int ifr_refe, int ifr) const;
    double dist_in_cb(Harmo h1, Harmo h2) const; ///< Returns the distance between two Harmo in terms of cb
    int ifr(double frq) const; ///< Returns projection of frq on velues of fr
    int ifr_fast(double f) const; ///< Rougher ifr

    std::string output_file_name; ///< File name prefix used for logfile and scorefile
    std::ofstream* logfile; 
    std::ofstream* scorefile;

    void scorefile_init();
    void logfile_init();    
private:
    int INVE_NET;
    double NET;
};



/// Contains the information contained in a pure sine signal
struct Harmo
{
    Harmo(int int_freq=0, double freq=0, double deci=0); ///< Initializes an empty Harmo
    ~Harmo();
        
    int f; ///< Discrete frequency
    double fr; ///< Double frequency
    double db; ///< power in db

    static bool equa(const Harmo& a, const Harmo& b); ///< Returns wether or not Harmo a and b are equal
};


/// Contains a vector of Harmo
struct Pitch
{
    Pitch(int n_harmo_pitch=0); ///< Initializes an empty pitch
    ~Pitch();
    
    std::vector<Harmo> harmo; ///< The Harmo content of the Pitch
};

/// Contains a vector of Pitch
struct Chord
{
    Chord(int n_pitch_chord=0, int n_harmo_pict=0); ///< Initializes an empty chord
    ~Chord();
        
    std::vector<std::vector<Harmo>> harmo; ///< The harmo content of the Chord
    
    Pitch pitch(int i); ///< Returns the ith pitch of the Chord
};



/// Contains a vector of Chord
struct ChordS
{
    ChordS(int n_chord=0, int n_pitch_chord=0, int n_harmo_pict=0); ///< Initializes an empty ChordS
    ~ChordS();
        
    std::vector<std::vector<std::vector<Harmo>>> harmo; ///< The Harmo content of the ChordS
   
    Chord chord(int i); ///< Returns the ith Chord of the ChordS
    Pitch pitch(int i,int j); ///< Returns the jth pitch of the ith Chord of the ChordS
};


/// Contains a vector of ChordS
struct Score 
{
    int T_FIN;
    int N_CHOR;
    int N_PITC_CHOR;
    int N_HARM_PITC;
    
    Score(int t_fin=0, int n_chord=0, int n_pitch_chord=0, int n_harmo_pict=0); ///< Initializes an empty Score
    ~Score();
    
    std::vector<std::vector<std::vector<std::vector<Harmo>>>> harmo; ///< The Harmo content of the Score
    std::vector<double> time; ///< The transition times of the Score
  
    ChordS t(int i); ///< ChordS at ith transition time
    Chord chord(int i,int j); ///< jth Chord or ChordS at ith transition time
    Pitch pitch(int i,int j, int k);///< kth Pitch of jth Chord or ChordS at ith transition time

    void print(const Generic& cg, int t_net=1); ///< Writes the t_net transition multiples Score in CSound .sco file
};

/// Contains a equally tempered scale of n_octaves octaves, n_pitch_octaves Pitch per octave with harmo_ratios Harmo ratios 
struct Scale
{
    Scale(double b_p, int n_o, int n_p_o, std::vector<double> r);
    ~Scale();
        
    int base_pitch; ///< db level of main Harmo
    int n_octave; ///< Number of octaves
    int n_pitch_octave; ///< Number of Pitch per octave
    std::vector<double> harmo_ratios; ///< Freq ratios of the Harmo of every Pitch

    Chord chord(const Generic& cg) const; ///< Returns the Scale as a Chord
};

#endif
